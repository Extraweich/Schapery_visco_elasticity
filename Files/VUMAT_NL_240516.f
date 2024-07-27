c ======================================================================
c User Subroutine UMAT for non-linear Schapery model with 3 Kelvin-Voigt elements
c By Nicolas Christ (Karlsruhe Institute of Technology) - May 2024
c ======================================================================
!DEC$ FREEFORM

MODULE FUNCTIONS
CONTAINS

    FUNCTION HYDROSTATIC(EPS, ntens) RESULT (HYDRO)
	  ! Calculate the hydrostatic strain/stress without changing factors (covariant/contravariant)
	  IMPLICIT NONE
	  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: EPS
	  INTEGER, INTENT(IN) :: ntens
	  DOUBLE PRECISION, DIMENSION(ntens) :: HYDRO
	  
	  DOUBLE PRECISION :: TRACE
	  INTEGER :: I
	  
	  TRACE = EPS(1)+EPS(2)+EPS(3)
	  DO I=1, 3
	    HYDRO(I) = 1.D0/3.D0*TRACE
	    HYDRO(I+3) = 0.D0
	  ENDDO
    END FUNCTION
	
    FUNCTION DEVIATORIC(EPS, HYDRO, ntens) RESULT (DEV)
	  ! Calculate the deviatoric strain/stress without changing factors (covariant/contravariant)
	  IMPLICIT NONE
	  DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: EPS, HYDRO
	  INTEGER, INTENT(IN) :: ntens
	  DOUBLE PRECISION, DIMENSION(ntens) :: DEV
	  
	  INTEGER :: I
	  
	  DO I=1, 3
	    DEV(I) = EPS(I)-HYDRO(I)
	    DEV(I+3) = EPS(I+3)
      ENDDO
    END FUNCTION
    
    FUNCTION GET_IND_STRESS(SIGMA) RESULT (IND_STRESS)
      ! Calculate indicator stress IND_STRESS for a given stress state SIGMA
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(6), INTENT(IN) :: SIGMA
      DOUBLE PRECISION :: IND_STRESS
      
      IND_STRESS = sqrt(SIGMA(1)**2.D0+SIGMA(2)**2.D0+SIGMA(3)**2.D0 &
                        -SIGMA(1)*SIGMA(2)-SIGMA(1)*SIGMA(3)-SIGMA(2)*SIGMA(3) &
                        +3.D0*(SIGMA(4)**2.D0+SIGMA(5)**2.D0+SIGMA(6)**2.D0))
						
	  IF (IND_STRESS .ne. IND_STRESS) THEN
	    print *, 'Probably underflow, SIGMA:', SIGMA
		IND_STRESS = 0.D0
	  END IF
                        
    END FUNCTION
    
    FUNCTION GET_G0(IND_STRESS) RESULT (G0)
      ! Calculate G0 for a given indicator stress
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: IND_STRESS
      DOUBLE PRECISION :: G0
      
      ! G0 = 6E-5*IND_STRESS**2 + 0.00036D0*IND_STRESS + 1.D0
	  G0 = 1.D0
	  
	  ! IF (G0 > 1.D0) THEN
	    ! print *, 'G0', G0
      ! END IF
      
    END FUNCTION
    
    FUNCTION GET_G1G2(IND_STRESS) RESULT (G1G2)
      ! Calculate G0 for a given indicator stress
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: IND_STRESS
      DOUBLE PRECISION :: G1G2
      
	  G1G2 = 9.016174D0 + (0.9808576D0 - 9.016174D0)/(1.D0 + (IND_STRESS/34.69898D0)**4.223144D0)
      ! G1G2 = 1.D0
      
    END FUNCTION
    
    FUNCTION GET_A_SIG(IND_STRESS) RESULT (A_SIG)
      ! Calculate G0 for a given indicator stress
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: IND_STRESS
      DOUBLE PRECISION :: A_SIG
      
      A_SIG = 2.286324D0 + (0.9821686D0 - 2.286324D0)/(1.D0 + (IND_STRESS/22.15665D0)**4.349705D0)
	  ! A_SIG = 1.D0
      
    END FUNCTION
    
    FUNCTION GET_DPSI(dtime, A_SIG) RESULT (DPSI)
      ! Calculate the reduced time increment 
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: dtime
      DOUBLE PRECISION :: A_SIG
      DOUBLE PRECISION :: DPSI
      
      DPSI = dtime/A_SIG
      
    END FUNCTION
    
    FUNCTION GET_RES_NORM(DSTRAN_TRIAL, dstran, ntens, DEBUG) RESULT (STRAN_RES)
      ! Calculate the strain residuum norm
      IMPLICIT NONE
      DOUBLE PRECISION, DIMENSION(:), INTENT(IN) :: DSTRAN_TRIAL, dstran
      INTEGER, INTENT(IN) :: ntens
      LOGICAL, INTENT(IN) :: DEBUG
      INTEGER :: I
      DOUBLE PRECISION :: NORM_STRAN, NORM_DIFF
      DOUBLE PRECISION, DIMENSION(ntens) :: STRAN_DIFF
      DOUBLE PRECISION :: STRAN_RES
      
      STRAN_DIFF = DSTRAN_TRIAL - dstran
      
      NORM_DIFF = 0.D0
      NORM_STRAN = 0.D0
      DO I=1, ntens
        NORM_DIFF = NORM_DIFF + STRAN_DIFF(I)**2
	    NORM_STRAN = NORM_STRAN + dstran(I)**2
	  ENDDO
      
      NORM_DIFF = sqrt(NORM_DIFF)
      NORM_STRAN = sqrt(NORM_STRAN)
	  
	  ! IF (ABS(NORM_STRAN) < 1E-8) THEN
	    ! print *, 'NORM_STRAN:', NORM_STRAN
	  ! END IF
	  
	  
	  ! print *, 'dstran:', dstran
	  ! print *, 'DSTRAN_TRIAL:', DSTRAN_TRIAL
      
      STRAN_RES = NORM_DIFF/NORM_STRAN
      
      IF (STRAN_RES .ne. STRAN_RES) THEN
        IF (DEBUG == .TRUE.) THEN
	      print *, 'Probably underflow, STRAN_RES:', STRAN_RES
        END IF
		STRAN_RES = 0.D0
	  END IF
      
    END FUNCTION
	
	FUNCTION GET_G(G0_NL, G1G2_NL, DPSI, G0, G1, G2, G3, TAU1, TAU2, TAU3) RESULT (G)
	  ! Return the effective viscous shear modulus
	  IMPLICIT NONE
	  DOUBLE PRECISION, INTENT(IN) :: G0_NL, G1G2_NL, DPSI, G0, G1, G2, G3, TAU1, TAU2, TAU3
	  DOUBLE PRECISION :: G, DECAY1, DECAY2, DECAY3
	  
	  DECAY1 = 1.D0-exp(-DPSI/TAU1)
      DECAY2 = 1.D0-exp(-DPSI/TAU2)
      DECAY3 = 1.D0-exp(-DPSI/TAU3)
	  
	  G = 1.D0/(G0_NL/G0+G1G2_NL*(1.D0/G1*(1.D0-TAU1/DPSI*DECAY1)+1.D0/G2*(1.D0-TAU2/DPSI*DECAY2)+1.D0/G3*(1.D0-TAU3/DPSI*DECAY3)))
	
	END FUNCTION
	
	FUNCTION GET_K(G0_NL, G1G2_NL, DPSI, K0, K1, K2, K3, TAU1, TAU2, TAU3) RESULT (K)
	  ! Return the effective viscous shear modulus
	  IMPLICIT NONE
	  DOUBLE PRECISION, INTENT(IN) :: G0_NL, G1G2_NL, DPSI, K0, K1, K2, K3, TAU1, TAU2, TAU3
	  DOUBLE PRECISION :: K, DECAY1, DECAY2, DECAY3
	  
	  DECAY1 = 1.D0-exp(-DPSI/TAU1)
      DECAY2 = 1.D0-exp(-DPSI/TAU2)
      DECAY3 = 1.D0-exp(-DPSI/TAU3)
	  
	  K = 1.D0/(G0_NL/K0+G1G2_NL*(1.D0/K1*(1.D0-TAU1/DPSI*DECAY1)+1.D0/K2*(1.D0-TAU2/DPSI*DECAY2)+1.D0/K3*(1.D0-TAU3/DPSI*DECAY3)))
	
	END FUNCTION
	
	FUNCTION GET_DEPS(stress, DSIGMA, STRAN_INH1, STRAN_INH2, STRAN_INH3, &
	  dtime, ntens, G0, G1, G2, G3, K0, K1, K2, K3, TAU1, TAU2, TAU3, M_PRIME, M_CIRC) RESULT (DEPS)
	  ! Calculate strain increment based on stress increment
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: ntens
	  DOUBLE PRECISION, INTENT(IN), DIMENSION(ntens) :: stress, DSIGMA, STRAN_INH1, STRAN_INH2, STRAN_INH3
	  DOUBLE PRECISION, INTENT(IN) :: dtime, G0, G1, G2, G3, K0, K1, K2, K3, TAU1, TAU2, TAU3
	  DOUBLE PRECISION, INTENT(IN), DIMENSION(6,6) :: M_PRIME, M_CIRC
	  DOUBLE PRECISION, DIMENSION(ntens) :: DEPS
	  DOUBLE PRECISION, DIMENSION(ntens) :: STRESS_NEW, HISTORY_STRAIN
	  DOUBLE PRECISION :: IND_STRESS, G0_NL, G1G2_NL, A_SIG_NL, DPSI, DECAY1, DECAY2, DECAY3, K, G
	  
	  STRESS_NEW = stress + DSIGMA
	  IND_STRESS = GET_IND_STRESS(STRESS_NEW) 
	  
	  G0_NL = GET_G0(IND_STRESS)
	  G1G2_NL = GET_G1G2(IND_STRESS)
	  A_SIG_NL = GET_A_SIG(IND_STRESS)
	  
	  DPSI = GET_DPSI(dtime, A_SIG_NL)
	  
	  G = GET_G(G0_NL, G1G2_NL, DPSI, G0, G1, G2, G3, TAU1, TAU2, TAU3)
	  K = GET_K(G0_NL, G1G2_NL, DPSI, K0, K1, K2, K3, TAU1, TAU2, TAU3)
	  
	  DECAY1 = 1.D0-exp(-DPSI/TAU1)
      DECAY2 = 1.D0-exp(-DPSI/TAU2)
      DECAY3 = 1.D0-exp(-DPSI/TAU3)
	  
	  HISTORY_STRAIN = (STRAN_INH1*DECAY1 + STRAN_INH2*DECAY2 + STRAN_INH3*DECAY3)
	  
	  DEPS = matmul((1.D0/(2.D0*G)*M_PRIME + 1.D0/(3.D0*K)*M_CIRC), DSIGMA) + HISTORY_STRAIN
	  
    END FUNCTION
	
	FUNCTION GET_JACOBIAN(stress, DSIGMA, STRAN_INH1, STRAN_INH2, STRAN_INH3, &
	  dtime, ntens, G0, G1, G2, G3, K0, K1, K2, K3, TAU1, TAU2, TAU3, M_PRIME, M_CIRC, h) RESULT (JACOBIAN)
	  ! Calculate Jacobian (numerically)
	  IMPLICIT NONE
	  INTEGER, INTENT(IN) :: ntens
	  DOUBLE PRECISION, INTENT(IN), DIMENSION(ntens) :: stress, DSIGMA, STRAN_INH1, STRAN_INH2, STRAN_INH3
	  DOUBLE PRECISION, INTENT(IN) :: dtime, G0, G1, G2, G3, K0, K1, K2, K3, TAU1, TAU2, TAU3
	  DOUBLE PRECISION, INTENT(IN), DIMENSION(6,6) :: M_PRIME, M_CIRC
	  DOUBLE PRECISION, DIMENSION(6,6) :: JACOBIAN
	  DOUBLE PRECISION, DIMENSION(ntens) :: DH
	  DOUBLE PRECISION :: h
	  INTEGER :: I, J
	  
	  DO J=1, ntens
	    DO I=1, ntens
		  IF (I == J) THEN
	        DH(I) = h
		  ELSE
		    DH(I) = 0.D0
	      END IF
        ENDDO
	    JACOBIAN(:,J) = 1.D0/(2.D0*h)*(GET_DEPS(stress, DSIGMA+DH, STRAN_INH1, STRAN_INH2, STRAN_INH3, &
	                    dtime, ntens, G0, G1, G2, G3, K0, K1, K2, K3, TAU1, TAU2, TAU3, M_PRIME, M_CIRC) &
	                    - GET_DEPS(stress, DSIGMA-DH, STRAN_INH1, STRAN_INH2, STRAN_INH3, &
	                    dtime, ntens, G0, G1, G2, G3, K0, K1, K2, K3, TAU1, TAU2, TAU3, M_PRIME, M_CIRC))
	  ENDDO
      
      write(*,*) 'Jacobian: ', JACOBIAN
      
    END FUNCTION

END MODULE FUNCTIONS

subroutine vumat(nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,&
    stepTime, totalTime, dt, cmname, coordMp, charLength,&
    props, density, strainInc, relSpinInc,&
    tempOld, stretchOld, defgradOld, fieldOld,&
    stressOld, stateOld, enerInternOld, enerInelasOld,&
    tempNew, stretchNew, defgradNew, fieldNew,&
    stressNew, stateNew, enerInternNew, enerInelasNew)

    USE FUNCTIONS
    include 'vaba_param.inc'

    dimension props(nprops), density(nblock), coordMp(nblock,*),&
        charLength(nblock), strainInc(nblock,ndir+nshr),&
        relSpinInc(nblock,nshr), tempOld(nblock),&
        stretchOld(nblock,ndir+nshr),&
        defgradOld(nblock,ndir+nshr+nshr),&
        fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),&
        stateOld(nblock,nstatev), enerInternOld(nblock),&
        enerInelasOld(nblock), tempNew(nblock),&
        stretchNew(nblock,ndir+nshr),&
        defgradNew(nblock,ndir+nshr+nshr),&
        fieldNew(nblock,nfieldv),&
        stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),&
        enerInternNew(nblock), enerInelasNew(nblock)

    character*80 cmname
    
    integer													:: i,j, n, info, ntens
	
	INTEGER                                                 :: ITER
    INTEGER, PARAMETER	                                    :: ITER_MAX = 50
	DOUBLE PRECISION, PARAMETER                             :: h = 1E-6
	DOUBLE PRECISION, PARAMETER                             :: RES_MIN = 1E-4
	
	DOUBLE PRECISION :: E0, NU, E1, TAU1, E2, TAU2, E3, TAU3, G0, K0, G1, K1, G2, K2, G3, K3, G, K 
	DOUBLE PRECISION :: DECAY1, DECAY2, DECAY3, A_SIG, DPSI, STRAN_RES, G0_NL, G1G2_NL, A_SIG_NL
	DOUBLE PRECISION, DIMENSION(6) :: DSTRESS, DSIGMA_TRIAL, DSTRAN_TRIAL, DDSTRESS, R
	DOUBLE PRECISION, DIMENSION(6) :: STRAN_INH1, STRAN_INH2, STRAN_INH3, HISTORY_STRAIN, STRAIN_DIFF
	DOUBLE PRECISION, DIMENSION(6,6) :: M_PRIME, M_CIRC, JACOBIAN
	
	DOUBLE PRECISION, DIMENSION(3,3) :: A, A_INV
	DOUBLE PRECISION, DIMENSION(3) :: b, x
    
	DOUBLE PRECISION :: IND_STRESS
	
	LOGICAL, PARAMETER :: DEBUG  = .FALSE.
       
    
    M_PRIME = reshape( 1.D0/6.D0*(/4.D0, -2.D0, -2.D0, 0.D0, 0.D0, 0.D0, &
                                      -2.D0, 4.D0, -2.D0, 0.D0, 0.D0, 0.D0, &
                                      -2.D0, -2.D0, 4.D0, 0.D0, 0.D0, 0.D0, & 
                                      0.D0, 0.D0, 0.D0, 6.D0, 0.D0, 0.D0, &
                                      0.D0, 0.D0, 0.D0, 0.D0, 6.D0, 0.D0, &
                                      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 6.D0/), (/6, 6/) )
									  
	M_CIRC = reshape( 1.D0/3.D0*(/1.D0, 1.D0, 1.D0, 0.D0, 0.D0, 0.D0, &
                                      1.D0, 1.D0, 1.D0, 0.D0, 0.D0, 0.D0, &
                                      1.D0, 1.D0, 1.D0, 0.D0, 0.D0, 0.D0, & 
                                      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
                                      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
                                      0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0/), (/6, 6/) )
    
! inputs
    
    ntens = ndir + nshr
    
    E0 = PROPS(1)
	NU = PROPS(2)
	E1 = PROPS(3)
	TAU1 = PROPS(4)
    E2 = PROPS(5)
	TAU2 = PROPS(6)
    E3 = PROPS(7)
	TAU3 = PROPS(8)
    
    G0 = E0/(2.D0*(1.D0+NU))
	K0 = E0/(3.D0*(1.D0-2.D0*NU))
	G1 = E1/(2.D0*(1.D0+NU))
	K1 = E1/(3.D0*(1.D0-2.D0*NU))
    G2 = E2/(2.D0*(1.D0+NU))
	K2 = E2/(3.D0*(1.D0-2.D0*NU))
	G3 = E3/(2.D0*(1.D0+NU))
	K3 = E3/(3.D0*(1.D0-2.D0*NU))
    
    
    
    
! stress increment evaluation for each element
    do i = 1, nblock
         
      ! Build state variables
	  STRAN_INH1(1) = stateOld(i, 1)
	  STRAN_INH1(2) = stateOld(i, 2)
	  STRAN_INH1(3) = stateOld(i, 3)
	  STRAN_INH1(4) = stateOld(i, 4)
	  STRAN_INH1(5) = stateOld(i, 5)
	  STRAN_INH1(6) = stateOld(i, 6)

	  STRAN_INH2(1) = stateOld(i, 7)
	  STRAN_INH2(2) = stateOld(i, 8)
	  STRAN_INH2(3) = stateOld(i, 9)
	  STRAN_INH2(4) = stateOld(i, 10)
	  STRAN_INH2(5) = stateOld(i, 11)
	  STRAN_INH2(6) = stateOld(i, 12)

	  STRAN_INH3(1) = stateOld(i, 13)
	  STRAN_INH3(2) = stateOld(i, 14)
	  STRAN_INH3(3) = stateOld(i, 15)
	  STRAN_INH3(4) = stateOld(i, 16)
	  STRAN_INH3(5) = stateOld(i, 17)
	  STRAN_INH3(6) = stateOld(i, 18)
      
      ! Old nonlinear parameters
      IF (totalTime == 0.D0) THEN
        G0_NL = 1.D0
        G1G2_NL = 1.D0
        A_SIG_NL = 1.D0
      ELSE
        G0_NL = stateOld(i, 19)
        G1G2_NL = stateOld(i, 20)
        A_SIG_NL = stateOld(i, 21)
        
        IF (G0_NL == 0.D0) THEN
          G0_NL = 1.D0
        END IF
        IF (G1G2_NL == 0.D0) THEN
          G1G2_NL = 1.D0
        END IF
        IF (A_SIG_NL == 0.D0) THEN
          A_SIG_NL = 1.D0
        END IF
      END IF
      
      ! Get new reduced time increment and time dependent material parameters
      DPSI = GET_DPSI(dt, A_SIG_NL)
      
      DECAY1 = 1.D0-exp(-DPSI/TAU1)
      DECAY2 = 1.D0-exp(-DPSI/TAU2)
      DECAY3 = 1.D0-exp(-DPSI/TAU3)
      
      G = GET_G(G0_NL, G1G2_NL, DPSI, G0, G1, G2, G3, TAU1, TAU2, TAU3)
	  K = GET_K(G0_NL, G1G2_NL, DPSI, K0, K1, K2, K3, TAU1, TAU2, TAU3)
      
      HISTORY_STRAIN = (STRAN_INH1*DECAY1 + STRAN_INH2*DECAY2 + STRAN_INH3*DECAY3)
      STRAIN_DIFF = strainInc(i,:) - HISTORY_STRAIN
      
      ! Calculate trial stress based on nonlinear parameters from previous time step
      DSIGMA_TRIAL = matmul((2.D0*G*M_PRIME + 3.D0*K*M_CIRC), STRAIN_DIFF)
      
      IF (DEBUG == .TRUE.) THEN
	    print *, 'time', totalTime
        print *, 'i', i
	    print *, 'A_SIG_NL', A_SIG_NL
	    print *, 'dtime', dt
	    print *, 'DPSI', DPSI
	    print *, 'dstran', strainInc(i,:)
	    print *, 'HISTORY_STRAIN', HISTORY_STRAIN
        print *, 'G', G
        print *, 'K', K
	    print *, 'STRAIN_DIFF', STRAIN_DIFF
	    print *, 'stress', stressOld(i,:)
	    print *, 'DSIGMA_TRIAL', DSIGMA_TRIAL
	    print *, 'IND_STRESS', GET_IND_STRESS(stressOld(i,:)+DSIGMA_TRIAL)
	  END IF
	
      
      ! Calculate trial strain to check if trial stress is sufficient
      DSTRAN_TRIAL = GET_DEPS(stressOld(i,:), DSIGMA_TRIAL, STRAN_INH1, STRAN_INH2, STRAN_INH3, &
	    dt, ntens, G0, G1, G2, G3, K0, K1, K2, K3, TAU1, TAU2, TAU3, M_PRIME, M_CIRC)
      
      ! Calculate residuum
      STRAN_RES = GET_RES_NORM(DSTRAN_TRIAL, strainInc(i,:), ntens, DEBUG)
      
      ! ---- Newton-Raphson ----
	  ITER = 0
	  DSTRESS = DSIGMA_TRIAL
	  n = ntens ! make parameter later
      DO WHILE(ITER < ITER_MAX .AND. STRAN_RES > RES_MIN)
	    IF (DEBUG == .TRUE.) THEN
	      IF (ITER == 0) THEN
	        print *, 'Going into Newton-Raphson...'
          END IF
	    END IF
	  
	    ! Calculate Jacobian
	    JACOBIAN = GET_JACOBIAN(stressOld(i,:), DSTRESS, STRAN_INH1, STRAN_INH2, STRAN_INH3, &
	                     dt, ntens, G0, G1, G2, G3, K0, K1, K2, K3, TAU1, TAU2, TAU3, &
	  				     M_PRIME, M_CIRC, h)
	    
	    ! Solve for DDSTRESS
	    R = DSTRAN_TRIAL - strainInc(i,:)
	    call linsolv(n,JACOBIAN,-R,DDSTRESS,info)
	    IF (DEBUG == .TRUE.) THEN
	      if (info.ne. 0) then 
	        print *, 'ITER:', ITER
	        print *, 'JACOBIAN:', JACOBIAN
	  	  print *, 'stress:', stressOld(i,:)
	  	  print *, 'DSTRESS:', DSTRESS
	        CALL XIT()
	      end if
	    END IF
	   
	    ! Update DSTRESS
	    DSTRESS = DSTRESS + DDSTRESS
	   
	    ! Update trial strain
	    DSTRAN_TRIAL = GET_DEPS(stressOld(i,:), DSTRESS, STRAN_INH1, STRAN_INH2, STRAN_INH3, &
	    dt, ntens, G0, G1, G2, G3, K0, K1, K2, K3, TAU1, TAU2, TAU3, M_PRIME, M_CIRC)
	  
	    STRAN_RES = GET_RES_NORM(DSTRAN_TRIAL, strainInc(i,:), ntens, DEBUG)
	    ITER = ITER + 1
	  END DO
      
      IF (ITER == ITER_MAX) THEN
	    print *, 'Newton-Raphson did not converge (time, noel)', totalTime
	    IF (DEBUG == .TRUE.) THEN
	      print *, 'STRAN_RES', STRAN_RES
	      print *, 'DSTRESS', DSTRESS
        END IF
	  END IF
	
	  IF (DEBUG == .TRUE.) THEN
	    print *, ''
	  END IF
      
      ! Update stress
      stressNew(i, :) = stressOld(i,:) + DSTRESS
	  IND_STRESS = GET_IND_STRESS(stressNew(i, :))
	  
	  G0_NL = GET_G0(IND_STRESS)
	  G1G2_NL = GET_G1G2(IND_STRESS)
	  A_SIG_NL = GET_A_SIG(IND_STRESS)
	  DPSI = GET_DPSI(dt, A_SIG_NL)
      
      DECAY1 = 1.D0-exp(-DPSI/TAU1)
      DECAY2 = 1.D0-exp(-DPSI/TAU2)
      DECAY3 = 1.D0-exp(-DPSI/TAU3)
      
      ! Update state variables
	  STRAN_INH1 = STRAN_INH1 + STRAN_INH1*(-1.D0)*DECAY1 + G1G2_NL*TAU1/DPSI*DECAY1*matmul((1.D0/(2.D0*G1)*M_PRIME + 1.D0/(3.D0*K1)*M_CIRC), DSTRESS)
	  STRAN_INH2 = STRAN_INH2 + STRAN_INH2*(-1.D0)*DECAY2 + G1G2_NL*TAU2/DPSI*DECAY2*matmul((1.D0/(2.D0*G2)*M_PRIME + 1.D0/(3.D0*K2)*M_CIRC), DSTRESS)
	  STRAN_INH3 = STRAN_INH3 + STRAN_INH3*(-1.D0)*DECAY3 + G1G2_NL*TAU3/DPSI*DECAY3*matmul((1.D0/(2.D0*G3)*M_PRIME + 1.D0/(3.D0*K3)*M_CIRC), DSTRESS)	
	  
	  ! Write state variables to STATEV
	  stateNew(i, 1) = STRAN_INH1(1)
	  stateNew(i, 2) = STRAN_INH1(2)
	  stateNew(i, 3) = STRAN_INH1(3)
	  stateNew(i, 4) = STRAN_INH1(4)
	  stateNew(i, 5) = STRAN_INH1(5)
	  stateNew(i, 6) = STRAN_INH1(6)
 
	  stateNew(i, 7) = STRAN_INH2(1)
	  stateNew(i, 8) = STRAN_INH2(2)
	  stateNew(i, 9) = STRAN_INH2(3)
	  stateNew(i, 10) = STRAN_INH2(4)
	  stateNew(i, 11) = STRAN_INH2(5)
	  stateNew(i, 12) = STRAN_INH2(6)

	  stateNew(i, 13) = STRAN_INH3(1)
	  stateNew(i, 14) = STRAN_INH3(2)
	  stateNew(i, 15) = STRAN_INH3(3)
	  stateNew(i, 16) = STRAN_INH3(4)
	  stateNew(i, 17) = STRAN_INH3(5)
	  stateNew(i, 18) = STRAN_INH3(6)

      stateNew(i, 19) = G0_NL
      stateNew(i, 20) = G1G2_NL
      stateNew(i, 21) = A_SIG_NL

      stateNew(i, 22) = STRAN_RES
	  stateNew(i, 23) = IND_STRESS

    end do
return
end subroutine vumat

subroutine invert (n,A,Ainv) 
	implicit none 
	integer n,info
	integer ipiv(n)
	integer nsize 
	parameter (nsize=9) ! matrix size 
	integer LWORK_MKL
	parameter (LWORK_MKL=64*nsize)
	real*8 WORK_MKL (LWORK_MKL)
	real*8 A(n,n), Ainv (n,n)
	Ainv = A
	! Intially, Ainv = A.
	! After linear solve, Ainv = inv(A)
	CALL dgetrf (n,n,Ainv,n, ipiv,info)
	CALL dgetri (n,Ainv,n, ipiv, WORK_MKL, LWORK_MKL, info)
end subroutine invert

subroutine linsolv(n,Ain,b,x, info)
	implicit none 
	integer n,info,nrhs
	integer ipiv(n)
	real*8 A(n, n), Ain (n,n), b(n), x(n)
	parameter (nrhs=1)
	A = Ain
	x = b
	
	! Initially x = b
	! After linear solve, x=inv(A)*b
	call dgesv(n, nrhs, A, n, ipiv, x, n, info)
	
	if (info.ne. 0) then 
	  write(*,*) '***ERROR INVERTING LOCAL JACOBIAN***', info
      write(*,*) 'n: ', n
      write(*,*) 'Jacobian: ', Ain
	end if
end subroutine linsolv
