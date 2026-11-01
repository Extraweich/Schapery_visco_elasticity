!DEC$ FREEFORM
!=======================================================================
!  NLVE_SCHAPERY - nonlinear viscoelastic material model after Schapery
!
!  Generalized Kelvin-Voigt model (spring + 3 Kelvin-Voigt elements in
!  series) with the stress-dependent nonlinearizing functions g0, g1*g2
!  and a_sigma of the Schapery integral. Implicit stress update by a local
!  Newton-Raphson iteration on the strain residual (total-strain form).
!
!  Author : Nicolas Christ (Karlsruhe Institute of Technology)
!  Theory : PhD thesis N. Christ, chapter "Matrix", Eqs. (63ext)-(74ext)
!  License: MIT (see LICENSE). Developed during the author's PhD project at KIT.
!           If you use this code in published work, please cite the thesis
!           (see README.md / CITATION.cff).
!
!  ONE source file for both solvers - pass it directly to Abaqus:
!     abaqus job=<job> user=nlve_schapery.f              (Standard: UMAT)
!     abaqus job=<job> user=nlve_schapery.f              (Explicit: VUMAT)
!  Abaqus/Explicit must run in double precision (default: double=explicit).
!  The file uses no Abaqus utility routines and no LAPACK, so the same
!  object links into Abaqus/Standard and Abaqus/Explicit.
!
!  Supported elements: 3D continuum, plane strain, axisymmetric
!     (UMAT: NDI=3, NSHR=3 or 1; VUMAT: NDIR=3, NSHR=3 or 1).
!  NOT supported: plane stress, shells, beams. Abaqus/Explicit cannot tell
!  plane stress from plane strain inside a VUMAT - do not use CPS/S elements.
!
!  Implicit -> explicit: the state variables are stored identically by
!  UMAT and VUMAT, so a Standard job can be continued in Explicit with
!  *IMPORT, STATE=YES (and vice versa).
!
!  Units: any consistent set; the parameter functions take the von Mises
!  stress in the same stress unit as the moduli.
!
!-----------------------------------------------------------------------
!  *USER MATERIAL, CONSTANTS=23   (UMAT: add UNSYMM, see DDSDDE below)
!  PROPS( 1) E0      instantaneous modulus
!  PROPS( 2) NU      Poisson's ratio (same for all springs)
!  PROPS( 3) E1      PROPS( 4) TAU1     Kelvin-Voigt element 1
!  PROPS( 5) E2      PROPS( 6) TAU2     Kelvin-Voigt element 2
!  PROPS( 7) E3      PROPS( 8) TAU3     Kelvin-Voigt element 3
!  PROPS( 9) type of g0         PROPS(10:13) its parameters p1..p4
!  PROPS(14) type of g1*g2      PROPS(15:18) its parameters p1..p4
!  PROPS(19) type of a_sigma    PROPS(20:23) its parameters p1..p4
!
!  Function types f(s), s = von Mises stress:
!   0  linear     f = 1                    (p1..p4 ignored)
!   1  sigmoid    f = p4 + (p1-p4)/(1 + (s/p3)**p2)    [a,b,c,d = p1..p4]
!   2  gaussian   f = p1 + p2*exp(-(s-p3)**2/(2*p4**2))
!  All types 0 -> linear viscoelasticity (verification against *VISCOELASTIC).
!
!-----------------------------------------------------------------------
!  *DEPVAR = 29   (same layout for UMAT and VUMAT, 3D and axisymmetric)
!  Tensors in Abaqus/Standard notation: order 11,22,33,12,13,23 and
!  engineering shear strains. In 2D the components 13 and 23 stay zero.
!   1- 6  inherited strain of Kelvin-Voigt element 1
!   7-12  inherited strain of Kelvin-Voigt element 2
!  13-18  inherited strain of Kelvin-Voigt element 3
!  19     g0         20  g1*g2        21  a_sigma   (at end of increment)
!  22     relative strain residual of the local Newton iteration
!  23     indicator (von Mises) stress
!  24-29  total strain (needed by the VUMAT and for *IMPORT)
!  SDVINI (*INITIAL CONDITIONS, TYPE=SOLUTION, USER) is optional: zero
!  values of SDV 19-21 are treated as 1.
!
!-----------------------------------------------------------------------
!  DDSDDE (UMAT) is the consistent tangent, i.e. the inverse of the local
!  Jacobian d(eps)/d(sigma) at the converged stress. It is not symmetric:
!  use *USER MATERIAL, UNSYMM for full quadratic convergence. Without UNSYMM
!  Abaqus uses its symmetric part, which still converges well.
!
!  Calls with a zero time increment (e.g. initial accelerations in
!  *DYNAMIC) return the instantaneous response and leave the state
!  variables unchanged. In the VUMAT data-check call (stepTime = totalTime
!  = 0) the initial elastic response with E0 is returned.
!
!  Input errors (wrong DEPVAR, CONSTANTS, element type, invalid PROPS)
!  stop the analysis with ERROR STOP and a message in the .log file.
!=======================================================================

MODULE NLVE_SCHAPERY
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: DP = KIND(1.D0)
  INTEGER, PARAMETER, PUBLIC :: NLVE_NPROPS = 23
  INTEGER, PARAMETER, PUBLIC :: NLVE_NSDV   = 29

  ! ---- numerical settings of the local Newton-Raphson iteration
  INTEGER,  PARAMETER :: ITER_MAX = 50       ! maximum number of iterations
  REAL(DP), PARAMETER :: RES_TOL  = 1.D-8    ! tolerance of relative strain residual
  REAL(DP), PARAMETER :: EPS_REF  = 1.D-6    ! strain floor -> abs. tol. 1e-14 for dstran=0
  REAL(DP), PARAMETER :: H_REL    = 1.D-6    ! relative step of the numerical Jacobian
  INTEGER,  PARAMETER :: MAX_WARN = 20       ! printed non-convergence warnings

  ! ---- projection matrices, Voigt notation 11,22,33,12,13,23, engineering shear strain
  REAL(DP), PARAMETER :: C23 = 2.D0/3.D0, C13 = 1.D0/3.D0
  ! spherical projector (identical for stresses and strains)
  REAL(DP), PARAMETER :: MC(6,6) = RESHAPE([ C13, C13, C13, 0.D0, 0.D0, 0.D0, &
                                              C13, C13, C13, 0.D0, 0.D0, 0.D0, &
                                              C13, C13, C13, 0.D0, 0.D0, 0.D0, &
                                              0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
                                              0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, &
                                              0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0], [6, 6])
  ! deviatoric projector strain -> stress (shear factor 1/2)
  REAL(DP), PARAMETER :: MPE(6,6) = RESHAPE([ C23, -C13, -C13, 0.D0, 0.D0, 0.D0, &
                                              -C13,  C23, -C13, 0.D0, 0.D0, 0.D0, &
                                              -C13, -C13,  C23, 0.D0, 0.D0, 0.D0, &
                                              0.D0, 0.D0, 0.D0, 0.5D0, 0.D0, 0.D0, &
                                              0.D0, 0.D0, 0.D0, 0.D0, 0.5D0, 0.D0, &
                                              0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.5D0], [6, 6])
  ! deviatoric projector stress -> strain (shear factor 2)
  REAL(DP), PARAMETER :: MPS(6,6) = RESHAPE([ C23, -C13, -C13, 0.D0, 0.D0, 0.D0, &
                                              -C13,  C23, -C13, 0.D0, 0.D0, 0.D0, &
                                              -C13, -C13,  C23, 0.D0, 0.D0, 0.D0, &
                                              0.D0, 0.D0, 0.D0, 2.D0, 0.D0, 0.D0, &
                                              0.D0, 0.D0, 0.D0, 0.D0, 2.D0, 0.D0, &
                                              0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 2.D0], [6, 6])
  ! deviatoric projector without change of basis (stress -> stress)
  REAL(DP), PARAMETER :: MP(6,6) = RESHAPE([ C23, -C13, -C13, 0.D0, 0.D0, 0.D0, &
                                             -C13,  C23, -C13, 0.D0, 0.D0, 0.D0, &
                                             -C13, -C13,  C23, 0.D0, 0.D0, 0.D0, &
                                             0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, &
                                             0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, &
                                             0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0], [6, 6])

  TYPE, PUBLIC :: NLVE_MAT
    REAL(DP) :: E(0:3), NU, TAU(3), G(0:3), K(0:3)
    INTEGER  :: FTYPE(3)          ! 1: g0, 2: g1*g2, 3: a_sigma
    REAL(DP) :: FPAR(4, 3)
  END TYPE

  INTEGER, SAVE :: NWARN = 0

  PUBLIC :: NLVE_SET_MATERIAL, NLVE_UPDATE, NLVE_ELASTIC, NLVE_ROT_STRAIN, NLVE_FATAL, NLVE_WARN

CONTAINS

  ! ===================================================================
  !  error handling (no XIT / XPLB_EXIT so that one object links for both solvers)
  ! ===================================================================
  SUBROUTINE NLVE_FATAL(MSG)
    CHARACTER(LEN=*), INTENT(IN) :: MSG
    WRITE(*, '(A)') ' '
    WRITE(*, '(A)') ' *** NLVE_SCHAPERY ERROR: '//TRIM(MSG)
    WRITE(*, '(A)') ' '
    ERROR STOP 'NLVE_SCHAPERY: input error, see message above'
  END SUBROUTINE

  SUBROUTINE NLVE_WARN(MSG, RES, TIME, IELEM)
    CHARACTER(LEN=*), INTENT(IN) :: MSG
    REAL(DP), INTENT(IN) :: RES, TIME
    INTEGER,  INTENT(IN) :: IELEM
    NWARN = NWARN + 1
    IF (NWARN <= MAX_WARN) THEN
      WRITE(*, '(A,A,ES10.3,A,ES12.5,A,I9)') ' NLVE_SCHAPERY WARNING: ', TRIM(MSG), RES, &
            '  time =', TIME, '  element/block point =', IELEM
      IF (NWARN == MAX_WARN) WRITE(*, '(A)') ' NLVE_SCHAPERY: further warnings suppressed'
    END IF
  END SUBROUTINE

  ! ===================================================================
  !  material parameters
  ! ===================================================================
  SUBROUTINE NLVE_SET_MATERIAL(PROPS, NPROPS, MAT, CHECK)
    INTEGER,  INTENT(IN)  :: NPROPS
    REAL(DP), INTENT(IN)  :: PROPS(NPROPS)
    TYPE(NLVE_MAT), INTENT(OUT) :: MAT
    LOGICAL,  INTENT(IN)  :: CHECK
    INTEGER :: I
    CHARACTER(LEN=120) :: TXT

    IF (NPROPS < NLVE_NPROPS) THEN
      WRITE(TXT, '(A,I3,A,I3)') 'CONSTANTS must be ', NLVE_NPROPS, ' but is ', NPROPS
      CALL NLVE_FATAL(TXT)
    END IF
    MAT%E(0) = PROPS(1);  MAT%NU = PROPS(2)
    MAT%E(1) = PROPS(3);  MAT%TAU(1) = PROPS(4)
    MAT%E(2) = PROPS(5);  MAT%TAU(2) = PROPS(6)
    MAT%E(3) = PROPS(7);  MAT%TAU(3) = PROPS(8)
    DO I = 1, 3
      MAT%FTYPE(I)   = NINT(PROPS(9 + 5*(I-1)))
      MAT%FPAR(:, I) = PROPS(10 + 5*(I-1) : 13 + 5*(I-1))
    END DO

    IF (CHECK) THEN
      IF (ANY(MAT%E <= 0.D0))   CALL NLVE_FATAL('moduli E0..E3 (PROPS 1,3,5,7) must be > 0')
      IF (ANY(MAT%TAU <= 0.D0)) CALL NLVE_FATAL('relaxation times (PROPS 4,6,8) must be > 0')
      IF (MAT%NU <= -1.D0 .OR. MAT%NU >= 0.5D0) CALL NLVE_FATAL('Poisson ratio (PROPS 2) must be in (-1, 0.5)')
      DO I = 1, 3
        SELECT CASE (MAT%FTYPE(I))
        CASE (0)
        CASE (1)
          IF (MAT%FPAR(3, I) <= 0.D0) CALL NLVE_FATAL('sigmoid parameter p3 (c) must be > 0')
        CASE (2)
          IF (MAT%FPAR(4, I) == 0.D0) CALL NLVE_FATAL('gaussian parameter p4 must be /= 0')
        CASE DEFAULT
          CALL NLVE_FATAL('function type (PROPS 9,14,19) must be 0, 1 or 2')
        END SELECT
      END DO
    END IF

    MAT%G = MAT%E / (2.D0*(1.D0 + MAT%NU))
    MAT%K = MAT%E / (3.D0*(1.D0 - 2.D0*MAT%NU))
  END SUBROUTINE

  ! nonlinearizing function of type ITYPE with parameters P at von Mises stress S
  PURE REAL(DP) FUNCTION NLFUN(ITYPE, P, S)
    INTEGER,  INTENT(IN) :: ITYPE
    REAL(DP), INTENT(IN) :: P(4), S
    SELECT CASE (ITYPE)
    CASE (1)
      NLFUN = P(4) + (P(1) - P(4)) / (1.D0 + (S/P(3))**P(2))
    CASE (2)
      NLFUN = P(1) + P(2)*EXP(-(S - P(3))**2 / (2.D0*P(4)**2))
    CASE DEFAULT
      NLFUN = 1.D0
    END SELECT
  END FUNCTION

  ! g0, g1*g2, a_sigma for the von Mises stress S
  PURE SUBROUTINE NL_PARAMS(MAT, S, G0V, G12V, ASV)
    TYPE(NLVE_MAT), INTENT(IN) :: MAT
    REAL(DP), INTENT(IN)  :: S
    REAL(DP), INTENT(OUT) :: G0V, G12V, ASV
    G0V  = MAX(NLFUN(MAT%FTYPE(1), MAT%FPAR(:, 1), S), 0.D0)
    G12V = MAX(NLFUN(MAT%FTYPE(2), MAT%FPAR(:, 2), S), 0.D0)
    ASV  = MAX(NLFUN(MAT%FTYPE(3), MAT%FPAR(:, 3), S), 1.D-12)
  END SUBROUTINE

  ! von Mises stress (stress vector, tensor components)
  PURE REAL(DP) FUNCTION VMISES(S)
    REAL(DP), INTENT(IN) :: S(6)
    VMISES = SQRT(MAX(0.D0, S(1)**2 + S(2)**2 + S(3)**2 - S(1)*S(2) - S(1)*S(3) - S(2)*S(3) &
                            + 3.D0*(S(4)**2 + S(5)**2 + S(6)**2)))
  END FUNCTION

  ! (1-exp(-x))/x and 1-(1-exp(-x))/x without cancellation for small x
  PURE REAL(DP) FUNCTION PHI(X)
    REAL(DP), INTENT(IN) :: X
    IF (X < 1.D-4) THEN
      PHI = 1.D0 - X/2.D0 + X*X/6.D0 - X*X*X/24.D0
    ELSE
      PHI = (1.D0 - EXP(-X)) / X
    END IF
  END FUNCTION

  PURE REAL(DP) FUNCTION EXPM(X)          ! exp(-x) without underflow
    REAL(DP), INTENT(IN) :: X
    IF (X > 600.D0) THEN
      EXPM = 0.D0
    ELSE
      EXPM = EXP(-X)
    END IF
  END FUNCTION

  PURE REAL(DP) FUNCTION OMPHI(X)
    REAL(DP), INTENT(IN) :: X
    IF (X < 1.D-4) THEN
      OMPHI = X/2.D0 - X*X/6.D0 + X*X*X/24.D0
    ELSE
      OMPHI = 1.D0 - (1.D0 - EXP(-X)) / X
    END IF
  END FUNCTION

  ! effective moduli G~, K~ (Eq. 68ext), factors D', D° (Eq. 70ext), decay factors
  PURE SUBROUTINE COEFFS(MAT, G0V, G12V, ASV, G12_OLD, DT, GT, KT, DDEV, DSPH, EXPF)
    TYPE(NLVE_MAT), INTENT(IN) :: MAT
    REAL(DP), INTENT(IN)  :: G0V, G12V, ASV, G12_OLD, DT
    REAL(DP), INTENT(OUT) :: GT, KT, DDEV, DSPH, EXPF(3)
    REAL(DP) :: X(3), PH(3), OM(3)
    INTEGER  :: I
    DO I = 1, 3
      X(I)    = (DT/ASV) / MAT%TAU(I)
      EXPF(I) = EXPM(X(I))
      PH(I)   = PHI(X(I))
      OM(I)   = OMPHI(X(I))
    END DO
    GT   = 1.D0 / (G0V/MAT%G(0) + G12V*SUM(OM/MAT%G(1:3)))
    KT   = 1.D0 / (G0V/MAT%K(0) + G12V*SUM(OM/MAT%K(1:3)))
    DDEV = G12_OLD*SUM(PH/(2.D0*MAT%G(1:3)))
    DSPH = G12_OLD*SUM(PH/(3.D0*MAT%K(1:3)))
  END SUBROUTINE

  ! total strain for a given stress at the end of the increment (Eq. 74ext)
  PURE FUNCTION EPS_OF_SIG(MAT, SIG, SIG_OLD, HIST, G12_OLD, DT) RESULT(EPS)
    TYPE(NLVE_MAT), INTENT(IN) :: MAT
    REAL(DP), INTENT(IN) :: SIG(6), SIG_OLD(6), HIST(6, 3), G12_OLD, DT
    REAL(DP) :: EPS(6)
    REAL(DP) :: G0V, G12V, ASV, GT, KT, DDEV, DSPH, EXPF(3)
    CALL NL_PARAMS(MAT, VMISES(SIG), G0V, G12V, ASV)
    CALL COEFFS(MAT, G0V, G12V, ASV, G12_OLD, DT, GT, KT, DDEV, DSPH, EXPF)
    EPS = MATMUL(MPS/(2.D0*GT) + MC/(3.D0*KT), SIG) + MATMUL(DDEV*MPS + DSPH*MC, SIG_OLD) &
          - (EXPF(1)*HIST(:, 1) + EXPF(2)*HIST(:, 2) + EXPF(3)*HIST(:, 3))
  END FUNCTION

  ! numerical Jacobian d(eps)/d(sigma), central differences
  PURE FUNCTION JACOBIAN(MAT, SIG, SIG_OLD, HIST, G12_OLD, DT) RESULT(JAC)
    TYPE(NLVE_MAT), INTENT(IN) :: MAT
    REAL(DP), INTENT(IN) :: SIG(6), SIG_OLD(6), HIST(6, 3), G12_OLD, DT
    REAL(DP) :: JAC(6, 6), DS(6), H
    INTEGER  :: J
    H = H_REL*MAX(MAXVAL(ABS(SIG)), 1.D-4*MAT%E(0))
    DO J = 1, 6
      DS = 0.D0;  DS(J) = H
      JAC(:, J) = (EPS_OF_SIG(MAT, SIG + DS, SIG_OLD, HIST, G12_OLD, DT) &
                 - EPS_OF_SIG(MAT, SIG - DS, SIG_OLD, HIST, G12_OLD, DT)) / (2.D0*H)
    END DO
  END FUNCTION

  ! solve A x = b (6x6, Gauss elimination with partial pivoting)
  PURE SUBROUTINE SOLVE6(A_IN, B, X, OK)
    REAL(DP), INTENT(IN)  :: A_IN(6, 6), B(6)
    REAL(DP), INTENT(OUT) :: X(6)
    LOGICAL,  INTENT(OUT) :: OK
    REAL(DP) :: A(6, 6), F, TMP(6), T
    INTEGER  :: I, J, P
    A = A_IN;  X = B;  OK = .FALSE.
    DO I = 1, 6
      P = I - 1 + MAXLOC(ABS(A(I:6, I)), 1)
      IF (ABS(A(P, I)) <= 1.D-300 .OR. ABS(A(P, I)) /= ABS(A(P, I))) RETURN
      IF (P /= I) THEN
        TMP = A(I, :);  A(I, :) = A(P, :);  A(P, :) = TMP
        T = X(I);  X(I) = X(P);  X(P) = T
      END IF
      DO J = I + 1, 6
        F = A(J, I) / A(I, I)
        A(J, I:6) = A(J, I:6) - F*A(I, I:6)
        X(J) = X(J) - F*X(I)
      END DO
    END DO
    DO I = 6, 1, -1
      X(I) = (X(I) - DOT_PRODUCT(A(I, I+1:6), X(I+1:6))) / A(I, I)
    END DO
    OK = .TRUE.
  END SUBROUTINE

  ! rotate an engineering-shear strain vector (11,22,33,12,13,23): T' = R T R^T
  ! (replaces ROTSIG, which is not available in Abaqus/Explicit)
  PURE FUNCTION NLVE_ROT_STRAIN(V, R) RESULT(W)
    REAL(DP), INTENT(IN) :: V(6), R(3, 3)
    REAL(DP) :: W(6), T(3, 3)
    T(1,1) = V(1);        T(2,2) = V(2);        T(3,3) = V(3)
    T(1,2) = 0.5D0*V(4);  T(2,1) = T(1,2)
    T(1,3) = 0.5D0*V(5);  T(3,1) = T(1,3)
    T(2,3) = 0.5D0*V(6);  T(3,2) = T(2,3)
    T = MATMUL(R, MATMUL(T, TRANSPOSE(R)))
    W = [T(1,1), T(2,2), T(3,3), 2.D0*T(1,2), 2.D0*T(1,3), 2.D0*T(2,3)]
  END FUNCTION

  ! initial elastic response with E0 (Abaqus/Explicit data check)
  PURE FUNCTION NLVE_ELASTIC(MAT, DEPS) RESULT(DSIG)
    TYPE(NLVE_MAT), INTENT(IN) :: MAT
    REAL(DP), INTENT(IN) :: DEPS(6)
    REAL(DP) :: DSIG(6)
    DSIG = MATMUL(2.D0*MAT%G(0)*MPE + 3.D0*MAT%K(0)*MC, DEPS)
  END FUNCTION

  ! ===================================================================
  !  stress update at one material point
  !    SIG_OLD, EPS_OLD, DEPS : stress, total strain, strain increment
  !                             (6x1, 11,22,33,12,13,23, engineering shear)
  !    SDV_OLD / SDV_NEW      : state variables (layout see header)
  !    CTAN                   : consistent tangent (only if WANT_TAN)
  !    STATUS                 : 0 converged, 1 not converged, 2 singular Jacobian
  ! ===================================================================
  SUBROUTINE NLVE_UPDATE(MAT, DT, SIG_OLD, EPS_OLD, DEPS, SDV_OLD, SIG_NEW, SDV_NEW, &
                         WANT_TAN, CTAN, STATUS)
    TYPE(NLVE_MAT), INTENT(IN) :: MAT
    REAL(DP), INTENT(IN)  :: DT, SIG_OLD(6), EPS_OLD(6), DEPS(6), SDV_OLD(NLVE_NSDV)
    REAL(DP), INTENT(OUT) :: SIG_NEW(6), SDV_NEW(NLVE_NSDV), CTAN(6, 6)
    LOGICAL,  INTENT(IN)  :: WANT_TAN
    INTEGER,  INTENT(OUT) :: STATUS

    REAL(DP) :: HIST(6, 3), EPS_NEW(6), G0_OLD, G12_OLD, AS_OLD
    REAL(DP) :: GT, KT, DDEV, DSPH, EXPF(3), SIG(6), R(6), DSIG(6), JAC(6, 6), SIG_T(6), R_T(6)
    REAL(DP) :: RES, RES_T, DNORM, LAM, G0V, G12V, ASV, X, E6(6)
    INTEGER  :: I, ITER, LS
    LOGICAL  :: OK

    STATUS = 0
    HIST(:, 1) = SDV_OLD(1:6)
    HIST(:, 2) = SDV_OLD(7:12)
    HIST(:, 3) = SDV_OLD(13:18)
    ! nonlinear parameters of the previous increment (0 = not initialized -> 1)
    G0_OLD  = SDV_OLD(19);  IF (.NOT. (G0_OLD  > 0.D0)) G0_OLD  = 1.D0
    G12_OLD = SDV_OLD(20);  IF (.NOT. (G12_OLD > 0.D0)) G12_OLD = 1.D0
    AS_OLD  = SDV_OLD(21);  IF (.NOT. (AS_OLD  > 0.D0)) AS_OLD  = 1.D0
    EPS_NEW = EPS_OLD + DEPS
    DNORM   = MAX(NORM2(DEPS), EPS_REF)

    ! ---- trial stress with the nonlinear parameters of the previous increment (Eq. 73ext)
    CALL COEFFS(MAT, G0_OLD, G12_OLD, AS_OLD, G12_OLD, DT, GT, KT, DDEV, DSPH, EXPF)
    SIG = MATMUL(2.D0*GT*MPE + 3.D0*KT*MC, &
                 EPS_NEW + EXPF(1)*HIST(:, 1) + EXPF(2)*HIST(:, 2) + EXPF(3)*HIST(:, 3)) &
        - MATMUL(2.D0*GT*DDEV*MP + 3.D0*KT*DSPH*MC, SIG_OLD)

    ! ---- local Newton-Raphson on the strain residual (Eqs. 74-81)
    R   = EPS_OF_SIG(MAT, SIG, SIG_OLD, HIST, G12_OLD, DT) - EPS_NEW
    RES = NORM2(R) / DNORM
    ITER = 0
    DO WHILE (RES > RES_TOL .AND. ITER < ITER_MAX)
      ITER = ITER + 1
      JAC = JACOBIAN(MAT, SIG, SIG_OLD, HIST, G12_OLD, DT)
      CALL SOLVE6(JAC, -R, DSIG, OK)
      IF (.NOT. OK) THEN
        STATUS = 2
        EXIT
      END IF
      ! simple backtracking: accept the full step unless it increases the residual
      LAM = 1.D0
      DO LS = 0, 6
        SIG_T = SIG + LAM*DSIG
        R_T   = EPS_OF_SIG(MAT, SIG_T, SIG_OLD, HIST, G12_OLD, DT) - EPS_NEW
        RES_T = NORM2(R_T) / DNORM
        IF (RES_T < RES .OR. LS == 6) EXIT
        LAM = 0.5D0*LAM
      END DO
      SIG = SIG_T;  R = R_T;  RES = RES_T
    END DO
    IF (.NOT. (RES <= RES_TOL)) THEN       ! also catches NaN
      IF (STATUS == 0) STATUS = 1
    END IF
    SIG_NEW = SIG

    ! ---- state variables
    SDV_NEW = SDV_OLD
    IF (DT > 0.D0) THEN
      CALL NL_PARAMS(MAT, VMISES(SIG_NEW), G0V, G12V, ASV)
      DO I = 1, 3
        X  = (DT/ASV) / MAT%TAU(I)
        E6 = MATMUL(MPS/(2.D0*MAT%G(I)) + MC/(3.D0*MAT%K(I)), G12V*SIG_NEW - G12_OLD*SIG_OLD)
        SDV_NEW(6*I-5:6*I) = EXPM(X)*HIST(:, I) + PHI(X)*E6          ! Eq. (72ext)
      END DO
      SDV_NEW(19) = G0V
      SDV_NEW(20) = G12V
      SDV_NEW(21) = ASV
      SDV_NEW(22) = RES
      SDV_NEW(23) = VMISES(SIG_NEW)
      SDV_NEW(24:29) = EPS_NEW
    END IF

    ! ---- consistent tangent d(sigma)/d(eps) = inverse of the local Jacobian
    IF (WANT_TAN) THEN
      JAC = JACOBIAN(MAT, SIG_NEW, SIG_OLD, HIST, G12_OLD, DT)
      DO I = 1, 6
        E6 = 0.D0;  E6(I) = 1.D0
        CALL SOLVE6(JAC, E6, CTAN(:, I), OK)
        IF (.NOT. OK) EXIT
      END DO
      IF (.NOT. OK) THEN         ! fallback: secant-type tangent of the trial step
        CALL COEFFS(MAT, G0_OLD, G12_OLD, AS_OLD, G12_OLD, DT, GT, KT, DDEV, DSPH, EXPF)
        CTAN = 2.D0*GT*MPE + 3.D0*KT*MC
      END IF
    ELSE
      CTAN = 0.D0
    END IF
  END SUBROUTINE NLVE_UPDATE

END MODULE NLVE_SCHAPERY


!=======================================================================
!  Abaqus/Standard
!=======================================================================
SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, DRPLDE, DRPLDT, &
                STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP, PREDEF, DPRED, CMNAME, &
                NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS, COORDS, DROT, PNEWDT, &
                CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC)
  USE NLVE_SCHAPERY
  IMPLICIT NONE
  CHARACTER(LEN=80) :: CMNAME
  INTEGER :: NDI, NSHR, NTENS, NSTATV, NPROPS, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
  DOUBLE PRECISION :: STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS), SSE, SPD, SCD, RPL
  DOUBLE PRECISION :: DDSDDT(NTENS), DRPLDE(NTENS), DRPLDT, STRAN(NTENS), DSTRAN(NTENS)
  DOUBLE PRECISION :: TIME(2), DTIME, TEMP, DTEMP, PREDEF(*), DPRED(*), PROPS(NPROPS)
  DOUBLE PRECISION :: COORDS(3), DROT(3, 3), PNEWDT, CELENT, DFGRD0(3, 3), DFGRD1(3, 3)

  TYPE(NLVE_MAT)   :: MAT
  DOUBLE PRECISION :: SIG_OLD(6), EPS_OLD(6), DEPS(6), SIG_NEW(6), SDV_OLD(NLVE_NSDV)
  DOUBLE PRECISION :: SDV_NEW(NLVE_NSDV), CTAN(6, 6)
  INTEGER :: I, STATUS
  CHARACTER(LEN=120) :: TXT

  ! ---- element type and DEPVAR
  IF (NDI /= 3 .OR. (NSHR /= 3 .AND. NSHR /= 1)) THEN
    WRITE(TXT, '(A,2I3,A)') 'unsupported element (NDI, NSHR =', NDI, NSHR, &
                            '): only 3D, plane strain, axisymmetric'
    CALL NLVE_FATAL(TXT)
  END IF
  IF (NSTATV < NLVE_NSDV) THEN
    WRITE(TXT, '(A,I3,A,I3)') 'DEPVAR must be at least ', NLVE_NSDV, ' but is ', NSTATV
    CALL NLVE_FATAL(TXT)
  END IF
  CALL NLVE_SET_MATERIAL(PROPS, NPROPS, MAT, .TRUE.)

  ! ---- Abaqus components -> internal 6x1 vectors (same order and engineering shear;
  !      for NTENS = 4 the components 13 and 23 are zero)
  SIG_OLD = 0.D0;  EPS_OLD = 0.D0;  DEPS = 0.D0
  SIG_OLD(1:NTENS) = STRESS
  EPS_OLD(1:NTENS) = STRAN
  DEPS(1:NTENS)    = DSTRAN
  SDV_OLD = STATEV(1:NLVE_NSDV)

  ! ---- NLGEOM: Abaqus rotates STRESS and STRAN, tensor state variables are rotated here
  DO I = 1, 3
    SDV_OLD(6*I-5:6*I) = NLVE_ROT_STRAIN(SDV_OLD(6*I-5:6*I), DROT)
  END DO

  CALL NLVE_UPDATE(MAT, DTIME, SIG_OLD, EPS_OLD, DEPS, SDV_OLD, SIG_NEW, SDV_NEW, &
                   .TRUE., CTAN, STATUS)

  IF (STATUS == 1) THEN
    CALL NLVE_WARN('UMAT local Newton not converged, residual =', SDV_NEW(22), TIME(2), NOEL)
    PNEWDT = MIN(PNEWDT, 0.5D0)
  ELSE IF (STATUS == 2) THEN
    CALL NLVE_WARN('UMAT singular local Jacobian, residual =', SDV_NEW(22), TIME(2), NOEL)
    PNEWDT = MIN(PNEWDT, 0.25D0)
  END IF

  STRESS = SIG_NEW(1:NTENS)
  STATEV(1:NLVE_NSDV) = SDV_NEW
  DDSDDE = CTAN(1:NTENS, 1:NTENS)
END SUBROUTINE UMAT


SUBROUTINE SDVINI(STATEV, COORDS, NSTATV, NCRDS, NOEL, NPT, LAYER, KSPT)
  ! optional: *INITIAL CONDITIONS, TYPE=SOLUTION, USER
  IMPLICIT NONE
  INTEGER :: NSTATV, NCRDS, NOEL, NPT, LAYER, KSPT
  DOUBLE PRECISION :: STATEV(NSTATV), COORDS(NCRDS)
  STATEV = 0.D0
  IF (NSTATV >= 21) STATEV(19:21) = 1.D0       ! g0, g1*g2, a_sigma
END SUBROUTINE SDVINI


!=======================================================================
!  Abaqus/Explicit
!  Explicit order 11,22,33,12,23,31 and tensor shear strains; the state
!  variables are stored in the Standard notation of the header.
!=======================================================================
SUBROUTINE VUMAT(NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL, &
                 STEPTIME, TOTALTIME, DT, CMNAME, COORDMP, CHARLENGTH, &
                 PROPS, DENSITY, STRAININC, RELSPININC, &
                 TEMPOLD, STRETCHOLD, DEFGRADOLD, FIELDOLD, &
                 STRESSOLD, STATEOLD, ENERINTERNOLD, ENERINELASOLD, &
                 TEMPNEW, STRETCHNEW, DEFGRADNEW, FIELDNEW, &
                 STRESSNEW, STATENEW, ENERINTERNNEW, ENERINELASNEW)
  USE NLVE_SCHAPERY
  IMPLICIT NONE
  INTEGER :: NBLOCK, NDIR, NSHR, NSTATEV, NFIELDV, NPROPS, LANNEAL
  CHARACTER(LEN=80) :: CMNAME
  DOUBLE PRECISION :: STEPTIME, TOTALTIME, DT
  DOUBLE PRECISION :: PROPS(NPROPS), DENSITY(NBLOCK), COORDMP(NBLOCK, *), CHARLENGTH(NBLOCK)
  DOUBLE PRECISION :: STRAININC(NBLOCK, NDIR+NSHR), RELSPININC(NBLOCK, NSHR)
  DOUBLE PRECISION :: TEMPOLD(NBLOCK), STRETCHOLD(NBLOCK, NDIR+NSHR)
  DOUBLE PRECISION :: DEFGRADOLD(NBLOCK, NDIR+NSHR+NSHR), FIELDOLD(NBLOCK, NFIELDV)
  DOUBLE PRECISION :: STRESSOLD(NBLOCK, NDIR+NSHR), STATEOLD(NBLOCK, NSTATEV)
  DOUBLE PRECISION :: ENERINTERNOLD(NBLOCK), ENERINELASOLD(NBLOCK)
  DOUBLE PRECISION :: TEMPNEW(NBLOCK), STRETCHNEW(NBLOCK, NDIR+NSHR)
  DOUBLE PRECISION :: DEFGRADNEW(NBLOCK, NDIR+NSHR+NSHR), FIELDNEW(NBLOCK, NFIELDV)
  DOUBLE PRECISION :: STRESSNEW(NBLOCK, NDIR+NSHR), STATENEW(NBLOCK, NSTATEV)
  DOUBLE PRECISION :: ENERINTERNNEW(NBLOCK), ENERINELASNEW(NBLOCK)

  TYPE(NLVE_MAT)   :: MAT
  DOUBLE PRECISION :: SIG_OLD(6), EPS_OLD(6), DEPS(6), SIG_NEW(6), SDV_OLD(NLVE_NSDV)
  DOUBLE PRECISION :: SDV_NEW(NLVE_NSDV), CTAN(6, 6)
  INTEGER :: KM, NTENS, STATUS
  LOGICAL :: DATACHECK
  CHARACTER(LEN=120) :: TXT

  NTENS = NDIR + NSHR
  IF (NDIR /= 3 .OR. (NSHR /= 3 .AND. NSHR /= 1)) THEN
    WRITE(TXT, '(A,2I3,A)') 'unsupported element (NDIR, NSHR =', NDIR, NSHR, &
                            '): only 3D, plane strain, axisymmetric'
    CALL NLVE_FATAL(TXT)
  END IF
  IF (NSTATEV < NLVE_NSDV) THEN
    WRITE(TXT, '(A,I3,A,I3)') 'DEPVAR must be at least ', NLVE_NSDV, ' but is ', NSTATEV
    CALL NLVE_FATAL(TXT)
  END IF
  ! data check call (fictitious strains, stepTime = totalTime = 0): initial elastic response
  DATACHECK = (STEPTIME == 0.D0 .AND. TOTALTIME == 0.D0)
  CALL NLVE_SET_MATERIAL(PROPS, NPROPS, MAT, .NOT. DATACHECK)

  DO KM = 1, NBLOCK
    ! ---- Explicit -> internal (Standard order, engineering shear)
    SIG_OLD = 0.D0;  DEPS = 0.D0
    SIG_OLD(1:4) = STRESSOLD(KM, 1:4)
    DEPS(1:3)    = STRAININC(KM, 1:3)
    DEPS(4)      = 2.D0*STRAININC(KM, 4)
    IF (NSHR == 3) THEN
      SIG_OLD(5) = STRESSOLD(KM, 6);   DEPS(5) = 2.D0*STRAININC(KM, 6)    ! 13 <- 31
      SIG_OLD(6) = STRESSOLD(KM, 5);   DEPS(6) = 2.D0*STRAININC(KM, 5)    ! 23 <- 23
    END IF
    SDV_OLD = STATEOLD(KM, 1:NLVE_NSDV)
    EPS_OLD = SDV_OLD(24:29)

    IF (DATACHECK) THEN
      SIG_NEW = SIG_OLD + NLVE_ELASTIC(MAT, DEPS)
      SDV_NEW = SDV_OLD
    ELSE
      CALL NLVE_UPDATE(MAT, DT, SIG_OLD, EPS_OLD, DEPS, SDV_OLD, SIG_NEW, SDV_NEW, &
                       .FALSE., CTAN, STATUS)
      IF (STATUS /= 0) CALL NLVE_WARN('VUMAT local Newton not converged, residual =', &
                                      SDV_NEW(22), TOTALTIME, KM)
    END IF

    ! ---- internal -> Explicit
    STRESSNEW(KM, 1:4) = SIG_NEW(1:4)
    IF (NSHR == 3) THEN
      STRESSNEW(KM, 5) = SIG_NEW(6)        ! 23
      STRESSNEW(KM, 6) = SIG_NEW(5)        ! 31
    END IF
    STATENEW(KM, 1:NLVE_NSDV) = SDV_NEW
    IF (NSTATEV > NLVE_NSDV) STATENEW(KM, NLVE_NSDV+1:NSTATEV) = STATEOLD(KM, NLVE_NSDV+1:NSTATEV)

    ! ---- specific internal energy (stress power), needed for a meaningful energy balance
    ENERINTERNNEW(KM) = ENERINTERNOLD(KM) + 0.5D0*DOT_PRODUCT(SIG_OLD + SIG_NEW, DEPS) / DENSITY(KM)
    ENERINELASNEW(KM) = ENERINELASOLD(KM)
  END DO
END SUBROUTINE VUMAT
