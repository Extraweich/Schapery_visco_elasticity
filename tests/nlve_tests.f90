! =====================================================================
!  Single-material-point tests for nlve_schapery.f without Abaqus
!
!  Build and run (gfortran, no other libraries):
!     gfortran -O2 -ffree-form src/nlve_schapery.f tests/nlve_tests.f90 -o nlve_tests
!     ./nlve_tests
!  (./nlve_tests error   checks the input-error message and exit code)
! =====================================================================
module tst
  use nlve_schapery, only: NLVE_NSDV, NLVE_NPROPS
  implicit none
  integer, parameter :: dp = kind(1.d0), nsdv = NLVE_NSDV, np = NLVE_NPROPS
  integer :: nfail = 0, nit_total = 0, ninc_total = 0
  integer :: tangent_mode = 0          ! 0: DDSDDE as returned, 1: symmetric part
  integer :: seed = 12345
contains

  ! ------------------------------------------------------------------ checks
  subroutine check(name, val, tol)
    character(*), intent(in) :: name
    real(dp), intent(in) :: val, tol
    character(4) :: tag
    if (val <= tol) then
      tag = 'ok'
    else
      tag = 'FAIL'; nfail = nfail + 1
    end if
    write(*, '(4x,a4,2x,a,t78,es10.3,a,es8.1,a)') tag, name, val, '  (tol ', tol, ')'
  end subroutine

  real(dp) function rnd()            ! deterministic uniform [0,1)
    seed = mod(1103515245*seed + 12345, 2147483647)
    if (seed < 0) seed = seed + 2147483647
    rnd = real(seed, dp) / 2147483647.d0
  end function

  ! ------------------------------------------------------------------ materials
  function mat_diss() result(p)       ! dissertation, PA6 23 C / 10 % r.H. (nu assumed)
    real(dp) :: p(np)
    p = 0
    p(1:8)   = [2600.d0, 0.35d0, 11676.d0, 932.d0, 7509.d0, 21258.d0, 5434.d0, 201109.d0]
    p(9:13)  = [1.d0, 1.00d0, 3.49d0, 104.23d0, 19.57d0]
    p(14:18) = [1.d0, 0.98d0, 4.61d0,  25.19d0,  3.72d0]
    p(19:23) = [1.d0, 1.01d0, 4.22d0,  22.42d0,  0.26d0]
  end function
  function mat_example() result(p)    ! minimal example (.inp) with the 250309 functions
    real(dp) :: p(np)
    p = 0
    p(1:8)   = [1000.d0, 0.38d0, 500.d0, 10.d0, 200.d0, 50.d0, 999999.d0, 1000.d0]
    p(9:13)  = [1.d0, 0.9948d0, 5.2292d0, 19.389d0, 1.3612d0]
    p(14:18) = [1.d0, 0.9918d0, 5.3329d0, 30.573d0, 7.587d0]
    p(19:23) = [2.d0, 1.d0, 1.49d0, 28.349d0, 8.154d0]
  end function

  ! independent implementation of the parameter functions
  real(dp) function fnl(p, k, s)
    real(dp), intent(in) :: p(np), s
    integer, intent(in) :: k            ! 1 g0, 2 g1g2, 3 a_sigma
    real(dp) :: q(4)
    integer :: t
    t = nint(p(9 + 5*(k-1))); q = p(10 + 5*(k-1):13 + 5*(k-1))
    select case (t)
    case (1); fnl = q(4) + (q(1) - q(4)) / (1 + (s/q(3))**q(2))
    case (2); fnl = q(1) + q(2)*exp(-(s - q(3))**2/(2*q(4)**2))
    case default; fnl = 1
    end select
  end function

  ! closed-form creep compliance for a step load (diss. Eq. 52), shear=.true.: gamma/tau
  real(dp) function creep_J(p, t, sbar, shear)
    real(dp), intent(in) :: p(np), t, sbar
    logical, intent(in) :: shear
    real(dp) :: M(0:3), psi, f
    integer :: i
    M = [p(1), p(3), p(5), p(7)]
    f = 1; if (shear) f = 2*(1 + p(2))
    psi = t / fnl(p, 3, sbar)
    creep_J = fnl(p, 1, sbar) * f / M(0)
    do i = 1, 3
      creep_J = creep_J + fnl(p, 2, sbar) * f / M(i) * (1 - exp(-psi/p(2*i+2)))
    end do
  end function

  ! ------------------------------------------------------------------ conversions
  ! internal: Standard order 11,22,33,12,13,23 (+ engineering shear strain)
  ! explicit: 11,22,33,12,23,31 (tensor shear strain)
  function s_int2ex(s, nshr) result(x)
    real(dp), intent(in) :: s(6); integer, intent(in) :: nshr
    real(dp) :: x(3+nshr)
    if (nshr == 3) then; x = [s(1), s(2), s(3), s(4), s(6), s(5)]; else; x = s(1:4); end if
  end function
  function s_ex2int(x, nshr) result(s)
    integer, intent(in) :: nshr; real(dp), intent(in) :: x(3+nshr)
    real(dp) :: s(6)
    s = 0
    if (nshr == 3) then; s = [x(1), x(2), x(3), x(4), x(6), x(5)]; else; s(1:4) = x; end if
  end function
  function e_int2ex(e, nshr) result(x)
    real(dp), intent(in) :: e(6); integer, intent(in) :: nshr
    real(dp) :: x(3+nshr)
    if (nshr == 3) then; x = [e(1), e(2), e(3), e(4)/2, e(6)/2, e(5)/2]; else; x = [e(1:3), e(4)/2]; end if
  end function

  ! ------------------------------------------------------------------ solver calls
  subroutine umat_call(ntens, p, stress0, statev0, stran, dstran, dt, ttot, stress, statev, C, pnew, drot_in)
    integer, intent(in) :: ntens
    real(dp), intent(in) :: p(np), stress0(ntens), statev0(nsdv), stran(ntens), dstran(ntens), dt, ttot
    real(dp), intent(out) :: stress(ntens), statev(nsdv), C(ntens, ntens), pnew
    real(dp), intent(in), optional :: drot_in(3, 3)
    real(dp) :: sse, spd, scd, rpl, drpldt, ddsddt(ntens), drplde(ntens), predef(1), dpred(1)
    real(dp) :: time(2), temp, dtemp, coords(3), drot(3, 3), celent, dfgrd0(3, 3), dfgrd1(3, 3)
    character(80) :: cmname
    integer :: i
    stress = stress0; statev = statev0; pnew = 1
    sse = 0; spd = 0; scd = 0; rpl = 0; drpldt = 0; ddsddt = 0; drplde = 0; predef = 0; dpred = 0
    time = [ttot, ttot]; temp = 23; dtemp = 0; coords = 0; celent = 1; cmname = 'PA6'
    drot = 0; dfgrd0 = 0; dfgrd1 = 0
    do i = 1, 3; drot(i, i) = 1; dfgrd0(i, i) = 1; dfgrd1(i, i) = 1; end do
    if (present(drot_in)) drot = drot_in
    call UMAT(stress, statev, C, sse, spd, scd, rpl, ddsddt, drplde, drpldt, stran, dstran, time, dt, &
              temp, dtemp, predef, dpred, cmname, 3, ntens-3, ntens, nsdv, p, np, coords, drot, pnew, &
              celent, dfgrd0, dfgrd1, 1, 1, 0, 0, 1, 1)
  end subroutine

  ! VUMAT for a block of nb points (explicit conventions in and out)
  subroutine vumat_call(nb, nshr, p, dt, steptime, totaltime, sold, stold, dstr, snew, stnew, eint)
    integer, intent(in) :: nb, nshr
    real(dp), intent(in) :: p(np), dt, steptime, totaltime, sold(nb, 3+nshr), stold(nb, nsdv), dstr(nb, 3+nshr)
    real(dp), intent(out) :: snew(nb, 3+nshr), stnew(nb, nsdv), eint(nb)
    real(dp) :: density(nb), coordmp(nb, 3), charlength(nb), relspininc(nb, nshr), tempold(nb)
    real(dp) :: stretchold(nb, 3+nshr), defgradold(nb, 3+2*nshr), fieldold(nb, 1), enerinternold(nb)
    real(dp) :: enerinelasold(nb), tempnew(nb), stretchnew(nb, 3+nshr), defgradnew(nb, 3+2*nshr)
    real(dp) :: fieldnew(nb, 1), enerinelasnew(nb), sold_c(nb, 3+nshr), stold_c(nb, nsdv), dstr_c(nb, 3+nshr)
    character(80) :: cmname
    density = 1.2d-9; coordmp = 0; charlength = 1; relspininc = 0; tempold = 23; tempnew = 23
    stretchold = 0; defgradold = 0; fieldold = 0; enerinternold = 0; enerinelasold = 0
    stretchnew = 0; defgradnew = 0; fieldnew = 0; cmname = 'PA6'
    sold_c = sold; stold_c = stold; dstr_c = dstr
    snew = -999; stnew = -999                               ! must be overwritten
    call VUMAT(nb, 3, nshr, nsdv, 1, np, 0, steptime, totaltime, dt, cmname, coordmp, charlength, &
               p, density, dstr_c, relspininc, tempold, stretchold, defgradold, fieldold, &
               sold_c, stold_c, enerinternold, enerinelasold, tempnew, stretchnew, defgradnew, fieldnew, &
               snew, stnew, eint, enerinelasnew)
  end subroutine

  ! stress-controlled UMAT increment (global Newton with DDSDDE)
  subroutine stress_inc(ntens, p, target, stress, statev, stran, ttot, dt, tol)
    integer, intent(in) :: ntens
    real(dp), intent(in) :: p(np), target(ntens), dt, tol
    real(dp), intent(inout) :: stress(ntens), statev(nsdv), stran(ntens), ttot
    real(dp) :: de(ntens), s(ntens), sv(nsdv), C(ntens, ntens), r(ntens), pnew, A(6, 6), b(6), x(6)
    integer :: it
    logical :: ok
    de = 0
    do it = 1, 200
      call umat_call(ntens, p, stress, statev, stran, de, dt, ttot, s, sv, C, pnew)
      r = s - target
      if (maxval(abs(r)) < tol*max(1.d0, maxval(abs(target)))) exit
      nit_total = nit_total + 1
      if (tangent_mode == 1) C = 0.5d0*(C + transpose(C))
      A = 0; b = 0; A(1:ntens, 1:ntens) = C; b(1:ntens) = -r
      if (ntens < 6) then; A(5, 5) = 1; A(6, 6) = 1; end if
      call gauss6(A, b, x, ok)
      de = de + x(1:ntens)
    end do
    ninc_total = ninc_total + 1
    stress = s; statev = sv; stran = stran + de; ttot = ttot + dt
  end subroutine

  subroutine gauss6(A_in, b, x, ok)
    real(dp), intent(in) :: A_in(6, 6), b(6)
    real(dp), intent(out) :: x(6)
    logical, intent(out) :: ok
    real(dp) :: A(6, 6), f, tmp(6), t
    integer :: i, j, k
    A = A_in; x = b; ok = .false.
    do i = 1, 6
      k = i - 1 + maxloc(abs(A(i:6, i)), 1)
      if (abs(A(k, i)) < 1d-300) return
      tmp = A(i, :); A(i, :) = A(k, :); A(k, :) = tmp
      t = x(i); x(i) = x(k); x(k) = t
      do j = i + 1, 6
        f = A(j, i)/A(i, i); A(j, :) = A(j, :) - f*A(i, :); x(j) = x(j) - f*x(i)
      end do
    end do
    do i = 6, 1, -1
      x(i) = (x(i) - dot_product(A(i, i+1:6), x(i+1:6)))/A(i, i)
    end do
    ok = .true.
  end subroutine

  subroutine init_sdv(sv)
    real(dp), intent(out) :: sv(nsdv)
    sv = 0                                     ! Abaqus default (no SDVINI)
  end subroutine

  ! random strain increment (engineering, internal order); 2D: 13, 23 = 0
  function rnd_deps(scale, nshr) result(de)
    real(dp), intent(in) :: scale; integer, intent(in) :: nshr
    real(dp) :: de(6)
    integer :: i
    do i = 1, 6; de(i) = scale*(2*rnd() - 1); end do
    de(1) = de(1) + 0.6d0*scale                 ! drift: keep the point loaded
    if (nshr == 1) de(5:6) = 0
  end function

  real(dp) function rnd_dt(lo, hi)
    real(dp), intent(in) :: lo, hi
    rnd_dt = lo*(hi/lo)**rnd()
  end function
end module tst

! =====================================================================
program nlve_tests
  use tst
  implicit none
  real(dp) :: p(np), pa(np), pb(np), stress(6), statev(nsdv), stran(6), ttot, target(6), Jn, Ja, err, lat
  real(dp) :: s1(6), sv1(nsdv), s2(6), sv2(nsdv), C(6, 6), Cn(6, 6), sp(6), sm(6), de(6), pnew, dt
  real(dp) :: e3(6), e4(6), sig3(6), sig4(6), sv3(nsdv), sv4(nsdv), s4(4), C4(4, 4), sdv_err, sig_err
  real(dp) :: xs(1, 6), xsv(1, nsdv), xe(1, 6), xn(1, 6), xsvn(1, nsdv), ei(1), R(3, 3), ang, hh
  real(dp) :: bs(4, 6), bsv(4, nsdv), bde(4, 6), bsn(4, 6), bsvn(4, nsdv), bei(4), ref(4, 6)
  real(dp) :: x4(1, 4), xsv4(1, nsdv), xe4(1, 4), xn4(1, 4), loads(4), smax
  integer :: i, k, n, il, nshr, im, ib
  character(len=16) :: arg
  character(len=8) :: mname(2) = ['diss    ', 'example ']

  pa = mat_diss(); pb = mat_example()

  call get_command_argument(1, arg)
  if (arg == 'error') then               ! plane stress element -> ERROR STOP with message
    stress = 0; stran = 0; de = 0; statev = 0
    call UMAT(stress, statev, C, err, err, err, err, s1, s2, err, stran, de, [0.d0, 0.d0], 1.d0, &
              err, err, s1, s1, 'PA6', 2, 1, 3, nsdv, pa, np, s1, R, pnew, err, R, R, 1, 1, 0, 0, 1, 1)
    stop 'not reached'
  end if

  ! ------------------------------------------------------------------ [1] creep
  print '(a)', '[1] uniaxial creep vs. closed-form Schapery solution (Eq. 52), UMAT 3D and axisymmetric'
  do im = 1, 2
    if (im == 1) then; p = pa; loads = [10, 20, 30, 40]; else; p = pb; loads = [5, 15, 25, 35]; end if
    do n = 4, 6, 2                                       ! ntens = 4 (axisymmetric), 6 (3D)
      err = 0; lat = 0
      do il = 1, 4
        call init_sdv(statev); stress = 0; stran = 0; ttot = 0
        target = 0; target(2) = loads(il)                ! load along 2 (axis of symmetry)
        call stress_inc(n, p, target(1:n), stress(1:n), statev, stran(1:n), ttot, 1.d-3, 1.d-11)
        dt = 1.d-2*p(4)
        do while (ttot < 50*p(6))
          call stress_inc(n, p, target(1:n), stress(1:n), statev, stran(1:n), ttot, dt, 1.d-11)
          Jn = stran(2)/loads(il); Ja = creep_J(p, ttot - 0.5d-3, loads(il), .false.)
          err = max(err, abs(Jn - Ja)/Ja)
          lat = max(lat, abs(stran(1)/stran(2) + p(2)), abs(stran(3)/stran(2) + p(2)))
          dt = min(1.2d0*dt, 2*p(6))
        end do
      end do
      call check('material '//trim(mname(im))//', ntens='//achar(48+n)//': max rel. error J(t), 4 load levels', &
                 err, 1.d-7)
      call check('material '//trim(mname(im))//', ntens='//achar(48+n)//': max |eps_lat/eps_axial + nu|', &
                 lat, 1.d-10)
    end do
  end do

  print '(a)', '[2] pure shear creep (engineering shear factors), UMAT 3D'
  do im = 1, 2
    if (im == 1) then; p = pa; else; p = pb; end if
    err = 0
    do k = 4, 6                                          ! tau12, tau13, tau23
      call init_sdv(statev); stress = 0; stran = 0; ttot = 0
      target = 0; target(k) = 20.d0/sqrt(3.d0)
      call stress_inc(6, p, target, stress, statev, stran, ttot, 1.d-3, 1.d-11)
      dt = 1.d-2*p(4)
      do while (ttot < 50*p(6))
        call stress_inc(6, p, target, stress, statev, stran, ttot, dt, 1.d-11)
        Jn = stran(k)/target(k); Ja = creep_J(p, ttot - 0.5d-3, 20.d0, .true.)
        err = max(err, abs(Jn - Ja)/Ja)
        dt = min(1.2d0*dt, 2*p(6))
      end do
    end do
    call check('material '//trim(mname(im))//': max rel. error gamma/tau (12, 13, 23)', err, 1.d-7)
  end do

  ! ------------------------------------------------------------------ [3]-[5] equivalence
  print '(a)', '[3] UMAT vs VUMAT on random multiaxial strain paths (300 increments)'
  do im = 1, 2
    if (im == 1) then; p = pa; else; p = pb; end if
    do nshr = 3, 1, -2
      n = 3 + nshr
      call init_sdv(sv3); sig3 = 0; e3 = 0; ttot = 0
      xsv(1, :) = 0; xs = 0; sig_err = 0; sdv_err = 0; smax = 0
      do k = 1, 300
        de = rnd_deps(merge(1.5d-3, 3.d-3, im == 1), nshr)
        dt = rnd_dt(1.d-3*p(4), 0.5d0*p(6))
        ! UMAT
        call umat_call(n, p, sig3(1:n), sv3, e3(1:n), de(1:n), dt, ttot, sig4(1:n), sv4, C(1:n, 1:n), pnew)
        ! VUMAT
        xe(1, 1:n) = e_int2ex(de, nshr)
        call vumat_call(1, nshr, p, dt, ttot + dt, ttot + dt, xs(:, 1:n), xsv, xe(:, 1:n), xn(:, 1:n), xsvn, ei)
        s1 = s_ex2int(xn(1, 1:n), nshr)
        sig_err = max(sig_err, maxval(abs(s1(1:n) - sig4(1:n))))
        sdv_err = max(sdv_err, maxval(abs(xsvn(1, :) - sv4)/max(abs(sv4), 1.d-3)))
        smax = max(smax, maxval(abs(sig4(1:n))))
        sig3(1:n) = sig4(1:n); sv3 = sv4; e3(1:n) = e3(1:n) + de(1:n); ttot = ttot + dt
        xs(:, 1:n) = xn(:, 1:n); xsv = xsvn
      end do
      call check('material '//trim(mname(im))//', nshr='//achar(48+nshr)//': max |sigma_UMAT - sigma_VUMAT| / max|sigma|', &
                 sig_err/smax, 1.d-10)
      call check('material '//trim(mname(im))//', nshr='//achar(48+nshr)//': max rel. difference of all 29 SDVs', &
                 sdv_err, 1.d-8)
    end do
  end do

  print '(a)', '[4] UMAT axisymmetric (NTENS=4) vs UMAT 3D with eps13 = eps23 = 0'
  do im = 1, 2
    if (im == 1) then; p = pa; else; p = pb; end if
    call init_sdv(sv3); call init_sdv(sv4); sig3 = 0; sig4 = 0; e3 = 0; ttot = 0; sig_err = 0; sdv_err = 0
    do k = 1, 300
      de = rnd_deps(merge(1.5d-3, 3.d-3, im == 1), 1)
      dt = rnd_dt(1.d-3*p(4), 0.5d0*p(6))
      call umat_call(6, p, sig3, sv3, e3, de, dt, ttot, s1, sv1, C, pnew)
      call umat_call(4, p, sig4(1:4), sv4, e3(1:4), de(1:4), dt, ttot, s4, sv2, C4, pnew)
      sig_err = max(sig_err, maxval(abs(s4 - s1(1:4))), maxval(abs(s1(5:6))))
      sdv_err = max(sdv_err, maxval(abs(sv2 - sv1)))
      sig3 = s1; sv3 = sv1; sig4(1:4) = s4; sv4 = sv2; e3 = e3 + de; ttot = ttot + dt
    end do
    call check('material '//trim(mname(im))//': max |stress difference| (incl. sigma13, sigma23 of 3D)', sig_err, 1.d-9)
    call check('material '//trim(mname(im))//': max |SDV difference|', sdv_err, 1.d-12)
  end do

  print '(a)', '[5] VUMAT with nblock = 4 vs. four single-point calls'
  p = pb
  bsv = 0; bs = 0; err = 0
  do k = 1, 100
    dt = rnd_dt(1.d-2, 5.d0)
    do ib = 1, 4; bde(ib, :) = e_int2ex(rnd_deps(3.d-3, 3), 3); end do
    call vumat_call(4, 3, p, dt, 1.d0 + k, 1.d0 + k, bs, bsv, bde, bsn, bsvn, bei)
    do ib = 1, 4
      call vumat_call(1, 3, p, dt, 1.d0 + k, 1.d0 + k, bs(ib:ib, :), bsv(ib:ib, :), bde(ib:ib, :), xn, xsvn, ei)
      err = max(err, maxval(abs(xn(1, :) - bsn(ib, :))), maxval(abs(xsvn(1, :) - bsvn(ib, :))))
    end do
    bs = bsn; bsv = bsvn
  end do
  call check('max difference stress/SDV block vs single', err, 0.d0)

  ! ------------------------------------------------------------------ [6] import
  print '(a)', '[6] implicit -> explicit: UMAT for 150 increments, then VUMAT (vs. UMAT throughout)'
  do nshr = 3, 1, -2
    p = pb; n = 3 + nshr
    call init_sdv(sv3); sig3 = 0; e3 = 0; ttot = 0; sig_err = 0; sdv_err = 0; smax = 0
    do k = 1, 300
      de = rnd_deps(3.d-3, nshr)
      dt = rnd_dt(1.d-2, 5.d0)
      call umat_call(n, p, sig3(1:n), sv3, e3(1:n), de(1:n), dt, ttot, sig4(1:n), sv4, C(1:n, 1:n), pnew)
      if (k == 150) then                   ! *IMPORT: stress and SDVs handed over
        xs(1, 1:n) = s_int2ex(sig4, nshr); xsv(1, :) = sv4
      else if (k > 150) then
        xe(1, 1:n) = e_int2ex(de, nshr)
        call vumat_call(1, nshr, p, dt, ttot + dt, ttot + dt, xs(:, 1:n), xsv, xe(:, 1:n), xn(:, 1:n), xsvn, ei)
        s1 = s_ex2int(xn(1, 1:n), nshr)
        sig_err = max(sig_err, maxval(abs(s1(1:n) - sig4(1:n))))
        sdv_err = max(sdv_err, maxval(abs(xsvn(1, :) - sv4)/max(abs(sv4), 1.d-3)))
        xs(:, 1:n) = xn(:, 1:n); xsv = xsvn
      end if
      smax = max(smax, maxval(abs(sig4(1:n))))
      sig3(1:n) = sig4(1:n); sv3 = sv4; e3(1:n) = e3(1:n) + de(1:n); ttot = ttot + dt
    end do
    call check('nshr='//achar(48+nshr)//': max |sigma difference| / max|sigma| after import', sig_err/smax, 1.d-10)
    call check('nshr='//achar(48+nshr)//': max rel. SDV difference after import', sdv_err, 1.d-8)
  end do

  ! ------------------------------------------------------------------ [7] NLGEOM
  print '(a)', '[7] NLGEOM: rigid rotation (DROT) during a creep increment, UMAT 3D'
  p = pa
  call init_sdv(statev); stress = 0; stran = 0; ttot = 0
  target = 0; target(1) = 30
  call stress_inc(6, p, target, stress, statev, stran, ttot, 1.d-3, 1.d-11)
  dt = 1
  do while (ttot < 3600); call stress_inc(6, p, target, stress, statev, stran, ttot, dt, 1.d-11); dt = 1.3d0*dt; end do
  de = 0
  call umat_call(6, p, stress, statev, stran, de, 600.d0, ttot, s1, sv1, C, pnew)
  ang = 0.7d0                                            ! rotation about axis (1,2,2)/3
  R = rot_axis([1.d0, 2.d0, 2.d0]/3.d0, ang)
  call umat_call(6, p, rot_s(stress, R), statev, rot_e(stran, R), de, 600.d0, ttot, s2, sv2, C, pnew, R)
  call check('||sigma(rotated run) - R sigma R^T|| / ||sigma||', norm2(s2 - rot_s(s1, R))/norm2(s1), 1.d-12)

  ! ------------------------------------------------------------------ [8]-[9] edge cases
  print '(a)', '[8] zero time increment and VUMAT data-check call'
  p = pb
  call init_sdv(sv3); sig3 = 0; e3 = 0
  de = 0; de(1) = 0.01d0
  call umat_call(6, p, sig3, sv3, e3, de, 0.d0, 0.d0, s1, sv1, C, pnew)
  call check('UMAT dtime=0: |sigma11 - g0-scaled instantaneous value| (finite, no NaN)', &
             abs(s1(1) - inst_uniax_strain(p, 0.01d0)), 1.d-6)
  call check('UMAT dtime=0: SDVs unchanged', maxval(abs(sv1 - sv3)), 0.d0)
  xs = 0; xsv = 0; xsv(1, 19:21) = 7; xe = 0; xe(1, 1) = 1.d-3
  call vumat_call(1, 3, p, 0.d0, 0.d0, 0.d0, xs, xsv, xe, xn, xsvn, ei)
  call check('VUMAT data check: |sigma11 - (K0 + 4/3 G0) eps11|', &
             abs(xn(1, 1) - (p(1)/(3*(1 - 2*p(2))) + 4.d0/3*p(1)/(2*(1 + p(2))))*1.d-3), 1.d-9)
  call check('VUMAT data check: SDVs unchanged', maxval(abs(xsvn - xsv)), 0.d0)

  print '(a)', '[9] relaxation hold with dstran = 0 exactly'
  p = pb
  call init_sdv(sv3); sig3 = 0; e3 = 0; ttot = 0
  de = 0; de(1) = 0.02d0
  call umat_call(6, p, sig3, sv3, e3, de, 1.d-3, ttot, s1, sv1, C, pnew)
  sig3 = s1; sv3 = sv1; e3 = e3 + de; ttot = 1.d-3; de = 0; err = 0; smax = s1(1)
  do k = 1, 20
    call umat_call(6, p, sig3, sv3, e3, de, 5.d0, ttot, s1, sv1, C, pnew)
    err = max(err, sv1(22), 1 - pnew)
    sig3 = s1; sv3 = sv1; ttot = ttot + 5
  end do
  call check('max(residual SDV22, 1 - PNEWDT) over 20 hold increments', err, 1.d-6)
  call check('stress relaxes during the hold (1 if not)', merge(0.d0, 1.d0, s1(1) < smax), 0.d0)

  ! ------------------------------------------------------------------ [10] tangent
  print '(a)', '[10] DDSDDE vs. numerical tangent; global Newton iterations per increment (creep)'
  p = pa
  call init_sdv(statev); stress = 0; stran = 0; ttot = 0
  target = 0; target(1) = 40
  call stress_inc(6, p, target, stress, statev, stran, ttot, 1.d-3, 1.d-11)
  dt = 1
  do while (ttot < 3600); call stress_inc(6, p, target, stress, statev, stran, ttot, dt, 1.d-11); dt = 1.3d0*dt; end do
  de = 0; de(1) = 2.d-4; de(2:3) = -0.7d-4; hh = 1.d-7
  do i = 1, 6
    de(i) = de(i) + hh;   call umat_call(6, p, stress, statev, stran, de, 600.d0, ttot, sp, sv1, C, pnew)
    de(i) = de(i) - 2*hh; call umat_call(6, p, stress, statev, stran, de, 600.d0, ttot, sm, sv1, C, pnew)
    de(i) = de(i) + hh;   Cn(:, i) = (sp - sm)/(2*hh)
  end do
  call umat_call(6, p, stress, statev, stran, de, 600.d0, ttot, s1, sv1, C, pnew)
  call check('||DDSDDE - numerical tangent|| / ||numerical tangent||', norm2(C - Cn)/norm2(Cn), 1.d-5)
  do k = 0, 1
    tangent_mode = k; nit_total = 0; ninc_total = 0
    do il = 1, 4
      call init_sdv(statev); stress = 0; stran = 0; ttot = 0
      target = 0; target(1) = 10*il
      call stress_inc(6, p, target, stress, statev, stran, ttot, 1.d-3, 1.d-6)
      dt = 1
      do while (ttot < 576000); call stress_inc(6, p, target, stress, statev, stran, ttot, dt, 1.d-6)
        dt = min(1.15d0*dt, 2.d4); end do
    end do
    print '(4x,a,f5.2)', merge('full DDSDDE (UNSYMM):         ', 'symmetric part (no UNSYMM):   ', k == 0), &
          real(nit_total, dp)/ninc_total
  end do
  tangent_mode = 0

  ! ------------------------------------------------------------------ [11] explicit time steps
  print '(a)', '[11] relaxation with explicit-size time increments (convergence in dt)'
  p = pb
  do k = 4, 6
    dt = 10.d0**(-k)
    xs = 0; xsv = 0; xe = 0; xe(1, 1) = 0.01d0
    call vumat_call(1, 3, p, dt, dt, dt, xs, xsv, xe, xn, xsvn, ei)
    xs = xn; xsv = xsvn; xe = 0
    do i = 1, nint(0.1d0/dt)
      call vumat_call(1, 3, p, dt, i*dt, i*dt, xs, xsv, xe, xn, xsvn, ei)
      xs = xn; xsv = xsvn
    end do
    if (k == 4) then; sm(1) = xs(1, 1); else; sp(k-4) = xs(1, 1); end if
  end do
  call check('|sigma11(dt=1e-5) - sigma11(dt=1e-4)| / sigma11', abs(sp(1) - sm(1))/sm(1), 1.d-4)
  call check('|sigma11(dt=1e-6) - sigma11(dt=1e-5)| / sigma11 (should shrink)', abs(sp(2) - sp(1))/sm(1), 1.d-5)

  print '(a)', ''
  if (nfail == 0) then
    print '(a)', 'ALL TESTS PASSED'
  else
    print '(a,i3,a)', 'FAILED: ', nfail, ' test(s)'
    error stop 1
  end if

contains

  function rot_axis(a, phi) result(Q)
    real(dp), intent(in) :: a(3), phi
    real(dp) :: Q(3, 3), K(3, 3)
    integer :: i
    K = reshape([0.d0, a(3), -a(2), -a(3), 0.d0, a(1), a(2), -a(1), 0.d0], [3, 3])
    Q = sin(phi)*K + (1 - cos(phi))*matmul(K, K)
    do i = 1, 3; Q(i, i) = Q(i, i) + 1; end do
  end function
  function rot_s(v, Q) result(w)          ! stress vector
    real(dp), intent(in) :: v(6), Q(3, 3)
    real(dp) :: w(6), T(3, 3)
    T = reshape([v(1), v(4), v(5), v(4), v(2), v(6), v(5), v(6), v(3)], [3, 3])
    T = matmul(Q, matmul(T, transpose(Q)))
    w = [T(1, 1), T(2, 2), T(3, 3), T(1, 2), T(1, 3), T(2, 3)]
  end function
  function rot_e(v, Q) result(w)          ! engineering strain vector
    real(dp), intent(in) :: v(6), Q(3, 3)
    real(dp) :: w(6), u(6)
    u = v; u(4:6) = 0.5d0*u(4:6)
    w = rot_s(u, Q); w(4:6) = 2*w(4:6)
  end function

  ! sigma11 for an instantaneous uniaxial-strain jump (dt = 0): sigma = C0/g0(sbar) eps, solved for sbar
  real(dp) function inst_uniax_strain(p, e) result(s11)
    real(dp), intent(in) :: p(np), e
    real(dp) :: G0, K0, sb, g
    integer :: it
    G0 = p(1)/(2*(1 + p(2))); K0 = p(1)/(3*(1 - 2*p(2)))
    sb = 2*G0*e
    do it = 1, 200
      g = fnl(p, 1, sb); sb = 2*G0*e/g
    end do
    s11 = (K0 + 4.d0/3*G0)/g*e             ! g0 scales bulk and shear compliance
  end function
end program
