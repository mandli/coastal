! ==============================================================================
!  Riemann solver for the one-dimensional shallow water equations
!  
!  Uses an fwave HLL type of approach and includes bathymetry terms.  The
!  common block must contain values for the acceleration due to gravity `g` and
!  the tolerance used to decide whether a cell is dry or not `dry_tolerance`.
! ==============================================================================
subroutine rp1(maxmx, num_eqn, num_waves, num_aux, num_ghost, mx, ql, qr,      &
               auxl, auxr, fwave, s, amdq, apdq)

    implicit none

    ! Input
    integer, intent(in) :: maxmx, num_eqn, num_waves, num_ghost, num_aux, mx
    real(kind=8), intent(in out) :: ql(num_eqn,1-num_ghost:maxmx+num_ghost)
    real(kind=8), intent(in out) :: qr(num_eqn,1-num_ghost:maxmx+num_ghost)
    real(kind=8), intent(in out) :: auxl(num_aux,1-num_ghost:maxmx+num_ghost)
    real(kind=8), intent(in out) :: auxr(num_aux,1-num_ghost:maxmx+num_ghost)

    ! Output
    real(kind=8), intent(in out) :: s(num_waves, 1-num_ghost:maxmx+num_ghost)
    real(kind=8), intent(in out) :: fwave(num_eqn, num_waves, 1-num_ghost:maxmx+num_ghost)
    real(kind=8), intent(in out) :: amdq(num_eqn, 1-num_ghost:maxmx+num_ghost)
    real(kind=8), intent(in out) :: apdq(num_eqn, 1-num_ghost:maxmx+num_ghost)

    ! Locals
    integer :: i, mw
    real(kind=8) :: h_l, h_r, hu_l, hu_r, u_l, u_r, b_l, b_r, phi_l, phi_r
    real(kind=8) :: delta(2), beta(2)


    ! Common block
    integer :: mcapa
    real(kind=8) :: g, dry_tolerance, earth_radius, deg2rad
    common /cparam/ g, dry_tolerance, earth_radius, deg2rad, mcapa

    fwave = 0.d0
    s = 0.d0
    amdq = 0.d0
    apdq = 0.d0

    do i=2-num_ghost,mx+num_ghost
        ! Extract state
        h_l = qr(1, i-1)
        h_r = ql(1, i)
        hu_l = qr(2, i-1)
        hu_r = ql(2, i)
        b_l = auxr(1, i-1)
        b_r = auxl(1, i)

        if (h_l < dry_tolerance) then
            u_l = 0.d0
        else
            u_l = hu_l / h_l
        end if

        if (h_r < dry_tolerance) then
            u_r = 0.d0
        else
            u_r = hu_r / h_r
        end if
        
        phi_l = h_l * u_l**2 + 0.5d0 * g * h_l**2
        phi_r = h_r * u_r**2 + 0.5d0 * g * h_r**2

        ! Speeds
        s(1, i) = u_l - sqrt(g * h_l)
        s(2, i) = u_r + sqrt(g * h_r)

        ! Solve for wave strengths
        delta(1) = hu_r - hu_l
        delta(2) = phi_r - phi_l + g * 0.5d0 * (h_r + h_l) * (b_r - b_l)

        beta(1) = (s(2, i) * delta(1) - delta(2)) / (s(2, i) - s(1, i))
        beta(2) = (delta(2) - s(1, i) * delta(1)) / (s(2, i) - s(1, i))

        ! Construct waves
        fwave(1, 1, i) = beta(1)
        fwave(2, 1, i) = beta(1) * s(1, i)

        fwave(1, 2, i) = beta(2)
        fwave(2, 2, i) = beta(2) * s(2, i)


        ! Calculate accumulated flux
        do mw=1,num_waves
            if (s(mw, i) < 0.d0) then
                amdq(:, i) = amdq(:, i) + fwave(:, mw, i)
            else if (s(mw, i) > 0.d0) then
                apdq(:, i) = apdq(:, i) + fwave(:, mw, i)
            else
                amdq(:, i) = amdq(:, i) + 0.5d0 * fwave(:, mw, i)
                apdq(:, i) = apdq(:, i) + 0.5d0 * fwave(:, mw, i)
            end if
        end do
    enddo

end subroutine rp1