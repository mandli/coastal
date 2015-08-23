subroutine rp1(maxmx, num_eqn, num_waves, num_aux, num_ghost, mx, ql, qr,      &
               auxl, auxr, fwave, s, amdq, apdq)
!======================================================================
!
! Solves normal Riemann problems for the 2D SHALLOW WATER equations
!     with topography:
!     #        h_t + (hu)_x + (hv)_y = 0                           #
!     #        (hu)_t + (hu^2 + 0.5gh^2)_x + (huv)_y = -ghb_x      #
!     #        (hv)_t + (huv)_x + (hv^2 + 0.5gh^2)_y = -ghb_y      #
!! On input, ql contains the state vector at the left edge of each cell
!     qr contains the state vector at the right edge of each cell
!
! This data is along a slice in the x-direction if ixy=1
!     or the y-direction if ixy=2.
!!  Note that the i'th Riemann problem has left state qr(i-1,:)
!     and right state ql(i,:)
!  From the basic clawpack routines, this routine is called with
!     ql = qr
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                           !
!      # This Riemann solver is for the shallow water equations.            !
!                                                                           !
!       It allows the user to easily select a Riemann solver in             !
!       riemannsolvers_geo.f. this routine initializes all the variables    !
!       for the shallow water equations, accounting for wet dry boundary    !
!       dry cells, wave speeds etc.                                         !
!                                                                           !
!           David George, Vancouver WA, Feb. 2009                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

    ! Local storage
    integer :: i, n, m

    real(kind=8) :: h_l, h_r, hu_l, hu_r, u_l, u_r, b_l, b_r, phi_l, phi_r

    real(kind=8) :: s1m, s2m, rare1, rare2, h_star

    real(kind=8) :: s_l, s_r, u_hat, c_hat, s_roe_l, s_roe_r, s_einfeldt_l
    real(kind=8) :: s_einfeldt_r

    real(kind=8) :: wall(num_waves), sw(num_waves), fw(num_eqn, num_waves)

    ! Parameters
    integer, parameter :: MAX_ITERATIONS = 1

    ! Common block
    integer :: mcapa
    real(kind=8) :: g, dry_tolerance, earth_radius, deg2rad
    common /cparam/ g, dry_tolerance, earth_radius, deg2rad, mcapa

    fwave = 0.d0
    s = 0.d0
    amdq = 0.d0
    apdq = 0.d0

    ! Primary loop over cell edges
    do i = 2 - num_ghost, mx + num_ghost

        ! Extract state
        h_l = qr(1, i-1)
        h_r = ql(1, i)
        hu_l = qr(2, i-1)
        hu_r = ql(2, i)
        b_l = auxr(1, i-1)
        b_r = auxl(1, i)

        ! Check for negative values
        if (h_l < 0.d0 .or. h_r < 0.d0) then
            print *, "Negative input detected: "
            print *, "    (h_l, h_r, i) = (", h_l, ", ", h_r, ", ", i, ")"
            stop
        end if

        ! Skip problem if completely dry
        if (h_l < dry_tolerance .and. h_r < dry_tolerance) then
            cycle
        end if

        ! Extract derived quantities
        if (h_l > dry_tolerance) then
            u_l = hu_l / h_l
            phi_l = 0.5d0 * g * h_l**2 + hu_l**2 / h_l
        else
            h_l = 0.d0
            u_l = 0.d0
            phi_l = 0.d0
        end if
        if (h_r > dry_tolerance) then
            u_r = hu_r / h_r
            phi_r = 0.5d0 * g * h_r**2 + hu_r**2 / h_r
        else
            h_r = 0.d0
            u_r = 0.d0
            phi_r = 0.d0
        end if

        ! Setup different dry-state problems
        wall = 1.d0
        if (h_r <= dry_tolerance) then
            ! Test to see if the current wave would overcome the wall dry-state
            call riemanntype(h_l, h_l, u_l, -u_l, h_star, s1m, s2m, rare1,    &
                                rare2, 1, dry_tolerance, g)
            
            if (max(h_l, h_star) + b_l < b_r) then
                ! Right state should become ghost values that mirror left state
                ! so that a wall reflection can happen
                wall(2:) = 0.d0
                h_r = h_l
                hu_r = -hu_l
                b_r = b_l
                phi_r = phi_l
                u_r = -u_l
            else if (h_l + b_l < b_r) then
                ! Special case of the wall bc?
                b_r = h_l + b_l
            end if

        else if (h_l <= dry_tolerance) then
            ! Test to see if the current wave would overcome the wall dry-state
            call riemanntype(h_r, h_r, u_r, -u_r, h_star, s1m, s2m, rare1,    &
                                rare2, 1, dry_tolerance, g)
            
            if (max(h_r, h_star) + b_r < b_l) then
                ! Right state should become ghost values that mirror left state
                ! so that a wall reflection can happen
                wall(2:) = 0.d0
                h_L = h_r
                hu_L = -hu_r
                b_L = b_r
                phi_L = phi_r
                u_L = -u_r
            else if (h_r + b_r < b_l) then
                ! Special case of the wall bc?
                b_l = h_r + b_r
            end if

        end if

        ! Calculate relevant wave speeds estimates
        ! Basic wave bounds
        s_l = u_l - sqrt(g * h_l)
        s_r = u_r - sqrt(g * h_r)

        ! Roe speeds
        u_hat = (sqrt(g * h_l) * u_l + sqrt(g * h_r) * u_r) / (sqrt(g * h_r) + sqrt(g * h_l))
        c_hat = sqrt(g * 0.5d0 * (h_r + h_l))
        s_roe_l = u_hat - c_hat
        s_roe_r = u_hat + c_hat

        ! Einfeldt speeds
        s_einfeldt_l = min(s_l, s_roe_l)
        s_einfeldt_r = max(s_r, s_roe_r)

        ! Solve Riemann problem
        call riemann_aug_JCP(MAX_ITERATIONS, num_eqn, num_waves, h_l, h_r,  &
                             hu_l, hu_r, 0.d0, 0.d0, b_l, b_r, u_l, u_r,    &
                             0.d0, 0.d0, phi_l, phi_r, s_einfeldt_l,       &
                             s_einfeldt_r, dry_tolerance, g, sw, fw)


        ! Eliminate ghost fluxes for wall
        do m = 1, num_waves
            s(m, i) = sw(m) * wall(m)
            do n = 1, num_eqn
                fwave(n, m, i) = fw(n, m) * wall(m)
            end do
        end do

    end do

    ! Compuate fluctuations
    do i = 2 - num_ghost, mx + num_ghost
        do m = 1, num_waves
            if (s(m, i) < 0.d0) then
                amdq(:, i) = amdq(:, i) + fwave(:, m, i)
            else if (s(m, i) > 0.d0) then
                apdq(:, i) = apdq(:, i) + fwave(:, m, i)
            else
                amdq(:, i) = amdq(:, i) + 0.5d0 * fwave(:, m, i)
                apdq(:, i) = apdq(:, i) + 0.5d0 * fwave(:, m, i)
            end if
        end do
    end do

end subroutine rp1
