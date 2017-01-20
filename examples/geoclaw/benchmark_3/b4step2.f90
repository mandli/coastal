subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)

    use geoclaw_module
    use topo_module
    use dz_module

    implicit none

    ! Input
    integer(kind=4), intent(in) :: mbc, mx, my, meqn, maux
    real(kind=8), intent(in) :: xlower, ylower, dx, dy, t, dt

    real(kind=8), intent(in out) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(in out) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer(kind=4) :: i, j, i1, i2
    character(len=25) :: bfname
    real(kind=8) :: bt, bx, runup1, runup2, x, axi, bxi, eta, x2, xc, xi, z
    real(kind=8) :: zeta

    ! parameters
    real(kind=8), parameter :: tol = 1d-6


    ! Common block
    real(kind=8) :: ybar1, ybar2
    common /comrunup/ ybar1, ybar2

    if (dx <= dxfine) then
        ! Only do this on finest grid
        do i = 1, mx
            if (q(1, i, 1) < 0.001d0) then
                i1 = 1
            end if
            if (q(1, i, 1) < 0.0005d0) then
                i2 = i
            end if
        end do
        runup1 = tanth * (xlower + i1 * dx)
        runup2 = tanth * (xlower + i2 * dx)
        write(32, "(3e20.10)") t, runup1, runup2
    end if

    do j = 1 - mbc, my + mbc
        eta = ylower + (j - 0.5d0) * dy
        do i = 1 - mbc, mx + mbc
            x = xlower + (i - 0.5d0) * dx

            ! Slope with no mass
            z = -x * tanth
            if ((abs(x - xc) < xlength) .and. (eta < width)) then
                ! Region of landslide mass
                x4zeroin = x
                eta4zeroin = eta
                axi = xic - length
                bxi = xic + length
                ! TODO: not sure where this function is coming from
                xi = zeroin(axi, bxi, fxi, tol)
                zeta = zeta_fcn(xi - xic, eta)
                x2 = xi * costh + zeta * sinth
                z = -xi * sinth + zeta * costh
                if (abs(x - x2) > tol) then
                    stop "ERROR"
                end if
            end if
            aux(1, i, j) = max(z, -1.5d0)
        end do
    end do

    ! Check for NaNs
    call check4nans(meqn, mbc, mx, my, q, t, 1)


    forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc,q(1,i,j) < dry_tolerance)
        q(1,i,j) = max(q(1,i,j),0.d0)
        q(2:3,i,j) = 0.d0
    end forall

end subroutine b4step2