! ============================================
subroutine b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,t,dt,maux,aux)
! ============================================
!
! # called before each call to step
! # use to set time-dependent aux arrays or perform other tasks.
!
! This particular routine sets negative values of q(1,i,j) to zero,
! as well as the corresponding q(m,i,j) for m=1,meqn.
! This is for problems where q(1,i,j) is a depth.
! This should occur only because of rounding error.
!
! Also calls movetopo if topography might be moving.

    use geoclaw_module, only: dry_tolerance
    use geoclaw_module, only: g => grav
    use topo_module, only: num_dtopo,topotime
    use topo_module, only: aux_finalized
    use topo_module, only: xlowdtopo,xhidtopo,ylowdtopo,yhidtopo

    use amr_module, only: xlowdomain => xlower
    use amr_module, only: ylowdomain => ylower
    use amr_module, only: xhidomain => xupper
    use amr_module, only: yhidomain => yupper
    use amr_module, only: xperdom,yperdom,spheredom,NEEDS_TO_BE_SET

    use storm_module, only: set_storm_fields

    implicit none

    ! Subroutine arguments
    integer, intent(in) :: meqn
    integer, intent(inout) :: mbc,mx,my,maux
    real(kind=8), intent(inout) :: xlower, ylower, dx, dy, t, dt
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Local storage
    integer :: index,i,j,k,dummy,jbar1, jbar2, i1, i2, istatus
    real(kind=8) :: h, u, v, x, y, bt, bx, runup1, runup2, xb1
    character(len=80) :: bfname

    ! Common block storage
    real(kind=8) :: ybar1, ybar2
    common /comrunup/ ybar1, ybar2

    if (dx < 0.02d0) then
        jbar1 = int((ybar1 - ylower) / dy + 1)
        jbar2 = int((ybar2 - ylower) / dy + 1)
        do i = 1, mx
            if (q(1, i, jbar1) < 0.001d0) then
                i1 = i
            endif
            if (q(1, i, jbar2) < 0.001d0) then
                i2 = i
            end if
        end do
        runup1 = 0.5d0 * (xlower + i1 * dx)
        runup2 = 0.5d0 * (xlower + i2 * dx)
        write(32, "(3e20.10)") t, runup1, runup2
    end if

    bfname = "../block2.tx"

    if (t < 3.d0) then
        open(unit=50, file=bfname, status="old", form='formatted')
        do k = 1, 601
            read(50, fmt=*, iostat=istatus) bt, bx
            if (bt >= t) then
                xb1 = 0.01d0 * bx + 0.2d0
                close(50)
                exit
            end if
        end do
    else
        xb1 = 4.045933071 + 0.2d0
    end if

    do i = 1, mx
        x = xlower + (i - 0.5d0) * dx
        do j = 1, my
            y = ylower + (j - 0.5d0) * dy
            if (abs(y) <= 0.305d0) then
                if ((x >= xb1) .and. (x <= xb1 + 0.91d0)) then
                    aux(1, i, j) = -0.5d0 * xb1
                else
                    aux(1, i, j) = -0.5d0 * x
                end if
            else
                aux(1, i, j) = -0.5d0 * x
            end if
        end do
    end do

    ! Check for NaNs in the solution
    call check4nans(meqn,mbc,mx,my,q,t,1)

    ! check for h < 0 and reset to zero
    ! check for h < drytolerance
    ! set hu = hv = 0 in all these cells
    forall(i=1-mbc:mx+mbc, j=1-mbc:my+mbc,q(1,i,j) < dry_tolerance)
        q(1,i,j) = max(q(1,i,j),0.d0)
        q(2:3,i,j) = 0.d0
    end forall


    if (aux_finalized < 2) then
        ! topo arrays might have been updated by dtopo more recently than
        ! aux arrays were set unless at least 1 step taken on all levels
        aux(1,:,:) = NEEDS_TO_BE_SET ! new system checks this val before setting
        call setaux(mbc,mx,my,xlower,ylower,dx,dy,maux,aux)
    endif

    ! Set wind and pressure aux variables for this grid
!     call set_storm_fields(maux,mbc,mx,my,xlower,ylower,dx,dy,t,aux)

end subroutine b4step2
