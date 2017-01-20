subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)

    use geoclaw_module, only: sea_level
    use dz_module

    implicit none

    ! Input
    integer(kind=4), intent(in) :: meqn, mbc, mx, my, maux
    real(kind=8), intent(in) :: xlower, ylower, dx, dy
       
    real(kind=8), intent(in out) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(in out) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

    ! Locals
    integer :: i, j

    dxfine = min(dxfine, dx)
    write(35, *) 'dxfine: ',dxfine

    ! Topo is set in b4step2:
    call b4step2(mbc,mx,my,meqn,q,xlower,ylower,dx,dy,0.d0,dt,maux,aux)

    do i = 1 - mbc, mx + mbc
        do j = 1 - mbc, my + mbc

            q(1, i, j) = max(0.d0, sea_level - aux(1, i, j))
            q(2, i, j) = 0.d0
            q(3, i, j) = 0.d0
        end do
    end do

end subroutine qinit

