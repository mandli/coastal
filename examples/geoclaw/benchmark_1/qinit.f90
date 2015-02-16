
subroutine qinit(meqn,mbc,mx,my,xlower,ylower,dx,dy,q,maux,aux)
    
    use qinit_module, only: qinit_type,add_perturbation
    use geoclaw_module, only: sea_level
    
    implicit none
    
    ! Subroutine arguments
    integer, intent(in) :: meqn,mbc,mx,my,maux
    real(kind=8), intent(in) :: xlower,ylower,dx,dy
    real(kind=8), intent(inout) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)
    
    ! Locals
    integer :: i, j
    real(kind=8) :: x, y, eta, h0

    ! Parameters
    real(kind=8) :: hmax, gamma, xs

    hmax = 0.0185d0
    gamma = sqrt(0.75d0 * hmax)
    xs = 38.d0
    
    ! Set soliton IC
    do i = 1 - mbc, mx + mbc
        x = xlower + (i - 0.5d0) * dx
        do j = 1 - mbc, my + mbc
            y = ylower + (j - 0.5d0) * dy
            eta = 4.d0 * hmax * (1.d0 / (exp(gamma * (x - xs)) + exp(-gamma * (x - xs))))**2
            q(1, i, j) = max(0.d0, eta - aux(1, i, j))
            h0 = max(0.d0, -aux(1, i, j))
            q(2, i, j) = q(1, i, j) * 2.d0 * (sqrt(h0) - sqrt(q(1, i, j)))
            q(3, i, j) = 0.d0
        end do
    end do
    
end subroutine qinit
