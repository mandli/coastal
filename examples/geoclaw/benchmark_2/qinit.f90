
subroutine qinit(meqn, mbc, mx, my, xlower, ylower, dx, dy, q, maux, aux)

    implicit none

    ! Input
    integer(kind=4), intent(in) :: meqn, mbc, mx, my, maux
    real(kind=8), intent(in) :: xlower, ylower, dx, dy

    ! In/Out
    real(kind=8), intent(in out) :: q(meqn,1-mbc:mx+mbc,1-mbc:my+mbc)
    real(kind=8), intent(in out) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)


    ! Locals
    real(kind=8) :: x, y, eta, h_0
    integer :: i, j

    real(kind=8), parameter :: h_max = 0.3d0
    real(kind=8), parameter :: gamma = sqrt(0.75d0 * h_max)
    real(kind=8), parameter :: x_s = 25.d0
    
    do i = 1 - mbc, mx + mbc
        x = xlower + (i - 0.5d0) * dx

        do j = 1 - mbc, my + mbc
            y = ylower + (j - 0.5d0) * dy

            eta = 4.0d0 * h_max * (1.d0 / (exp(gamma * (x - x_s))    &
                                         + exp(-gamma * (x - x_s))))**2
    
            q(1, i, j) = max(0.d0, eta - aux(1, i, j))
            h_0 = max(0.d0, -aux(1, i, j))
            q(2, i, j) = q(1, i, j) * 2.d0 * (sqrt(h_0) - sqrt(q(1, i, j)))
            q(3, i, j) = 0.d0

        end do
            
    end do

end subroutine qinit

