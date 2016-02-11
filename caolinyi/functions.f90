subroutine functionF(irow, xc, grid, PP, qx, qdx, N, A, alpha, beta, mu, b, &
                 kappa, delta, xcr, fx, fdx)
    implicit none
    integer(4)      :: i, j
    integer(4)      :: irow, N    
    real(8) 	    :: xc(N)
    real(8) 	    :: A, alpha, beta, mu, b, kappa, delta, xcr
    real(8)         :: grid(N), PP(N,N), q_x(N), q_dx(N)
    real(8)         :: fxi, fdxi
    real(8)         :: qx(N), qdx(N)
    real(8)         :: fx(N), fdx(N)
    real(8)         :: tmp0
    
    tmp0 = 0

    do i = 1, N
        call functionQ(xc(i), A, alpha, xcr, qx(i), qdx(i))
    end do
        
    do j = 1, N
        tmp0 = tmp0 + PP(irow, j) *((1-mu)*(grid(j)-b) - kappa*mu*grid(j)  &             ! I add a z_j term here
               * exp(xc(j)) + (1-delta)*kappa / qx(j))    
    end do
 
    fxi  = beta * qx(irow) * tmp0 - kappa
    fdxi = beta * qdx(irow) * tmp0 - beta * qx(irow) * PP(irow, irow)  &
         *(kappa*mu*grid(irow)*exp(xc(irow)) - (1-delta)*kappa/qx(irow)**2*qdx(irow))    ! I add a z_i term here
    
    fx(irow) = fxi
    fdx(irow) = fdxi

end subroutine functionF
    
    
    
    
subroutine functionQ(x, A, alpha, xcr, qx, qdx)
    implicit none    
    real(8) 	                :: x, A, alpha
    real(8) 	                :: xcr
    real(8) 	                :: qx, qdx
    
    if(x <= xcr) then
        qx  = 1
        qdx = 0
    else
        qx  = A * exp(-alpha * x)
        qdx = - alpha * A * exp(-alpha * x)
    end if
    
end subroutine functionQ
    
    
subroutine nrm2(x, N, x1)
    implicit none
    integer(4)   :: i
    integer(4)   :: N
    real(8)      :: x(N)
    real(8)      :: x1      ! 2-nd norm
    x1 = 0.0
    do i = 1, N
        x1 = x1 + x(i)**2
    end do    
    x1 = sqrt(x1)    
    end subroutine nrm2
    
    subroutine init_random_seed()
    integer :: i, n
    integer, dimension(:), allocatable :: seed
    
    call random_seed(size = n)
    allocate(seed(n))
    
! the code will get the same random number 
    seed = 1990010219890801 + 37 * (/ (i, i = 1, n) /)
    call random_seed(put = seed)
    
    deallocate(seed)
 end subroutine init_random_seed    