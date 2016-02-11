program main
    use mkl95_lapack            !for solving the eigen vector of PP
    implicit none
	integer(4)	                :: N 
    integer(4)	                :: irow, i, j, k
    integer(4)	                :: iter1
	real(8) 	                :: nu,rho,sigma
	real(8), allocatable		:: PP(:,:)
    real(8), allocatable		:: PP0(:,:)
	real(8), allocatable		:: grid(:)
	real(8), allocatable		:: theta(:), xc(:), xp(:), dx(:)
    real(8) 	                :: A, alpha, delta, mu, kappa, beta, b
    real(8) 	                :: xcr
	real(8) 	                :: tol 
    real(8) 	                :: normx
    real(8), allocatable		:: fx(:), fdx(:)
    real(8), allocatable		:: qx(:), qdx(:)
    real(8), allocatable		:: eigr(:), eigi(:), egvector(:,:)
    real(8), allocatable		:: PDF(:)    ! stationary distribution
    real(8), allocatable		:: CDF(:)
    real(8)                     :: y, grid_0, theta_0
	integer(4)	                :: T
    real(8), allocatable		:: Y_(:), grid_(:), theta_(:)
    real(8)                     :: avg_z,va_z, sd_z, avg_theta,va_theta, sd_theta
    
    open(unit=20,file='Output_(0).txt',status='unknown')   
    
    N     = 50
    nu    = 1.1 
    rho   = 0.9895
    sigma = 0.0034
    
    delta = 0.0081
    alpha = 0.72
    A     = 0.158
    mu    = 0.72 
    b     = 0.4   
    !mu    = 0.05
    !b     = 0.95
    kappa = 0.6
    beta  = 0.999
    
    tol   = 1.0E-9
    iter1 = 0
    
    xcr = log(A)/alpha
        
    allocate(PP(N,N),grid(N),theta(N))
    allocate(xc(N),xp(N),dx(N))
    allocate(qx(N),qdx(N))
    allocate(fx(N),fdx(N))
    allocate(eigr(N),eigi(N),egvector(N,N),PDF(N),CDF(N))
    
    call rouwenhorst(N,nu,rho,sigma,grid,PP)
        
    xp = 0.0
    xc = xp
    do  ! 
        iter1  = iter1 + 1
        do irow = 1, N
            do   ! solve the irow-th equation
                ! compute the fxi and f'xi
                call functionF(irow, xc, grid, PP, qx, qdx, N, A, alpha, beta, mu, b, &
                     kappa, delta, xcr, fx, fdx)
                ! update the current value
                xc(irow) = xc(irow) - fx(irow) / fdx(irow)
                ! convergence creteria
                if(abs(fx(irow)/fdx(irow)) < tol .or. abs(fx(irow)) < tol) exit
            end do
        end do                
        
        dx = xc - xp
        call nrm2(dx, N, normx)
        if(normx < tol) exit
        xp = xc !updating

        write(*,"(4x,'Iteration:',i7)") iter1
        write(*,"(4x,'Convergence:',e15.5)") normx
        write(*,*)
    end do
     
    theta = exp(xc) 
     
    !   simulation
    !   initial step
    PP0 = PP
    call geev(PP0,eigr,eigi,egvector)
    PDF = egvector(:,1)/sum(egvector(:,1))
    do i = 1, N
        CDF(i) = sum(PDF(1:i))
    end do
    
    T = 100000
    allocate(grid_(T+1),Y_(T+1),theta_(T+1))      !T+1 draws
    call RANDOM_SEED()
    
    call init_random_seed()    
    do i = 1, T+1
        call RANDOM_NUMBER(Y_(i))
    end do
    y = Y_(1)
    
    if(y <=CDF(1)) irow = 1
    do j = 2, N
        if(y > CDF(j-1) .and. y <=CDF(j)) then
            irow  = j
            grid_0 = grid(irow)
            theta_0 = theta(irow)
            exit
        end if
    end do
    grid_(1) = grid_0
    theta_(1) = theta_0
    
    !  iteration   
    do i = 1, T
        PDF = PP(irow,:)/sum(PP(irow,:))
        do j = 1, N
            CDF(j) = sum(PDF(1:j))
        end do
        
        y = Y_(i+1)
        if(y <=CDF(1)) irow = 1
        do k = 2, N
            if(y > CDF(k-1) .and. y <=CDF(k)) then
                irow  = k
                grid_(i+1) = grid(irow)
                theta_(i+1) = theta(irow)
                exit
            end if
        end do 
!        write(20,"(2i5,3f15.5)") i, irow, y, grid_(i+1), theta_(i+1)
    end do
 
!   calculate sample standard deviation, remember we have T+1 draws
    avg_z = sum(grid_) / (T+1)
    va_z = 0
    do i = 1, T+1
        va_z = va_z + (grid_(i)-avg_z)**2
    end do
    va_z = va_z / T
    sd_z = sqrt (va_z)
    
    avg_theta = sum(theta_) / (T+1)
    va_theta = 0
    do i = 1, T+1
        va_theta = va_theta + (theta_(i)-avg_theta)**2
    end do
    va_theta = va_theta / T
    sd_theta = sqrt (va_theta)
    
    
    !write(20,*) "Transition Matrix"
    !do i = 1, N
        !write(*,"(1x,i4,<N>E15.5)") i, PP(i,:)
        !write(20,"(1x,i4,<N>E15.5)") i, PP(i,:)
    !end do
    write(*,*)
    write(20,*)
    write(*,*) "       grid           theta"
    write(20,*) "       grid           theta"
    do i = 1, N
        write(*,"(1x,i4,2F15.5)") i, grid(i), theta(i)
        write(20,"(1x,i4,2F15.5)") i, grid(i), theta(i)
    end do
    
    write(*,*)
    write(*,"(1x,a,F15.5)") "avg_z     =", avg_z
    write(*,"(1x,a,F15.5)") "sd_z      =", sd_z
    write(*,"(1x,a,F15.5)") "avg_theta =", avg_theta
    write(*,"(1x,a,F15.5)") "sd_theta  =", sd_theta
    write(20,*)
    write(20,"(1x,a,F15.5)") "avg_z     =", avg_z
    write(20,"(1x,a,F15.5)") "sd_z      =", sd_z
    write(20,"(1x,a,F15.5)") "avg_theta =", avg_theta
    write(20,"(1x,a,F15.5)") "sd_theta  =", sd_theta
    
    close(20)
    
    
    
end program main
