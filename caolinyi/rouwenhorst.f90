

subroutine rouwenhorst(N,nu,rho,sigma,grid,PP)

!Purpose:    Finds a Markov chain whose sample paths approximate those of
!            the AR(1) process
!                z(t+1) = (1-rho)*nu + rho * z(t) + eps(t+1)
!            where eps are normal with stddev sigma
!
!Format:     [Z, PI] = rouwenhorst(N,nu,rho,sigma)
!
!Input:      N       scalar, number of nodes for Z
!            nu      scalar, unconditional mean of process
!            rho     scalar
!            sigma   scalar, std. dev. of epsilons
!
!Output:     Z       N*1 vector, nodes for Z
!            PI      N*N matrix, transition probabilities
!Adapted from Martin Floden's Matlab code by David Wiczer

	integer, intent(in)	:: N
	real(8), intent(in)	:: nu,rho,sigma
	real(8), dimension(N,N)	, intent(out)	:: PP
	real(8), dimension(N)	, intent(out)	:: grid
    real(8), dimension(N)                	:: num
	real(8), dimension(N,N)	:: PPZ, PZP, ZPP	
	real(8)	:: sigmaz, p, sd
	integer :: i
    
   
    PP	= 0.0
	PPZ	= 0.0
	PZP	= 0.0
	ZPP	= 0.0	
    
	sigmaz	= sigma / sqrt(1.0-rho*rho)
	p	= (1.0+rho)/2.0
	PP(1,1:2)	= (/ p	 , 1.0-p/)
	PP(2,1:2)	= (/1.0-p,  p 	/)
    
	if (N>=3) then
	do i= 3,N
		PPZ	= 0.0
		PZP	= 0.0
		ZPP	= 0.0
		PPZ(1:i-1,2:i)	= PP(1:i-1,1:i-1)
		PZP(2:i,1:i-1)	= PP(1:i-1,1:i-1)
		ZPP(2:i,2:i)	= PP(1:i-1,1:i-1)
		PP(1:i,1:i) 	= p*PP(1:i,1:i) + (1.0 - p)*PPZ(1:i,1:i) + (1.0 - p)*PZP(1:i,1:i) + p*ZPP(1:i,1:i)
		PP(2:i-1,1:i)	= PP(2:i-1,1:i)/2
	enddo
	endif
	sd	= sqrt(real(N)-1.0)*sigmaz
	num	= (/ (i, i=0,N-1) /)
	grid	= num *(2.0*sd / (real(N) - 1.0)) - sd + nu
end subroutine rouwenhorst
