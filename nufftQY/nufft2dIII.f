c	subroutine nufft2dIII(nj,k,x,nk,iflag,ns,rt,tol,U,V,kksub,r)
c	implicit none
c
c	integer :: ns,rt,iflag,nj,j,r,xsub(nj),i,nk
c	real  :: tol
c	real*8 pi,x(nj),k(nk)
c	parameter (pi=3.141592653589793238462643383279502884197d0)
c	complex*16 fftconst,U(nk,r),V(nj,r)
c	complex*16 M(nk,nj)
c
c	fftconst = iflag*dcmplx(0,1)/ns*2*pi
c
c	do i = 1,nk
c	   do j = 1,nj
c	      M(i,j) = exp(fftconst*(x(j,1)*(k(i,1)-floor(k(i,1)+0.5))+
c	&     x(j,2)*(k(i,2)-floor(k(i,2)+0.5))))
c	   enddo
c	enddo
c
c	call lowrankfac(M,tol,rt,rt,U,V)
c
c	ksub = mod(floor(k+0.5),ns)+1
c	do i = 1,nk
c	   kksub(i) = ksub(i,2)*ns-ns+ksub(i,1)
c	enddo
c
c	r = size(V,2)
c
c	end subroutine

	subroutine nufft2dIIIapp(nj,nk,plan,c,U,V,U1,
     &  V1,xxsub,kksub,ns,kflag,r,r1,S)
	integer  r,i,j,k,nj,ns,kflag,num,nk
	integer mm,r1
	integer xxsub(nj),kksub(nk)
	complex*16 M(r,nj),S(nk),c(nj),U(r,nk),V(r,nj)
	complex*16 U1(ns*ns,r1),V1(r1,nj),M1(r,ns*ns)


	M=0

	do i = 1,nj
	   do k = 1,r
	      M(k,i) = V(k,i)*c(i)
	   enddo
	enddo

	do i = 1,r
	call nufft2dIapp(nj,plan,M(i,:),U1,V1,xxsub,ns,kflag,r1,M1(i,:))
	enddo

	S=sum(U*M1(:,kksub),1)
	end subroutine

