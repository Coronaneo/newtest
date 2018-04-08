	subroutine nufft1dIII(nj,k,x,nk,iflag,ns,rt,tol,U,V,xsub,ksub,r)
	implicit none

	integer :: ns,rt,iflag,nj,j,r,xsub(nj),i,nk,ksub(nk)
	real  :: tol
	real*8 pi,x(nj),k(nk)
	parameter (pi=3.141592653589793238462643383279502884197d0)
	complex*16 fftconst,U(nk,r),V(nj,r)
	complex*16 M(nk,nj),w1,w2

	fftconst = iflag*dcmplx(0,1)/ns*2*pi

	do i = 1,nk
	   do j = 1,nj
	      M(i,j) = exp(fftconst*(k(i)*x(j)-
     &        floor(k(i)+0.5)*floor(x(j)+0.5)))
	   enddo
	enddo

	!call lowrankfac(M,tol,rt,rt,U,V)

	xsub = mod(floor(x+0.5),ns)+1
	ksub = mod(floor(k+0.5),ns)+1

	r = size(V,2)

	end subroutine

	subroutine nufft1dIIIapp(nj,nk,plan,c,U,V,xsub,ksub,ns,r,S)
	integer  r,i,j,k,nj,ns,num,nk
	integer mm
	integer xsub(nj),ksub(nk)
	complex*16 M(r,nj),N(r,nj),S(nk),c(nj),U(r,nk),V(r,nj),VV(r,nj)
	complex*16 NN(nk,r),Idx(nj,nj),Idk(nk,nk),w1,w2
	double complex in1, out1
	real*16  time_begin,time_end,countrage,countmax
	dimension in1(nj), out1(nj)
	integer*8 :: plan
        w1=dcmplx(1,0)
        w2=dcmplx(0,0)
	M=0
        !VV=V*spread(c,1,r)
	do i = 1,nj
	   do k = 1,r
	      M(k,xsub(i)) = M(k,xsub(i))+V(k,i)*c(i)
	   enddo
	enddo
        
c        do k = 1,r
c           VV(:,k) = V(:,k)*c
c        enddo
c        print *,'V(1,1:5)=',V(1,1:5)
c        do i = 1,nj
c          M(xsub(i),:) = M(xsub(i),:)+VV(i,:)
c        enddo

c        print *,'M(1,1:5)=',M(66,1:5)
	do i = 1,r
	   in1 = M(i,:)
	   call dfftw_execute_dft(plan, in1, out1)
	   N(i,:) = out1
	enddo
        !print *,'N(1,1:5)=',N(66,1:5)
	
c	do i = 1,nk
c	   NN(i,:) = N(ksub(i),:)
c	enddo
        !print *,'NN(1,1:5)',NN(66,1:5)
        !NN=matmul(Idk,N)
        !call zgemm('n','n',nj,r,nj,w1,Idk,nj,N,nj,w2,NN,nj)
        !S = sum(U*NN,2)
        do i = 1,nk
           S(i)=dot_product(U(:,i),N(:,ksub(i)))
        enddo
	!S = sum(U*N(:,ksub),1)

        return
	end subroutine
