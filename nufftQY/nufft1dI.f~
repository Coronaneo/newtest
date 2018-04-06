	subroutine nufft1dI(nj,x,iflag,ns,rt,tol,U,V,xsub,r)
	implicit none

	integer :: ns,rt,iflag,nj,k,j,r,xsub(nj)
	real  :: tol
	real*8 pi,x(nj)
	parameter (pi=3.141592653589793238462643383279502884197d0)
	complex*16 fftconst,U(ns,r),V(nj,r)
	complex*16 M(ns,nj)


        

	fftconst = iflag*dcmplx(0.0,1.0)/ns*2*pi;
	do k = 1,ns
	   do j = 1,nj
	      M(k,j) = exp(fftconst*(k-1)*(x(j)-floor(x(j)+0.5)))
	   enddo
	enddo

	!call lowrankfac(M,tol,rt,rt,U,V)

	r = size(V,2)

	xsub = mod(floor(x+0.5),ns)+1

	return

	end subroutine



        subroutine nufft1dIapp(nj,c,U,V,xsub,ns,iflag,r,S,plan)
	implicit none
	integer  r,i,j,k,nj,ns,iflag,num
        integer mm
	integer xsub(nj)
	complex*16 M(nj,r),N(ns,r),S(ns),c(nj),U(ns,r),V(nj,r)
        complex*16,allocatable :: NN(:,:)
	double complex in1, out1
        real*16  time_begin,time_end,countrage,countmax
	dimension in1(nj), out1(ns)
	integer*8 :: plan
        integer FFTW_FORWARD,FFTW_MEASURE
        parameter (FFTW_FORWARD=-1)
        parameter (FFTW_MEASURE=0)


	M=0
        num=100

        !M=conjg(V)*reshape(c,1,r)

        !do j = 1,num
       ! M=0
	do i = 1,nj
	   do k = 1,r
	      M(xsub(i),k) = M(xsub(i),k)+V(i,k)*c(i)
	   enddo
	enddo
       ! enddo


        !do i = 1,nj
         !  X(xsub(i),i) = 1
        !enddo
        !print *,sum(X)
        !M=matmul(X,M)
	!if (iflag < 0) 



	
        !do j = 1,num
	  do i = 1,r
	     in1 = conjg(M(:,i))

             call dfftw_execute_dft(plan, in1, out1)
	     N(:,i) = out1
	  enddo
        !enddo
        
        
        

 
        !call dfftw_destroy_plan(plan)
	!end
        mm=floor(ns/2.0+0.6)
        allocate(NN(mm,r))
        NN=N(1:mm,:)
        N(1:ns-mm,:)=N(mm+1:ns,:)
        N(ns-mm+1:ns,:)=NN

        N=U*N
	S = sum(N,2)

        
	end
