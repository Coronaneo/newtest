	program test
	implicit none
        integer r
        integer*8 nj,ms,nk
        parameter (r=21,ms=131072)
        parameter (nj=ms,nk=ms)
        integer*8 i,iflag,ier,num
        integer*8,allocatable :: xsub(:),ksub(:)
        real*16 begin1,end1
        integer*8  time_begin,time_end,countrage,countmax
        real*16,allocatable :: U1(:,:),V1(:,:),U2(:,:),V2(:,:),x(:)
        real*16 time1,time2
        real*16 arr(4)
        real*8 pi,eps,error
        parameter (pi=3.141592653589793238462643383279502884197d0)
        complex*16,allocatable :: c(:),S(:),U(:,:),V(:,:),fk(:)
        real*8,allocatable :: x1(:),k(:)
        double complex in1, out1
        dimension in1(ms), out1(ms)
	integer*8 :: plan
        integer FFTW_FORWARD,FFTW_MEASURE
        parameter (FFTW_FORWARD=-1)
        parameter (FFTW_MEASURE=0)
    
        character*8 date
        character*10 time
        character*5 zone 
        integer*4 values1(8),values2(8)

        allocate(U1(r,nk))
        allocate(U2(r,nk)) 
        allocate(V1(r,nj))
        allocate(V2(r,nj))
        allocate(x1(nj))
        allocate(S(nk))      
        allocate(U(r,nk)) 
        allocate(V(r,nj))
        allocate(c(nj))
        allocate(fk(nk))
        allocate(xsub(nj))
        allocate(x(nj))
        allocate(k(nk))
        allocate(ksub(nk))
        
        arr(1)=3600
        arr(2)=60
        arr(3)=1
        arr(4)=0.001

        iflag=-1
        eps=1d-4
        num=100
        open(unit = 10,file = 'Ur3.txt')
        read(10,*) U1
        open(unit = 20,file = 'Vr3.txt')
        read(20,*) V1
        open(unit = 10,file = 'Ui3.txt')
        read(10,*) U2
        open(unit = 10,file = 'Vi3.txt')
        read(10,*) V2

        call dfftw_plan_dft_1d(plan,nj,in1,out1,FFTW_FORWARD,0)
        print *,'start 1D type 3 testing:'
        print *,'nj  =',nj,'nk  =',nk,'ms  =',ms
        print *,'rank  r = ',r
        print *,'eps            =',eps
      
        U=dcmplx(U1,U2)
        U=conjg(U)
        V=dcmplx(V1,V2)

        !print *,'V(1,1:5)=',V(1,1:5)
        !print *,'U(1,1:5)=',U(1,1:5)
        do i = 1,nj
           x(i) = dsin(-pi*i/nj)/2*nj
        enddo
        xsub=mod(floor(x+0.5)+ms,ms)+1
        !print *,'xsub(1:5)=',xsub(65:74)
        do i = 1,nk
           k(i) = 48*dcos(i*pi/ms)
        enddo
        ksub=mod(floor(k+0.5)+ms,ms)+1
        !print *,'ksub(1:5)=',ksub(65:74)
        do i = 1,nj
           c(i) = dcmplx( dcos(pi*i/nj), -dsin(pi*i/nj))
        enddo
        !print *,'c =',c(-64:-59)
        do i = 1,nj
           x1(i) = pi * dsin(-pi*i/nj)
        enddo
        !print *,'ok'


        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft1dIIIapp(nj,nk,plan,c,U,V,xsub,ksub,ms,r,S)
        enddo
        call date_and_time(date,time,zone,values2)
        !print *,'S(1:5)=',S(65:75)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our         = ',time1/num

        call date_and_time(date,time,zone,values1)
        do i=1,num
           call nufft1d3f90(nj,x1,c,iflag,eps,nk,k,fk,ier)
        enddo
        call date_and_time(date,time,zone,values2)
        time2=sum((values2(5:8)-values1(5:8))*arr)
     
        !print *,'ier=',ier
        !print *,'fk(1:5)=',fk(1:5)
        print *,' T_nyu         = ',time2/num
        print *,' T_our/T_nyu   = ',time1/time2
c        do i = 1,128
c           print *,i,re(i)-fk(i)
c        enddo
c        error=sqrt(real(sum((S-re)*conjg(S-re))/
c     &  sum(re*conjg(re))))
c        print *,' relative error1= ',error
c        error=sqrt(real(sum((fk-re)*conjg(fk-re))/
c     &  sum(re*conjg(re))))
c        print *,' relative error2= ',error
        error=sqrt(real(sum((S-fk)*conjg(S-fk))/
     &  sum(S*conjg(S))))
        print *,' relative error= ',error
        call dfftw_destroy_plan(plan)
        

	end program
