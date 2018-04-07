	program test
	implicit none
        integer nj,ns,r,mm,kflag
        parameter (r=12,ns=131072,kflag=-1)
        parameter (nj=ns)
        integer i,iflag,ier,num
        integer,allocatable :: xsub(:)
        real*16 begin1,end1
        integer*8  time_begin,time_end,countrage,countmax
        real*16,allocatable :: U1(:,:),V1(:,:),U2(:,:),V2(:,:),x(:)
        real*16 pi,time1,time2
        real*16 arr(4)
        parameter (pi=3.141592653589793238462643383279502884197d0)
        complex*16,allocatable :: U(:,:),V(:,:),c(:),S(:)
        complex*16,allocatable :: fk(:)
        complex*16,allocatable :: UU(:,:)
        real*8 x1(nj),eps,error
        double complex in1, out1
        dimension in1(ns), out1(ns)
	integer*8 :: plan
        integer FFTW_FORWARD,FFTW_MEASURE
        parameter (FFTW_FORWARD=-1)
        parameter (FFTW_MEASURE=0)
    
        character*8 date
        character*10 time
        character*5 zone 
        integer*4 values1(8),values2(8)

        allocate(U1(r,ns))
        allocate(U2(r,ns)) 
        allocate(V1(r,nj))
        allocate(V2(r,nj))
        allocate(x1(nj))
        allocate(S(ns))      
        allocate(U(r,ns)) 
        allocate(V(r,nj))
        allocate(c(nj))
        allocate(fk(ns))
        allocate(xsub(nj))
        allocate(x(nj))
        
        arr(1)=3600
        arr(2)=60
        arr(3)=1
        arr(4)=0.001


        iflag=-1
        eps=1E-4
        num=100
        open(unit = 10,file = 'Ur1.txt')
        read(10,*) U1
        open(unit = 20,file = 'Vr1.txt')
        read(20,*) V1
        open(unit = 10,file = 'Ui1.txt')
        read(10,*) U2
        open(unit = 10,file = 'Vi1.txt')
        read(10,*) V2
        print *,'ok'
        call dfftw_plan_dft_1d(plan,nj,in1,out1,FFTW_FORWARD,0)
        print *,'start 1D type 1 testing:'
        print *,'nj  =',nj,'ns  =',ns
        print *,'rank r = ',r
        print *,'eps             =',eps
      
        U=dcmplx(U1,U2)
        if (kflag .lt. 0) then
           mm=floor(ns/2.0+0.6)
           allocate(UU(r,mm))
           UU=U(:,1:mm)
           U(:,1:mm)=U(:,mm+1:ns)
           U(:,mm+1:ns)=UU
        endif
        V=dcmplx(V1,V2)
        V=conjg(V)
        !print *,V(2,:)
        !print *,U(1,:)
        do i = 1,nj
           x(i) = i*pi/8
        enddo
        xsub=mod(floor(x+0.5),ns)+1
        do i = 1,nj
           c(i) = exp(-dcmplx(0,1)*i/ns)
        enddo
        !print *,'c(1:5)=',c(1:5)
        do i = 1,nj
           x1(i) = i*pi*2*pi/(8*nj)
        enddo

        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft1dIapp(nj,plan,c,U,V,xsub,ns,kflag,r,S)
        enddo
        call date_and_time(date,time,zone,values2)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our         = ',time1/num

        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft1d1f90(nj,x1,c,iflag,eps,ns,fk,ier)
        enddo
        call date_and_time(date,time,zone,values2)
        time2=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_nyu         = ',time2/num
        print *,' T_our/T_nyu   = ',time1/time2
        error=sqrt(real(sum((S-fk*nj)*conjg(S-fk*nj))/
     &  sum(S*conjg(S))))
        print *,' relative error= ',error
        call dfftw_destroy_plan(plan)
        

	end program
