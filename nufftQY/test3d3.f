	program test
	implicit none
        
        integer r,r1
        parameter (r=74,r1=560)
        integer ms
        parameter (ms=48)
        integer nj
        parameter (nj=ms*ms*ms)
        integer nk
        parameter (nk=ms*ms*ms)
        integer i,iflag,ier,num,j
        integer n1,n2,mt,k1,k2,mm,k3
        integer,allocatable :: xsub(:,:),ksub(:,:),xxsub(:),kksub(:)
        integer*8  time_begin,time_end,countrage,countmax
        real*8,allocatable :: U1(:,:),V1(:,:),U2(:,:),V2(:,:)
        real*8 arr(4),time1,time2,pi,eps,error
        real*8,allocatable :: UU1(:,:),VV1(:,:),UU2(:,:)
        real*8,allocatable :: xj(:),yj(:),zj(:),VV2(:,:)
        parameter (pi=3.141592653589793238462643383279502884197d0)
        complex*16,allocatable :: U(:,:),V(:,:),cj(:),S(:)
        complex*16,allocatable :: fk(:,:,:),fk1(:),UU(:,:),VV(:,:)
        complex*16,allocatable :: NN(:,:,:)
        real*8,allocatable :: x(:,:),k(:,:),sk(:),tk(:),uk(:)
        double complex in1, out1
        dimension in1(ms,ms,ms), out1(ms,ms,ms)
	integer*16 :: plan
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
        allocate(xj(nj))
        allocate(yj(nj))
        allocate(zj(nj))
        allocate(S(nk))
        allocate(fk1(nj))
        allocate(U(r,nk))
        allocate(V(r,nj))
        allocate(cj(nj))
        allocate(fk(ms,ms,ms))
        allocate(UU(ms*ms*ms,r1))
        allocate(xsub(nj,3))
        allocate(xxsub(nj))
        allocate(x(nj,3))
        allocate(ksub(nk,3))
        allocate(kksub(nk))
        allocate(UU1(r1,ms*ms*ms))
        allocate(VV1(r1,nj))
        allocate(UU2(r1,ms*ms*ms))
        allocate(VV2(r1,nj))
        allocate(VV(r1,nj))
        allocate(k(nk,3))
        allocate(sk(nk))
        allocate(tk(nk))
        allocate(uk(nk))
        arr(1)=3600
        arr(2)=60
        arr(3)=1
        arr(4)=0.001


        iflag=-1
        eps=1E-4
        num=1
        
        open(unit = 10,file = 'U3r3.txt')
        read(10,*) U1
        open(unit = 20,file = 'V3r3.txt')
        read(20,*) V1
        open(unit = 10,file = 'U3i3.txt')
        read(10,*) U2
        open(unit = 10,file = 'V3i3.txt')
        read(10,*) V2
        U=dcmplx(U1,U2)
        V=dcmplx(V1,V2)
        V=conjg(V)
        !U=conjg(U)
        !print *,V(2,:)
        !print *,U(1,:)

        open(unit = 10,file = 'UU3r3.txt')
        read(10,*) UU1
        open(unit = 10,file = 'VV3r3.txt')
        read(10,*) VV1
        open(unit = 10,file = 'UU3i3.txt')
        read(10,*) UU2
        open(unit = 10,file = 'VV3i3.txt')
        read(10,*) VV2
        UU=transpose(dcmplx(UU1,UU2))
        VV=dcmplx(VV1,VV2)
        VV=conjg(VV)
        !U=conjg(U)
        !print *,V(2,:)
        !print *,U(1,:)


        call dfftw_plan_dft_3d(plan,ms,ms,ms,in1,out1,FFTW_FORWARD,0)

        print *,'start 3D type 3 testing:'
        print *,'nj  =',nj,'nk  =',nk,'ms  =',ms
        print *,'rank of the first factorization r = ',r
        print *,'rank of the second factorization  r1 = ',r1
        print *,'eps             =',eps
        !re=dcmplx(re1,re2)        

        
        do k3 = -ms/2,(ms-1)/2
	 do k1 = -ms/2, (ms-1)/2
	   do k2 = -ms/2, (ms-1)/2
	      j =  (k1+ms/2+1) + (k2+ms/2)*ms + (k3+ms/2)*ms*ms
	      xj(j) = pi*dcos(-pi*k1/ms)
	      yj(j) = pi*dcos(-pi*k2/ms)
	      zj(j) = pi*dcos(-pi*k3/ms)
	      cj(j) = dcmplx(dsin(pi*j/ms),dcos(pi*j/ms))
	   enddo
	 enddo
        enddo
        x(:,1)=xj
        x(:,2)=yj
        x(:,3)=zj
        xsub=mod(floor(x/2/pi*ms+0.5),ms)+1
        do i = 1,nj
           xxsub(i)=xsub(i,3)*ms*ms-ms*ms+xsub(i,2)*ms-ms+xsub(i,1)
        enddo
        !print *,'xxsub=',size(xxsub)


	do i = 1,nk
	   sk(i) = 10*(cos(i*pi/nk-pi/2));
           tk(i) = 8*(sin(i*pi/nk));
           uk(i) = 10*(cos(i*pi/nk-pi/2));
	enddo
        k(:,1)=sk
        k(:,2)=tk
        k(:,3)=uk
        ksub=mod(floor(k+0.5),ms)+1
        do i = 1,nk
           kksub(i)=ksub(i,3)*ms*ms-ms*ms+ksub(i,2)*ms-ms+ksub(i,1)
        enddo


        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft3dIIIapp(nj,nk,plan,
     &  cj,U,V,UU,VV,xxsub,kksub,ms,1,r,r1,S)
        enddo
        call date_and_time(date,time,zone,values2)
        !print *,'S(1:5)=',S(1:5)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our         = ',time1/num

        call date_and_time(date,time,zone,values1)        
        do i=1,num
        call nufft3d3f90(nj,xj,yj,zj,cj,iflag,eps,nk,sk,tk,uk,fk,ier)
        enddo
        call date_and_time(date,time,zone,values2)
        time2=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_nyu         = ',time2/num
        print *,' T_our/T_nyu   = ',time1/time2
        fk1=reshape(fk,(/nj/))
        error=sqrt(real(sum((S-fk1)*conjg(S-fk1))/
     &  sum(fk1*conjg(fk1))))
        print *,' relative error= ',error
        
        call dfftw_destroy_plan(plan)
        

	end program
