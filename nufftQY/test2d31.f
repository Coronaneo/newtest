	program test
	implicit none
        
        integer r,r1
        parameter (r=42,r1=95)
        integer ms
        parameter (ms=64)
        integer nj
        parameter (nj=ms*ms)
        integer nk
        parameter (nk=ms*ms)
        integer i,iflag,xsub(nj,2),ier,num,j,xxsub(nj)
        integer n1,n2,mt,k1,k2,mm,ksub(nk,2),kksub(nk)
        integer*8  time_begin,time_end,countrage,countmax
        real*16 U1(r,nk),V1(r,nj),U2(r,nk),V2(r,nj)
        real*16 re1(nj),re2(nj),time1,time2
        real*16 arr(4),UU1(r1,ms*ms),VV1(r1,nj),UU2(r1,ms*ms),VV2(r1,nj)
        real*8 pi,xj(nj),yj(nj)
        parameter (pi=3.141592653589793238462643383279502884197d0)
        complex*16 U(r,nk),V(r,nj),cj(nj),S(nk),re(nj)
        complex*16 fk(ms,ms),fk1(nj),UU(ms*ms,r1),VV(r1,nj)
        complex*16,allocatable :: NN(:,:,:)
        real*8 x(nj,2),eps,error,k(nk,2),sk(nk),tk(nk)
        double complex in1, out1
        dimension in1(ms,ms), out1(ms,ms)
	integer*16 :: plan
        integer FFTW_FORWARD,FFTW_MEASURE
        parameter (FFTW_FORWARD=-1)
        parameter (FFTW_MEASURE=0)
    
        character*8 date
        character*10 time
        character*5 zone 
        integer*4 values1(8),values2(8)
        
        arr(1)=3600
        arr(2)=60
        arr(3)=1
        arr(4)=0.001


        iflag=-1
        eps=1E-12
        num=1
        
        open(unit = 10,file = 'U2r3.txt')
        read(10,*) U1
        open(unit = 20,file = 'V2r3.txt')
        read(20,*) V1
        open(unit = 10,file = 'U2i3.txt')
        read(10,*) U2
        open(unit = 10,file = 'V2i3.txt')
        read(10,*) V2
        U=dcmplx(U1,U2)
        V=dcmplx(V1,V2)
        V=conjg(V)
        !U=conjg(U)
        !print *,V(2,:)
        !print *,U(1,:)

        open(unit = 10,file = 'UU2r3.txt')
        read(10,*) UU1
        open(unit = 10,file = 'VV2r3.txt')
        read(10,*) VV1
        open(unit = 10,file = 'UU2i3.txt')
        read(10,*) UU2
        open(unit = 10,file = 'VV2i3.txt')
        read(10,*) VV2
        UU=transpose(dcmplx(UU1,UU2))
        VV=dcmplx(VV1,VV2)
        VV=conjg(VV)
        !U=conjg(U)
        !print *,V(2,:)
        !print *,U(1,:)


        call dfftw_plan_dft_2d(plan,ms,ms,in1,out1,FFTW_FORWARD,0)

        print *,'start 2D type 3 testing:'
        print *,'nj  =',nj,'nk  =',nk,'ms  =',ms
        print *,'rank of the first factorization r = ',r
        print *,'rank of the second factorization  r1 = ',r1
        print *,'eps             =',eps
        !re=dcmplx(re1,re2)        

        
        do k1 = -ms/2, (ms-1)/2
           do k2 = -ms/2, (ms-1)/2
              j = (k2+ms/2+1) + (k1+ms/2)*ms
              xj(j) = pi*dcos(-pi*k1/ms)
              yj(j) = pi*dcos(-pi*k2/ms)
              cj(j) = dcmplx(dsin(pi*j/ms),dcos(pi*j/ms))
           enddo
        enddo
        x(:,1)=xj
        x(:,2)=yj
        xsub=mod(floor(x/2/pi*ms+0.5),ms)+1
        do i = 1,nj
           xxsub(i)=xsub(i,2)*ms-ms+xsub(i,1)
        enddo
        !print *,'xxsub=',size(xxsub)


	do i = 1,nk
	   sk(i) = 48*(dcos(i*pi/nk-pi/2))
	   tk(i) = 32*(dsin(i*pi/nk))
	enddo
        k(:,1)=sk
        k(:,2)=tk
        ksub=mod(floor(k+0.5),ms)+1
        do i = 1,nk
           kksub(i)=ksub(i,2)*ms-ms+ksub(i,1)
        enddo


        call date_and_time(date,time,zone,values1)
        do i=1,num
        call nufft2dIIIapp(nj,nk,plan,
     &  cj,U,V,UU,VV,xxsub,kksub,ms,1,r,r1,S)
        enddo
        call date_and_time(date,time,zone,values2)
        !print *,'S(1:5)=',S(1:5)
        time1=sum((values2(5:8)-values1(5:8))*arr)
        print *,' T_our         = ',time1/num

        call date_and_time(date,time,zone,values1)        
        do i=1,num
        call nufft2d3f90(nj,xj,yj,cj,iflag,eps,nk,sk,tk,fk,ier)
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
