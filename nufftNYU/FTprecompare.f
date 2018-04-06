cc Copyright (C) 2004-2009: Leslie Greengard and June-Yub Lee 
cc Contact: greengard@cims.nyu.edu
cc 
cc This software is being released under a FreeBSD license
cc (see license.txt in this directory). 
c
      program testfft
      implicit none
c
c --- local variables
c
      integer i,ier,iflag,j,k1,mx,ms,nj,ip
      parameter (mx=100000)
      parameter (nj=65536)
      parameter (ms=65536)
      integer*8  time_begin,time_end,time_begin1,time_end1
      integer*8 countrage,countmax
      integer  num
      parameter (num=500)
      integer r
      parameter (r=12)
      real*8 xj(mx), sk(mx)
      real*8 err,eps,pi,ratio
      parameter (pi=3.141592653589793238462643383279502884197d0)
      complex*16 cj(mx),cj0(mx),cj1(mx)
      complex*16 fk0(mx),fk1(mx)
      real*16,allocatable :: M1(:,:)
      real*16,allocatable :: M2(:,:)
      real*16 ::  a(nj),b(nj)
      complex*16,allocatable :: M(:,:)
      double complex in, out
      dimension in(nj), out(nj)
      integer*16 plan
      integer FFTW_FORWARD,FFTW_MEASURE
      parameter (FFTW_FORWARD=-1)
      parameter (FFTW_MEASURE=0)
c   
c     --------------------------------------------------
c     create some test data
c     --------------------------------------------------
     
      allocate(M1(r,nj))
      allocate(M2(r,nj))
      allocate(M(r,nj))
      do k1 = -nj/2, (nj-1)/2
         j = k1+nj/2+1
         xj(j) = j/nj*2*pi*pi/4
         cj(j) = dcmplx( dcos(pi*j/nj/pi), dsin(pi*j/nj/pi))
      enddo
      call dfftw_plan_dft_1d(plan,nj,in,out,FFTW_FORWARD,FFTW_MEASURE)

c   
c     --------------------------------------------------
c     start tests
c     --------------------------------------------------
c
      iflag = -1
      print*,' Start 1D comparing: ', ' nj =',nj, ' ms =',ms
      do i = 3,3
         if (i.eq.1) eps=1d-6
         if (i.eq.2) eps=1d-8
         if (i.eq.3) eps=1d-12
         if (i.eq.4) eps=1d-16
c extended/quad precision tests
         if (i.eq.5) eps=1d-20
         if (i.eq.6) eps=1d-24
         if (i.eq.7) eps=1d-28
         if (i.eq.8) eps=1d-32
	 print*,' '
  	 print*,' Requested precision eps =',eps
	 print*,' '
c
c     -----------------------
c     call 1D Type1 method
c     -----------------------
c
         print *,' rank r = ',r
         print *,' '
c         call dirft1d1(nj,xj,cj,iflag, ms,fk0)
         print *,'type 1 :'
	 open(unit=10,file='M1rh.txt')
	 read(10,*) M1
	 open(unit=20,file='M1ih.txt')
	 read(20,*) M2
	 a=M1(1,:)
	 b=M2(1,:)
	 in=dcmplx(a,-b)
	 call system_clock(time_begin,countrage,countmax)
	 do ip=1,num
	  call dfftw_execute_dft(plan, in, out)
	 enddo

	 call system_clock(time_end,countrage,countmax)
         
	 print *,' T_fft         = ',(time_end-time_begin)/num
	 print *,' r*T_fft       = ',r*(time_end-time_begin)/num

            
         call system_clock(time_begin1,countrage,countmax)
         do ip=1,num
            call nufft1d1f90(nj,xj,cj,iflag,eps, ms,fk1,ier)
         enddo
         call system_clock(time_end1,countrage,countmax)
         print *,' T_g           = ',(time_end1-time_begin1)/num
         ratio=1.000*r*(time_end-time_begin)/(time_end1-time_begin1)
         print *,' T_fft / T_g   = ',ratio/r
         print *,' r*T_fft / T_g = ',ratio
c         call errcomp(fk0,fk1,ms,err)
c         print *,' ier           = ',ier
c         print *,' type 1 error  = ',err
c         print *,' '
c
c     -----------------------
c     call 1D Type2 method
c     -----------------------
c
c         call dirft1d2(nj,xj,cj0,iflag, ms,fk0,ier)
         print *,'type 2 :'
	 open(unit=10,file='M2rh.txt')
	 read(10,*) M1
	 open(unit=20,file='M2ih.txt')
	 read(20,*) M2
	 a=M1(1,:)
	 b=M2(1,:)
	 in=dcmplx(a,-b)
	 call system_clock(time_begin,countrage,countmax)
	 do ip=1,num
	  call dfftw_execute_dft(plan, in, out)
	 enddo
	 call system_clock(time_end,countrage,countmax)
	 print *,' T_fft         = ',(time_end-time_begin)/num
	 print *,' r*T_fft       = ',r*(time_end-time_begin)/num

         call system_clock(time_begin1,countrage,countmax)
         do ip=1,num
           call nufft1d2f90(nj,xj,cj1,iflag, eps, ms,fk0,ier)
         enddo
         call system_clock(time_end1,countrage,countmax)
         print *,' T_g           = ',(time_end1-time_begin1)/num
         ratio=1.000*r*(time_end-time_begin)/(time_end1-time_begin1)
         print *,' T_fft / T_g   = ',ratio/r
         print *,' r*T_fft / T_g = ',ratio
c         call errcomp(cj0,cj1,nj,err)
c         print *,' ier           = ',ier
c         print *,' type 2 error  = ',err
c         print *,' '
c
c     -----------------------
c     call 1D Type3 method
c     -----------------------
         do k1 = 1, ms
            sk(k1) = 48*dcos(k1*pi/ms)
         enddo
c         call dirft1d3(nj,xj,cj,iflag, ms,sk,fk0)
         print *,'type 3 :'
	 open(unit=10,file='M3rh.txt')
	 read(10,*) M1
	 open(unit=20,file='M3ih.txt')
	 read(20,*) M2
	 a=M1(1,:)
	 b=M2(1,:)
	 in=dcmplx(a,-b)
	 call system_clock(time_begin,countrage,countmax)
	 do ip=1,num
	  call dfftw_execute_dft(plan, in, out)
	 enddo
	 call system_clock(time_end,countrage,countmax)
	 print *,' T_fft         = ',(time_end-time_begin)/num
	 print *,' r*T_fft       = ',r*(time_end-time_begin)/num

         call system_clock(time_begin1,countrage,countmax)
         do ip=1,num
           call nufft1d3f90(nj,xj,cj,iflag,eps, ms,sk,fk1,ier)
         enddo
         call system_clock(time_end1,countrage,countmax)
         print *,' T_g           = ',(time_end1-time_begin1)/num
         ratio=1.000*r*(time_end-time_begin)/(time_end1-time_begin1)
         print *,' T_fft / T_g   = ',ratio/r
         print *,' r*T_fft / T_g = ',ratio
         
c         call errcomp(cj0,cj1,nj,err)
c         print *,' ier           = ',ier
c         print *,' type 3 error  = ',err
c         print *,' '
      enddo
      stop
      end
c
c
c
c
c
      subroutine errcomp(fk0,fk1,n,err)
      implicit none
      integer k,n
      complex*16 fk0(n), fk1(n)
      real *8 salg,ealg,err
c
      ealg = 0d0
      salg = 0d0
      do k = 1, n
         ealg = ealg + cdabs(fk1(k)-fk0(k))**2
         salg = salg + cdabs(fk0(k))**2
      enddo
      err =sqrt(ealg/salg)
      return
      end
