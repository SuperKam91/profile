module mikes_fft

!   use fft
   
contains
   
! -----*-----------------------------------------------------------------
   
      subroutine makefft(nx,ny,map,job,iform,work)
      implicit none
      integer  nx,ny,job,iform
      real     work(2*nx)
      complex  map(nx,ny) 

      integer  ndim,nn(2)

! ... initialise variables
      ndim=2
      nn(1)=nx
      nn(2)=ny

! ... calculate FFT
      call cheqboard(map,nx,ny)
      call fourt(map,nn,ndim,job,iform,work) 
      call cheqboard(map,nx,ny)
      call normfft(map,nx,ny,job)
      return
   end subroutine makefft

!-----*-----------------------------------------------------------------

      subroutine cheqboard(array,nx,ny)
      implicit none
      integer    nx,ny
      complex    array(nx,ny)
      
      integer    i,j,isign

      do i=1,nx
        do j=1,ny
          isign=(-1)**(i+j)
          array(i,j)=array(i,j)*float(isign)
        end do
      end do      

      return
   end subroutine cheqboard

! -----*-----------------------------------------------------------------

      subroutine normfft(array,nx,ny,job)
      implicit none
      integer  nx,ny,job
      complex  array(nx,ny)

      integer  i,j
      real     pi,factor

      pi=3.1415926535

      if (job.eq.1) then
        factor=1.                     ! forward transform (+ve exponential)
      else if (job.eq.-1) then
        factor=(1./real(nx))**2       ! backward transform (-ve exponential)
      end if

      do i=1,nx
        do j=1,ny
          array(i,j)=array(i,j)*factor
        end do
      end do

      return
   end subroutine normfft


! -----*-----------------------------------------------------------------

      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
   

end FUNCTION ran1

! -----*-----------------------------------------------------------------


end module mikes_fft
