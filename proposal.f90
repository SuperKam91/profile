
double precision function proposal(idum,n,ndim,cube,width)
   implicit none
   
   integer  idum,n,ndim,nc
   double precision cube(ndim),width(ndim)
   
   double precision ran,gdev
!      external ran,gdev
!      proposal=ran(idum) ! uniform proposal distribution
   proposal=cube(n)+width(n)*gdev(idum) ! Gaussian prop. dist.
   
end function proposal

! -----------------------------------------------------------------

   function proposal2(idum,n,ndim,cube,width)
      
      use kind_def
      
      implicit none
      
      real(dp) :: proposal2
      integer :: idum,n,ndim,nc
      real(dp) :: cube(ndim),width(ndim)
      
      real(dp) :: gdev
! !         external gdev
      
! !      proposal=ran(idum) ! uniform proposal distribution
      proposal2=cube(n)+(width(n)/sqrt(real(ndim)))*gdev(idum) ! Gaussian 
      

   end function proposal2
       
! ------------------------------------------------------------------

       FUNCTION gdev(idum)

          use kind_def

          implicit none

       INTEGER :: idum
       real(dp) :: gdev
 !CU    USES ran
       INTEGER :: iset
       real(dp) :: fac,gset,rsq,v1,v2,ran
       SAVE iset,gset
       DATA iset/0/

 !      external ran
       if (iset.eq.0) then
 1       v1=2d0*ran(idum)-1d0
         v2=2d0*ran(idum)-1d0
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=dsqrt(-2d0*dlog(rsq)/rsq)
         gset=v1*fac
         gdev=v2*fac
         iset=1
       else
         gdev=gset
         iset=0
       endif

    END FUNCTION gdev

! --------------------------------------------------------------------

! assigns a value between 1 and n to a discrete integer random variable 
! according to the probabilities p(n)

       integer function irv(idum,n,p)
       implicit none

       integer  idum,n
       double precision p(n)

       integer  i
       double precision dum,tmp,ran,gdev
 !      external ran,gdev
       dum=ran(idum)

       i=0
       tmp=0d0
  10   i=i+1
       tmp=tmp+p(i)
       if (dum.gt.tmp) goto 10
       irv=i


    end function irv

! ---------------------------------------------------------

     FUNCTION ran(idum)

          use kind_def

          implicit none

!      INTEGER :: idum,IA,IM,IQ,IR,NTAB,NDIV
!      real(dp) ::  ran,AM,EPS,RNMX
       INTEGER :: idum
       real(dp) ::  ran
       integer, PARAMETER :: IA=16807
       integer, parameter :: IM=2147483647
       real(dp), parameter :: AM=1d0/IM
       integer, parameter :: IQ=127773
       integer, parameter :: IR=2836
       integer, parameter :: NTAB=32
       integer, parameter :: NDIV=1+(IM-1)/NTAB
       real(dp), parameter :: EPS=1.2d-7
       real(dp), parameter :: RNMX=1d0-EPS
       INTEGER :: j,k,iv(NTAB),iy
       SAVE iv,iy
       DATA iv /NTAB*0/, iy /0/

       real(dp) :: gdev
! !      external gdev
! write(*,*) 'got to ran'
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
       ran=min(AM*iy,RNMX)


 end FUNCTION ran

!-------------------------------------------------------------

  double precision function dmin(x,y)
      implicit none
      double precision x,y

      double precision gdev
!      external gdev

      if (x.le.y) then
         dmin=x
      else
         dmin=y
      end if


   end function dmin

!-------------------------------------------------------------
