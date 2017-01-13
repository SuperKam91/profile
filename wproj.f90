
module pswf_globals

   ! stuff for PSWF visibility convolution:
   integer, parameter :: support = 3
   integer, parameter :: oversample = 128
   integer, parameter :: csize = 2*support+1  ! 7
   integer, parameter :: ccentre = support+1    ! 4
   integer, parameter :: nplanes = oversample*oversample
   integer            :: wstatus
   integer            :: qnx, qny
   integer            :: wsupport
   integer            :: wcentre
   complex, allocatable  :: convfunc(:,:,:)
   
end module

! --------------------------------------------------------------------------

module w_proj

contains

   subroutine wproj(l,m,wu,wre,wim)

! evaluates G(l,m,w), equ. 5 from Cornwell, Golap & Bhatnagar, Astronomical data analysis software and systems XIV, ASP Conference series, Vol. 347, 2005.
      
      use kind_def
      use sz_globals
      
      real(dp),intent(in)  :: l,m,wu
      real(dp),intent(out) :: wre,wim
      real(dp) :: exponent
      
      
      exponent = -1d0*2*pi*(wu*(sqrt(real(1-l**2-m**2)) - 1))
      wre = cos(exponent)
      wim = sin(exponent)
      
      
   end subroutine wproj
   
! --------------------------------------------------------------------------
   
   subroutine pswf(vnu,value)
      
! evaluates the PROLATE SPHEROIDAL WAVEFUNCTION with m=6, alpha = 1 from Schwab, Indirect Imaging. Code adapted from Tim Cornwell's C++ SphFuncVisGridder developed for CONRAD for ASKAP.
      
      use kind_def
      use sz_globals
  
      implicit none
    
      integer              :: k,part
      real(dp)             :: top,bot,delnusq,nuend
      real(dp)             :: factor
      real(dp),intent(in)  :: vnu
      real(dp),intent(out) :: value 
      integer, parameter   :: np = 4
      integer, parameter   :: nq = 2
      real(dp)             :: p(2,5)
      real(dp)             :: q(2,3) 
      
      p(1,1) = 8.203343e-2
      p(1,2) = -3.644705e-1
      p(1,3) = 6.278660e-1
      p(1,4) = -5.335581e-1
      p(1,5) = 2.312756e-1
      p(2,1) = 4.028559e-3
      p(2,2) = -3.697768e-2
      p(2,3) = 1.021332e-1
      p(2,4) = -1.201436e-1
      p(2,5) = 6.412774e-2

      q(1,1) = 1.0000000
      q(1,2) = 8.212018e-1
      q(1,3) = 2.078043e-1
      q(2,1) = 1.0000000
      q(2,2) = 9.599102e-1
      q(2,3) = 2.918724e-1

      value = 0.
      
      if ((vnu.ge.0.).and.(vnu.lt.0.75)) then
         part = 1
         nuend = 0.75
      elseif ((vnu.ge.0.75).and.(vnu.le.1.)) then
         part = 2
         nuend = 1.0
      else
         value = 0.
         return
      end if
      
      top = p(part,1)
      bot = q(part,1)
      delnusq = vnu**2 - nuend**2
      do k = 1,np
         factor = delnusq**k
         top = top+p(part,k+1)*factor
      end do
      do k = 1,nq
         factor = delnusq**k
         bot = bot+q(part,k+1)*factor
      end do
      
      if (bot.ne.0.) then
         value = top/bot
      else
         value = 0.
      end if
      
      if (value.lt.0.) value = 0.
      
      
   end subroutine pswf
   
   
! -------------------------------------------------------------------------

subroutine sinckernel(l,kernel)

! taken from Mapper

   implicit none

   real,intent(in)  :: l
   real             :: x
   real, parameter  :: sigma1 = 0.4
   real, parameter  :: sigma2 = 4.0
   real, intent(out) :: kernel

   
   x = l/sigma1
   kernel = 0.

   if (x.gt.0.) then
      kernel = ((sin(x)/x)*exp(-l*l/sigma2))
   else
      kernel = 1.
   end if
  

end subroutine sinckernel

! ------------------------------------------------------------------------

end module w_proj
