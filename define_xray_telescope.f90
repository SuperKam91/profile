module define_xray_telescope

   use kind_def
   use sz_globals
   use maths

   implicit none

   public :: make_rosat_beam
   public :: make_chandra_beam
   public :: make_xmm_beam
   public :: rosat_beam
   public :: chandra_beam
   public :: xmm_beam
   public :: get_xray_telescope
   public :: make_xray_telescope

contains

!***************************************************************************

   subroutine make_rosat_beam

      implicit none
      real(kind=dp) :: temp1, temp2, r, value, centre
      integer :: i, j, k, status
      character*80 :: plot_device

!  set up beam
      centre = maxsize/2+1      
      do j = 1,maxsize
         do i = 1, maxsize
            temp1 = i-centre
            temp2 = j-centre
            r = sqrt(temp1**2+temp2**2)
            call rosat_beam(r*cellsize,value)
            ros_beam_re(i,j) = value
         end do
      end do
      ros_beam_im = 0.0

      rosat_beam_sum = sum(ros_beam_re)
      rosat_beam_sum = 1.0/(maxsize**2*rosat_beam_sum)

! FFT Rosat beam
      call do_2d_fft(ros_beam_re,ros_beam_im,maxsize,maxsize,1)

   end subroutine make_rosat_beam

!*************************************************************************

! 09/03/04 AMS
! This subroutine makes the chandra beam and FFTs it

   subroutine make_chandra_beam

      implicit none
      real(kind=dp) :: temp1, temp2, r, value, centre
      integer :: i, j, k, status
      character*80 :: plot_device

!  set up beam

      centre = maxsize/2+1      
      do j = 1,maxsize
         do i = 1, maxsize
            temp1 = i-centre
            temp2 = j-centre
            r = sqrt(temp1**2+temp2**2)
            call chandra_beam(r*cellsize,value)
            ros_beam_re(i,j) = value
         end do
      end do
      ros_beam_im = 0.0

      rosat_beam_sum = sum(ros_beam_re)
      rosat_beam_sum = 1.0/(maxsize**2*rosat_beam_sum)

! FFT Chandra beam
      call do_2d_fft(ros_beam_re,ros_beam_im,maxsize,maxsize,1)

   end subroutine make_chandra_beam

! **************************************************************************

   subroutine make_xmm_beam(alpha, r0)

      implicit none
      real(dp), intent(inout) :: alpha, r0
      real(kind=dp) :: temp1, temp2, r, value, centre
      integer :: i, j, k, status
      character*80 :: plot_device

!  set up beam

      centre = maxsize/2+1      
      do j = 1,maxsize
         do i = 1, maxsize
            temp1 = i-centre
            temp2 = j-centre
            r = sqrt(temp1**2+temp2**2)
            call xmm_beam(r*cellsize,value, r0, alpha)
            ros_beam_re(i,j) = value
         end do
      end do
      ros_beam_im = 0.0

      rosat_beam_sum = sum(ros_beam_re)
      rosat_beam_sum = 1.0/(maxsize**2*rosat_beam_sum)

! FFT XMM beam
      call do_2d_fft(ros_beam_re,ros_beam_im,maxsize,maxsize,1)

   end subroutine make_xmm_beam

! **************************************************************************

!  returns value of approximate rosat beam at radius r arcsec
   subroutine rosat_beam(r,taper)

! NB This is an approximation - a rather better to be inserted at some point 
! in the future as given in ROSAT user guide.

      implicit none
      real(kind=dp), intent(in) :: r
      real(kind=dp), intent(out) :: taper
      real(kind=dp) :: r2

      r2 = r*2.0

!   taper = 0.78*exp(-r*r)

      taper = 0.9*exp(-r2*r2*7.701635e-4)+0.1*exp(-r2*r2*2.43547e-4)

      if (taper.lt.1e-3) taper = 0.0

   end subroutine rosat_beam

! ********************************************************************

! Added 09/03/04 AMS

! This is also an approximation

   subroutine chandra_beam(r,taper)

      implicit none
      real(kind=dp), intent(in) :: r
      real(kind=dp), intent(out) :: taper
      real(kind=dp) :: r2

      r2 = r*2.0

      taper = 0.78*exp(-r2*r2)

      if (taper.lt.1e-3) taper = 0.0

   end subroutine chandra_beam

! ********************************************************************

! Added 11/11/04 AMS

! This is also an approximation. The parameters of the King model 
! are energy dependent and those used here are only averages over
! the full energy range of each instrument.
! Also, the PSF is still being fitted. There is an intention to
! model it as a King profile plus a Gaussian to account for central
! excess.

   subroutine xmm_beam(r, taper, r0, alpha)

      implicit none
      real(dp), intent(in) :: r
      real(dp), intent(out) :: taper
      character*80 :: inst
      integer :: status
      real(dp), intent(in) :: alpha, r0
      real(dp) :: norm

!   norm = 0.015
      norm = 633

      taper = (1+(r/(cellsize*r0))**2)**(-1.0*alpha)
      taper = norm*taper

      if (taper.lt.1e-3) taper = 0.0

   end subroutine xmm_beam

! ***************************************************************
   
   subroutine get_xray_telescope

      implicit none

      integer :: status

      write(*,*) 'Possible X-ray telescope:'
      write(*,*) 'C  - Chandra'
      write(*,*) 'P  - PSPC on ROSAT'
      write(*,*) 'H  - HRI on ROSAT'
      write(*,*) 'X1 - mos1 on XMM'
      write(*,*) 'X2 - mos2 on XMM'
      write(*,*) 'XP - pn on XMM'
      write(*,*)
      call io_getc('Select an X-ray telescope','*',x_tel,status)

      call make_xray_telescope

      call io_getd('Exposure Time (s):','*',x_exp,status)
      call io_getd('X-ray emission constant:','*',const,status)

   end subroutine get_xray_telescope

! ***************************************************************

   subroutine make_xray_telescope

      implicit none
      real(dp) :: alpha, r0      

      select case(x_tel)
      case('c','C')
         call make_chandra_beam
         x_beam = .true.
      case('p','P')
         call make_rosat_beam
         x_beam = .true.
      case('x1','X1')
         alpha = 1.35
         r0 = 3.65
         call make_xmm_beam(alpha, r0)
         x_beam = .true.
      case('x2','X2')
         alpha = 1.4
         r0 = 3.65
         call make_xmm_beam(alpha, r0)
         x_beam = .true.
      case('xp','XP')
         alpha = 1.55
         r0 = 5.25
         call make_xmm_beam(alpha, r0)
         x_beam = .true.
      case default
         x_beam = .false.
      end select

   end subroutine make_xray_telescope

! ***************************************************************

end module define_xray_telescope
