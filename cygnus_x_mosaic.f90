
! this program calculates the central positions of 23 fields used to
! make a mosaic over the cygnus X region. it calculates it first in arcmin
! and then changes coordinates to RA and dec

! it is written into the Makefile

program cygX_mosaic


   implicit none

   integer, parameter :: mapsize = 512 
   real, dimension(mapsize,mapsize) :: sensi  !sensitivity
   real, dimension(2,23) :: position  !duh
   real, dimension(2,23) :: centre
   real, parameter :: pi = 3.14159265
   real :: sigma, gaussian,distance,fwhm, count_deg,count,spacing
   integer ::  i, x, y, n_beam, status, just, np
   real :: maxlev, minlev, x1, x2, y1, y2, deg2rad
   real :: amin2deg ! arcmin to degrees
   real :: amin2RAmin ! arcmin to minutes of RA
   integer, allocatable :: RA_h(:),RA_min(:)
   integer, allocatable :: dec_deg(:), dec_min(:)
   real, allocatable :: RA_sec(:), dec_sec(:)
   real, allocatable :: RA_min_pt(:), dec_min_pt(:)

   real, dimension(6) :: tr
   character*80 plot_device
   integer  pgopen
   external pgopen

   open(2, file = 'mosaic_grid.txt')

   fwhm = 1.2 ! degrees
   fwhm = fwhm*60 ! arcmins
   spacing = 0.8*fwhm ! arcmins
   sigma = fwhm/(2*sqrt(2*log(2.0))) ! arcmins
   n_beam = 23 ! number of fields
   deg2rad = pi/180.
   amin2deg = 1/60.
   amin2RAmin = 1/15.

   do x = 1,mapsize
      do y = 1,mapsize
         sensi(x,y) = 0.0
      end do
   end do

! define centre as (x/y, circle no.)
! there are 27 fields in this mosaic...


   centre(1,1) = mapsize/2
   centre(2,1) = mapsize/2

   do i = 2,7
      centre(1,i) = mapsize/2+spacing*Sin(real((i-2)*60*deg2rad))
!      write(*,*) i, real((i-2)*60), sin(real((i-2)*60*deg2rad))
      centre(2,i) = mapsize/2+spacing*Cos(real((i-2)*60*deg2rad))
!      write(*,*) centre(1,i), centre(2,i)
   end do

   do i = 8,13
      centre(1,i) = mapsize/2+2*cos(30.*deg2rad)*spacing*sin(real(((i-8)*60)+30)*deg2rad) 
      centre(2,i) = mapsize/2+2*cos(30.*deg2rad)*spacing*cos(real(((i-8)*60)+30)*deg2rad)
   end do

   do i = 14,19
      centre(1,i) = mapsize/2+2*spacing*sin(real((i-14)*60)*deg2rad)
      centre(2,i) = mapsize/2+2*spacing*cos(real((i-14)*60)*deg2rad)
   end do


   centre(1,20) = mapsize/2+2*(spacing/cos(40.89*deg2rad))*sin(40.89*deg2rad)
   centre(2,20) = mapsize/2+2*(spacing/cos(40.89*deg2rad))*cos(40.89*deg2rad)

   centre(1,21) = mapsize/2+2*(spacing/cos(40.89*deg2rad))*sin((180-40.89)*deg2rad)
   centre(2,21) = mapsize/2+2*(spacing/cos(40.89*deg2rad))*cos((180-40.89)*deg2rad)

   centre(1,22) = mapsize/2+2*(spacing/cos(40.89*deg2rad))*sin((180+40.89)*deg2rad)
   centre(2,22) = mapsize/2+2*(spacing/cos(40.89*deg2rad))*cos((180+40.89)*deg2rad)

   centre(1,23) = mapsize/2+2*(spacing/cos(40.89*deg2rad))*sin(-40.89*deg2rad)
   centre(2,23) = mapsize/2+2*(spacing/cos(40.89*deg2rad))*cos(-40.89*deg2rad)


   do np = 1, n_beam
      do x = 1,mapsize
         do y = 1,mapsize

            position(1,np) = x
            position(2,np) = y

! find the distance from the centre
            distance = ((position(1,np)-centre(1,np))**2+&
                        (position(2,np)-centre(2,np))**2)**0.5

! calculate value of beam at this position
            gaussian = exp((-(distance**2))/(2*sigma**2))

! add in contribution from this pointing to sensitivity map
            sensi(x,y) = (sensi(x,y)**2 + gaussian**2)**0.5 

         end do
      end do
   end do

! count how many pixels have sensitivities above 0.5
   count = 0

   do x = 1,mapsize
      do y = 1,mapsize
         if (sensi(x,y) > 0.5) then
            count = count + 1
         else
            count = count
         endif

      end do
   end do

   count_deg = count/3600.0

   write (*,*) 'coverage in square degrees = ', count_deg

   call io_getc('Plot device:','/xwindow',plot_device,status)
   call pgbegin(0,plot_device,1,1)

   maxlev = maxval(sensi)
   write(*,*) 'Maximum level is ',maxlev
   minlev = 0.0
   just = 1

   tr(1) = 0
   tr(4) = 0
   tr(2) = 1.0
   tr(6) = 1.0

   x1 = 1.0
   x2 = mapsize
   y1 = 1.0
   y2 = mapsize

   call pgenv(x1,x2,y1,y2,just,0)

   call pggray(sensi,mapsize,mapsize,1,mapsize,1,mapsize,&
        maxlev,minlev,tr)

   call pgsci(2)
   call pgsfs(2)
   do np = 1, n_beam
      call pgpt1(centre(1,np),centre(2,np), -4)
      call pgsls(4)
      call pgcirc(centre(1,np),centre(2,np),fwhm/2.0)
   end do

   call pgend

!   call io_getc('Plot device:','/xwindow',plot_device,status)
   i = pgopen('cygX_coverage.ps/PS')
!   call pgbegin(0,plot_device,1,1)

   maxlev = maxval(sensi)
   write(*,*) 'Maximum level is ',maxlev
   minlev = 0.0
   just = 1

   tr(1) = 0
   tr(4) = 0
   tr(2) = 1.0
   tr(6) = 1.0

   x1 = 1.0
   x2 = mapsize
   y1 = 1.0
   y2 = mapsize

   call pgenv(x1,x2,y1,y2,just,0)

   call pggray(sensi,mapsize,mapsize,1,mapsize,1,mapsize,&
        maxlev,minlev,tr)

   call pgsci(2)
   call pgsfs(2)
   do np = 1, n_beam
      call pgpt1(centre(1,np),centre(2,np), -4)
      call pgsls(4)
      call pgcirc(centre(1,np),centre(2,np),fwhm/2.0)
   end do
   call pgclos()
   i = pgopen('/xwindow')
!   call pgend

! convert to RA and dec

   centre(1,1) = 30 ! minutes of RA
   centre(2,1) = 41 ! degrees

   do i = 2,7
      centre(1,i) = 30 + amin2RAmin*spacing*Sin(real((i-2)*60*deg2rad)) 
      centre(2,i) = 41 + amin2deg*spacing*Cos(real((i-2)*60*deg2rad))
      write(*,*) centre(1,i), centre(2,i)
   end do

   do i = 8,13
      centre(1,i) = 30 + amin2RAmin*2*cos(30.*deg2rad)*spacing*sin(real(((i-8)*60)+30)*deg2rad) 
      centre(2,i) = 41 + amin2deg*2*cos(30.*deg2rad)*spacing*cos(real(((i-8)*60)+30)*deg2rad)
   end do

   do i = 14,19
      centre(1,i) = 30 + amin2RAmin*2*spacing*sin(real((i-14)*60)*deg2rad)
      centre(2,i) = 41 + amin2deg*2*spacing*cos(real((i-14)*60)*deg2rad)
   end do


   centre(1,20) = 30 + amin2RAmin*2*(spacing/cos(40.89*deg2rad))*sin(40.89*deg2rad)
   centre(2,20) = 41 + amin2deg*2*(spacing/cos(40.89*deg2rad))*cos(40.89*deg2rad)

   centre(1,21) = 30 + amin2RAmin*2*(spacing/cos(40.89*deg2rad))*sin((180-40.89)*deg2rad)
   centre(2,21) = 41 + amin2deg*2*(spacing/cos(40.89*deg2rad))*cos((180-40.89)*deg2rad)

   centre(1,22) = 30 + amin2RAmin*2*(spacing/cos(40.89*deg2rad))*sin((180+40.89)*deg2rad)
   centre(2,22) = 41 + amin2deg*2*(spacing/cos(40.89*deg2rad))*cos((180+40.89)*deg2rad)

   centre(1,23) = 30 + amin2RAmin*2*(spacing/cos(40.89*deg2rad))*sin(-40.89*deg2rad)
   centre(2,23) = 41 + amin2deg*2*(spacing/cos(40.89*deg2rad))*cos(-40.89*deg2rad)

   call io_getc('Plot device:','/xwindow',plot_device,status)
   call pgbegin(0,plot_device,1,1)

!   maxlev = maxval(sensi)
!   write(*,*) 'Maximum level is ',maxlev
!   minlev = 0.0

   just = 1

   tr = 0
   tr(2) = 1
   tr(6) = 1

   x1 = 10
   x2 = 50
   y1 = 37
   y2 = 45

   call pgenv(x1,x2,y1,y2,just,0)

!   call pggray(sensi,mapsize,mapsize,1,mapsize,1,mapsize,&
!        maxlev,minlev,tr)

   call pgsci(2)
   call pgsfs(2)
   do np = 1, n_beam
      call pgpt1(centre(1,np),centre(2,np), -4)
!      call pgsls(4)
!      call pgcirc(centre(1,np),centre(2,np),fwhm/2.0)
   end do

   call pgend

   allocate(RA_h(np))
   allocate(RA_min_pt(np))
   allocate(RA_min(np))
   allocate(RA_sec(np))
   allocate(dec_deg(np))
   allocate(dec_min_pt(np))
   allocate(dec_min(np))
   allocate(dec_sec(np))

   do i = 1, np
      RA_h(i) = 20
      RA_min_pt(i) = centre(1,i)
      RA_min(i) = int(RA_min_pt(i))
      RA_sec(i) = (RA_min_pt(i) - RA_min(i))*60

      dec_deg(i) = int(centre(2,i))
      dec_min_pt(i) = (centre(2,i) - dec_deg(i))*60
      dec_min(i) = int(dec_min_pt(i))
      dec_sec(i) = (dec_min_pt(i) - dec_min(i))*60
   end do

   do i = 1, np
      write(2,*)'SE_CYGX_',i,RA_h(i),RA_min(i),RA_sec(i), dec_deg(i),dec_min(i),dec_sec(i)
   end do


end program cygX_mosaic
