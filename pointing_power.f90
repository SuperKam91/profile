
subroutine pointing_power

   use kind_def
   use sz_globals
   use physics
   use kgraph
   use file_handling
   use map_handling

   implicit none

   integer            :: i,j,offset,ncell
   real               :: gau, mappower, mappower2
   real               :: norm1, norm2, tr(6), radius
!   real(dp), allocatable  :: temp_sky(:,:)
   character*80       :: chr1
   integer  pgopen
   external pgopen

   ncell = maxsize
   
   allocate(temp_sky(maxsize,maxsize))
   temp_sky = szsky(:,:,nchan)

! multiply map by Gaussian at center:

   do i = 1, maxsize
      do j = 1, maxsize
         gau = exp(-1.*(abs(i-(maxsize/2))**2+abs(j-(maxsize/2))**2)**(0.5)/(2*31**2))
         temp_sky(i,j) = temp_sky(i,j)*gau
      end do
   end do
   
   norm1 = temp_sky(maxsize/2,maxsize/2)

   mappower = sum(temp_sky**2)
   write(*,*) 'Total power in map:', mappower
   write(*,*) 'Displaying sky'
   call display_map(temp_sky, ncell, graph_adj,cellsize,cellsize,0D0,0D0)


   call io_geti('Pointing offset (amin):','10',offset,status)
   write(*,*) '                     ---'
   write(*,*) 'Pointing will be offset to the NORTH - direction'
   write(*,*) '  should not matter because beam is symmetric.'
   write(*,*) '                     ---'
   
!convert arcmin to pixels:
   offset = offset*60./cellsize

! get back the original sky:
   temp_sky = szsky(:,:,nchan)

! multiply by new offset Gaussian:
   do i = 1, maxsize
      do j = 1, maxsize
         gau = exp(-1.*(abs(i-(maxsize/2))**2+abs(j-((maxsize/2)+offset))**2)**(0.5)/(2*31**2))
         temp_sky(i,j) = temp_sky(i,j)*gau
      end do
   end do
   
   norm2 = (norm1/temp_sky(maxsize/2,maxsize/2))
   temp_sky = norm2*temp_sky

   mappower2 = sum(temp_sky**2)
   write(*,*) 'Total power in map:', mappower2
   write(*,*) 'Ratio of powers:', mappower2/mappower
   write(*,*) 'Displaying sky'
   call display_map(temp_sky, ncell, graph_adj,cellsize,cellsize,0D0,0D0)


   deallocate(temp_sky)

   call io_getc('Do you want to create a false pointing offset in your map?','y',chr1,status)
   if (chr1=='y') then
      do i = 1, maxsize
         do j = 1, maxsize
            gau = exp(-1.*(abs(i-(maxsize/2))**2+abs(j-((maxsize/2)+offset))**2)**(0.5)/(2*31**2))
            szsky(i,j,nchan) = szsky(i,j,nchan)*gau
         end do
      end do
      szsky(:,:,nchan) = norm2*szsky(:,:,nchan) ! to account for offset calibrator
      norm2 = cellsize*sec2rad
      norm2 = norm2*norm2
      norm2 = norm2*2.0*k_b*obsfreq**2/const_c2*1.0D26
      szsky(:,:,nchan) = norm2*szsky(:,:,nchan) ! convert from temperature to Jy
   end if

   write(*,*) 'Displaying sz sky'
   call display_map(szsky(:,:,nchan), ncell, graph_adj,cellsize,cellsize,0D0,0D0)


! Plotting routine for thesis [AMS]:

!   tr = 0.
!   tr(2) = 1.
!   tr(6) = 1.
!   radius = 1.2*60.*60./cellsize
!   i = pgopen('pointing_offset.ps/CPS')
!   call pgenv(1.,real(maxsize),1.,real(maxsize),1,-1)
!   call pggray(real(szsky),512,512,1,512,1,512,maxval(real(szsky)),minval(real(szsky)),tr)
!   call pgsfs(2)
!   call pgcirc(real(maxsize/2),real(maxsize/2),radius)
!   call pgpt1(real(maxsize/2),real(maxsize/2),2)
!   call pgsls(2)
!   call pgcirc(real(maxsize/2),real((maxsize/2)+offset),radius)
!   call pgpt1(real(maxsize/2),(real(maxsize/2)+offset),2)
!   call pgsls(1)
!   call pgclos()
!      


end subroutine pointing_power
