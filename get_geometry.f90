! Changes position of primary beam and position of centre of cluster
! and the size and gridding of the map arrays
! 25-02-2004 hacked to include chandra beam option - AMS
! 06/08/09   removed choice of x-ray telescope [KG] 

subroutine get_geometry

   use kind_def
   use sz_globals
   use physics
   use define_telescope

   implicit none

   integer :: old_maxsize,i
   real(kind=dp) :: old_cellsize
   character*1 :: chr1,ch1
   real(dp) :: alpha, r0
   character*80 :: inst

   old_maxsize = maxsize
   old_cellsize = cellsize

   write(*,*) 'cellsize ',cellsize

   write(*,*) 'All displacements relative to pointing centre'
   write(*,*) 'Positive displacement shifts to the west in RA'
   write(*,*) 'Positive displacement shifts to the north in Dec'
   call io_getd('Displacement of Pointing Centre in RA (arc-seconds)&
   &:','*',pb_x,status)
   call io_getd('Displacement of Pointing Centre in Dec (arc-seconds)&
   &:','*',pb_y,status)
   call io_getd('Displacement of Phase Centre in &
                 &RA (arc-seconds):','*',ph_x,status)
   call io_getd('Displacement of Phase Centre in & 
                 &Dec (arc-seconds):','*',ph_y,status)
   if ((ph_x.ne.0.0).or.(ph_y.ne.0.0)) write(*,*) &
        'Are you sure? Generally not useful to move phase centre....'
   call io_geti('Size of sky arrays:','*',maxsize,status)
   call io_getd('Cellsize for model /arcsec:','*',cellsize,status)
   
! Warning if field of view to small
   if (maxsize*cellsize*sec2rad/2.0.lt.&
          1.0/(dish_diameter*obsfreq/const_c)) then
      write(*,*) 
      write(*,*) &
 '** WARNING ** Field of view too small - increase cellsize or number of cells'
      write(*,*)
   end if

! Check whether it is necessary to reallocate all global arrays
   if (old_maxsize.ne.maxsize) then
      call do_deallocation
      call do_allocation
      write(*,*) 'All maps now blanked'
   end if

end subroutine get_geometry




