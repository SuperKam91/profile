subroutine mcg_x_fit

   use kind_def
   use sz_globals
   use mcg_xray
   use kgraph
   use map_handling
   use file_handling

   implicit none

   integer :: num_unknowns
   real(kind=dp), dimension(:), allocatable :: parms
   logical, dimension(:), allocatable :: large_chng
   real(kind=dp), dimension(maxsize,maxsize) :: nat_wt 
   character*1 :: chr1
   integer :: num_steps
   real(kind=dp) :: countlevel, countmax, countmin
   integer :: i,j

! Switch off messages writen while making skies
   verbose = .false.

! Ensure that model sky is always recalculated from scratch
   overwrite = .true.

! Switch off plots
   do_plot = .false.

! Determine whether or not to fit to beta
   call io_getc('Allow beta to vary in fit (y/n):','n',chr1,status)
   if ((chr1=='n').or.(chr1=='N')) then
      vary_beta = .false.
   else
      vary_beta = .true.
   end if

! Determine whether or not to fit to position
   call io_getc('Allow position to vary in fit (y/n):','n',chr1,status)
   if ((chr1=='n').or.(chr1=='N')) then
      vary_posn = .false.
   else
      vary_posn = .true.
   end if

! Determine which areas of the Rosat map to fit to
   call io_getc('Fit to entire Rosat map (y/n):','n',chr1,status)

   if ((chr1=='n').or.(chr1=='N'))then

! Define region by hand or use previously saved version?
      call io_getc('Read in region from file (y/n):','y',chr1,status)

      if ((chr1=='y').or.(chr1=='Y')) then

! Load region
         call read_flagpixel(ok_data,maxsize)

         call io_getc('Modify region (y/n):','n',chr1,status)
         if ((chr1=='y').or.(chr1=='Y')) then

! Display rosat map
            write(*,*) 'Adjust grey scale of map (press q when finished)'
            call display_map(rosat_map,maxsize,graph_adj,&
                             cellsize,cellsize,0D0,0D0)
            call mod_map_region
! Save this region?
            call io_getc('Write this region to file (y/n):','y',&
                 & chr1,status)
            if ((chr1=='y').or.(chr1=='Y')) then
               call write_flagpixel(ok_data,maxsize)
            end if
         end if

      else

! Display rosat map
         write(*,*) 'Adjust grey scale of map (press q when finished)'
         call display_map(rosat_map,maxsize,graph_adj,&
                          cellsize,cellsize,0D0,0D0)

!          write(*,*) 'debugging here...'

! Define area by hand, or by the numbers?
         call io_getc('Fit down to a count level, or draw an area, or both (c/d/b):','d',chr1,status)
         select case (chr1)
         case ('d','D')
            call def_map_region
         case ('c','C')
            write(*,*) 'NB : USE WITH CARE - MAY BIAS FITTING PROCESS...'
            countmax = maxval(rosat_map)
            countmin = minval(rosat_map)
            write (*,*) ' Range of counts in map is from',countmin,' to ',countmax 
            call io_getd('Count level:','*',countlevel,status)
            do i=1, maxsize
               do j=1, maxsize 
                  if (rosat_map(i,j) < countlevel) then
                     ok_data(i,j) = .false.
                  else 
                     ok_data(i,j) = .true.
                  end if
               end do
            end do
         case ('b','B')
            write(*,*) 'NB : USE WITH CARE - MAY BIAS FITTING PROCESS...'
            call def_map_region 
            countmax = maxval(rosat_map)
            countmin = minval(rosat_map)
            write (*,*) ' Range of counts in map is from',countmin,' to ',countmax 
            call io_getd('Count level:','*',countlevel,status)
            do i=1, maxsize
               do j=1, maxsize 
                  if (rosat_map(i,j) < countlevel) then
                     ok_data(i,j) = .false.
!                   else                        Can't do this here, otherwise
!                      ok_data(i,j) = .true.  excluded regions will be included
                  end if
               end do
            end do
!             call def_map_region 
         end select


! Determine which area of the map to use for fitting
!          call def_map_region

! Save this region?
         call io_getc('Write this region to file (y/n):','y',&
              & chr1,status)
         if ((chr1=='y').or.(chr1=='Y')) then
            call write_flagpixel(ok_data,maxsize)
         end if
      end if
   else

! Use entire map
      ok_data = .true.
   end if
! Determine which areas of the Rosat map to fit to
!   call io_getc('Fit to entire Rosat map (y/n):','n',chr1,status)

!   if ((chr1=='n').or.(chr1=='N'))then

! Display rosat map
!      write(*,*) 'Adjust grey scale of map (press q when finished)'
!      call display_map(rosat_map,maxsize,graph_adj,cellsize,cellsize,0.0,0.0)

! Determine which area of the map to use for fitting
!      call def_map_region

!   else
!      ok_data = .true.
!   end if

! Ask whether or not to adaptively smooth the model x-ray map
   call io_getc('Apply smoothing to model x-ray map (y/n):','n',chr1,status)
   if ((chr1=='y').or.(chr1=='Y')) then
      smooth_x_model = .true.
      call io_getd('min conv size? (arcsec)','5',mingauss,status)
      call io_getd('max conv size? (arcsec)','60',maxgauss,status)
      mingauss = mingauss/cellsize
      maxgauss = maxgauss/cellsize
   else
      smooth_x_model = .false.
   end if

! Adjust the weighting given to every point on the map
   write(*,*) 'Available weightings : p(oisson); u(niform)'
   call io_getc('What weighting for map points (p/u)','u',chr1,status)
   weighting : select case(chr1)

! Poisson weights
   case('p')
      nat_wt = rosat_map-back
      where (nat_wt==0.0)
         nat_wt = 0.0
      elsewhere
         nat_wt = 1.0/nat_wt
      end where


! Uniform weights
   case default
      nat_wt = 1.0
   end select weighting

! Call fitting routine dependent upon what type of model is being fitted to
   if (ellipsoid) then
      write(*,*) 'Fitting ellipsoidal model'
      if (vary_beta) then
         num_unknowns = 7

! Allocate array for the parameters
         allocate (parms(num_unknowns))
         allocate (large_chng(num_unknowns))

! Define initial guess
         parms(1) = n0
         parms(2) = beta   
         parms(3) = theta_c(1)
         parms(4) = theta_c(2)
         parms(5) = el_alpha
         parms(6) = sz_x
         parms(7) = sz_y

! Make sure that only the position can be changed by a large amount in one
! iteration
         large_chng = .false.
         large_chng(6) = .true.
         large_chng(7) = .true.
      else
         num_unknowns = 6

! Allocate array for the parameters
         allocate (parms(num_unknowns))
         allocate (large_chng(num_unknowns))

! Define initial guess
         parms(1) = n0
         parms(2) = theta_c(1)
         parms(3) = theta_c(2)
         parms(4) = el_alpha
         parms(5) = sz_x
         parms(6) = sz_y

! Make sure that only the position can be changed by a large amount in one
! iteration
         large_chng = .false.
         large_chng(5) = .true.
         large_chng(6) = .true.
      end if

! Find best fit solution
      call x_root_mcg(maxsize,maxsize,rosat_map,ok_data,nat_wt,parms,&
           & large_chng,num_steps,num_unknowns)

! Redefine parameters as the best fit values
      if (vary_beta) then
         n0 = parms(1)
         beta = parms(2)      
         theta_c(1) = parms(3)
         theta_c(2) = parms(4) 
         el_alpha = parms(5)
         sz_x = parms(6)
         sz_y = parms(7)
      else
         n0 = parms(1)
         theta_c(1) = parms(2)
         theta_c(2) = parms(3) 
         el_alpha = parms(4)
         sz_x = parms(5)
         sz_y = parms(6)
      end if

      write(*,*) 'Found fit to rosat map after ',num_steps, ' iterations'
      write(*,*) 'central density     ', n0
      write(*,*) 'beta                ', beta
      write(*,*) 'theta 1             ', theta_c(1)
      write(*,*) 'theta 2             ', theta_c(2)
      write(*,*) 'inclination angle   ', el_alpha/deg2rad
      write(*,*) 'Offset in RA        ', sz_x
      write(*,*) 'Offset in Dec       ', sz_y
   else

! Do spherical model fit
      write(*,*) 'Fitting spherical model'
      if (vary_beta) then
         num_unknowns = 5

! Allocate arrays for parameter
         allocate (parms(num_unknowns))
         allocate (large_chng(num_unknowns))

! Define initial guess
         parms(1) = n0
         parms(2) = beta   
         parms(3) = theta_cx
         parms(4) = sz_x
         parms(5) = sz_y

! Make sure that only the position can be changed by a large amount in one
! iteration
         large_chng = .false.
         large_chng(4) = .true.
         large_chng(5) = .true.
      else
         num_unknowns = 4

! Allocate arrays for parameter
         allocate (parms(num_unknowns))
         allocate (large_chng(num_unknowns))

! Define initial guess
         parms(1) = n0
         parms(2) = theta_cx
         parms(3) = sz_x
         parms(4) = sz_y

! Make sure that only the position can be changed by a large amount in one
! iteration
         large_chng = .false.
         large_chng(3) = .true.
         large_chng(4) = .true.
      end if

! Find best fit solution
      call x_root_mcg(maxsize,maxsize,rosat_map,ok_data,nat_wt,parms,&
           & large_chng,num_steps,num_unknowns)

! Redefine parameters as the best fit values
      if (vary_beta) then
         n0 = parms(1)
         beta = parms(2)
         theta_cx = parms(3)
         sz_x = parms(4)
         sz_y = parms(5)
      else
         n0 = parms(1)
         theta_cx = parms(2)
         sz_x = parms(3)
         sz_y = parms(4)
      end if

      write(*,*) 'Found fit to rosat map after ',num_steps, ' iterations'
      write(*,*) 'central density     ', n0
      write(*,*) 'beta                ', beta
      write(*,*) 'theta               ', theta_cx
      write(*,*) 'Offset in RA        ', sz_x
      write(*,*) 'Offset in Dec       ', sz_y
   end if

! Switch back on messages written while making skies
   verbose = .true.

   deallocate(parms)
   deallocate(large_chng)

end subroutine mcg_x_fit



