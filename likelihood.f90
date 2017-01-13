module likelihood

   use kind_def
   use structures
   use maths
   use stuff
   use minimum
   use sz_globals
   use kgraph
   use map_handling
   use file_handling

   implicit none

   public :: plot_likelihood

contains

! *****************************************************************************

   subroutine plot_likelihood
! Changed 20/3/01 to do marginalisation properly KG

      integer :: num_unknowns
      real(kind=dp), dimension(:), allocatable :: parms, n_arr, b_arr, th_arr
      real(kind=dp), dimension(:,:,:), allocatable :: likeli 
      real(kind=dp), dimension(:,:), allocatable :: prob2, sz_flux
      real(kind=dp) :: n_min, n_max, th_min, th_max, b_min, b_max
      real(kind=dp) :: n_inc, th_inc, b_inc, probmax, ellip,dummy
      real(kind=dp) :: b_cent, th_cent
      real(kind=dp) :: mean_re, re_pt, im_pt, psi, inc
      real(kind=dp) :: max_like, u, v, rad
      integer :: npoints
      character*1 :: chr1
      integer :: status, i, j, k, which_par, cube_s, cube_3
      type(data) :: exp_data
      type(data), target :: model_data, temp_data
      type(param) :: par
      type(model) :: model_stuff, temp_model

      npoints = 20
      rad = 870.0

! Switch off messages writen while making skies
      verbose = .false.

! Ensure that model sky is always recalculated from scratch
      overwrite = .true.

! Switch off plots
      do_plot = .false.

! Ask for background level to use
      call io_getd('Background level:','*',back,status)

! Determine which areas of the Rosat map to fit to
      call io_getc('Fit to entire Rosat map (y/n):','n',chr1,status)

      if ((chr1=='n').or.(chr1=='N'))then

! Define region by hand or use previously saved version?
         call io_getc('Read in region from file (y/n):','y',chr1,status)

         if ((chr1=='y').or.(chr1=='Y')) then

! Load region
            call read_flagpixel(ok_data,maxsize)
         else

! Define region by hand

! Display rosat map
            write(*,*) 'Adjust grey scale of map (press q when finished)'
            call display_map(rosat_map,maxsize,graph_adj,&
                             cellsize,cellsize,0D0,0D0)

! Determine which area of the map to use for fitting
            call def_map_region

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

! Count up number of good pixels to be used in fitting
      num_dat_pt = count(ok_data)
      write(*,*) 'Using ',num_dat_pt,' pixels'

! Allocate space to arrays
      allocate(exp_data%values(num_dat_pt))
      allocate(model_data%values(num_dat_pt))
      model_stuff%m_data=>model_data
      allocate(temp_data%values(num_dat_pt))
      temp_model%m_data=>temp_data
      call create_data(maxsize,maxsize,rosat_map,exp_data)

! Do not do adaptive smoothing on the model x-ray map
      smooth_x_model = .false.

! Allow for the possibility to vary beta
      vary_beta = .true.

! Write out assumptions
      if (.not. ellipsoid) then
         ellipsoid = .true.
         theta_c = theta_cx
         write(*,*) 'Converting to ellipsoidal model'
      end if

      ellip = theta_c(1)/theta_c(2)
      write(*,*) 'Assumptions:'
      write(*,*) 'Ellipticity       ',ellip
      write(*,*) 'Inclination angle ',el_alpha
      write(*,*) 'RA offset         ',sz_x
      write(*,*) 'Dec offset        ',sz_y
      write(*,*) 'These parameters are marginalised over'
      write(*,*) 

      num_unknowns = 7

! Allocate array for the parameters
      allocate (parms(num_unknowns))
      

! Ask for size of likelihood cube to calculate
      cube_s = 10
      cube_3 = 20
      call io_geti('Size of likelihood cube','*',cube_s,status)
      call io_geti('Dimension in n0 direction','*',cube_3,status)

      allocate(likeli(cube_s,cube_s,cube_3))
      allocate(n_arr(cube_3))
      allocate(b_arr(cube_s))
      allocate(th_arr(cube_s))
      allocate(prob2(cube_s,cube_s))
      allocate(sz_flux(cube_s,cube_s))

      call io_getd('Minimum value of density:','8.0D-3',n_min,status)
      call io_getd('Maximum value of density:','12D-3',n_max,status)
      call io_getd('Minimum value of beta:','0.6',b_min,status)
      call io_getd('Maximum value of beta:','0.7',b_max,status)
      call io_getd('Minimum value of theta:','50.0',th_min,status)
      call io_getd('Maximum value of theta:','70.',th_max,status)
      call io_getd('Baseline length for SZ:','*',rad,status)
      
      n_inc = (n_max-n_min)/(cube_3-1)  
      b_inc = (b_max-b_min)/(cube_s-1)  
      th_inc = (th_max-th_min)/(cube_s-1) 
      b_cent = (b_max+b_min)/2.0
      th_cent = (th_max+th_min)/2.0

      do i = 1, cube_3
         n_arr(i) = n_min+(i-1)*n_inc
      end do
      do i = 1, cube_s
         b_arr(i) = b_min+(i-1)*b_inc
         th_arr(i) = th_min+(i-1)*th_inc
      end do
  
      call create_param(num_unknowns,parms,par)
      parms(2) = b_arr(1)
      parms(3) = th_arr(1)
      parms(4) = th_arr(1)*ellip
      parms(1) = n_arr(1)
      parms(5) = el_alpha
      parms(6) = sz_x
      parms(7) = sz_y
      call set_model_param(par,temp_model)
      dummy = create_and_calculate(exp_data,temp_model)

! Calculate likelihood cube

! Loop over beta
      do i = 1, cube_s
         write(*,*) 'Calculating row ',i
         parms(2) = b_arr(i)

! Loop over theta
         do j = 1, cube_s
            parms(3) = th_arr(j)
            parms(4) = th_arr(j)*ellip

! Calculate for first n0 value
            parms(1) = n_arr(1)
            call set_model_param(par,model_stuff)
            likeli(i,j,1) = create_and_calculate(exp_data,model_stuff)

! Find mean SZ Real for this model
            call make_aperture
            mean_re = 0.0
            psi = 0.0
            inc = pi/npoints
            do k = 1, npoints
               u = rad*sin(psi)
               v = rad*cos(psi)
               call extract_visibility&
                    & (re,im,maxsize,u,v,cellsize,re_pt,im_pt)
               mean_re = mean_re+re_pt
               psi = psi+inc
            end do
            mean_re = mean_re/npoints

            max_like = likeli(i,j,1)
            n_max = n_arr(1)

! Loop over n0
            do k = 2, cube_3

               parms(1) = n_arr(k)
               temp_model%m_data%values = &
                   (model_stuff%m_data%values-back)*&
                   (n_arr(k)/n_arr(1))**2+back

               call calculate_likely(exp_data,temp_model)

               likeli(i,j,k) = temp_model%likely
               if (likeli(i,j,k).gt.max_like) then
                  max_like = likeli(i,j,k)
                  n_max = n_arr(k)
               end if

!               write(*,*) parms(2),parms(3),parms(1),likeli(i,j,k)

            end do
            sz_flux(i,j) = mean_re*n_max/n_arr(1)

!            do k = 1,cube_s
!               parms(2) = b_arr(i)
!               parms(3) = th_arr(j)
!               parms(4) = th_arr(j)*ellip
!               parms(1) = n_arr(k)
!               parms(5) = el_alpha
!               parms(6) = sz_x
!               parms(7) = sz_y
!               call set_model_param(par,model_stuff)

!               dummy = create_and_calculate(exp_data,model_stuff)
!               write(*,*) likeli(i,j,k),dummy
!            end do
         end do
      end do

! Convert these log probabilities to likelihood
      probmax = maxval(likeli)
      likeli = likeli-probmax
      likeli = exp(likeli)

! Plotting loop
      plot_loop : do
         
         call io_getc('Command ','p',chr1,status)
         
         if (chr1.eq.'q') exit plot_loop
         
         if (chr1.eq.'p') then

! Do marginalisation
            prob2 = sum(likeli,dim=3)
            write(*,*) 'Marginalised'

! Normalising to the maximum likelihood
            probmax = maxval(prob2)
            if (probmax.ne.0.0) prob2 = prob2/probmax

            write(*,*) 'plotting'
            call pgpage
            call display_map(prob2,cube_s,graph_adj,&
                             b_inc,th_inc,b_cent,th_cent)
            call superpose_cntr(sz_flux,cube_s,&
                                b_inc,th_inc,b_cent,th_cent)
         end if
      end do plot_loop

      deallocate(likeli)
      deallocate(parms)
      deallocate(n_arr)
      deallocate(b_arr)
      deallocate(th_arr)
      deallocate(prob2)      
      deallocate(exp_data%values)
      deallocate(model_data%values)
      deallocate(temp_data%values)
      deallocate(sz_flux)

   end subroutine plot_likelihood

! ****************************************************************************

end module likelihood
