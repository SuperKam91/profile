
! *** HACKED 24/03/04 AMS *** 


module bayesian

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

   public :: bayesian_x_fit

contains

! *************************************************************************

   subroutine bayesian_x_fit

      integer :: num_unknowns,num_iter1, num_iter2, num_iter
      real(kind=dp), dimension(:), allocatable :: parms,simplex_y
      character*1 :: chr1
      integer :: status,i,j,k
      type(data) :: exp_data
      type(data), target :: model_data
      type(param) :: par
      type(model) :: model_stuff
      real(kind=dp) :: temp
      real(kind=dp), dimension(:,:), allocatable :: simplex_p
      real(kind=dp) :: countlevel, countmax, countmin

! for anna hack
      real :: b_min, b_max, th_min, th_max, b_cent, th_cent 
!      real :: b_min_s, b_max_s, th_min_s, th_max_s
      real(kind = dp) :: th_inc, b_inc
      integer :: cube_s, just, axis, da, numcont, grid
      real(kind=dp), dimension(:), allocatable :: b_arr, th_arr
      real(kind=dp), dimension(:,:), allocatable :: w, temppost, postg
      real, dimension(:,:), allocatable :: w_s
      real(kind=dp), dimension(num_dat_pt) :: prob
      real(kind=dp) :: dummy
      integer :: allo_stat
      type(data) :: data_stuff
      real :: peakpost2, xmin, xmax, ymin, ymax
      real, dimension(6) :: tr
      character*80 :: plot_device
      integer :: pgopen
      external :: pgopen
      logical flag
      real(kind=dp) :: ellip, delta, totalvol, max, clev, bit
      real(kind=dp) :: likely, fraction, volume, min_g68, min_g95, min_g99
      real(kind=dp), dimension(:), allocatable :: cg


! Switch off messages writen while making skies
      verbose = .false.

! Ensure that model sky is always recalculated from scratch
      overwrite = .true.

! Switch off plots
      do_plot = .false.

! A warning about a known problem
      write(*,*) 'NB : ROSAT map must have actual photon counts rather than',&
                 'count rates because we are using Poisson statistics'
      write(*,*)

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


! Ask for background level to use
      call io_getd('Background level:','*',back,status)

      loud_fitting = .false.

      call io_getc('Loud fitting (for debugging)','y',chr1,status)
      if ((chr1=='y').or.(chr1=='Y')) then
         loud_fitting = .true.
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

! Count up number of good pixels to be used in fitting
      num_dat_pt = count(ok_data)
      write(*,*) 'Using ',num_dat_pt,' pixels in fitting'

! Allocate space to arrays
      allocate(exp_data%values(num_dat_pt),STAT=allo_stat)
      if (allo_stat.ne.0) then 
         write(*,*) 'Error allocating exp_data'
      else
         write(*,*) 'No problem allocating exp_data'
      end if
      allocate(model_data%values(num_dat_pt),STAT=allo_stat)
      if (allo_stat.ne.0) then 
         write(*,*) 'Error allocating model_data'
      else
         write(*,*) 'No problem allocating model_data'
      end if
      model_stuff%m_data=>model_data
      call create_data(maxsize,maxsize,rosat_map,exp_data)

! Do not do adaptive smoothing on the model x-ray map
      smooth_x_model = .false.

! *****************************************************************************
! *****************************************************************************
! HACK:
      if (ellipsoid) then

         num_unknowns = 7
         ellip = theta_c(1)/theta_c(2)

! how many points on likelihood plot
         call io_geti('size of likelihood grid?','25',cube_s,status)

         allocate(b_arr(cube_s))
         allocate(th_arr(cube_s))
         allocate(w(cube_s,cube_s))
         allocate(temppost(cube_s,cube_s))
         allocate(postg(cube_s,cube_s))
         allocate(w_s(cube_s,cube_s))
         allocate(parms(num_unknowns))

! get initial values
         call io_getr('Minimum value of beta:','0.6',b_min,status)
         call io_getr('Maximum value of beta:','0.7',b_max,status)
         call io_getr('Minimum value of theta:','50.0',th_min,status)
         call io_getr('Maximum value of theta:','70.',th_max,status)

! define increments
         b_inc = (b_max-b_min)/(cube_s-1)  
         th_inc = (th_max-th_min)/(cube_s-1) 
! DEBUG
!write(*,*) 'found increments'

         write(*,*) 'plotting likelihood...this might take a while.'

         tr(1) = 0.0
         tr(2) = 1.0
         tr(3) = 0.0
         tr(4) = 0.0
         tr(5) = 0.0
         tr(6) = 1.0

         call pgask(flag)
write(*,*) 'pgask'

         xmin = 0.5
         xmax = 0.8
         ymin = 30.0
         ymax = 80.0
         just = 0
         axis = 1
         
         write(*,*) xmin,xmax,ymin,ymax,just,axis

         call pgenv(xmin,xmax,ymin,ymax,just,axis)
write(*,*) 'pgenv'
         call pglab('(beta)','(theta)',&
              'posterior distribution of beta and theta')
write(*,*) 'pglabel'



!goto 10

! loop over beta        
         do i = 1, cube_s
            b_arr(i) = b_min+(i-1)*b_inc
! loop over theta            
            do j = 1, cube_s
               th_arr(j) = th_min+(j-1)*th_inc
               
!               call create_param(num_unknowns,parms,par)

               parms(1) = n0
               parms(2) = b_arr(i)             
               parms(3) = th_arr(j)              
               parms(4) = th_arr(j)*ellip           
               parms(5) = el_alpha
               parms(6) = sz_x              
               parms(7) = sz_y
! DEBUG
               write(*,*) 'theta:', th_arr(j)
               write(*,*) 'put in parms'

               call create_param(num_unknowns,parms,par)         
! DEBUG
               write(*,*) 'created parameters'

               call set_model_param(par,model_stuff)
! DEBUG
               write(*,*) 'set parameters'

               call create_model(model_stuff)
! DEBUG
               write(*,*) 'created model'
!               write(*,*) model_stuff%m_data%values(num_dat_pt)

               do k=1,num_dat_pt

! Don't pay too much attention to points where the model predicts virtually 
! no signal. NOT TOO SURE THAT THIS IS A VALID THING TO DO. 
                  if(model_stuff%m_data%values(k).le.0.1)then                     
                     prob(k)=0.0                     
                  else
!        prob(i) = given%values(i)*log(calculated%m_data%values(i))&
!             & -1.0*calculated%m_data%values(i)-&
!             & log(1.0*factorial(given%values(i)))
                  dummy=1.0*nint(exp_data%values(k))
! DEBUG
!                  write(*,*) 'done dummy'

                     prob(k) = exp_data%values(k)*log(model_stuff%m_data%values(k))&
                          & -1.0*model_stuff%m_data%values(k)-&
                          & log_factorial(dummy)
                  end if
               end do
!               model_stuff%likely=sum(prob)
! DEBUG
write(*,*) 'worked out log probabilities'

               w(i,j) = sum(prob)
               write(*,*) i,j,w(i,j)
            end do
         end do

10 write(*,*)'jump'

         do j = 1,cube_s
            do k = 1,cube_s
               w(j,k) = exp(w(j,k))
            end do
         end do

! DEBUG
write(*,*) 'exponentiated'         

         w(1,1) = 1.5

         peakpost2 = maxval(w)
         w_s =sngl(w)

         write(*,*) 'done peakpost'

         b_cent = (b_max - b_min)/2
         th_cent = (th_max - th_min)/2
write(*,*) 'debug1'
!         call display_map(w, cube_s, graph_adj, b_inc, th_inc, b_cent, th_cent)

 write(*,*) 'pglabel'

 write(*,*) w_s

 write(*,*) tr


 write(*,*) 

         call pggray(w_s,cube_s,cube_s,1,cube_s,1,cube_s,peakpost2,0.0,tr)

write(*,*) 'pggray'
         call pgtext(0.,0.,'99,95 and 68 percent contours')
!         read(*,*)
!         call pgask(flag)
!         read(*,*)
        
         temppost = w
         delta = 0.0d0
         da = 500
         totalvol = sum(temppost)
         max = maxval(temppost)
         
         clev = 0.68
         bit = 0.08
         postg = 0.0d0
         grid = cube_s

         allocate(cg(3))

         do i = 1, da                 
            temppost = w
            do j = 1,grid
               do k = 1,grid              
                  if (temppost(j,k) < (max - delta*max)) then
                     temppost(j,k) = 0
                  endif
               end do
            end do
            volume = sum(temppost) 
            fraction = volume/totalvol 
            delta = delta + (1/dble(da))                        
            if ((0.675 < fraction) .AND. (fraction < 0.685)) then
               do j = 1, grid
                  do k = 1, grid
                     if ((temppost(j,k)/max > (1 - delta) - bit) .AND. (temppost(j,k)/max < (1 - delta) + bit)) then
!                        write(2,*) m_grid(j), c_grid(k)
                        min_g68 = 1 - delta
                     endif
                  enddo
               enddo
            endif
            if (fraction .GE. clev) exit 
         enddo
         
         write(*,*) 'confidence level = ', fraction
         
         clev = 0.95
         bit = 0.01
         do i = 1, da                 
            temppost = w
            do j = 1,grid
               do k = 1,grid              
                  if (temppost(j,k) < (max - delta*max)) then
                     temppost(j,k) = 0
                  endif
               end do
            end do
            volume = sum(temppost) 
            fraction = volume/totalvol 
            delta = delta + (1/dble(da))                        
            if ((0.945 < fraction) .AND. (fraction < 0.955)) then
               do j = 1, grid
                  do k = 1, grid
                     if ((temppost(j,k)/max > (1 - delta) - bit) .AND. (temppost(j,k)/max < (1 - delta) + bit)) then
!                        write(2,*) m_grid(j), c_grid(k)
                        min_g95 = 1 - delta
                     endif
                  enddo
               enddo
            endif
            if (fraction .GE. clev) exit 
         enddo
         
         write(*,*) 'confidence level = ', fraction
         
         
         clev = 0.99
         bit = 0.003
         do i = 1, da                 
            temppost = w
            do j = 1,grid
               do k = 1,grid              
                  if (temppost(j,k) < (max - delta*max)) then
                     temppost(j,k) = 0
                  endif
               end do
            end do
            volume = sum(temppost) 
            fraction = volume/totalvol 
            delta = delta + (1/dble(da))                        
            if ((0.985 < fraction) .AND. (fraction < 0.995)) then
               do j = 1, grid
                  do k = 1, grid
                     if ((temppost(j,k)/max > (1 - delta) - bit) .AND. (temppost(j,k)/max < (1 - delta) + bit)) then
!                        write(2,*) m_grid(j), c_grid(k)
                        min_g99 = 1 - delta
                     endif
                  enddo
               enddo
            endif
            if (fraction .GE. clev) exit 
         enddo
         
         write(*,*) 'confidence level = ', fraction
         write(*,*)
         
         numcont = 3
         postg = w
         cg(1) = min_g68*max
         cg(2) = min_g95*max
         cg(3) = min_g99*max
         
         call pgcont(postg,grid,grid,1,grid,1,grid,cg,numcont,tr)
         call pgask(flag)
         
         
      end if

! ****************************************************************************
! ****************************************************************************

      goto 20


! Call fitting routine dependent upon what type of model is being fitted to
      if (ellipsoid) then
         write(*,*) 'Fitting ellipsoidal model'
         if (vary_beta) then
            num_unknowns = 7

! Allocate array for the parameters
            allocate (parms(num_unknowns))


! Define initial guess
            parms(1) = n0         ! set this as fixed
            parms(2) = beta   
            parms(3) = theta_c(1)
            parms(4) = theta_c(2)
            parms(5) = el_alpha   ! fixed
            parms(6) = sz_x       ! fixed
            parms(7) = sz_y       ! fixed
         else
            num_unknowns = 6

! Allocate array for the parameters
            allocate (parms(num_unknowns))

! Define initial guess
            parms(1) = n0
            parms(2) = theta_c(1)
            parms(3) = theta_c(2)
            parms(4) = el_alpha
            parms(5) = sz_x
            parms(6) = sz_y
         end if

         call create_param(num_unknowns,parms,par)
!DEBUG
!       write (*,*) 'Created params'


         call set_model_param(par,model_stuff)
!DEBUG
!       write (*,*) 'set model parameters'
         allocate(simplex_p(num_param+1,num_param))
         allocate(simplex_y(num_param+1))

!DEBUG
!       write (*,*) 'allocated space'
         call setup_simplex(exp_data,model_stuff,simplex_p,simplex_y)
!DEBUG
!       write (*,*) 'set up simplex'
         call amoeba(exp_data,model_stuff,simplex_p,simplex_y,num_iter1)   ! in minimum.f90
!DEBUG
!       write (*,*) 'called amoeba'
         call amoeba(exp_data,model_stuff,simplex_p,simplex_y,num_iter2)

! Redefine parameters as the best fit values

         parms=simplex_p(1,:)

         if (vary_beta) then

! new parameter best values
            n0 = parms(1)
            beta = parms(2)      
            theta_c(1) = parms(3)
            theta_c(2) = parms(4) 
            el_alpha = parms(5)
            sz_x = parms(6)
            sz_y = parms(7)
         else

! new parameter best values
            n0 = parms(1)
            theta_c(1) = parms(2)
            theta_c(2) = parms(3) 
            el_alpha = parms(4)
            sz_x = parms(5)
            sz_y = parms(6)
         end if

         num_iter=num_iter1+num_iter2

         write(*,*) 'Found fit to rosat map after ',num_iter, ' iterations'
         write(*,*) 'central density     ', n0
         write(*,*) 'beta                ', beta
         write(*,*) 'theta 1             ', theta_c(1)
         write(*,*) 'theta 2             ', theta_c(2)
         write(*,*) 'inclination angle   ', el_alpha/deg2rad
         temp = modulo((el_alpha/deg2rad),360.0d0)
         write(*,*) 'incln angle, mod360 ', temp
         write(*,*) 'Offset in RA        ', sz_x
         write(*,*) 'Offset in Dec       ', sz_y
         theta_c(3) = sqrt(theta_c(1)*theta_c(2))
         write(*,*) 'set theta 3 to GM   ', theta_c(3)
         ang(1) = el_alpha
      else

! Do spherical model fit
         write(*,*) 'Fitting spherical model'
         if (vary_beta) then
            num_unknowns = 5

! Allocate arrays for parameter
            allocate (parms(num_unknowns))

! Define inital guess
            parms(1)=n0
            parms(2)=beta
            parms(3)=theta_cx
            parms(4)=sz_x
            parms(5)=sz_y

         else
            num_unknowns = 4

! Allocate arrays for parameter
            allocate (parms(num_unknowns))

! Define initial guess
            parms(1) = n0
            parms(2) = theta_cx
            parms(3) = sz_x
            parms(4) = sz_y
         end if
         call create_param(num_unknowns,parms,par)
         call set_model_param(par,model_stuff)
         allocate(simplex_p(num_param+1,num_param))
         allocate(simplex_y(num_param+1))
         call setup_simplex(exp_data,model_stuff,simplex_p,simplex_y)
         call amoeba(exp_data,model_stuff,simplex_p,simplex_y,num_iter1)
         call amoeba(exp_data,model_stuff,simplex_p,simplex_y,num_iter2)

         parms=simplex_p(1,:)

         if (vary_beta) then
! Redefine parameters as the best fit values
            n0 = parms(1)
            beta = parms(2)
            theta_cx = parms(3)
            sz_x = parms(4)
            sz_y = parms(5)
         else
! Redefine parameters as the best fit values
            n0 = parms(1)
            theta_cx = parms(2)
            sz_x = parms(3)
            sz_y = parms(4)
         end if

         num_iter=num_iter1+num_iter2
         write(*,*) 'Found fit to rosat map after ',num_iter, ' iterations'
         write(*,*) 'central density     ', n0
         write(*,*) 'beta                ', beta
         write(*,*) 'theta               ', theta_cx
         write(*,*) 'Offset in RA        ', sz_x
         write(*,*) 'Offset in Dec       ', sz_y
      end if

! Switch back on messages written while making skies
20      verbose = .true.

      deallocate(parms)
!      deallocate(simplex_p)
!      deallocate(simplex_y)
      deallocate(exp_data%values)
      deallocate(model_data%values)
      deallocate(b_arr)
      deallocate(th_arr)
      deallocate(w)
      deallocate(temppost)
      deallocate(postg)
      deallocate(w_s)

   end subroutine bayesian_x_fit

! *************************************************************************

end module bayesian



