module batchbayesian

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

   public :: multi_bayesian_x_fit

!! Guess what... its a hacked around procedure for doing multiple bayesian fits
!!  to a particular image at different count levels. Won't do anything too
!!  different to the bayesian_x_fit routine....
!!
!!  Hacked around by Will Grainger (surprise!)
!!
!! Will _only_ do elliptical fits, with beta varying. 
!! Doesn't force a map display
!! REQUIRES the fitting to be silent!

contains

! *************************************************************************

  subroutine multi_bayesian_x_fit

    integer :: num_unknowns,num_iter1, num_iter2, num_iter
    real(kind=dp), dimension(:), allocatable :: parms,simplex_y
    character*1 :: chr1
    integer :: status,i,j
    type(data) :: exp_data
    type(data), target :: model_data
    type(param) :: par
    type(model) :: model_stuff
    real(kind=dp) :: temp
    real(kind=dp), dimension(:,:), allocatable :: simplex_p
    real(kind=dp) :: countlevel, countmax, countmin

    real(kind=dp), dimension(:), allocatable :: countlevels
    integer noofruns, loopcounter

    ! Ask if want to see the map first....
    call io_getc('Display the Rosat map (y/n):','y',chr1,status)
    if ((chr1=='y').or.(chr1=='Y')) then
       call display_map(rosat_map,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
    end if

    ! Ask for background level to use
    call io_getd('Background level:','*',back,status)

    ! define the area by counts.....
    countmax = maxval(rosat_map)
    countmin = minval(rosat_map)
    write (*,*) ' Range of counts in map is from',countmin,' to ',countmax 

    call io_geti('Number of fits to do:','10',noofruns,status)

    allocate(countlevels(noofruns))

    do i=1,noofruns
       call io_getd('Count level:','*',countlevels(i),status)
    end do

    do loopcounter=1,noofruns
       ! Switch off messages writen while making skies
       verbose = .false.
       ! Ensure that model sky is always recalculated from scratch
       overwrite = .true.

       ! Switch off plots
       do_plot = .false.

       ! We _shall vary beta
       vary_beta = .true.

       do i=1, maxsize
          do j=1, maxsize 
             if (rosat_map(i,j) < countlevels(loopcounter)) then
                ok_data(i,j) = .false.
             else 
                ok_data(i,j) = .true.
             end if
          end do
       end do

       ! Count up number of good pixels to be used in fitting
       num_dat_pt = count(ok_data)
       write(*,*) 'Using ',num_dat_pt,' pixels in fitting'

       ! Allocate space to arrays
       allocate(exp_data%values(num_dat_pt))
       allocate(model_data%values(num_dat_pt))
       model_stuff%m_data=>model_data
       call create_data(maxsize,maxsize,rosat_map,exp_data)

       ! Do not do adaptive smoothing on the model x-ray map
       smooth_x_model = .false.

       ! Call fitting routine dependent upon what type of model is being fitted to
       num_unknowns = 7

       ! Allocate array for the parameters
       allocate (parms(num_unknowns))

       ! Define initial guess
       parms(1) = n0
       parms(2) = beta   
       parms(3) = theta_c(1)
       parms(4) = theta_c(2)
       parms(5) = el_alpha
       parms(6) = sz_x
       parms(7) = sz_y

       call create_param(num_unknowns,parms,par)
       call set_model_param(par,model_stuff)
       allocate(simplex_p(num_param+1,num_param))
       allocate(simplex_y(num_param+1))

       call setup_simplex(exp_data,model_stuff,simplex_p,simplex_y)
       call amoeba(exp_data,model_stuff,simplex_p,simplex_y,num_iter1)   ! in minimum.f90
       call amoeba(exp_data,model_stuff,simplex_p,simplex_y,num_iter2)

       ! Redefine parameters as the best fit values

       parms=simplex_p(1,:)

       ! new parameter best values
       n0 = parms(1)
       beta = parms(2)      
       theta_c(1) = parms(3)
       theta_c(2) = parms(4) 
       el_alpha = parms(5)
       sz_x = parms(6)
       sz_y = parms(7)

       num_iter=num_iter1+num_iter2

       write(*,*) 'Found fit to rosat map after ',num_iter, ' iterations'
       write(*,*)
       write(*,*) 'Fit down to count   ', countlevels(loopcounter)
       write(*,*) 'central density     ', n0
       write(*,*) 'beta                ', beta
       write(*,*) 'theta 1             ', theta_c(1)
       write(*,*) 'theta 2             ', theta_c(2)
       write(*,*) 'inclination angle   ', el_alpha/deg2rad
       temp = modulo((el_alpha/deg2rad),360.0d0)
       write(*,*) 'incln angle, mod360 ', temp
       write(*,*) 'Offset in RA        ', sz_x
       write(*,*) 'Offset in Dec       ', sz_y
       theta_c(3) = sqrt(abs(theta_c(1)*theta_c(2)))
       write(*,*) 'set theta 3 to GM   ', theta_c(3)
       ang(1) = el_alpha

       ! Switch back on messages written while making skies
       verbose = .true.

       deallocate(parms)
       deallocate(simplex_p)
       deallocate(simplex_y)
       deallocate(exp_data%values)
       deallocate(model_data%values)

    end do

  end subroutine multi_bayesian_x_fit

! *************************************************************************
end module batchbayesian



