module mcg_xray

! Module for Method of Conjugate Gradients

   use kind_def

   implicit none  

! Globals which specify the required accuraccy of the fit
   integer, parameter :: max_steps = 50
   real(kind=dp), parameter :: accuracy = 0.01D0

! Global which specifies by what factor a parameter can change in one iteration
   real(kind=dp), parameter :: max_chng = 0.25D0

! Global which determines the size of the small step in the determination of
! the derivative
   real(kind=dp), parameter :: eps = 0.0001D0

   interface x_root_mcg
      module procedure x_rootfinder
   end interface

contains

!*************************************************************************

   subroutine x_rootfinder(num_observables_x,num_observables_y,rec_data,&
        & ok_data,nat_wt,parms,large_chng,num_steps,num_unknowns)
! Subroutine which calls conjgr and differences and stops fitting when accuracy
! is sufficient

      integer, intent(in) :: num_observables_x, num_observables_y, num_unknowns
      real(kind=dp), intent(in), &
           dimension(num_observables_x,num_observables_y) :: rec_data, nat_wt
      logical, intent(in), &
           dimension(num_observables_x,num_observables_y) :: ok_data
      real(kind=dp), intent(inout), dimension(num_unknowns) :: parms
      logical, intent(in), dimension(num_unknowns) :: large_chng
      real(kind=dp), &
           dimension(num_observables_x,num_observables_y,num_unknowns) :: der
      real(kind=dp), &
           dimension(num_observables_x,num_observables_y) :: diff  
      real(kind=dp), dimension(num_unknowns) :: dpm, chng
      integer, intent(out) :: num_steps
      integer :: i
      real(kind=dp) :: step_length, chng_factor

      num_steps = 0

! Main loop in which differences (and derivatives) between model and data 
! are calculated and then model parameters are adjusted accordingly
      finder: do
         num_steps = num_steps+1

! Find difference between model and data, and local derivatives 
         call differences(num_observables_x,num_observables_y,rec_data,&
        & ok_data,parms,diff,der,nat_wt,num_unknowns)

!Check whether model fits data exactly ie all components of diff are 0
         if (sum(abs(diff))==0.0) then
            write(*,*) 'Initial guess correct'
            exit finder
         end if

! Find optimum change in model parameters
! Does this by calculating dpm assuming that
! (Model)(Model parms) + (Derivative of model)(dpm) = (recorded data)
         call conjgr(num_observables_x,num_observables_y,diff,der,dpm,&
        & num_unknowns)

! Update model parameters checking first that none of the changes are
! too much greater than the parameters as they are now
         where (parms.ne.0.0)
            chng = dpm/parms
         elsewhere
            chng = 0.0
         end where
         if (any((abs(chng).gt.max_chng).and.(.not.large_chng))) then
            write(*,*) 'Model far from optimum - stepping slowly'
            chng_factor = max_chng/maxval(abs(chng),mask=(.not.large_chng))
            parms = parms+dpm*chng_factor
         else
            parms = parms+dpm
         end if

! Step length is the sum of the absolute values of the changes in the 
! parameters 
         step_length = sum(abs(dpm))

! Check for completion or if no solution has yet been found in the maximum
! number of iterations
         if ((num_steps==max_steps).or.(step_length<accuracy)) exit finder

! Display what the current values of the parameters are
         write(*,*) '-----------------------------------------------------'
         write(*,*) 'Iteration ',num_steps, ' step length ',step_length
         write(*,*)
         do i = 1, num_unknowns
            write(*,*) parms(i)
         end do
         write(*,*) '-----------------------------------------------------'
      end do finder
      call differences(num_observables_x,num_observables_y,rec_data,&
        & ok_data,parms,diff,der,nat_wt,num_unknowns)

      if (step_length>1.0) then
         write(*,*) 'Did not converge'
      end if

   end subroutine x_rootfinder

!*************************************************************************

   subroutine differences(num_observables_x,num_observables_y,rec_data,&
        & ok_data,parms,diff,der,nat_wt,num_unknowns)

! Returns the differences, diff, between measured values, rec_data, 
! and estimated values, est_data, for the parameter set parms. Also
! calculates the differentials in the array der at this point.

      implicit none

      integer, intent(in) :: num_observables_x, num_observables_y, num_unknowns
      real(kind=dp), intent(in),&
           dimension(num_observables_x,num_observables_y) :: rec_data,nat_wt
      logical, intent(in), &
           dimension(num_observables_x,num_observables_y) :: ok_data
      real(kind=dp), intent(inout), dimension(num_unknowns) :: parms
      real(kind=dp), dimension(num_observables_x,num_observables_y) :: &
           est_data, dummy_est_data
      real(kind=dp), intent(out), &
           dimension(num_observables_x,num_observables_y) :: diff  
      real(kind=dp), intent(out), &
           dimension(num_observables_x,num_observables_y,num_unknowns) :: der
      real(kind=dp), dimension(num_unknowns) :: dummy_parms
      real(kind=dp) :: a, delta_dummy_parms_j
      integer :: i, j, k

      diff = 0.0
      der = 0.0

! Generate model data from the parameters and find weighted difference from 
! recorded  data
      call model_fn(parms,num_unknowns,est_data)

! Only use unflagged observables
      where (ok_data)     
         diff = (rec_data-est_data)*nat_wt
      elsewhere 
         diff = 0.0
      end where

! Find derivative of the model wrt each of the parameters at this point by
! calculating the model at a small distance (1+eps) times the parameter in
! question, subtracting the value at the point, and dividing by the small 
! distance. The derivatives are then also weighted.
      do j = 1, num_unknowns
         dummy_parms = parms
         if (dummy_parms(j)==0.0) then
            delta_dummy_parms_j = eps
         else
            delta_dummy_parms_j = dummy_parms(j)*eps
         end if
         dummy_parms(j) = dummy_parms(j)+delta_dummy_parms_j
         call model_fn(dummy_parms,num_unknowns,dummy_est_data)

! Loop over all observables
         do i = 1, num_observables_x
            do k = 1,num_observables_y

! Only use unflagged observables
               if (ok_data(i,k)) then     
                  der(i,k,j) = ((dummy_est_data(i,k)-est_data(i,k))&
                       *nat_wt(i,k))/delta_dummy_parms_j
               end if
            end do
         end do
      end do

   end subroutine differences

!*************************************************************************

   subroutine conjgr(num_observables_x,num_observables_y,diff,der,dpm,&
        & num_unknowns)

!Solves for dpm=diff.INV(der) by method of conjugate gradients

      integer, intent(in) :: num_observables_x, num_observables_y, num_unknowns
      real(kind=dp), dimension(num_unknowns) :: lws, tws, sws
      real(kind=dp), dimension(num_unknowns,num_unknowns) :: mws
      real(kind=dp) :: u2, w2, u, w
      logical :: error_stat
      real(kind=dp), intent(in),&
           dimension(num_observables_x,num_observables_y) :: diff
      real(kind=dp), intent(out), dimension(num_unknowns) :: dpm
      real(kind=dp), intent(in), &
           dimension(num_observables_x,num_observables_y,num_unknowns) :: der
      integer :: tot_observables, ji

      tot_observables = num_observables_x*num_observables_y

! Set error check
      error_stat = .false. 
      dpm = 0.0


      lws = matmul(transpose(&
           reshape(der,shape=(/(tot_observables),num_unknowns/))),&
           reshape(diff,shape=(/(tot_observables)/)))

      mws = matmul(transpose(&
           reshape(der,shape=(/(tot_observables),num_unknowns/))),&
           reshape(der,shape=(/(tot_observables),num_unknowns/)))

      sws = lws

      loop_2 : do ji = 1, num_unknowns
         if (error_stat) exit loop_2

         tws = matmul(mws,sws)
         u = sum(lws*sws)
         w = sum(tws*sws)
         if (w==0.0) then
            error_stat = .true.
            write(*,*) 'w = 0 in conjgr'
            exit loop_2
         end if
         u = u/w
         lws = lws-u*tws
         dpm = dpm+sws*u

         u2 = sum(lws*tws)
         w2 = sum(sws*tws)
         if (w2==0.0) then
            error_stat = .true.
            write(*,*) 'w2 = 0 in conjgr'
            exit loop_2
         end if
         u2 = -u2/w2
         sws = lws+u2*sws

      end do loop_2

! Check if things went horribly wrong
      if (error_stat) then
         dpm = 0.0
         write(*,*) 'Error in conjgr'
      end if

   end subroutine conjgr

!*************************************************************************

end module mcg_xray
