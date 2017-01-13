subroutine model_fn(parms,num_unknowns,est_data)

   use kind_def
   use sz_globals
   use map_handling

   implicit none
   integer, intent(in) :: num_unknowns
   real(kind=dp), intent(in), dimension(num_unknowns) :: parms
   real(kind=dp), dimension(maxsize,maxsize), intent(out) :: est_data
   real(kind=dp), dimension(:,:), allocatable :: smooth_map

   if (ellipsoid) then

! new parameter setup for elliptical model
      if (vary_beta) then

! Allowing beta to vary
         if (num_unknowns.ne.7) then
            write(*,*) 'num_unknowns incorrect in model_fn', num_unknowns
         end if
         n0 = parms(1)
         beta = parms(2)      
         theta_c(1) = parms(3)
         theta_c(2) = parms(4)
         theta_c(3) = sqrt(theta_c(1)*theta_c(2))
         el_alpha = parms(5)
         if (vary_posn) then
            sz_x = parms(6)
            sz_y = parms(7)
         end if
      else

! Not allowing beta to vary
         if (num_unknowns.ne.6) then
            write(*,*) 'num_unknowns incorrect in model_fn', num_unknowns
         end if
         n0 = parms(1)
         theta_c(1) = parms(2)
         theta_c(2) = parms(3)
         theta_c(3) = sqrt(theta_c(1)*theta_c(2))
         el_alpha = parms(4)
         if (vary_posn) then
            sz_x = parms(5)
            sz_y = parms(6)
         end if
      end if
   else

! new parameter setup for spherical model
      if (vary_beta) then

! Allowing beta to vary
         if (num_unknowns.ne.5) then
            write(*,*) 'num_unknowns incorrect in model_fn'
         end if
         n0 = parms(1)
         beta = parms(2)
         theta_cx = parms(3)
         if (vary_posn) then
            sz_x = parms(4)
            sz_y = parms(5)
         end if
      else

! Not allowing beta to vary
         n0 = parms(1)
         theta_cx = parms(2)
         if (vary_posn) then
            sz_x = parms(3)
            sz_y = parms(4)
         end if
      end if
   end if

! Use these parameters to make an x-ray sky and observe with telescope
   call make_cluster
   call make_model_xray_map

! Smooth model x-ray map if necessary
   if (smooth_x_model) then

! Allocate space to smooth map
      allocate(smooth_map(maxsize,maxsize))

! Produce smoothed model map
      call adaptive_smooth(xmap,smooth_map,maxsize)     

! Copy smoothed map to model map
      xmap = smooth_map

! Deallocate space for smoothed map
      deallocate(smooth_map)
   end if

! Write model x-ray map to the estimated data
   est_data = xmap

end subroutine model_fn


