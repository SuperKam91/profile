! Routine which produces 2D T model with a halo
! Moved from make_skies 28/08/09 [KG] 
subroutine make_model_temp_halo(T)

   use kind_def
   use sz_globals
   use kgraph

   implicit none

   integer :: i,j, i1, j1, which_i
   real(kind=dp), dimension(maxsize,maxsize) :: T
   real(kind=dp) :: rad

! If Using Halo model then has T_central until Tcorerad then goes linearly 
! to T_halo at 1.5*Tcorerad
   do j = 1, maxsize
      do i = 1, maxsize
         i1 = i-1
         j1 = j-1
         rad = sqrt(float(i1**2+j1**2))*cellsize
         if   (rad.lt.Tcorerad) then
            T(i,j) = 1.
         else if (rad.lt.Tcorerad*1.5) then
            T(i,j) = (T_central + (T_halo-T_central)*&
                 (rad/Tcorerad -1.0)*2.0)/T_central
         else
            T(i,j) = T_halo/T_central
         end if
      end do
   end do

   if (do_plot) then
      call display_map(T,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
   end if

end subroutine make_model_temp_halo
