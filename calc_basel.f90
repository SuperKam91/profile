subroutine calc_basel

   use kind_def
   use sz_globals

   implicit none

   integer :: i, j
   real(kind=dp) :: obslambda

   obslambda = const_c/obsfreq

   do i = 1, n_antennas
      do j = 1, n_antennas
         basel_x(i,j) = (x_pos(i)-x_pos(j))/obslambda
         basel_y(i,j) = (y_pos(i)-y_pos(j))/obslambda
         basel_z(i,j) = (z_pos(i)-z_pos(j))/obslambda
      end do
   end do

end subroutine calc_basel
