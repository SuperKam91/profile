! Calculates u, v, w coordinates for a given HA and Dec at the central 
! telescope frequency obsfreq

subroutine calc_uvw(ha,dec,uu,vv,ww,shadow)

   use kind_def
   use sz_globals

   implicit none

   integer :: i, j, l
   real(kind=dp), dimension(n_antennas,n_antennas) :: uu,vv,ww
   integer, dimension(n_antennas,n_antennas) :: shadow
   real(kind=dp), intent(in) :: ha, dec

   do i = 1, n_antennas
      do j = 1, n_antennas
         if (i.ne.j) then
            uu(i,j) = basel_x(i,j)*cos(ha)+&
                      basel_y(i,j)*(-sin(ha)*sin(tel_lat))+&
                      basel_z(i,j)*(-sin(ha)*cos(tel_lat))

            vv(i,j) = basel_x(i,j)*sin(ha)*sin(dec)+&
                      basel_y(i,j)*(cos(dec)*cos(tel_lat)&
                                   +cos(ha)*sin(dec)*sin(tel_lat))+&
                      basel_z(i,j)*(-cos(dec)*sin(tel_lat)&
                                   +cos(ha)*sin(dec)*cos(tel_lat))

            ww(i,j) = basel_x(i,j)*(-sin(ha)*cos(dec))+&
                      basel_y(i,j)*(sin(dec)*cos(tel_lat)&
                                   -cos(ha)*cos(dec)*sin(tel_lat))+&
                      basel_z(i,j)*(-sin(dec)*sin(tel_lat)&
                                   -cos(ha)*cos(dec)*cos(tel_lat))
         else
            uu(i,j) = 0.0
            vv(i,j) = 0.0
            ww(i,j) = 0.0
            
         end if
      end do
   end do

! Check for shadowing
   shadow = 1
   do i = 1, n_antennas-1
      do j = i+1, n_antennas
         if (sqrt(uu(i,j)**2+vv(i,j)**2).lt.dish_diameter*obsfreq/const_c) then
            if (ww(i,j).lt.0.0) then
               do l = 1, i-1
                  shadow(l,i) = i+1
                  shadow(i,l) = i+1
               end do
               do l = i+1, n_antennas
                  shadow(l,i) = i+1
                  shadow(i,l) = i+1
               end do
            else
               do l = 1, j-1
                  shadow(l,j) = j+1
                  shadow(j,l) = j+1
               end do
               do l = j+1, n_antennas
                  shadow(l,j) = j+1
                  shadow(j,l) = j+1
               end do

            end if
         end if
      end do
   end do

! Ensure that shadow does not exceed 12
   do i = 1, n_antennas
      do j = 1, n_antennas
         shadow_alloc: do
            if (shadow(i,j).gt.12) then
               shadow(i,j) = shadow(i,j)-10
            else
               exit shadow_alloc
            end if
         end do shadow_alloc
      end do
   end do

end subroutine calc_uvw


