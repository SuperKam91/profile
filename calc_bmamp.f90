
! routine to find the maginatude of a beam at an offset position
! to be used in timestream simulations of point sources [AMS 22/2/08]

subroutine calc_bamp(ha_cen,ra_cen,dec_cen,ra,dec,bmag)

   use kind_def
   use sz_globals

   implicit none

   integer   :: i,j,k
   real(dp)  :: ha_cen, ra_cen, dec_cen
   real(dp)  :: ra,dec  
   real(dp), intent(out) :: bmag
   real(dp)  :: ha_new,LST,azim,elev,c_azim,c_elev
   real(dp)  :: D_az, D_el,az_inc,el_inc


!  Convert rads to hours:
   
   ha_cen = ha_cen/hr2rad
   ra_cen = ra_cen/hr2rad
   LST = ha_cen+ra_cen
   
   ha_new = LST - ra

!  Convert hours to rad:

   ha_cen = ha_cen*hr2rad
   ha_new = ha_new*hr2rad

!  Find required azimuth and elevation:

   call sla_de2h(ha_cen, dec_cen, tel_lat, c_azim, c_elev)
   call sla_de2h(ha_new, dec, tel_lat, azim, elev)

!  Return azimuth/elevation in the range -pi/2 to pi/2:

   if ((azim.gt.pi/2.).and.(azim.le.pi)) then
      azim = pi-azim
      elev = -1.*(pi/2. - elev)
   elseif (azim.lt.-1.*pi/2.) then
      azim = -1.*(pi+azim)
      elev = -1.*(pi/2. - elev)
   end if
   if ((c_azim.gt.pi/2.).and.(c_azim.le.pi)) then
      c_azim = pi-c_azim
      c_elev = -1.*(pi/2. - c_elev)
   elseif (c_azim.lt.-1.*pi/2.) then
      c_azim = -1.*(pi+c_azim)
      c_elev = -1.*(pi/2. - c_elev)
   end if
   
!  Find difference between centre and offset source [radians]:

   D_az = c_azim - azim
   D_el = c_elev - elev

! change to projected co-ordinates:

   D_az = (D_az + pi/2.)/az_inc
   D_el = (D_el + pi/2.)/el_inc
   i = int(D_az+sign(0.5d0,D_az))
   j = int(D_el+sign(0.5d0,D_el))

! find beam magnitude at offset point:

   bmag = azelbeam(i,j)
   

end subroutine calc_bamp
