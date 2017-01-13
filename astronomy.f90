module astronomy

   use kind_def
   use sz_globals

   implicit none

   public :: calc_sepn
   public :: sin_proj_forward
   public :: sin_proj_back

contains

!**************************************************************************

    subroutine calc_sepn(ra1,dec1,ra2,dec2,sepn)

!
!  Finds the separation between two points on the sky - answer returned 
!  in radians
!
!  History: 22/8/00 - original version [KG]
!           5/12/00 - It used to return NaN, when it should 0 (occasionally)
!                      Added line 7 [AS]

      use kind_def

      implicit none

      real(kind=dp), intent(in) ::  ra1, ra2, dec1, dec2
      real(kind=dp), intent(out) :: sepn
      real(kind=dp) :: cos_sepn

      cos_sepn = sin(dec1)*sin(dec2)+cos(dec1)*cos(dec2)*cos(ra1-ra2)
      if (cos_sepn.gt.1.00000000) cos_sepn=1.0 
      sepn = acos(cos_sepn)

   end subroutine calc_sepn

! **************************************************************************

   subroutine sin_proj_forward(ra1,dec1,ra2,dec2,insize,dx,dy)
!  
!  Calculates the orthographic or sin projection where the celestial sphere
!  is projected onto a tangent plane and the centre of the projection is 
!  infinitely distant from the tangent plane.
!  ra1, dec1 is the position on the sky of interest (in radians)
!  ra2, dec2 is the position of the map centre (in radians)
!  insize is the cellsize in arcseconds
!  dx and dy are the offsets in pixels of the sky position from the map centre

      real(kind=dp), intent(in) :: ra1, dec1, ra2, dec2, insize
      real(kind=dp), intent(out) :: dx, dy

      dx = -1./(insize*sec2rad)*cos(dec1)*sin(ra1-ra2)
      dy = 1./(insize*sec2rad)*(cos(dec2)*sin(dec1)-&
           sin(dec2)*cos(dec1)*cos(ra1-ra2))

   end subroutine sin_proj_forward

! **************************************************************************

   subroutine sin_proj_back(ra1,dec1,ra2,dec2,insize,dx,dy)
!  
!  Calculates the orthographic or sin projection where the celestial sphere
!  is projected onto a tangent plane and the centre of the projection is 
!  infinitely distant from the tangent plane.
!  ra1, dec1 is the position on the sky of interest (in radians)
!  ra2, dec2 is the position of the map centre (in radians)
!  insize is the cellsize in arcseconds
!  dx and dy are the offsets in pixels of the sky position from the map centre

      real(kind=dp), intent(in) :: ra2, dec2, insize, dx, dy
      real(kind=dp), intent(out) :: ra1, dec1
      real(kind=dp) :: DD, dx_r, dy_r

      dx_r = dx*insize*sec2rad
      dy_r = dy*insize*sec2rad
      DD = asin(sqrt(dx_r**2+dy_r**2))

      if ((dy_r.eq.0.0).and.(dx_r.eq.0.0d0)) then
         ra1 = ra2
         dec1 = dec2
      else
         ra1 = ra2-atan(dx_r*sin(DD)/&
            (sqrt(dx_r**2+dy_r**2)*cos(dec2)*cos(DD)-dy_r*sin(dec2)*sin(DD)))
         dec1 = asin(cos(DD)*sin(dec2)+(dy_r*sin(DD)*cos(dec2))/&
             (sqrt(dx_r**2+dy_r**2)))
      end if

   end subroutine sin_proj_back

! **************************************************************************

end module astronomy
