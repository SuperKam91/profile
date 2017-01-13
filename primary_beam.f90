module primary_beam
   
   use kind_def
   use sz_globals
   use maths

   implicit none

   public :: gridded_pb
   public :: gaussian_pb
   public :: polynomial_pb
   public :: scaled_vla_pb
   public :: scaled_rt_pb
   public :: scaled_cat_pb
   public :: survey_pb
   public :: vsa_ea_a_pb
   public :: elliptical_pb
   public :: SPH10_beam
   public :: SPH11_beam
   public :: SPH22_beam
   public :: SPH31_beam
   public :: SPH52_beam
   public :: szi_beam
   public :: display_primary_beam
   public :: get_taper
   public :: unity_beam

   private :: ryle_pb
   private :: cat_pb
   private :: vla_pb
   private :: cat_poly

contains

! **************************************************************************

! returns primary beam (power) at given r (arcsec) from a gridded data set. 

subroutine gridded_pb(del_ra,del_dec,taper)

   use kind_def
   use sz_globals
   use maths

   implicit none
   integer order, horder

! order is the order of the polynomial fit done to the nearby data points to 
! estimate the taper at the required position
   parameter (order = 3)
   parameter (horder = order/2)
   real(kind=dp), intent(in) :: del_ra, del_dec
   real(kind=dp) :: taper
   real(kind=dp) :: dy, x1, x2 ,y
   real(kind=dp), dimension(order) :: x1a, x2a
   real(kind=dp), dimension(order,order) :: ya
   real(kind=dp) :: xtmp, ytmp, terr
   integer :: i, j, k, l, m, n

! x and y in arcseconds - convert to pixels

   xtmp = del_ra/pb_spacing
   ytmp = del_dec/pb_spacing

! shift to centre of datpoint * datpoint scan

   xtmp = xtmp+pb_points/2+1.0
   ytmp = ytmp+pb_points/2+1.0
   i = int(xtmp)
   j = int(ytmp)

! put the appropriate values into the arrays required by polin2
   do k = 1,order
      x1a(k) = 1.0*(i-horder+k)
      x2a(k) = 1.0*(j-horder+k)
      do l = 1,order
      if ((i+k-horder.ge.1).and.(i+k-horder.le.pb_points).and.& 
          (j+l-horder.ge.1).and.(j+l-horder.le.pb_points)) then
             ya(k,l) = pb_data(i+k-horder,j+l-horder)
         else
             ya(k,l) = 0.0
         end if
      end do
   end do

   m = order
   n = order

   call polin2(x1a,x2a,ya,m,n,xtmp,ytmp,taper,terr)

! make primary beam positive
   taper = abs(taper)

end subroutine gridded_pb

! ****************************************************************

! returns primary beam (power) at given r (arcsec) for a Gaussian PB
! sigma_pb = FWHM/2(sqrt (2 ln 2))
subroutine gaussian_pb(r,taper,sigma_pb)

   use kind_def
   use sz_globals

   implicit none
   real(kind=dp), intent(in) :: r, sigma_pb
   real(kind=dp), intent(out) :: taper

   taper=exp(-r**2/(2.*sigma_pb**2))

end subroutine gaussian_pb

! **************************************************************************

! Returns primary beam (power) at a given r (arcsec) based on a polynomial
subroutine polynomial_pb(r,taper)

   use kind_def
   use sz_globals

   real(kind=dp), intent(in) :: r
   real(kind=dp), intent(out) :: taper
   integer :: i

   if (r.lt.poly_cut) then
      taper = 1.0
      do i = 1, n_poly_pb
         taper = taper+poly_pb(i)*r**i
      end do
   else
      taper = 0.0
   end if 

end subroutine polynomial_pb

! ****************************************************************

! returns RT primary beam (power) at given r (arcsec). Polynomial fit.
subroutine ryle_pb(r,taper)

   use kind_def
   use sz_globals

   implicit none
   real(kind=dp), intent(in) :: r
   real(kind=dp), intent(out) :: taper

!     NB r in arcsec!

!     primary beam modelled by polynomial until r=286'' then goes linearly 
!     to zero at 300'' (pb = 0.05 at r=286'')

   if (r.lt.300) then
      taper =  1.0-2.08e-4*r-2.61e-5*r**2+5.32e-8*r**3
      if (taper.lt.0.05) then
         taper = 0.05 * (300-r)/14
      end if
   else
      taper = 0.0
   end if

end subroutine ryle_pb

! ***************************************************************************

subroutine cat_poly(xf,res)

   use kind_def

   implicit none
   real(kind=dp), intent(out) :: res
   real(kind=dp), intent(in) :: xf

   res =  0.9987 &
        + 1.2e-3*xf &
        - 2.543e-3*xf**2 &
        + 8.49e-6*xf**3 &
        +3.262e-6*xf**4 &
        -7.954e-8*xf**5 &
        +5.573e-10*xf**6 

end subroutine cat_poly

! ***************************************************************************

! returns CAT primary beam (power) at given r (arcsec). Polynomial fit.
subroutine cat_pb(r,taper)

   use kind_def
   use sz_globals

   implicit none
   real(kind=dp), intent(in) :: r
   real(kind=dp), intent(out) :: taper
   real(kind=dp) :: rdeg, freq, xf, cutoff

   rdeg = r/3600.
   freq = 15.5
   xf = rdeg*freq
   call cat_poly(3.D0*freq,cutoff)
   if (rdeg.lt.3.) then
      call cat_poly(xf,taper)
   else if (rdeg.lt.3.5) then
      taper = cutoff * (3.5-rdeg)/0.5
   else
      taper = 0.0
   end if

end subroutine cat_pb

! *************************************************************************

subroutine vla_pb(r,taper)
! factor of 1e9 to convert obsfreq into GHz

   use kind_def
   use sz_globals

   implicit none
   real(kind=dp), intent(in) :: r
   real(kind=dp), intent(out) :: taper
   real(kind=dp) :: x

! Convert angle into arcminutes
   x = (r/60.0*(obsfreq/1.0D9))**2

   taper = 1/(0.99204+0.99569e-3*x+0.38146e-5*x**2-&
        & 0.53117e-8*x**3+0.3981e-11*x**4)

end subroutine vla_pb

! **************************************************************************

! returns primary beam (power) at given r (arcsec) by scaling from
! the VLA polynomial 
subroutine scaled_vla_pb(r,taper)

   use kind_def
   use sz_globals

   implicit none
   real(kind=dp), intent(in) :: r
   real(kind=dp), intent(out) :: taper
   real(kind=dp) :: r1

! Scale r for differences between VLA and desired telescope
   r1 = r/((vla_dish_diameter/dish_diameter))

   call vla_pb(r1,taper)

end subroutine scaled_vla_pb

! **************************************************************************

! returns primary beam (power) at given r (arcsec) by scaling from
! the RT polynomial 
subroutine scaled_rt_pb(r,taper)

   use kind_def
   use sz_globals

   implicit none
   real(kind=dp), intent(in) :: r
   real(kind=dp), intent(out) :: taper
   real(kind=dp) :: r1

! Scale r for differences between RT and desired telescope
   r1 = r/(((rt_freq)/obsfreq)*(rt_dish_diameter/dish_diameter))

   call ryle_pb(r1,taper)

end subroutine scaled_rt_pb

! *************************************************************************

! returns primary beam (power) at given r (arcsec) by scaling from
! the CAT polynomial 
subroutine scaled_cat_pb(r,taper)

   use kind_def
   use sz_globals

   implicit none
   real(kind=dp), intent(in) :: r
   real(kind=dp), intent(out) :: taper
   real(kind=dp) :: r1

! Scale rdeg appropriate for desired telescope
   r1 = r/(((cat_freq)/obsfreq)*(cat_dish_diameter/dish_diameter))

   call cat_pb(r1,taper)

end subroutine scaled_cat_pb

! **************************************************************************

! A quick and nasty beam for simulating surveys
subroutine survey_pb(r,taper)
   
   use kind_def
   use sz_globals
   
   implicit none
   real(kind=dp), intent(in) :: r
   real(kind=dp), intent(out) :: taper

   if (r.lt.pb_sbeam) then
      taper = 1.0
   else 
      taper=exp(-(r-pb_sbeam)**2/(2*pb_rolloff**2))
   end if

end subroutine survey_pb

! **********************************************************************

subroutine vsa_ea_a_pb(sensi)

   use kind_def
   use sz_globals

   implicit none
  
   real(kind=dp), intent(out) :: sensi(maxsize,maxsize)
   real(kind=dp) :: r(19),rdeg(19),freq,xf,cat_poly,cutoff,tmp_taper,scale
   integer :: i, j, doscale, mapsize, spacing, np, n_beam, x, y
   real(dp), dimension(2,19) :: position
   real(dp), dimension(2,19) :: centre, offset
   real(dp) :: sigma, gaussian, norm, distance

   n_beam = 19
   mapsize = maxsize
   spacing = int((0.7*7603)/cellsize)
   sigma = 7603/(2*sqrt(2*log(2.0)))

!   write(*,*) 'field spacing:',spacing



   centre(1,1) =  mapsize/2
   centre(2,1) =  mapsize/2
   centre(1,2) =  mapsize/2-spacing
   centre(2,2) =  mapsize/2
   centre(1,3) =  mapsize/2+spacing
   centre(2,3) =  mapsize/2
   centre(1,4) =  mapsize/2-spacing/2.0
   centre(2,4) =  mapsize/2-((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,5) =  mapsize/2+spacing/2.0
   centre(2,5) =  mapsize/2-((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,6) =  mapsize/2-spacing/2.0
   centre(2,6) =  mapsize/2+((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,7) =  mapsize/2+spacing/2.0
   centre(2,7) =  mapsize/2+((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,8) =  mapsize/2-(2*spacing)
   centre(2,8) =  mapsize/2
   centre(1,9) =  mapsize/2+(2*spacing)
   centre(2,9) =  mapsize/2
   centre(1,10) =  centre(1,3)+spacing/2.0
   centre(2,10) =  centre(2,3)+((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,11) =  centre(1,3)+spacing/2.0
   centre(2,11) =  centre(2,3)-((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,12) =  centre(1,2)-spacing/2.0
   centre(2,12) =  centre(2,2)+((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,13) =  centre(1,2)-spacing/2.0
   centre(2,13) =  centre(2,2)-((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,14) =  centre(1,6)-spacing/2.0
   centre(2,14) =  centre(2,6)+((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,15) =  centre(1,6)+spacing/2.0
   centre(2,15) =  centre(2,6)+((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,16) =  centre(1,6)+(3.0*spacing/2.0)
   centre(2,16) =  centre(2,6)+((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,17) =  centre(1,4)-spacing/2.0
   centre(2,17) =  centre(2,4)-((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,18) =  centre(1,4)+spacing/2.0
   centre(2,18) =  centre(2,4)-((spacing**2-(spacing/2.0)**2)**0.5)
   centre(1,19) =  centre(1,4)+(3.0*spacing/2.0)
   centre(2,19) =  centre(2,4)-((spacing**2-(spacing/2.0)**2)**0.5)  


   do np = 1, n_beam
      do x = 1,mapsize
         do y = 1,mapsize
            
            position(1,np) = x
            position(2,np) = y

! find the distance from the centre
            distance = ((position(1,np)-centre(1,np))**2+&
                        (position(2,np)-centre(2,np))**2)**0.5

! calculate value of beam at this position
            gaussian = exp((-(distance**2))/(2*sigma**2))

! add in contribution from this pointing to sensitivity map
            sensi(x,y) = (sensi(x,y)**2 + gaussian**2)**0.5 

         end do
      end do
   end do

   norm = maxval(sensi)
   sensi = sensi/norm


end subroutine vsa_ea_a_pb

! ----------------------------------------------------------------

subroutine elliptical_pb(x,y,taper)

   use kind_def
   use sz_globals

   implicit none

   real(kind=dp), intent(in)  :: x,y
   real(kind=dp), intent(out) :: taper
   real(kind=dp)  :: xx,yy

! perform rotation
   xx = x*cos(el_theta*pi/180.)+y*sin(el_theta*pi/180.)
   yy = y*cos(el_theta*pi/180.)-x*sin(el_theta*pi/180.)

   taper = exp(-0.5*real((xx**2/el_sigma1**2)))
   taper = taper*exp(-0.5*real(yy**2/el_sigma2**2))
   
end subroutine elliptical_pb

!----------------------------------------------------------------

subroutine SPH10_beam(temp1,temp2,taper)

   use kind_def
   use sz_globals

   real(kind=dp), intent(in)  :: temp1,temp2
   real(kind=dp), intent(out) :: taper
   real(kind=dp)  :: theta,cosphi,sinphi,phi
   

   theta = dble(sqrt(real(temp1**2+temp2**2)))
   
! this is the real part:
   taper = cos(4.*theta*sec2rad)
   
end subroutine SPH10_beam

!----------------------------------------------------------------

subroutine SPH11_beam(temp1,temp2,taper)

   use kind_def
   use sz_globals

   real(kind=dp), intent(in)  :: temp1,temp2
   real(kind=dp), intent(out) :: taper
   real(kind=dp)  :: theta,cosphi,sinphi,phi
   

   theta = dble(sqrt(real(temp1**2+temp2**2)))
   if (phi==0.) then
      cosphi = 1.
      sinphi = 1.
   else
      cosphi = temp1/theta
      sinphi = temp2/theta     
   end if

! this is the real part:
   taper = sin(theta*sec2rad)*cosphi
   
end subroutine SPH11_beam

!----------------------------------------------------------------

subroutine SPH22_beam(temp1,temp2,taper)

   use kind_def
   use sz_globals

   real(kind=dp), intent(in)  :: temp1,temp2
   real(kind=dp), intent(out) :: taper
   real(kind=dp)  :: theta,cosphi,sinphi,cos2phi,phi
   

   theta = dble(sqrt(real(temp1**2+temp2**2)))
   if ((temp1.gt.0.).and.(temp2.gt.0.)) then
      phi = atan(temp2/temp1)
   elseif ((temp1.gt.0.).and.(temp2.lt.0.)) then
      phi = atan(-1.*temp2/temp1)
      phi = (2.*pi) - phi
   elseif ((temp1.lt.0.).and.(temp2.gt.0.)) then
      phi = atan(-1.*temp1/temp2)
      phi = phi + (pi/2.)
   elseif ((temp1.lt.0.).and.(temp2.lt.0.)) then
      phi = atan(temp2/temp1)
      phi = phi + pi
   elseif ((temp1.eq.0.).and.(temp2.gt.0.)) then
      phi = pi/2.
   elseif ((temp1.eq.0.).and.(temp2.lt.0.)) then
      phi = 3.*pi/2.
   elseif ((temp1.gt.0.).and.(temp2.eq.0.)) then
      phi = 0.
   elseif ((temp1.lt.0.).and.(temp2.eq.0.)) then
      phi = pi
   end if


! this is the real part:
   taper = sin(theta*sec2rad)**2*cos(2.*phi)
   
end subroutine SPH22_beam

!----------------------------------------------------------------

! l = 3, m = 1 ...

subroutine SPH31_beam(temp1,temp2,taper)

   use kind_def
   use sz_globals

   real(kind=dp), intent(in)  :: temp1,temp2
   real(kind=dp), intent(out) :: taper
   real(kind=dp)  :: theta,cosphi,sinphi,cos2phi, phi
   
   theta = dble(sqrt(real(temp1**2+temp2**2)))
   if ((temp1.gt.0.).and.(temp2.gt.0.)) then
      phi = atan(temp2/temp1)
   elseif ((temp1.gt.0.).and.(temp2.lt.0.)) then
      phi = atan(-1.*temp2/temp1)
      phi = (2.*pi) - phi
   elseif ((temp1.lt.0.).and.(temp2.gt.0.)) then
      phi = atan(-1.*temp1/temp2)
      phi = phi + (pi/2.)
   elseif ((temp1.lt.0.).and.(temp2.lt.0.)) then
      phi = atan(temp2/temp1)
      phi = phi + pi
   elseif ((temp1.eq.0.).and.(temp2.gt.0.)) then
      phi = pi/2.
   elseif ((temp1.eq.0.).and.(temp2.lt.0.)) then
      phi = 3.*pi/2.
   elseif ((temp1.gt.0.).and.(temp2.eq.0.)) then
      phi = 0.
   elseif ((temp1.lt.0.).and.(temp2.eq.0.)) then
      phi = pi
   end if

! this is the real part:
   taper = sin(theta*sec2rad)*cos(phi)*(5.*cos(theta*sec2rad)**2 - 1)
   
end subroutine SPH31_beam

!----------------------------------------------------------------

! l = 5, m = 2 ...

subroutine SPH52_beam(temp1,temp2,taper)

   use kind_def
   use sz_globals

   real(kind=dp), intent(in)  :: temp1,temp2
   real(kind=dp), intent(out) :: taper
   real(kind=dp)  :: theta,cosphi,sinphi,cos2phi, phi
   

   theta = dble(sqrt(real(temp1**2+temp2**2)))
   if ((temp1.gt.0.).and.(temp2.gt.0.)) then
      phi = atan(temp2/temp1)
   elseif ((temp1.gt.0.).and.(temp2.lt.0.)) then
      phi = atan(-1.*temp2/temp1)
      phi = (2.*pi) - phi
   elseif ((temp1.lt.0.).and.(temp2.gt.0.)) then
      phi = atan(-1.*temp1/temp2)
      phi = phi + (pi/2.)
   elseif ((temp1.lt.0.).and.(temp2.lt.0.)) then
      phi = atan(temp2/temp1)
      phi = phi + pi
   elseif ((temp1.eq.0.).and.(temp2.gt.0.)) then
      phi = pi/2.
   elseif ((temp1.eq.0.).and.(temp2.lt.0.)) then
      phi = 3.*pi/2.
   elseif ((temp1.gt.0.).and.(temp2.eq.0.)) then
      phi = 0.
   elseif ((temp1.lt.0.).and.(temp2.eq.0.)) then
      phi = pi
   end if

   phi = phi+(sph_theta*pi/180.)

! this is the real part:
   taper = sin(theta*sec2rad)**2*cos(2.*phi)*(3.*cos(theta*sec2rad)**3 - cos(theta*sec2rad))
   
end subroutine SPH52_beam

!----------------------------------------------------------------

subroutine szi_beam

   use kind_def
   use sz_globals
   use kgraph

   implicit none

   integer  :: i,j,k,l,m
   real(dp) :: taper,beammax
   real(dp) :: xspacing,yspacing,sigma2
   real(dp), allocatable :: centre(:,:)
   integer  pgopen
   external pgopen


!  FWHM of primary beam
   xspacing = 1.22*1./(dish_diameter)*(1./sec2rad)*1./cellsize
   xspacing = xspacing/2.
   yspacing = sqrt(3.)*xspacing
   
   allocate(centre(fparray,2))


   if (fparray==1) then
      centre(1,1) = 0.
      centre(1,2) = 0.
   elseif (fparray==7) then
      centre(1,1) = -1.*xspacing
      centre(1,2) = -1.*yspacing
      centre(2,1) = xspacing
      centre(2,2) = -1.*yspacing
      centre(3,1) = -2.*xspacing
      centre(3,2) = 0.
      centre(4,1) = 0.
      centre(4,2) = 0.
      centre(5,1) = 2*xspacing
      centre(5,2) = 0.
      centre(6,1) = -1.*xspacing
      centre(6,2) = yspacing
      centre(7,1) = xspacing
      centre(7,2) = yspacing
   elseif (fparray==19) then
      centre(1,1) = -2.*xspacing
      centre(1,2) = -2.*yspacing
      centre(2,1) = 0.
      centre(2,2) = -2.*yspacing
      centre(3,1) = 2.*xspacing
      centre(3,2) = -2.*yspacing
      centre(4,1) = -3.*xspacing
      centre(4,2) = -1.*yspacing
      centre(5,1) = -1.*xspacing
      centre(5,2) = -1.*yspacing
      centre(6,1) = xspacing
      centre(6,2) = -1.*yspacing
      centre(7,1) = 3.*xspacing
      centre(7,2) = -1.*yspacing
      centre(8,1) = -4.*xspacing
      centre(8,2) = 0.
      centre(9,1) = -2.*xspacing
      centre(9,2) = 0.
      centre(10,1)= 0.
      centre(10,2)= 0.
      centre(11,1)= 2.*xspacing
      centre(11,2)= 0.
      centre(12,1)= 4.*xspacing
      centre(12,2)= 0.
      centre(13,1)= -3.*xspacing
      centre(13,2)= yspacing
      centre(14,1)= -1.*xspacing
      centre(14,2)= yspacing
      centre(15,1)= xspacing
      centre(15,2)= yspacing
      centre(16,1)= 3.*xspacing
      centre(16,2)= yspacing
      centre(17,1)= -2.*xspacing
      centre(17,2)= 2.*yspacing
      centre(18,1)= 0.
      centre(18,2)= 2.*yspacing
      centre(19,1)= 2.*xspacing
      centre(19,2)= 2.*yspacing
   end if


   if (do_plot) then
!      i = pgopen('beammap.ps/CPS')
      call pgsfs(2)
      call pgenv(-1.*real(maxsize)/2.,real(maxsize)/2.,&
           -1.*real(maxsize)/2.,real(maxsize)/2.,1,1)
      i = 1
      call pgsci(3)
      call pgsls(1)
      call pgslw(5)
      call pgcirc(real(centre(i,1)),real(centre(i,2)),real(xspacing))
      call pgsls(2)
      call pgslw(2)
      call pgcirc(real(centre(i,1)+(xspacing)),real(centre(i,2)),real(xspacing))
      call pgcirc(real(centre(i,1)+(xspacing/2.)),&
           real(centre(i,2)-(yspacing/2.)),real(xspacing))
      call pgcirc(real(centre(i,1)-(xspacing/2.)),&
           real(centre(i,2)-(yspacing/2.)),real(xspacing))
      do i = 2,fparray
         call pgsci(1)
         call pgsls(1)
         call pgslw(5)
         call pgcirc(real(centre(i,1)),real(centre(i,2)),real(xspacing))
         call pgsci(2)
         call pgsls(2)
         call pgslw(2)
         call pgcirc(real(centre(i,1)+(xspacing)),real(centre(i,2)),real(xspacing))
         call pgcirc(real(centre(i,1)+(xspacing/2.)),&
              real(centre(i,2)-(yspacing/2.)),real(xspacing))
         call pgcirc(real(centre(i,1)-(xspacing/2.)),&
              real(centre(i,2)-(yspacing/2.)),real(xspacing))
      end do
      call pgsci(1)
      call pgsls(1)
!      call pgclos()
   end if

   sigma2 = (xspacing)**2/(2.*log(2.))
   taper = 0.
   szibeam = 0.
   do i = -maxsize/2,maxsize/2-1
      do j = -maxsize/2,maxsize/2-1
         do k = 1,fparray
            taper = exp(-1.*(real(i)-centre(k,1))**2/(2.*sigma2))&
                 *exp(-1.*(real(j)-centre(k,2))**2/(2.*sigma2))
            szibeam(i+maxsize/2+1,j+maxsize/2+1,1) = &
                 sqrt(szibeam(i+maxsize/2+1,j+maxsize/2+1,1)**2+taper**2)
            taper = exp(-1.*(real(i)-(xspacing)-centre(k,1))**2/(2.*sigma2))&
                 *exp(-1.*(real(j)-centre(k,2))**2/(2.*sigma2))
            szibeam(i+maxsize/2+1,j+maxsize/2+1,2) = &
                 sqrt(szibeam(i+maxsize/2+1,j+maxsize/2+1,2)**2+taper**2)
            taper = exp(-1.*(real(i)-(xspacing/2.)-centre(k,1))**2/&
                 (2.*sigma2))*exp(-1.*(real(j)+(yspacing/2.)-centre(k,2))**2/&
                 (2.*sigma2))
            szibeam(i+maxsize/2+1,j+maxsize/2+1,3) = &
                 sqrt(szibeam(i+maxsize/2+1,j+maxsize/2+1,3)**2+taper**2)
            taper = exp(-1.*(real(i)+(xspacing/2.)-centre(k,1))**2/&
                 (2.*sigma2))*exp(-1.*(real(j)+(yspacing/2.)-centre(k,2))**2/&
                 (2.*sigma2))
            szibeam(i+maxsize/2+1,j+maxsize/2+1,4) = &
                 sqrt(szibeam(i+maxsize/2+1,j+maxsize/2+1,4)**2+taper**2)
         end do
      end do
   end do  
   
   szioverlay = 0.
   if (overlay) then
      szioverlay = szibeam(:,:,1)+szibeam(:,:,2)&
                 +szibeam(:,:,3)+szibeam(:,:,4)
      beammax = maxval(real(szioverlay))
      szioverlay = szioverlay/beammax
      if (do_plot) then
         write(*,*) 'Displaying overlay:'
         call display_map(szioverlay,maxsize,graph_adj,&
              cellsize,cellsize,0D0,0D0)
      end if
   end if
         
   deallocate(centre)
   
end subroutine szi_beam
   
! ****************************************************************

subroutine get_taper(temp1,temp2,taper)
   
   use kind_def
   use sz_globals

   implicit none

   real(kind=dp), intent(in) :: temp1, temp2
   real(kind=dp) :: r, taper, sigma_pb

   r = sqrt(temp1**2+temp2**2)

! get primary beam
   apply_pb : select case(which_pb)
   case('r','R')
      call scaled_rt_pb(r,taper)               
   case('p','P')
      call polynomial_pb(r,taper)
   case('c','C')
      call scaled_cat_pb(r,taper)
   case('v','V')
      call scaled_vla_pb(r,taper)    
   case('s','S')
      call survey_pb(r,taper)
   case('e','E')
      call elliptical_pb(temp1,temp2,taper)
   case('h','H')
      if (which_sph_harm.eq.1) call SPH10_beam(temp1,temp2,taper)
      if (which_sph_harm.eq.2) call SPH11_beam(temp1,temp2,taper)
      if (which_sph_harm.eq.3) call SPH22_beam(temp1,temp2,taper)
      if (which_sph_harm.eq.4) call SPH31_beam(temp1,temp2,taper)
      if (which_sph_harm.eq.5) call SPH52_beam(temp1,temp2,taper)
   case('u','U')
      call unity_beam(taper)
   case default
      sigma_pb = fwhm(chan)/(2.D0*sqrt(2.D0*log(2.)))
      call gaussian_pb(r,taper,sigma_pb)
   end select apply_pb

end subroutine get_taper

! ****************************************************************

subroutine unity_beam(taper)

   use kind_def
   
   implicit none

   real(kind=dp), intent(out) :: taper

   taper = 1.0d0

end subroutine unity_beam

! ****************************************************************

subroutine display_primary_beam

   use kind_def
   use sz_globals
   use kgraph

   implicit none

   real(kind=dp), dimension(maxsize,maxsize) :: pb
   real(kind=dp) :: taper, temp1, temp2, sigma_pb
   integer :: i, j

   if (nchan.ne.1) then
      chan = nchan/2+1
      call io_geti('Which channel','*',chan,status)
   else
      chan = 1
   end if

! get primary beam
   do i = 1,maxsize
      do j = 1, maxsize            
         temp1 = (float(i-maxsize/2)-1.0)*cellsize-pb_x
         temp2 = (float(j-maxsize/2)-1.0)*cellsize-pb_y
         call get_taper(temp1,temp2,taper)
         pb(i,j) = taper
      end do
   end do

   call display_map(pb,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)

end subroutine display_primary_beam

! ****************************************************************

end module primary_beam




