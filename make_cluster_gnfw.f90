!---------------------------------------------------------------------
!     calculates projected sz profiles for GNFW model
!     Inputs: 
!     Either (for DM-GNFW model):
!     real MT200           total mass within r200
!     real fg200           gas fraction at r200
!     real c200            concentration parameter
!     Or (for observational GNFW model):
!     real Ytot            total Y, arcmin^-2
!     real theta_s         scaling size, arcmin
!     And (for both):
!     real c500, a, b, c   GNFW profile shape parameters
!     real x0, y0          centre of cluster in arcsec relative to the phase centre
!     real cellsize        cell size in arcsec
!     integer ncell        number of cells. 
!
!     Outputs:
!     real(*) sz_profile array of temperature decrement
!                          / 2nkT (Thompson cross section)
!                      T0  | -----------------------------   dl
!                          /            mc^2
!
!     all pixel values refer to value at centre of the pixel
!
!     the centre of the map is the centre of pixel maxsize/2+1,
!     maxsize/2+1.
!
subroutine make_cluster_gnfw

   use kind_def
   use sz_globals
   use maths
   use kgraph
   use physics
   use cosmology

   implicit none

   real(kind=dp), dimension(maxsize/2) :: y_profile, y_integrand
   real(kind=dp), dimension(maxsize/2) :: yarray, theta
   real(kind=dp) :: yintegral, yfunc, GNFW_YsphVolIntegrand, DM_GNFWsphVolInt, r500_fun
   real(kind=dp) :: Y_int, Y500_int, y0_int, Y500, y_coeff, theta500
   real(kind=dp) :: angfactor, rlimit1, pixsize, y0
   real(kind=dp) :: xx, yy, dx, dy, rr, logthetalim1
   real(kind=dp) :: a, b, c, c500, h
   real(kind=dp) :: rhocritz, Mg200, rhos_DM
   real(kind=dp) :: r500, Pei_GNFW, DM_GNFWgasVol, rhocrit, Y5R500, Y5R500_int
   real(kind=dp) :: p, f200, f500, x
   real(kind=dp), parameter :: a1=0.5116d0, a2=-0.4283d0, a3=-3.13d-3, a4=-3.52d-5
   real(kind=dp), parameter :: eps = 1.d-4
   real(kind=dp), parameter :: G=4.518d-48 !Gravitational constant in Mpc^3 M_sun^-1 s^-2
   real(kind=dp), parameter :: J2keV=6.24251d+15, mec2=511.d0
   integer :: i, j, halfsize, iunit
   integer :: centre
   logical, external :: io_yesno
   character*80 :: outfile
   external :: GNFW_YsphVolIntegrand, yintegral
   external :: DM_GNFWsphVolInt, r500_fun
     
   a=a_GNFW
   b=b_GNFW
   c=c_GNFW
   c500=c500_GNFW

   h=H0/100
   rhocrit=3.*(h*100.*m2Mpc*1000.)**2./(8.*pi*G) ! M_sol/MPc^3
   write(*,*) 'rhocrit = ', rhocrit, ' in M_sol/MPc^3'

   select case(model_type)
   case('n', 'N')
     ! Work out observational parameters from physical ones
     rhocritz = rhocrit/one_over_h(z)**2
     Mg200=MT200*fg200     !M_sun
     r200=((3.d0*MT200)/(4.d0*pi*200.d0*rhocritz))**(1.d0/3.d0)   !Mpc
     write(*,*) 'r200 = ', r200, ' MPc'
     rs_DM=r200/c200                   !Mpc
     rhos_DM=(200.d0/3.d0)*((r200/rs_DM)**3.d0)*&
		   ( rhocritz/(DLOG(1.d0 + c200) - ( 1.d0/(1.d0 + 1.d0/c200) )))   !M_sunMpc-3
     !r500=r200/1.5d0                 !Mpc. The coeff 1.5 was derived from testing the model over a wide
	                                            ! range of cluster masses see paper for the derivation
     !write(*,*) 'r500 = r200/1.5 = ', r500
     !r500 = zbrent(r500_fun, eps, r200, eps) ! solve the implicit function for r500 given r200, c200
     !write(*,*) 'Numerical solution for r500 is', r500
     ! Use Hu & Kravtsov 2003 fitting formula
     f200=(1.d0/c200)**3*(DLOG(1.d0+c200)-1.d0/(1+1.d0/c200))
     f500=500d0/200d0*f200
     p=a2+a3*DLOG(f500)+a4*(DLOG(f500))**2
     x=(a1*f500**(2.d0*p)+(3.d0/4.d0)**2)**(-1.d0/2.d0)+2.d0*f500
     r500=rs_DM/x
     write(*,*) 'Using analytical solution for r500:', r500
     rp_GNFW=r500/c500_GNFW    !Mpc
     call qtrap(DM_GNFWsphVolInt, 1.0d-4, r200, eps, DM_GNFWgasVol)
     write(*,*) DM_GNFWgasVol
     Pei_GNFW=((mu_n/mu_e)*(G*rhos_DM*rs_DM*rs_DM*rs_DM)*Mg200)/DM_GNFWgasVol
     Pei_GNFW= Pei_GNFW*(m_0*m2Mpc) *(J2keV)       !keVm-3
     theta_s=(rp_GNFW/D_theta)/min2rad
     Ytot=(4.d0*pi*sigma_T*Pei_GNFW*rp_GNFW)/(mec2*m2Mpc)*theta_s*theta_s/a_GNFW*&
	           Gammafun((3.d0-c_GNFW)/a_GNFW) *&
		   Gammafun((b_GNFW-3.d0)/a_GNFW)/Gammafun((b_GNFW-c_GNFW)/a_GNFW)
     if (verbose) then
       write(*,*)
       write(*,*) 'Critical density is', rhocritz, 'M_sol/Mpc^3 at redshift', z
       write(*,*) 'Dark matter concentration parameter c200 is', c200
       write(*,*) 'r200 is', r200, 'Mpc'
       write(*,*) 'Normalisation constant for NFW DM halo, rhos_DM is', rhos_DM, 'M_sol/Mpc^3'
       write(*,*) 'Normalisation constant for gas profile, Pei is', Pei_GNFW, 'keV/m^3'
       write(*,*) 'r500 is', r500, 'Mpc'
       write(*,*) 'theta_s is', theta_s, 'arcmin'
       write(*,*) 'Ytot is', Ytot, 'arcmin^2'
     end if
   end select
   theta500 = theta_s*c500
   halfsize = maxsize/2

   allocate(logyintegrand(halfsize))
   allocate(logtheta(halfsize))
   allocate(logyarray(halfsize))
   
!  Initialise arrays
   compton_y = 0.d0
   logyintegrand = 0.d0
   logyarray = 0.d0
   y_profile = 0.d0
   y_integrand = 0.d0
   yarray = 0.d0
   theta = 0.d0

!  Calculate the spherical integral of the electron pressure
   Y_int = Gammafun((3.0-c)/a)*Gammafun((b-3.0)/a)/Gammafun((-c+b)/a)/a
   call qtrap(GNFW_YsphVolIntegrand, 0d0, c500, eps, Y500_int)
   Y500 = Ytot/Y_int*Y500_int
   write(*,*) 'Y500 = Ytot/', Y_int/Y500_int
   call qtrap(GNFW_YsphVolIntegrand, 0d0, 5d0*c500, eps, Y5R500_int)
   Y5R500 = Ytot/Y_int*Y5R500_int
   write(*,*) 'Y5R500 = Ytot/', Y_int/Y5R500_int
   
!  Calculate the integral of the electron pressure over line of sight at the centre
   y0_int = Gammafun((1.0-c)/a)*Gammafun((b-1.0)/a)/Gammafun((-c+b)/a)/a
!  Use these to get the normalisation constants for y
   y_coeff = Ytot/Y_int/(theta_s**3)/4.0/pi
!  y value at the centre
   y0 = 2.0*y_coeff*theta_s*y0_int
    
!  Now set up y-integrand arrays
!  Make the limit a little bit bigger to avoid interpolation problems later
   logthetalim1 = log10(thetalimit) + (log10(thetalimit) - log10(thetamin))/dble(halfsize-2)
   DO i=1,halfsize
     logtheta(i) = log10(thetamin) + (logthetalim1 - log10(thetamin)) * (i-1)/dble(halfsize-1)
     theta(i) = 10.d0 ** logtheta(i)
     y_integrand(i)=y_coeff*((theta(i)/theta_s)**(-c))*&	     
	 ( (1.d0 + ((theta(i)/theta_s)**(a)))**((c - b)/a))
     logyintegrand(i)=funlog10(y_integrand(i))
   ENDDO

   if (verbose) then
      write(*,*)
      write(*,*) 'theta_500 is : ',theta500,' arcmin'
      write(*,*) 'Y_500 is     : ',Y500,' arcmin^2'
   end if
     
!  Project pressure and integrate over line of sight
   DO i=1,halfsize
     thetai=theta(i)
     rlimit1=sqrt( max(thetalimit*thetalimit - thetai*thetai, 0d0) )
     if( rlimit1 >0d0 ) call qtrap(yintegral,-rlimit1,rlimit1,eps,yarray(i))
     if(yarray(i) > 1d10) then
	write(*,*) 'i,theta(i),rlimit1,yarray(i)=',i,theta(i),rlimit1,yarray(i)		   
        write(*,*) 'QTRAP error.'
     endif
     logyarray(i)=funlog10(yarray(i))
   ENDDO

!  Write out a y(theta) profile
   !if (io_yesno('Write a y(theta) profile?', 'n', status)) then
   !      status = 0
   !      call io_nxtlun(iunit,status)
   !      call io_getc('Filename:','',outfile,status)
   !      open (iunit,file=outfile,form = 'formatted')
   !      write(iunit,*) '# theta  y'
   !      write(iunit,*) 0d0, y0
   !      do i=1, halfsize
   !         write(iunit,*) theta(i), yarray(i)
   !      end do
   !      close(iunit)
   !end if
   

!  Make map of Comptonisation parameter y
   angfactor=1.0/60.0
   do i=1,maxsize
     do j=1,maxsize
	xx=-float(halfsize+1)*cellsize+i*cellsize
	yy=-float(halfsize+1)*cellsize+j*cellsize
	dx=(xx-sz_x)
	dy=(yy-sz_y)
	if (ellipsoid) then  ! Elliptical geometry:
	  rr=(f_GNFW*cos(el_alpha)**2+sin(el_alpha)**2/f_GNFW)*dx**2&
            +(f_GNFW*sin(el_alpha)**2+cos(el_alpha)**2/f_GNFW)*dy**2&
            +2*cos(el_alpha)*sin(el_alpha)*(f_GNFW-1./f_GNFW)*dx*dy
          rr=sqrt(rr)*angfactor
	else                 ! Spherical geometry:
          rr=sqrt(dx*dx+dy*dy)*angfactor
	endif
	IF( abs(dx)<cellsize/2 .AND. abs(dy)<cellsize/2 )THEN
	  compton_y(i,j)=y0
	ELSE
	  compton_y(i,j)=compton_y(i,j)+yfunc(rr)
	ENDIF
     enddo
   enddo

   pixsize = cellsize**2
   if (do_plot) then
      write(*,*) '2d map of y'
      call display_map(compton_y,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
   end if

   if (do_plot) then
      call pgpage
      call display_profile(yarray/y0,halfsize,cellsize)
      write(*,*) 'y against radius'
   end if

! Now work out brightness temperatures of SZ effect at each channel
   do chan = 1,nchan
      if (overwrite) then
         szsky(:,:,chan) = compton_y*g_nu_x2(dble(nu(chan)))*T0
!         szsky(:,:,chan) = compton_y
      else
         szsky(:,:,chan) = szsky(:,:,chan)+ &
                           compton_y*g_nu_x2(dble(nu(chan)))*T0
      end if
   end do

   chan = nchan/2+1

   if (verbose) then
      write(*,*)
      write(*,*) 'Central y value                       ',maxval(compton_y)
      write(*,*) 'Central decrement (ie deltaT)         ',&
                 minval(szsky(:,:,chan)), ' K'
      write(*,*) 'Y_tot from y_map (arcmin^2)           ',sum(compton_y)*&
                 pixsize/3600.d0
   end if

   if (do_plot) then
      call pgpage
      write(*,*) 'Displaying SZ sky (values of deltaT)'
      call display_map(szsky(:,:,chan),maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
   end if
         
   deallocate(logyintegrand)
   deallocate(logtheta)
   deallocate(logyarray)
 
end subroutine make_cluster_gnfw

!================================================================================

FUNCTION yintegral(zz)
   
   use sz_globals
   use maths
	
   IMPLICIT NONE	
   REAL*8     :: yintegral
   real(dp)   :: r, zz, result, i
   integer    :: i1

   r = sqrt( thetai*thetai + zz*zz )
   IF( r > thetalimit ) then
     yintegral = 0.d0
   ELSE
     IF(r < thetamin ) THEN
       r = thetamin
     ENDIF
     i1 = locate_dp(logtheta, funlog10(r))
     IF (i1 == maxsize/2) THEN
       i = i1 + (funlog10(r) - logtheta(i1))/(logtheta(i1)-logtheta(i1-1))
     ELSE 
       i = i1 + (funlog10(r) - logtheta(i1))/(logtheta(i1+1)-logtheta(i1)) 
     ENDIF
     CALL interp_1d_db(logyintegrand,maxsize/2,i,result)
     yintegral = 10.0**result
   ENDIF
	
END FUNCTION yintegral

!========================================================

FUNCTION yfunc(r)

  use sz_globals
  use maths
	
  IMPLICIT NONE
	
  REAL*8     :: r, yfunc, result, i
  integer    :: i1

  if( r > thetalimit) then
    yfunc=0.
    return
  else if( r < thetamin ) then
    r = thetamin 
  endif
  i1 = locate_dp(logtheta, funlog10(r))
  IF (i1 == maxsize/2) THEN
     i = i1 + (funlog10(r) - logtheta(i1))/(logtheta(i1)-logtheta(i1-1))
  ELSE 
     i = i1 + (funlog10(r) - logtheta(i1))/(logtheta(i1+1)-logtheta(i1)) 
  ENDIF
  CALL interp_1d_db(logyarray,maxsize/2,i,result)
  yfunc=10.0**result
  return

END FUNCTION yfunc

!========================================================

function GNFW_YsphVolIntegrand(r)

  use sz_globals
	
  REAL*8     :: r, GNFW_YsphVolIntegrand

  GNFW_YsphVolIntegrand = r**(2-c_GNFW)*(1+r**a_GNFW)**((c_GNFW-b_GNFW)/a_GNFW)

end function GNFW_YsphVolIntegrand


!=========================================================================

function DM_GNFWsphVolInt(r)
	
  use sz_globals
	
  REAL *8               :: r, DM_GNFWsphVolInt

  if (r<rmin) then
    r=rmin
  endif

  if (r<=rmax) then
    DM_GNFWsphVolInt=((r*r*r)/((DLOG(1.0 + (r/rs_DM))) - (1.0/(1.0 +(rs_DM/r))))) *&
       ((r/rp_GNFW)**(-1.0*c_GNFW))*&
       ( (1.0 + (r/rp_GNFW)**(a_GNFW))**(-1.0*(a_GNFW +b_GNFW -c_GNFW)/a_GNFW))*&
       ((b_GNFW*((r/rp_GNFW)**(a_GNFW))) + c_GNFW)
  else
    DM_GNFWsphVolInt=0.d0		
  endif

end function DM_GNFWsphVolInt
  
!=========================================================================

function r500_fun(r500)

  use sz_globals

  REAL*8         :: r500, Eps, r500_fun

  Eps = (log(1+r500*c200/r200) - 1./(1+r200/r500/c200))/(log(1+c200)-1./(1+1./c200))
  r500_fun = r500 - r200*(2./5.*Eps)**(1./3.)
  !write(*,*) r200, c200, r500, r500_fun

end function r500_fun
