! Define cluster model
subroutine get_cluster_model

   use kind_def
   use sz_globals
   use file_handling
   use cosmology

   implicit none

   integer ::  field
   real(kind=dp) :: arcsec2kpc, arcmin2kpc, h
   real(kind=dp) :: a, b
   character*1 :: chr1, c200_type
   integer :: i
   integer :: dat_read
   logical :: do_kpc
   logical, external :: io_yesno

   from_file = .false.
   from_poly = .false.

   do_kpc = io_yesno('Specify dimensions in kpc rather than arcsec',&
            'no',status)
   call io_getd('Redshift:','*',z,status)
         
! Calculate the angular size distance and luminosity distance in metres
   call lumdist
   D_lum = D_lum*Mpc         ! D_lum in metres now....
   call angdist
   D_theta = D_theta*Mpc
   arcsec2kpc = D_theta/kpc*sec2rad
   arcmin2kpc = D_theta/kpc*min2rad
   write(*,*) 'Luminosity distance is : ',D_lum/Mpc,' Mpc'
   write(*,*) 'Angular distance is : ',D_theta/Mpc,' Mpc'
   write(*,*) '1 arcsec corresponds to ',arcsec2kpc,' kpc'

   write(*,*) 'Select model type'
   write(*,*) '(B)eta model - isothermal'
   write(*,*) '(D)ouble beta model'
   write(*,*) '(F)ile with n and T profiles'
   write(*,*) 'Beta model - (P)olynomial T model'
   write(*,*) 'Beta model - (X)-ray profile T model'
   write(*,*) 'Beta model - X-ray (M)ap T model'
   write(*,*) 'Beta model - temperature (H)alo'
   write(*,*) '(G)NFW model'
   write(*,*) 'DM-G(N)FW model'
   write(*,*)

   call io_getc('Model type','g',model_type,status)
   select case(model_type)
   case('d','D')
      call io_getd('Central number density /cm^-3:','*',n0,status)
      call io_getd('Central density for 2nd beta model:','*',nx_2,status) 
      call io_getd('Beta parameter:','*',beta,status)
      call io_getd('Second beta for 2beta model:','*',beta2,status) 

      ellipsoid = io_yesno('Ellipsoidal model (y/n):','*',status)
      if (ellipsoid) then
         write(*,*) 'Axes 1 and 2 on sky; axis 3 along los'
         if (do_kpc) then
            theta_c = theta_c*arcsec2kpc
            call io_getd('Core Radius Axis 1 (kpc):','*',theta_c(1),status)
            call io_getd('Core Radius Axis 2 (kpc):','*',theta_c(2),status)
            call io_getd('Core Radius Axis 3 (kpc):','*',theta_c(3),status)
            theta_c = theta_c/arcsec2kpc
            write(*,*) '2nd beta model assumed to have same ellipticity&
                       & as first'
            theta_cx2a = theta_cx2a*arcsec2kpc
            call io_getd('Major axis radius for 2nd beta model (kpc):','*',&
                          theta_cx2a,status) 
            theta_cx2a = theta_cx2a/arcsec2kpc
         else
            call io_getd('Core Radius Axis 1 ("):','*',theta_c(1),status)
            call io_getd('Core Radius Axis 2 ("):','*',theta_c(2),status)
            call io_getd('Core Radius Axis 3 ("):','*',theta_c(3),status)
            write(*,*) '2nd beta model assumed to have same ellipticity&
                       & as first'
            call io_getd('Major axis radius for 2nd beta model ("):','*',&
                          theta_cx2a,status) 
         end if
            
         if (theta_c(1).gt.theta_c(2)) then 
            theta_cx = theta_c(1)
         else
            theta_cx = theta_c(2)
         end if
 
!         write (*,*) 'Euler angles for ellipsoid'
!         do i=1,3
!            ang(i)=ang(i)/deg2rad
!         end do

!         call io_getd('Euler angle 1:','*',ang(1),status)
!         call io_getd('Euler angle 2:','*',ang(2),status)
!         call io_getd('Euler angle 3:','*',ang(3),status)
!         do i=1,3
!            ang(i)=ang(i)*deg2rad
!         end do
         call io_getd('Angle of ellipsoid on sky:','*',el_alpha,status)
!         el_alpha = ang(1) 
      else
         if (do_kpc) then
            theta_cx = theta_cx*arcsec2kpc
            call io_getd('Core Radius (kpc):','*',theta_cx,status)
            theta_cx = theta_cx/arcsec2kpc
            theta_cx2a = theta_cx2a*arcsec2kpc
            call io_getd('Core radius for 2nd beta model (kpc):','*',&
                         theta_cx2a,status) 
            theta_cx2a = theta_cx2a/arcsec2kpc
         else
            call io_getd('Core Radius ("):','*',theta_cx,status)
            call io_getd('Core radius for 2nd beta model ("):','*',&
                         theta_cx2a,status)
         end if
      end if

      write(*,*) 'Model applies Welch window between primary and secondary&
                  & cutoffs'
      if (do_kpc) then
         cutoff1 = cutoff1*theta_cx*arcsec2kpc
         call io_getd('Primary cutoff (kpc):','*',cutoff1,status)
         cutoff1 = cutoff1/(theta_cx*arcsec2kpc)
         cutoff2 = cutoff2*theta_cx*arcsec2kpc
         call io_getd('Secondary cutoff (kpc):','*',cutoff2,status)
         cutoff2 = cutoff2/(theta_cx*arcsec2kpc)
      else 
         call io_getd('Primary cutoff (units of core radii):','*',&
                      cutoff1,status)
         call io_getd('Secondary cutoff:','*',cutoff2,status)
      end if

   case('b','B','p','P','x','X','m','M')
      call io_getd('Central number density /cm^-3:','*',n0,status)
      call io_getd('Beta parameter:','*',beta,status)
      ellipsoid = io_yesno('Ellipsoidal model (y/n):','*',status)
      if (ellipsoid) then
         write(*,*) 'Axes 1 and 2 on sky; axis 3 along los'
         if (do_kpc) then
            theta_c = theta_c*arcsec2kpc
            call io_getd('Core Radius Axis 1 (kpc):','*',theta_c(1),status)
            call io_getd('Core Radius Axis 2 (kpc):','*',theta_c(2),status)
            call io_getd('Core Radius Axis 3 (kpc):','*',theta_c(3),status)
            theta_c = theta_c/arcsec2kpc
         else
            call io_getd('Core Radius Axis 1 ("):','*',theta_c(1),status)
            call io_getd('Core Radius Axis 2 ("):','*',theta_c(2),status)
            call io_getd('Core Radius Axis 3 ("):','*',theta_c(3),status)
         end if

!         write (*,*) 'Euler angles for ellipsoid'
!         do i=1,3
!            ang(i)=ang(i)/deg2rad
!         end do

!         call io_getd('Euler angle 1:','*',ang(1),status)
!         call io_getd('Euler angle 2:','*',ang(2),status)
!         call io_getd('Euler angle 3:','*',ang(3),status)
!         do i=1,3
!            ang(i)=ang(i)*deg2rad
!         end do

!         el_alpha = ang(1) 
         call io_getd('Angle of ellipsoid on sky:','*',el_alpha,status)
      else
         if (do_kpc) then
            theta_cx = theta_cx*arcsec2kpc
            call io_getd('Core Radius (kpc):','*',theta_cx,status)
            theta_cx = theta_cx/arcsec2kpc
         else
            call io_getd('Core Radius ("):','*',theta_cx,status)
         end if
      end if

      write(*,*) 'Model applies Welch window between primary and secondary&
                  & cutoffs'
      if (do_kpc) then
         cutoff1 = cutoff1*theta_cx*arcsec2kpc
         call io_getd('Primary cutoff (kpc):','*',cutoff1,status)
         cutoff1 = cutoff1/(theta_cx*arcsec2kpc)
         cutoff2 = cutoff2*theta_cx*arcsec2kpc
         call io_getd('Secondary cutoff (kpc):','*',cutoff2,status)
         cutoff2 = cutoff2/(theta_cx*arcsec2kpc)
      else 
         call io_getd('Primary cutoff (units of core radii):','*',&
                      cutoff1,status)
         call io_getd('Secondary cutoff:','*',cutoff2,status)
      end if
      hanning = io_yesno('Apply Hanning (as opposed to Welch) window',&
                         'yes',status)

   case('f','F')
      from_file= .true.     
      call io_geti('Number of n and T points:','5000',num_nt,status)
      call io_getd('step size:','0.01',ntstep,status)
      if (allocated(ndat)) deallocate(ndat)
      allocate(ndat(num_nt))
      if (allocated(Tdat)) deallocate(Tdat)
      allocate(Tdat(num_nt))
      call read_nt(ndat,Tdat,num_nt,dat_read)
      write(*,*) 'Read ',dat_read,' points.'
   case('n', 'N')
      MT200=MT200*1d-14
      call io_getd('Total mass at r200 (10^14 M_sol):', '*', MT200, status)
      MT200=MT200*1d14
      call io_getd('Gas mass fraction at r200:', '*', fg200, status)
      h=H0/100
      write(*,*) 'How to calculate c200:'
      write(*,*) '(N)eto et al 2007 for relaxed clusters'
      write(*,*) '(D)utton and Maccio 2014 for NFW halos'
      write(*,*) 'Enter (m)anually'
      write(*,*)
      call io_getc('c200 type','n',c200_type,status)
      c200_sel: select case(c200_type)
        case('n', 'N') c200_sel
          c200=5.26d0*(((MT200*h)/1.d14)**(-0.1d0))*(1.d0/(1.d0+z))
        case('d', 'D') c200_sel
          a=0.520d0+(0.905d0-0.520d0)*exp(-0.617d0*z**1.21)
          b=-0.101d0+0.026d0*z
          c200=10d0**(a+b*LOG10(MT200*h/1d12))
        case('m', 'M') c200_sel
          call io_getd('c200:', '*', c200, status)
      end select c200_sel
      write(*,*) 'c200 is', c200
      call io_getd('GNFW alpha parameter:', '*', a_GNFW, status)
      call io_getd('GNFW beta parameter:', '*', b_GNFW, status)
      call io_getd('GNFW gamma parameter:', '*', c_GNFW, status)
      call io_getd('GNFW concentration parameter c_500:', '*', c500_GNFW, status)
      ellipsoid = io_yesno('Ellipsoidal model (y/n):','*',status)
      if (ellipsoid) then
         call io_getd('Ratio (minor axis length)/(major axis length):', '*', f_GNFW, status)
         call io_getd('Angle on the sky (degrees, anticlockwise from x-axis):', '*', el_alpha, status)
         el_alpha=el_alpha*deg2rad
      end if
   case('g', 'G')
      call io_getd('GNFW alpha parameter:', '*', a_GNFW, status)
      call io_getd('GNFW beta parameter:', '*', b_GNFW, status)
      call io_getd('GNFW gamma parameter:', '*', c_GNFW, status)
      call io_getd('GNFW concentration parameter c_500:', '*', c500_GNFW, status)
      ellipsoid = io_yesno('Ellipsoidal model (y/n):','*',status)
      if (ellipsoid) then
         call io_getd('Ratio (minor axis length)/(major axis length):', '*', f_GNFW, status)
         call io_getd('Angle on the sky (degrees, anticlockwise from x-axis):', '*', el_alpha, status)
         el_alpha=el_alpha*deg2rad
      end if
      call io_getd('Y_tot (arcmin^2):', '*', Ytot, status)
      if (do_kpc) then
         theta_s = theta_s*arcmin2kpc
         call io_getd('theta_s (kpc):', '*', theta_s, status)
         theta_s = theta_s/arcmin2kpc
      else
         call io_getd('theta_s (arcmin):', '*', theta_s, status)
      end if

   end select

! Electron temperature for beta model 
   select case(model_type)
   case('d','D','b','B')
      Te = Te/keVtoK
      call io_getd('Electron Temperature (keV):','*',Te,status)
      Te = Te*keVtoK
   case('p','P')
      from_poly = .true.
      write(*,*) 'T=a + b * r + c *r^2 +d*r^3'
      call io_getd('Parameter 1','0.',poly_a,status)
      call io_getd('Parameter 2','0.',poly_b,status)
      call io_getd('Parameter 3','0.',poly_c,status)
      call io_getd('Parameter 4','0.',poly_d,status)
   case('h','H')
      T_central = T_central/keVtoK
      call io_getd('Central Temperature (keV):','*',T_central,status)
      T_central = T_central*keVtoK
      T_halo = T_halo/keVtoK
      call io_getd('Halo Temperature(keV):','*',T_halo,status)
      T_halo = T_halo*keVtoK
      if (do_kpc) then
         Tcorerad = Tcorerad*arcsec2kpc
         call io_getd('Halo Core radius (kpc):','*',Tcorerad,status)
         Tcorerad = Tcorerad/arcsec2kpc
      else
         call io_getd('Halo Core radius ("):','*',Tcorerad,status)
      end if
   end select

   write(*,*) 'All displacements relative to pointing centre'
   write(*,*) 'Positive displacement shifts to the west in RA'
   write(*,*) 'Positive displacement shifts to the north in Dec'
   call io_getd('Displacement of Cluster Centre in RA (arc-seconds)&
   &:','*',sz_x,status)
   call io_getd('Displacement of Cluster Centre in Dec (arc-seconds)&
   &:','*',sz_y,status)

   end subroutine get_cluster_model














