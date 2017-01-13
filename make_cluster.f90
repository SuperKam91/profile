!---------------------------------------------------------------------
!     calculates projected sz and x-ray profiles from relevant inputs
!     Inputs: 
!     real beta       beta parameter
!     real theta_cx   core radius in arcsec
!     real n0         central number density in cm^-3
!     real H0         Hubble constant in kms^-1 Mpc^-1
!     real z          redshift 
!     real cellsize   cell size in arcsec
!     integer ncell  number of cells. 
!
!     Outputs:
!     real(*) sz_profile array of temperature decrement
!                          / 2nkT (Thompson cross section)
!                      T0  | -----------------------------   dl
!                          /            mc^2
!
!     real(*) x_profile  array of x-ray counts
!                       / 
!                       | AS n^2 dl 
!                       /
!     where A = 3.279e-63 n^2 D^-2
!     and D is the luminosity distance 
!           D = (2c/H0){(1+z) - (1+z)^1/2} 
!     (assuming q0 = 1/2)
!     and S is the area (in cm2) of one pixel, using the angular size
!     relation
!           theta = (d/D) (1+z)^2
!
!     all pixel values refer to value at centre of the pixel
!
!     the centre of the map is the centre of pixel maxsize/2+1,
!     maxsize/2+1. The central z is maxsize/2
!
subroutine make_cluster

   use kind_def
   use sz_globals
   use maths
   use kgraph
   use cosmology
   use physics
   use Te_map_ams

   implicit none

   real(kind=dp), dimension(maxsize) :: x_profile, n_profile, y_profile
   real(kind=dp), dimension(maxsize) :: m_tot, m_shell, rho, v_tot, v_shell
   real(kind=dp), dimension(maxsize,maxsize) :: n, T, sqrt_T, teemap
   real(kind=dp) :: temp1, temp2, r, rot_temp1, rot_temp2
   real(kind=dp) :: y_value, x_value, n_value
   real(kind=dp) :: theta2, theta_cx2, z1,theta_cx2b
   real(kind=dp) ::  k_x, k_n, k_y, k_rho, k_m
   real(kind=dp) :: rhocrit, rhocrit_b, r_2500, r_500, r_200, v_200
   real(kind=dp) :: m_2500, m_500, m_200, frac, f_b
   real(kind=dp) :: pixsize, temp, dl
   integer :: i, j, i1, j1
   integer :: centre
   logical, external :: io_yesno

   z1 = 1.0D0+z

   if (verbose) then
      write(*,*)
      write(*,*) 'Luminosity distance is : ',D_lum/Mpc,' Mpc'
      write(*,*) 'Angular distance is    : ',D_theta/Mpc,' Mpc'
      write (*,*) '1 Mpc is               : ',1.0/(D_theta/Mpc*sec2rad),&
                  ' arcsec'
   end if
!         pixsize in 1000 km squares: convert to cm^2 later
   pixsize = 1.0D-6*cellsize*sec2rad*D_theta

!         length integral in cm
   dl = pixsize*1.0D8

   pixsize = pixsize**2
   
!     work out constant for xrays: the 1e-14 is 1e-30
!     from Alastair's number and 1e16 for 1000km^2 -> cm^2
!     ^^^^^^^^^^^^^^^ This is so that we don't run out of significant digits
!     too early (!). 
!
! The CONSTANT const
! is count rate per emission integral (cm^-3) n_e^2 dV over 4 pi at cluster 
! redshift
!
! Therefore the counts per pixel are the CONSTANT * D_lum^-2 * n_e^2 * pixel 
! area all integrated along the line of sight.
!
! D_lum in Mpc ; n_e per cm^-3 ; pixel size and dl in cm
!
! All temperature dependance and k correction is pulled into this constant
! and its all done properly (detector response, galatic hydrogen etc. etc.)
! so we don't waste time in this code.

   k_x = (Mpc/D_lum)
   k_x = k_x**2
   k_x = k_x*pixsize*const
   k_x = k_x*1.0D-14*dl*x_exp 

! work out constant for SZ (dl in metres now), 1e4 to go from cm^-2 to m^-2
! This gives a *brightness temperature* defined as I = 2kT/lambda^2 for the 
! SZ effect. So use R-J to convert back to an intensity. 
 
   chan = 1

!   k_sz = -g_nu_x2(dble(nu(chan)))*T0*k_b*Te
!   k_sz = k_sz*dl*1.0D4
!   k_sz = k_sz*sigma_T/(m_e*const_c2)

   k_y = 1.0D4*dl*sigma_T*k_b*Te/(m_e*const_c2)

! work out constant for n; want nsky to be total number of electrons in
! a rectangular column along the line of sight of the pixel size
! nsky intrinsically in cm^-3 so x 10^6 to get m^-3
   k_n = (cellsize*sec2rad*D_theta)**3
   k_n = k_n*1.0D6

! work out constant for mass density, rho. This is the multiplicative factor 
! to convert from the 2d map of electron density n to mass density, so 
! multiply by 1.0D6 to convert into m^-3 from cm^-3 and then multiply by the
! average mass per electron mu_e. Constant k_m is then the total mass in a
! cube of the cellsize at this density
   k_rho = 1.d6*mu_e
   k_m = k_rho*(cellsize*sec2rad*D_theta)**3

! Initialise n Arrays
   select case(model_type)
   case('f','F')
! Note also initialises T
      call make_model_from_file(T,n)
   case default 
      call make_model_beta(n)
   end select
   write(*,*)
   write(*,*) 'n map calculated'
   if (do_plot) then
      write(*,*) '2d map of n'
      call display_map(n,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
   end if

! Determine where r_2500, r_500 and r_200 lie and then estimate total
! baryon mass within these radii
! NB this assumes a constant gas mass fraction through the cluster

   rhocrit = rhocrit0/one_over_h(z)**2
   f_b = omegab/omegam
   rhocrit_b = rhocrit*f_b

   if (verbose) then
      write(*,*)
      write(*,*) 'Rho_crit_0       ',rhocrit0
      write(*,*) 'Rho_crit_z       ',rhocrit
      write(*,*) 'Baryon fraction  ',f_b
      write(*,*) 
! r^3 = V/(4/3pi); V=M/rho
      write(*,*) 'r_200  (M=10^15) ',(((1.d15*m_0)/(200.d0*rhocrit*&
                                     4.d0/3.d0*pi))**(1.d0/3.d0))/kpc,&
                 ' kpc'
      write(*,*) 'T_sing (M=10^15) ',1.d0/(k_b*keVtoK)*&
                                 (1.d15*m_0*G_c)**(2.d0/3.d0)*&
                                 (H0*1.d3/Mpc/one_over_h(z))**(2.d0/3.d0)*&
                                 (mu_n/2.d0)*(200.d0/2.d0)**(1.d0/3.d0),&
                 ' keV'

   end if

! Initialise variables
   r_2500 = 0.0d0
   r_500 = 0.0d0
   r_200 = 0.0d0
   m_2500 = 0.0d0
   m_500 = 0.0d0
   m_200 = 0.0d0
   v_200 = 0.0d0
   m_tot = 0.0d0
   m_shell = 0.0d0
   v_tot = 0.0d0
   v_shell = 0.0d0

! Size of a pixel in m
   pixsize = cellsize*sec2rad*D_theta
   do i = 1, maxsize
      if (i.eq.1) then
         m_shell(1) = 4.d0/3.d0*pi*0.5d0**3*n(i,1)*k_m
         m_tot(1) = m_shell(1)
         v_shell(1) = 3.d0*G_c*m_shell(1)**2/(5.d0*(pixsize/2.d0))/f_b
         v_tot(1) = v_shell(1)
      else
         m_shell(i) = 4./3.*pi*((i-0.5)**3-(i-1.5)**3)*n(i,1)*k_m
         m_tot(i) = m_tot(i-1)+m_shell(i)
         v_shell(i) = G_c*m_shell(i)*(m_tot(i-1)+m_tot(i))/2.d0/&
                      ((i-0.5)*pixsize)/f_b
         v_tot(i) = v_tot(i-1)+v_shell(i)
      end if
      rho(i) = m_tot(i)/(4./3.*pi*((i-0.5)**3*(cellsize*sec2rad*D_theta)**3))
   end do

   do i = 1, maxsize-1
   
      if (rho(i).gt.2500.*rhocrit_b) then
         if (rho(i+1).lt.2500.*rhocrit_b) then
            frac = ((rho(i)-2500.*rhocrit_b)/(rho(i)-rho(i+1)))
            r_2500 = i-1+frac**0.5
            m_2500 = m_2500+frac*m_shell(i)
         else
            m_2500 = m_2500+m_shell(i)
         end if
      end if
      if (rho(i).gt.500.*rhocrit_b) then
         if  (rho(i+1).lt.500.*rhocrit_b) then
            frac = ((rho(i)-500.*rhocrit_b)/(rho(i)-rho(i+1)))
            r_500 = i-1+frac**0.5
            m_500 = m_500+frac*m_shell(i)
         else
            m_500 = m_500+m_shell(i)
         end if
      end if
      if (rho(i).gt.200.*rhocrit_b) then
         if (rho(i+1).lt.200.*rhocrit_b) then
            frac = ((rho(i)-200.*rhocrit_b)/(rho(i)-rho(i+1)))
            r_200 = i-1+frac**0.5
            m_200 = m_200+frac*m_shell(i)
            v_200 = v_200+frac**0.5*v_shell(i)
         else
            m_200 = m_200+m_shell(i)
            v_200 = v_200+v_shell(i)
         end if
      end if
   end do

   if (verbose) then
! Masses and radii
      write(*,*) 
      write(*,*) 'R_2500  (kpc)   ',r_2500*cellsize*sec2rad*D_theta/kpc
      write(*,*) 'M_2500  (M_o)   ',m_2500/m_0
      write(*,*) 'R_500   (kpc)   ',r_500*cellsize*sec2rad*D_theta/kpc
      write(*,*) 'M_500   (M_o)   ',m_500/m_0
      write(*,*) 'R_200   (kpc)   ',r_200*cellsize*sec2rad*D_theta/kpc
      write(*,*) 'M_200   (M_o)   ',m_200/m_0
      write(*,*) 'R_core  (kpc)   ',theta_cx*sec2rad*D_theta/kpc 
      write(*,*) 'M_total (M_o)   ',m_tot(maxsize)/m_0
      write(*,*) 
      write(*,*) '5R_500  (kpc)   ',5.d0*r_500*cellsize*sec2rad*D_theta/kpc
      write(*,*) 'Cutoff  (kpc)   ',cutoff2*theta_cx*sec2rad*D_theta/kpc
      write(*,*) 

! Potentials and virial temperatures 
      write(*,*)
! (-ve) gravitation potential energy out to r_200
      write(*,*) 'V_200           (J)  ',v_200
      write(*,*) 'V_tot           (J)  ',v_tot(maxsize)
! gravitational potential energy of the baryons for a uniform sphere
! 3/5 GM^2/r (but one M is total and one is baryon) 
      write(*,*) '3/5 G M_b M_T/r (J)  ',3.d0*G_c*m_200**2/&
                          (f_b*5.d0*r_200*cellsize*sec2rad*D_theta)
! Virial temperature using 3/2 nkT = V/2. n is m_200/0.6m_p
      write(*,*) 'V_200/3nk       (K)  ',v_200*mu_n/(3.d0*m_200*k_b)
      write(*,*) 'V_200/3nk       (keV)',v_200*mu_n/(3.d0*m_200*k_b*keVtoK)
! Virial temperature for a singular isothermal sphere from Voit 2005
      write(*,*) 'GmuM/2rk        (K)  ',G_c*m_200/f_b*mu_n/&
                                   (2.d0*r_200*cellsize*sec2rad*D_theta*k_b)
      write(*,*) 'GmuM/2rk        (keV)',G_c*m_200/f_b*mu_n/&
                               (2.d0*r_200*cellsize*sec2rad*D_theta*k_b*keVtoK)
      write(*,*) 'Ratio                ',&
                                  v_200*mu_n/(3.d0*m_200*k_b*keVtoK)/&
                                  (G_c*m_200/f_b*mu_n/&
                              (2.d0*r_200*cellsize*sec2rad*D_theta*k_b*keVtoK))
!      write(*,*) 'T(M=10^15M_0)   ',G_c*m_p*0.6*(1.d15*m_0)**(2./3.)/&
!                                    (2.d0*k_b*keVtoK*&
!                                   (3.d0/(rhocrit0*4.d0*pi*200))**(1.d0/3.d0))
      write(*,*)
    end if

   if (do_plot) then
      call pgpage
      call display_profile(m_shell/m_0,maxsize,cellsize*sec2rad*D_theta/kpc)
      write(*,*) 'Mass in spherical shell against radius'
      call pgpage
      call display_profile(rho,maxsize,cellsize*sec2rad*D_theta/kpc)
      write(*,*) 'Mean density within sphere against radius'
      call pgpage
      call display_profile(n(:,1),maxsize,cellsize*sec2rad*D_theta/kpc)
      write(*,*) 'Local density against radius'
      call pgpage
      call display_profile(v_tot/1.d50,maxsize,cellsize*sec2rad*D_theta/kpc)
      write(*,*) 'GPE against radius'
   end if


! Initialise T Arrays (T normalised to T_central)
   select case(model_type)
   case('x','X')
      call make_Te_map(T)
   case('m','M')
      T = 1.
      call two_d_tempmap(teemap)
      teemap = teemap/Te
      if (io_yesno('Plot temperature map?','n',status)) then
         call display_map(teemap,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
      end if
      T = teemap
! Use T and sqrt_T arrays defined in get_cluster_model
   case('h','H')
      call make_model_temp_halo(T)
   case('p','P')
      call make_model_from_poly(T)
   case('b','B','d','D')
      T = 1.0
   end select
   sqrt_T = sqrt(T)
   write(*,*) 'T map calculated'
   if (do_plot) then
      write(*,*) '2d map of T'
      call display_map(T,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
   end if

! Do integration

! initialise profiles. x_profile is for X-rays; y_profile for Comptonisation 
! parameter; n_profile for mass
   y_profile = 0.
   x_profile = 0.
   n_profile = 0.

! integrate columns. Make sure that column 1 in y-direction is only 
! counted once.
   do i = 1, maxsize
      y_profile(i) = T(1,i)*n(1,i)
      n_profile(i) = n(1,i)

! Added by WFG as occasionally T(i,j) = 0 if not all the data gets loaded...
      if (T(1,i) > 0.0) then 
         x_profile(i) = T(1,i)**(-0.3)*n(1,i)**2
      else
         x_profile(i) = 0.0
      end if

      do j = 2,maxsize
         y_profile(i) = y_profile(i)+2.0*T(j,i)*n(j,i)

! Added by WFG as occasionally T(i,j) = 0 if not all the data gets loaded...
         if (T(j,i) > 0.0) then 
            x_profile(i) = x_profile(i)+2.0*(T(j,i)**(-0.3))*(n(j,i)**2)
         end if
         n_profile(i) = n_profile(i)+2.0*n(j,i)
      end do
   end do

   x_profile = x_profile*k_x
   y_profile = y_profile*k_y
   n_profile = n_profile*k_n

   centre = maxsize/2+1
   do j = 1,maxsize
      do i = 1, maxsize
         temp1 = float(i)-centre
         temp1 = temp1-(sz_x/cellsize)
         temp2 = float(j)-centre
         temp2 = temp2-(sz_y/cellsize)
         if (ellipsoid) then
            rot_temp1 = temp1*cos(el_alpha)-temp2*sin(el_alpha)
            rot_temp2 = temp2*cos(el_alpha)+temp1*sin(el_alpha)
            r = theta_cx*sqrt((rot_temp1/theta_c(1))**2+&
                 & (rot_temp2/theta_c(2))**2)+1.0
            call interp(y_profile,maxsize,r,y_value)
            call interp(x_profile,maxsize,r,x_value)
            call interp(n_profile,maxsize,r,n_value)
            y_value = y_value*theta_c(3)/theta_cx
            x_value = x_value*theta_c(3)/theta_cx
            n_value = n_value*theta_c(3)/theta_cx
         else
            r = sqrt(temp1**2+temp2**2)+1.0
            call interp(y_profile,maxsize,r,y_value)
            call interp(x_profile,maxsize,r,x_value)
            call interp(n_profile,maxsize,r,n_value)
         end if

         compton_y(i,j) = y_value
         xsky(i,j) = x_value
         nsky(i,j) = n_value
      end do
   end do

! Now work out brightness temperatures of SZ effect at each channel
   do chan = 1,nchan
      if (overwrite) then
         szsky(:,:,chan) = compton_y*g_nu_x2(dble(nu(chan)))*T0
      else
         szsky(:,:,chan) = szsky(:,:,chan)+ &
                           compton_y*g_nu_x2(dble(nu(chan)))*T0
      end if
   end do

   chan = nchan/2+1

   if (verbose) then
      write(*,*)
      write(*,*) 'Central y value                       ',maxval(compton_y)
      write(*,*) 'Central decrement (ie deltaT)        ',&
                 minval(szsky(:,:,chan)), ' K'
      write(*,*) 'Y_tot (arcmin^2)                      ',sum(compton_y)*&
                 cellsize**2/3600.d0
      write(*,*) 'Maximum of X-ray emission             ',maxval(xsky)
      write(*,*) 'Total number of electrons in cluster  ',&
                 sum(nsky)
      write(*,*) 'Total cluster baryon mass             ',&
                 sum(nsky)*mu_e/m_0,' M_o' 
      write(*,*)
   end if

   if (do_plot) then
      write(*,*) 'Displaying SZ sky (values of deltaT)'
      call display_map(szsky(:,:,chan),maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
   end if

   if (do_plot) then
      write(*,*) 'Displaying X-ray sky'
      call display_map(xsky,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
   end if

end subroutine make_cluster















