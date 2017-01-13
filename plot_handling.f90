module plot_handling

   use kind_def
   use sz_globals

   implicit none

   public :: plot_cut
   public :: plot_max_amplitude
   public :: plot_visibility
   public :: plot_spectra
   public :: plot_cosmology

contains

! **************************************************************************

   subroutine plot_cut(map1,mapsize)

      use kind_def
      use sz_globals
      use maths
      use kgraph
      use map_handling

      implicit none

      character*80 :: plot_device
      integer, intent(in) :: mapsize
      real(kind=dp), dimension(mapsize, mapsize), intent(in) :: map1
      real,dimension(4) :: map_line
      integer, parameter :: resol = 200
      real :: pix_x,pix_y,cent_x,cent_y
      real(kind=dp), dimension(resol) :: profile
      character*1 :: cursor
      real(kind=dp) :: xstep, ystep, stepsize, xpt, ypt
      integer :: i, indev
      logical, external :: io_yesno
      integer, external :: pgopen
      logical :: sep_plot

      sep_plot = io_yesno("Plot cut in seperate window","yes",status) 
      if (sep_plot) then
         call io_getc('Input plot device:','/xwindow',plot_device,status)
      end if

      call display_map(map1,mapsize,graph_adj,&
           cellsize,cellsize,0D0,0D0)         

!      pix_x = pix_xd
!      pix_y = pix_yd
!      cent_x = cent_xd
!      cent_y = cent_yd

      cursor = '  '

      write(*,*) 'Draw cut line'

      do while (cursor .ne. 'A')
         call pgband(0,1,0,0,map_line(1), map_line(2), cursor)
         call pgband(1,0,map_line(1), map_line(2), map_line(3),&
              map_line(4), cursor)
      end do

! Prepare to plot cut
      if (sep_plot) then
         indev = pgopen(plot_device)
         call pgslct(indev)
         call pgsci(plot_col)
         call pgsls(line_sty)
         call pgslw(line_width)
         call pgsch(char_hgt)
      else
         call pgpage
      end if

! change from world co-ords to array positions

      pix_x = cellsize
      pix_y = cellsize
      cent_x = 0.
      cent_y = 0.
      map_line(1) = (map_line(1)-cent_x)/pix_x+mapsize/2+1
      map_line(2) = (map_line(2)-cent_y)/pix_y+mapsize/2+1
      map_line(3) = (map_line(3)-cent_x)/pix_x+mapsize/2+1
      map_line(4) = (map_line(4)-cent_y)/pix_y+mapsize/2+1

      xstep = (map_line(3)-map_line(1))/(resol-1.)
      ystep = (map_line(4)-map_line(2))/(resol-1.)

      stepsize = sqrt(xstep**2+ystep**2)*cellsize
      do i = 1,resol
         xpt = map_line(1)+(i-1.)*xstep
         ypt = map_line(2)+(i-1.)*ystep
         call interp(map1,mapsize,mapsize,xpt,ypt,profile(i))
      end do

      call display_profile(profile,resol,stepsize)

      if (sep_plot) then
         call pgslct(indev)
         call pgclos
         call pgslct(pl_dev1)
      end if

   end subroutine plot_cut

! **************************************************************************

!  Finds and plots the maximum amplitude along circular tracks of the aperture
!  plane
   subroutine plot_max_amplitude(re_map,im_map,mapsize,cell)

      use kind_def
      use sz_globals
      use kgraph
      use cosmology

      implicit none
      integer, intent(in) :: mapsize
      real(kind=dp), dimension(mapsize,mapsize), intent(in) :: re_map, im_map
      real(kind=dp), intent(in) :: cell
      integer, parameter :: resol=440 !***changed from original program***
      integer, parameter :: npoints = 50
      real (kind=dp) :: upt, vpt, psi, inc, rad_inc
      real (kind=dp), dimension(npoints) :: re_pt, im_pt, amp
      real (kind=dp), dimension(resol) :: max_amp, mean_amp, rad, min_amp
      real (kind=dp) :: rzerofl, izerofl
      real, dimension(resol) :: xpoints, ypoints
      integer :: i, j, iunit
      character :: char1
      character*80 :: outfile


      telescope : select case(which_tel)
      case('c','C','m','M')
         rad_inc = 150./dble(resol)
      case('e','E')
         rad_inc = 300./dble(resol)
      case('t','T')
         rad_inc = 550./dble(resol)
      case('d','D')
         call io_getd('Maximum lambda:','10000.',rad_inc,status)
         rad_inc = rad_inc/dble(resol)
      case default
         rad_inc = 2750./dble(resol)
      end select telescope

      inc = pi/npoints
      do j = 1, resol
         rad(j) = j*rad_inc
         psi = 0.
         do i = 1, npoints
            upt = rad(j)*sin(psi)
            vpt = rad(j)*cos(psi)
            call extract_visibility&
                 & (re_map,im_map,mapsize,upt,vpt,cell,re_pt(i),im_pt(i))
            psi = psi+inc
         end do
         amp = sqrt(re_pt**2+im_pt**2)
         max_amp(j) = maxval(amp)
         mean_amp(j) = sum(amp)/npoints
         min_amp(j) = minval(amp)
      end do

      upt = 0.0
      vpt = 0.0
      call extract_visibility(re_map,im_map,mapsize,upt,vpt,cell,rzerofl,izerofl)
      write (*,*)'Zerospacing flux is ',rzerofl,izerofl

      status = 0
      char1 = 'n'

      call io_getc('Dump the data to a file?','n',char1,status)
      if ((char1=='y').or.(char1=='Y')) then 
! ***
         write(*,*) '**********************************************'
         write(*,*) 'this version of profile has been hacked so'
!         write(*,*) 'that it only writes a file of baseline against'
!         write(*,*) 'MAXIMUM amplitude'
         write(*,*) 'that it writes a file of baseline against'
         write(*,*) 'min, max and mean amplitude'
         write(*,*) '**********************************************'
! ***   
         status = 0
         call io_nxtlun(iunit,status)
         call io_getc('Filename:','',outfile,status)
         open (iunit,file=outfile,form = 'formatted')
         write(iunit,*) '# Baseline min_amplitude max_amplitude mean_amplitude'
         write(iunit,*) 0., abs(rzerofl), abs(rzerofl), abs(rzerofl)
         do i=1, resol
            write(iunit,*) rad(i), min_amp(i), max_amp(i), mean_amp(i)
         end do
         close(iunit)
      end if

      if (autoscale) then
         write(*,*) 'Displaying maximum amplitudes'
         call pgpage
         call display_profile(rad,max_amp,resol,.true.)
         call pglabel('Observing baseline (\gl)','Predicted flux (&
              &Jy)','Maximum amplitude')
         call pgpage
         write(*,*) 'Displaying mean amplitudes'
         call display_profile(rad,mean_amp,resol,.true.)
         call pglabel('Observing baseline (\gl)','Predicted flux (&
              &Jy)','Mean amplitude')
      else
         call pgpage
         call display_profile(rad,max_amp,resol,prof_x1,prof_x2,prof_y1,prof_y2)
         call pglabel('Observing baseline (\gl)','Predicted flux (&
              &Jy)','Maximum amplitude')
         call pgpage
         write(*,*) 'Displaying mean amplitudes'
         call display_profile(rad,mean_amp,resol,prof_x1,prof_x2,prof_y1,prof_y2)
         call pglabel('Observing baseline (\gl)','Predicted flux (&
              &Jy)','Mean amplitude')
         call pglabel('Observing baseline (\gl)','Predicted flux (&
              &Jy)','Flux against baseline for clusters at z = 0.2, 0.5')
      end if

   end subroutine plot_max_amplitude

! **************************************************************************

! plots an elliptical track in an aperture plane
   subroutine plot_visibility(re_map,im_map,mapsize,rad,cell,sindec)

      use kind_def
      use sz_globals
      use kgraph

      implicit none
      integer, intent(in) :: mapsize
      real(kind=dp), dimension(mapsize,mapsize), intent(in) :: re_map, im_map
      real(kind=dp), intent(in) :: rad,cell,sindec
      integer, parameter :: npoints = 200
      real(kind=dp) :: upt,vpt,psi,inc, temp, ndata
      real(kind=dp) :: max_amp,mean_real
      real(kind=dp), dimension(npoints) :: re_pt, im_pt, amp, ph
      integer :: i
      character(1) :: chr1

      psi = 0.
      ndata = 0
      mean_real = 0.0
      inc = pi/npoints
      do i = 1,npoints
         upt = rad*sin(psi) 
         vpt = rad*cos(psi)*sindec
         call extract_visibility(re_map,im_map,mapsize,upt,vpt,cell,&
              re_pt(i),im_pt(i))
         telescope : select case(which_tel)
         case ('c','C') 
            if (upt**2+vpt**2.lt.57.0**2) then
               re_pt(i) = 0.
               im_pt(i) = 0.
            else
               mean_real = mean_real+re_pt(i)
               ndata = ndata+1
            end if
            psi = psi+inc
         case ('r','R','p','P')
            if (upt**2+vpt**2.lt.640.0**2) then
               re_pt(i) = 0.
               im_pt(i) = 0.
            else
               mean_real = mean_real+re_pt(i)
               ndata = ndata+1
            end if
            psi = psi+inc
         case default
            mean_real = mean_real+re_pt(i)
            ndata = ndata+1
            psi = psi+inc
         end select telescope
      end do

      inc = inc*12.0/pi

      max_amp = maxval((re_pt**2+im_pt**2))
      max_amp = sqrt(max_amp)

      write(*,*) 'Maximum amplitude round track is ',max_amp
      write(*,*) 'Mean real part = ', mean_real/ndata

      call io_getc('Plot real/imag or amp/phase (r/a):','r',chr1,status)

      if (chr1.eq.'r') then
         write(*,*) 'Plotting real part of visibility'
         if (autoscale) then
            call display_profile(re_pt,npoints,inc)

            call pgpage

            write(*,*) 'Plotting imaginary part of visibility'
            call display_profile(im_pt,npoints,inc)
         else
            call display_profile(re_pt,npoints,inc,prof_x1,prof_x2,&
                                 prof_y1,prof_y2)

            call pgpage

            write(*,*) 'Plotting imaginary part of visibility'
            call display_profile(im_pt,npoints,inc,prof_x1,prof_x2,&
                                 prof_y1,prof_y2)
         end if
      else
         amp = sqrt(re_pt**2+im_pt**2)
         ph = atan2(im_pt,re_pt)
         write(*,*) 'Plotting amplitude of visibility'
         if (autoscale) then
            call display_profile(amp,npoints,inc)

            call pgpage

            write(*,*) 'Plotting phase of visibility'
            call display_profile(ph,npoints,inc)
         else
            write(*,*) 'Plotting amplitude of visibility'
            call display_profile(amp,npoints,inc,prof_x1,&
                                 prof_x2,prof_y1,prof_y2)

            call pgpage

            write(*,*) 'Plotting phase of visibility'
            call display_profile(ph,npoints,inc,prof_x1,&
                                 prof_x2,prof_y1,prof_y2)
         end if
      end if

   end subroutine plot_visibility

! **************************************************************************

   subroutine plot_spectra

      use kind_def
      use sz_globals
      use kgraph
      use physics

      implicit none

      character*80 :: title, x_axis, y_axis, ps1  
      character*1 :: which_spec
      integer, parameter :: resol = 700
      integer :: i, j, k
      real(kind=dp) :: nu1, nu1_start, nu1_inc 
      real(kind=dp) :: clust_y, x, tau, beta_z, theta_e
      real(kind=dp), dimension(resol) :: freq, I_nu1, T_B, B_nu1, Del_T, T_B2
      real(kind=dp), dimension(resol) :: delta_T, flu_spc, diff_0_20, diff_0_10
      real(kind=dp), dimension(resol) :: rel_I_nu1_5,rel_I_nu1_10, rel_I_nu1_20
      real(kind=dp), dimension(resol) :: kinetic_sz, diff_0_5, I_nu1_over_nu1
      real(kind=dp), dimension(resol) :: ratio
      real, dimension(resol) :: x_points, y_points

      write(*,*) 'Possible spectra for plotting :-'
      write(*,*) '(P)lanck spectrum'   
      write(*,*) '(F)luctuation spectrum'   
      write(*,*) '(T)hermodynamic temperature of thermal SZ effect'
      write(*,*) '(B)rightness temperature of thermal SZ effect'
      write(*,*) '(I)ntensity of thermal SZ effect'
      write(*,*) '(R)elativistic thermal SZ intensity spectrum'
      write(*,*) '(K)inetic SZ effect'
      write(*,*) '(D)ifferentiating SZ spectra '
      write(*,*)
      call io_getc('Select one','p',&
           which_spec,status)     

      plot_col = 1
      line_sty = 1
      line_width = 4
      char_hgt = 1.0

      call pgslct(pl_dev1)
      call pgsci(plot_col)
      call pgsls(line_sty)
      call pgslw(line_width)
      call pgsch(char_hgt)

      clust_y = 2.0d-4
      beta_z = 1000.0d3/const_c
      Te = 10*keVtoK
      theta_e = Te*k_b/(m_e*const_c2)
      tau = clust_y/theta_e

      nu1_start = 0.0d9
      nu1_inc = 1.0d9

      do i = 1, resol
         nu1 = nu1_start+i*nu1_inc
         x = const_h*nu1/(k_b*T0)
         freq(i) = nu1
         B_nu1(i) = planck_fn(T0,nu1)
         delta_T(i) = planck_fn(2.74d0,nu1)-planck_fn(2.7399d0,nu1)
         flu_spc(i) = (x*exp(x)/(exp(x)-1.0))*planck_fn(T0,nu1)*&
                      1.0d-4/T0
         flu_spc(i) = dB_by_dT(T0,nu1)
         kinetic_sz(i) = (x*exp(x)/(exp(x)-1.0))*planck_fn(T0,nu1)*&
                      tau*beta_z
         Del_T(i) = (x*coth(x/2.0)-4.0)*clust_y*T0
         I_nu1(i) = (B_nu1(i)*x*exp(x)/(exp(x)-1.0))*(x*coth(x/2.0)-4.0)*&
                    clust_y
         I_nu1_over_nu1(i) = I_nu1(i)/nu1
         Te = 5*keVtoK
         rel_I_nu1_5(i) = clust_y*rel_corr_sz(nu1)
         Te = 10*keVtoK
         rel_I_nu1_10(i) = clust_y*rel_corr_sz(nu1)
         ratio(i) = I_nu1(i)/rel_I_nu1_10(i)
         Te = 20*keVtoK
         rel_I_nu1_20(i) = clust_y*rel_corr_sz(nu1)
         diff_0_20(i) = (0.938*I_nu1(i)-rel_I_nu1_20(i))
         diff_0_10(i) = (0.968*I_nu1(i)-rel_I_nu1_10(i))
         diff_0_5(i) = (0.984*I_nu1(i)-rel_I_nu1_5(i))
         T_B(i) = I_nu1(i)/(2.0*k_b*nu1**2/const_c2)
      end do

      write(*,*) 'Plots for cluster with y-parameter ',clust_y
      write(*,*) 'Beta_z ',beta_z
      write(*,*) 'Theta_e', theta_e
      write(*,*) 'tau', tau

      select case(which_spec)
      case('p','P')

         call pgpage
         call display_profile(freq,B_nu1,resol,.true.)
         title = 'Planck spectrum for T=2.73 K'
         x_axis = 'Frequency (Hz)'
         y_axis = 'B(\gn) (W m\u-2\dHz\u-1\dsr\u-1\d)'
         call pglab(x_axis,y_axis,title)   

      case('x','X')

! This is a brute force method for calculating the spectrum of a small
! temperature fluctuation - use "f" instead
         call pgpage
         call display_profile(freq,delta_T,resol,.true.)
         title = 'Spectrum of intensity change for \gDT=100\gmK'
         x_axis = 'Frequency (Hz)'
         y_axis = 'I(\gn) = \gDT(\gdB/\gdT) (W m\u-2\dHz\u-1\dsr\u-1\d)'
         call pglab(x_axis,y_axis,title)  
         write(*,*) 'Maximum value ',maxval(delta_T)

      case('f','F')

         call pgpage
         call display_profile(freq,flu_spc,resol,.true.)
         title = '\gdB/\gdT(\gn) for \gDT=100\gmK'
         x_axis = 'Frequency (Hz)'
         y_axis = 'T(\gdB/\gdT)(\gn) (W m\u-2\dHz\u-1\dsr\u-1\d)'
         call pglab(x_axis,y_axis,title)   
         write(*,*) 'Maximum value ',maxval(flu_spc)

      case('t','T')

         call pgpage
         title = 'Thermodynamic temperature of SZ effect'
         x_axis = 'Frequency (Hz)'
         y_axis = '\gDT\dSZ\u = I\dSZ\u/(\gdB)/(\gdT) (K)'
         call display_profile(freq,Del_T,resol,.true.)
         call pglab(x_axis,y_axis,title)   

      case('i','I')

         call pgpage
         title = 'Surface brightness of SZ effect'
         x_axis = 'Frequency (Hz)'
         y_axis = 'I\dSZ\u(\gn)/\gn (W m\u-2\dHz\u-1\dsr\u-1\d)'
         call display_profile(freq,I_nu1_over_nu1,resol,.true.)
         call pglab(x_axis,y_axis,title)  

      case('r','R')

         call pgpage
         title = 'Relativistic correction to thermal SZ effect'
         x_axis = 'Frequency (Hz)'
         y_axis = 'I\dSZ\u(\gn) (W m\u-2\dHz\u-1\dsr\u-1\d)'
         call display_profile(freq,I_nu1,resol,.true.)
         call pglab(x_axis,y_axis,title)  
         ps1 = "0 keV"

         call pgptext(100.0e9,12.0e-22,0.0,0.0,ps1)
         x_points = freq
         y_points = rel_I_nu1_5
         call pgsci(2)
         call pgline(resol,x_points,y_points)
         ps1 = "5 keV"

         call pgptext(100.0e9,10.0e-22,0.0,0.0,ps1)
         y_points = rel_I_nu1_10
         call pgsci(3)
         call pgline(resol,x_points,y_points)
         ps1 = "10 keV"

         call pgptext(100.0e9,8.0e-22,0.0,0.0,ps1)
         y_points = rel_I_nu1_20
         call pgsci(5)
         call pgline(resol,x_points,y_points)
         ps1 = "20 keV"

         call pgptext(100.0e9,6.0e-22,0.0,0.0,ps1)
         call pgsci(1)

      case('c','C')

         call pgpage
         title = 'Relativistically corrected surface brightness of SZ effect'
         x_axis = 'Frequency (Hz)'
         y_axis = 'I\dSZ\u(\gn) (W m\u-2\dHz\u-1\dsr\u-1\d)'
         call display_profile(freq,I_nu1,resol,.true.)
         x_points = freq
         y_points = rel_I_nu1_20
         call pgline(resol,x_points,y_points)
         call pglab(x_axis,y_axis,title)   

      case('d','D')

         call pgpage
         title = 'Discriminating between SZ components'
         x_axis = 'Frequency (Hz)'
         y_axis = '\gDI\dSZ\u(\gn) (W m\u-2\dHz\u-1\dsr\u-1\d)'
         call display_profile(freq,diff_0_20,resol,.true.)
         call pglab(x_axis,y_axis,title)   
         ps1 = "0.938I\d\gn\u(0 keV)-I\d\gn\u(20 keV)"
         call pgptext(3.0e10,-3.7e-22,0.0,0.0,ps1)

         call pgsci(2)
         x_points = freq
         y_points = diff_0_10
         call pgline(resol,x_points,y_points)
         ps1 = "0.968I\d\gn\u(0 keV)-I\d\gn\u(10 keV)"
         call pgptext(3.0e10,-4.4e-22,0.0,0.0,ps1)

         call pgsci(3)
         x_points = freq
         y_points = I_nu1/10.0
         call pgline(resol,x_points,y_points)
         ps1 = "I\d\gn\u(0 keV)/10"
         call pgptext(3.0e10,-5.8e-22,0.0,0.0,ps1)

         call pgsci(6)
         x_points = freq
         y_points = diff_0_5
         call pgline(resol,x_points,y_points)
         ps1 = "0.984I\d\gn\u(0 keV)-I\d\gn\u(5 keV)"
         call pgptext(3.0e10,-5.1e-22,0.0,0.0,ps1)

         call pgsci(5)
         y_points = kinetic_sz
         call pgline(resol,x_points,y_points)
         ps1 = "I\dKinetic\u (=93\gmK)"
         call pgptext(3.0e10,-6.5e-22,0.0,0.0,ps1)
         call pgsci(1)

      case('b','B')

         call pgpage
         title = 'Brightness Temperature of SZ effect'
         x_axis = 'Frequency (Hz)'
         y_axis = 'T\dB\u = I\d\gn\u\gl\u2\d/2k\dB\u (K)'
         call display_profile(freq,T_B,resol,.true.)
!   x_points = freq
!   y_points = T_B2*0.9
         call pglab(x_axis,y_axis,title)   

      case('k','K')

         call pgpage
         title = 'Spectrum of Kinetic S-Z effect'
         x_axis = 'Frequency (Hz)'
         y_axis = 'I\dk-SZ\u(\gn) (W m\u-2\dHz\u-1\dsr\u-1\d)'
         call display_profile(freq,kinetic_sz,resol,0.0d0,7.0d11,-3.0d-22,5.0d-22)
         call pglab(x_axis,y_axis,title)   
         ps1 = "I\dKinetic\u"
         call pgptext(6.0e11,4.0e-22,0.0,0.0,ps1)
         call pgsci(2)
         x_points = freq
         y_points = I_nu1/10
         call pgline(resol,x_points,y_points)
         ps1 = "I\dThermal\u/10"
         call pgptext(6.0e11,3.5e-22,0.0,0.0,ps1)
         write(*,*) 'Maximum value ',maxval(kinetic_sz)

      case('o','O')

! This is not particularly useful or informative so not included in the options
! list
         call pgpage
         title = 'Ratio of non-rel SZ to rel SZ'
         x_axis = 'Frequency (Hz)'
         y_axis = 'I\dk-SZ\u(\gn) (W m\u-2\dHz\u-1\dsr\u-1\d)'
         call display_profile(freq,ratio,resol,0.0d0,7.0d11,0.0d0,1.5d0)
         call pglab(x_axis,y_axis,title)   

      end select

   end subroutine plot_spectra

! **************************************************************************

! Produce plots for cosmology
! 06/08/09 First attempt [kg]

   subroutine plot_cosmology

      use kind_def
      use sz_globals
      use kgraph
      use cosmology

      implicit none

      character*80 :: title, x_axis, y_axis, ps1  
      character*8 :: chr_wm
      integer, parameter :: resol = 60
      real(kind=dp), dimension(resol) :: y_points, x_points
      real, dimension(resol) :: x1, y1, y2
      real(kind=dp) :: d_z, old_z, D_th, old_q0
      integer :: i

      d_z = 0.05
      old_z = z
      old_q0 = q0
      do i = 1, resol
         z = i*d_z
         x_points(i) = z
         x1(i) = z

! Use E.Komatsu's routines instead of Will G's for chosen cosmology
         call angdist3
         y_points(i) = D_theta

! Compare with Will's calculation
!         call angdist
!         y1(i) = D_theta

! Plot Lambda = 0, Omega_M = 1 and 0 for comparison
         q0 = 0.0
         call angdist(z,H0,D_th)
         y1(i) = D_th
         q0 = 0.5
         call angdist(z,H0,D_th)
         y2(i) = D_th

      end do

      call pgpage
      call display_profile(x_points,y_points,resol,.true.)
      title = ""
      x_axis = 'Redshift'
      y_axis = 'Angular diameter distance (MPc)'
      call pglab(x_axis,y_axis,title)   
      call pgsci(2)
      call pgline(resol,x1,y1)
      call pgsci(3)
      call pgline(resol,x1,y2)
      call pgsci(1)

      write(*,*) 'Hubble constant /km s^-1 Mpc^-1:   ', H0
      write(*,*) 'Omega_M:                           ', omegam
      write(*,*) 'Omega_L:                           ', omegal

      z = old_z
      q0 = old_q0

   end subroutine plot_cosmology


! **************************************************************************

end module plot_handling







