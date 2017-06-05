program profile

   use kind_def
   use sz_globals
   use helptext
   use file_handling
   use map_handling
   use plot_handling
   use kgraph
   use fits_handling
   use bayesian
   use batchbayesian
   use cosmology
   use likelihood
   use physics 
   use maths
   use astronomy
   use define_telescope
   use define_xray_telescope
   use w_proj
   use primary_beam

   implicit none

   real(kind=dp), dimension(:,:), allocatable :: smooth_map
   real(kind=dp), dimension(:,:), allocatable :: off_beam, tempmap
   real(kind=dp), dimension(:,:), allocatable :: array1, ssp, r
   real(kind=dp), dimension(:), allocatable :: std, xbar
   real(kind=dp), dimension(2) :: fr

   complex(kind=dp), dimension(:,:), allocatable :: ctempmap

   real :: mid_beam, sigma, offset, offrad2
   real :: synth_beam, bmaj, bmin, tr(6),maxx, maxy
   real :: xin, yin, s_30, al, slim, flux, src, temp, maxf
   real :: rescale, xback, datmin, datmax, tfac

   character(len=fname_len) :: pb_filename, fits_file, kernel, plist
   character(len=fname_len) :: plot_device, filename, simfile, fitsname
   character(len=fname_len) ::  corrfile
   character*1 :: chr_epoch, chr1
   character*64 :: line

   real(kind=dp) :: e, ave, sumoffset, x, dummy_dp, likely_sum, likely_sum2
   real(kind=dp) :: chisq_sum, chisq_sum2, true_back, norm, xsum
   real(kind=dp) :: x_off, y_off, smoo_scale, fits_freq, upt, vpt
   real(kind=dp) :: bandwidth, f_step, f_start
   real(kind=dp) :: wmax, baselev, ra, dec
   real(kind=dp) :: cutoff, redat, imdat, uvcell
   real(kind=dp) :: sindec, incell
   real(kind=dp) :: rad, chisq, likely, pa, var, cent_fwhm
   real(kind=dp) :: sigma1, sigma2, lmax, psf_max
   real(kind=dp) :: ra1,dec1,ra2,dec2,dx,dy

   integer :: blc_x, blc_y, dim, adjust
   integer :: noofbins
   integer :: write_time(3)
   integer :: i, j, k, idum, n_iter
   integer :: ncell, num, mapsize
   integer :: command, Ecommand, noofp, ichan, nbin
   integer :: noofradii, npix, unit, l_dir
   integer :: nwplanes, oversample, ifail
   integer :: icmd, iout, idev
   integer :: errorcode

   logical :: stat_fl, ex, fill_single, exact_spec
   logical :: report_stat, do_poisson, do_regrd

   integer, external :: pgopen, chr_lenb
   logical, external :: io_yesno
   external G02BAF

! Write greeting
   write(*,*) 'Profile     Ver 3.1.10 - 20/06/13'
   write(*,*) '(type explain-profile for help)'
   write(*,*)

!  Initialise directories
   call init_directory

! Initialise variables

! Flags
   plot_open = .false.
   do_plot = .false.
   overwrite = .false.
   verbose = .true.
   graph_adj = .false.
   autoscale = .true.
   wis = .false.

   fits_file = '0016-2.fits'
   write_fits_beam_table = .false.

! Define initial telescope
   which_tel = 'd'
   l_dir  = chr_lenb(telescope_dir)
   filename = telescope_dir(1:l_dir)//'amidc_sa.tel'
   call read_telescope(filename,errorcode)
   if (errorcode.ne.0) then
      write(*,*) 'NB no telescope defined'
   else
      write(*,*) 'Default telescope parameters read from '//trim(filename)
   end if
   tel_name = 'SA'
   chan = nchan/2+1
   call calc_basel
   cent_fwhm = const_c/(obsfreq*dish_diameter*sec2rad)
   which_pb = 'a'
   do chan = 1, nchan
      fwhm(chan) = (101.1/(nu(chan)*1d-9) + 1.89)*2.*sqrt(2.*log(2.))*60.
   end do
   n_pol = 1

! Initial map sizes
   maxsize = 512
   ncell = maxsize/2

! Default map
   flag_nd = 2
   map_type = 'y'
   p_sky => compton_y

! cellsize in arcseconds
   cellsize = 15.0D0
   incell = 15.0D0

! observation parameters
   pb_x = 0.D0
   pb_y = 0.D0
   sz_x = 0.D0
   sz_y = 0.D0
   sz_z = 0.D0
   ph_x = 0.D0
   ph_y = 0.D0
   obs_raref = 12.D0
   obs_decref = 30.D0
   n_samp = 0
   ha_step = 200.d0

! model parameters
   model_type = 'b'
   n0 = 1.D-2
   beta = 6.5D-1
   cutoff1 = 10.0d0
   cutoff2 = 12.0d0
   hanning = .true.
   Te = 7.8D7
   T_central = 1.0D0
   T_halo = 1.0D0
   Tcorerad = 150.D0
   x_exp = 44530.D0
   const = 3.279D-33
   back = 0.0D0
   theta_cx = 60.0D0
   theta_c = 60.0D0
   el_alpha = 0.0D0

! New GNFW model for initial model
   model_type = 'g'
   Ytot = 3.d-3
   theta_s = 7.5d0
! 'Universal' GNFW parameters (Planck)
   a_GNFW = 1.0510d0
   b_GNFW = 5.4905d0
   c_GNFW = 0.3081d0
   c500_GNFW = 1.177d0
   f_GNFW = 1.0d0
! Einasto DM shape parameter
   aEin_DM = 0.20d0
!   temp_halo = .false.
   ellipsoid = .true.
   vary_beta = .true.
   from_file = .false.
   from_poly = .false.
!  DM-GNFW defaults
   MT200=1.d15
   fg200=0.12

! Initialise variables for graphics operations
   prof_x1 = 0.0
   prof_x2 = 2750.0
   prof_y1 = 0.0
   prof_y2 = 0.01
   plot_col = 1
   line_sty = 1
   line_width = 6
   char_hgt = 2.0
   plot_device = '/xwindow'

! Initialise variables for specific commands
   n_iter = 10
   rescale = 1.0
   smoo_scale = 300.0

! Define cosmology
   z = 1.71D-1
   H0 = 70.0
   omegam = 0.3
   omegal = 0.7
   omegab = 0.044
   w_de = -1.d0
   call lumdist
   D_lum = D_lum*Mpc         ! D_lum in metres now....
   call angdist
   D_theta = D_theta*Mpc
   rhocrit0 = (3.0*(H0*1.D3/Mpc)**2)/(8.0*pi*G_c)

   status = 0

! Initialise io library
   call io_initio

! Generate a random seed for the random number generator based upon the 
! current time and make an initial call to generator
   if (io_yesno('Initialise random number generator based on clock','yes',&
                status)) then
      call itime(write_time)
      idum = -(write_time(1)+1)*(write_time(2)+1)*(write_time(3)+1)
      dummy_dp = ran2(idum)
   else
      idum = 1
      call io_geti('Random seed to use:', '*', idum, status)
   end if

! Allocate space to all common arrays
   call do_allocation

! Define x-ray telescope
   x_tel = 'p'
   call make_xray_telescope

! Open \xwindow plot device
   if (io_yesno('Open plot device','yes',status)) then
      pl_dev1 = pgopen(plot_device)
      call pgslct(pl_dev1)
      call pgsci(plot_col)
      call pgsls(line_sty)
      call pgslw(line_width)
      call pgsch(char_hgt)
      plot_open = .true.
   else
      plot_open = .false.
   end if

! Set up table with angular diameter distances
   call setup_da

! Main loop; user enters a command and appropriate subroutine then called
   main_loop : do

! Reset chan since this often goes to nchan+1
      chan = 1

! Reset status
      status = 0

! Get next command
      call io_getcmd('[Profile]>', comms, ncomm, command, status)

! 'Quit' or 'Exit'
      if ((command.eq.1).or.(command.eq.58)) then
         call io_getc('Are you sure? (y/n)','n',chr1,status)
         if (chr1=='y') then 
            exit main_loop
         end if

! 'Read-fits-map'
      else if (command.eq.2) then

         call select_map

! Get filename
         prompt = 'FITS filename'
         which_dir = profile_data_dir
         rootname = '0016-2'
         extn = '.fits'
         call get_read_filename(fits_file)
         call read_fits_map(fits_file)

         if (flag_nd.eq.3) then
           xsum=sum(p3_sky(:,:,1))
         else
           xsum = sum(p_sky)
         endif
         write(*,*) 'Total flux in map is ',xsum

! 'Display-map'
      else if (command.eq.3) then
         if (plot_open) then
            call select_map
            if (flag_nd.eq.2) then
               call display_map(p_sky,maxsize,graph_adj,&
                    cellsize,cellsize,0D0,0D0)         
            elseif (flag_nd.eq.3) then
               do chan = 1, nchan
                  p_sky = p3_sky(:,:,chan)
                  call display_map(p_sky,maxsize,graph_adj,&
                       cellsize,cellsize,0D0,0D0)
               end do
            end if
         else
            write(*,*) 'No plot device open'
         end if

! 'Make-aperture'
      else if (command.eq.4) then
         gsig = 0.
         call io_getc('Do plots (y/n)','n',chr1,status)
         if ((chr1=='y').or.(chr1=='Y')) then
            if (plot_open) then
               do_plot = .true.
            else
               write(*,*) 'No plot device open'
               do_plot = .false.
            end if
         else
            do_plot = .false.
         end if
         call make_aperture

! 'Difference-maps'
      else if (command.eq.5) then
         report_stat = io_yesno("Report fit statistics (assuming Poisson)",&
              "yes",status)
         if (report_stat) then
            write(*,*) "Select data map"
         end if
         call select_map
         if (flag_nd.eq.2) then
            temp_sky = p_sky
         else
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if
         if (report_stat) then
            write(*,*) "Select model map"
         end if
         call select_map
         if (flag_nd.ne.2) then
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if
         call io_getc('Show only fit region (y/n):','n',chr1,status)
         if (chr1.eq.'y') then               
            call read_flagpixel(ok_data,maxsize)
         else
            ok_data = .true.
         end if
         call diff_map(temp_sky,p_sky,chisq,likely,report_stat)

! 'Copy-map'
      else if (command.eq.6) then
         write(*,*) 'Select map to copy'
         call select_map
         if (flag_nd.eq.2) then
            temp_sky = p_sky
         else
            write(*,*) 'Using channel 1'
            chan = 1
            temp_sky = p3_sky(:,:,chan)
         end if
         write(*,*) 'Select map to overwrite'
         call select_map
         if (flag_nd.eq.2) then
            p_sky = temp_sky
         elseif (flag_nd.eq.3) then
            do chan = 1, nchan
               p3_sky(:,:,chan) = temp_sky
            end do
         end if

! 'Make-cluster'
      else if (command.eq.7) then
         call angdist

         call io_getc('Do plots (y/n):','n',chr1,status)
         if ((chr1=='y').or.(chr1=='Y')) then
            if (plot_open) then
               do_plot = .true.
            else
               write(*,*) 'No plot device open'
               do_plot = .false.
            end if
         else
            do_plot = .false.
         end if

         call io_getc('Overwrite previous sky (y/n):','y',chr1,status)
         if ((chr1=='y').or.(chr1=='Y')) then
            overwrite = .true.
         else
            overwrite = .false.
         end if

         select case(model_type)
           case('g', 'G', 'n', 'N')
             call make_cluster_GNFW
           case default
             D_theta = D_theta * Mpc
             call make_cluster
         end select

         write(*,*)
         write(*,*) 'Skies made'

! 'Help'
      else if (command.eq.8) then
         do i = 1,ncomm+nsep
            if (index1(i).eq.0) then
               write(*,*)'--------------------'
            else if (index1(i).ne.99) then           
               call io_wrout(comms(index1(i)))
            end if
         end do

! 'Get-cluster-model'
      else if (command.eq.9) then
         call get_cluster_model

! 'plot-map-profile'
      else if (command.eq.10) then
         call select_map
         if (flag_nd.eq.2) then
            call make_map_profile(p_sky,maxsize,map_profile)
            if (autoscale) then
               call display_profile(map_profile,maxsize/2,cellsize)
            else
               call display_profile(map_profile,maxsize/2,cellsize,&
                    prof_x1, prof_x2, prof_y1, prof_y2)
            end if
         else
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if

! 'Histogram-map'
      else if (command.eq.11) then 
         call select_map
         if (flag_nd.eq.2) then
            datmin = minval(p_sky)
            datmax = maxval(p_sky)
            call io_geti('Number of bins','1000',nbin,status)
            call map_histogram(p_sky,maxsize,datmin,datmax,nbin)
         else
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if
 
! 'Make-model-xray-map'
      else if (command.eq.12) then 
         do_poisson = io_yesno('Make model x-ray map with poisson noise',&
                               'yes',status)
         call io_getd('Background level:','*',back,status)
         call make_model_xray_map(idum,do_poisson)
         write(*,*) 'Model x-ray map made'

! Flag-region
      else if (command.eq.13) then
         if (plot_open) then
            ok_data = .true.
            call select_map
            if (io_yesno('Read in a flagged region file','no',status))then
               call read_flagpixel(ok_data,maxsize)
               if (flag_nd.eq.2) then
                  where (.not. ok_data)
                     p_sky = 0.0
                  endwhere
               elseif (flag_nd.eq.3) then 
                  do chan = 1, nchan
                     p_sky => p3_sky(:,:,chan)
                     where (.not. ok_data)
                        p_sky = 0.0
                     endwhere
                  end do
               end if
            end if
            if (flag_nd.eq.2) then
               call remsource(p_sky,ok_data,maxsize)
               where (.not. ok_data)
                  p_sky = 0.0
               endwhere
            elseif (flag_nd.eq.3) then
               chan = 1
               p_sky => p3_sky(:,:,chan)
               call remsource(p_sky,ok_data,maxsize)
               do chan = 1, nchan
                  p_sky => p3_sky(:,:,chan)
                  where (.not. ok_data)
                     p_sky = 0.0
                  endwhere
               end do
            end if
            if (io_yesno('Write flagged region to file','no',status))then
               call write_flagpixel(ok_data,maxsize)
            end if
         else
            write(*,*) 'No plot device open'
         end if

! 'get-directory'
      else if (command.eq.14) then
         call get_directory

! 'Sum-maps' 
      else if (command.eq.15) then
         call select_map
         if (flag_nd.eq.2) then
            temp_sky = p_sky
         else
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if
         call select_map
         if (flag_nd.ne.2) then
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if
         if (plot_open) then
            call display_map(p_sky+temp_sky,maxsize,graph_adj,&
                 cellsize,cellsize,0D0,0D0)         
         end if
         if (io_yesno('Overwrite (2nd) sky with sum','no',status)) then
            p_sky = p_sky+temp_sky
         end if

! 'Open-plot-device'
      else if (command.eq.16) then
         call io_getc('Plot device:','/xwindow',plot_device,status)
         pl_dev1 = pgopen(plot_device)
         call pgslct(pl_dev1)
         call pgsci(plot_col)
         call pgsls(line_sty)
         call pgslw(line_width)
         call pgsch(char_hgt)
         plot_open = .true.
         if (plot_device.ne.'/xwindow') then
            graph_adj = .false.
         end if

! 'Close-plot-device'
      else if (command.eq.17) then
         call pgend
         plot_open = .false.

! 'Rescale-map'
      else if (command.eq.18) then
         call select_map
         call io_getr('Rescale factor','*',rescale,status)
         if (flag_nd.eq.2) then
            p_sky = p_sky*rescale
         elseif (flag_nd.eq.3) then
            p3_sky = p3_sky*rescale
         end if

! 'Smooth-map' - convolves with a gaussian
      else if (command.eq.19) then
         call io_getd('Smoothing scale (arcsec)','*',smoo_scale,status)
         do_plot = .false.
         call select_map
         if (flag_nd.eq.2) then
            call smoo_map(p_sky,maxsize,smoo_scale)
         elseif (flag_nd.eq.3) then 
            do chan = 1, nchan
               p_sky => p3_sky(:,:,chan)
               call smoo_map(p_sky,maxsize,smoo_scale)
            end do
         end if

! 'Plot-model-visibility'
      else if (command.eq.20) then
         if (plot_open) then
            call io_getd('Radius in wavelengths:','870',rad,status)
            telescope : select case(which_tel)
            case ('r','R','p','P')
               call io_getd('Declination:','46',sindec,status)
               sindec = sin(sindec*deg2rad)
            case default
               sindec = 1.0
            end select telescope
            call pgpage
            chan = 1
            if (nchan.gt.1) then
               call io_geti('Which channel','*',chan,status)
            end if
            call plot_visibility(re(:,:,chan),im(:,:,chan),maxsize,&
                 rad,cellsize,sindec)
         else
            write(*,*) 'No plot device open'
         end if

! 'Make-source-population
      else if (command.eq.21) then
         call make_source_population(idum)

! 'Plot-spectra'
      else if (command.eq.22) then
         call plot_spectra

! 'Get-cosmology'
      else if (command.eq.23) then
         call get_cosmology

! 'Estimate-map-background'
      else if (command.eq.24) then
         if (plot_open) then
            if (flag_nd.eq.2) then
               call select_map
               call estimate_map_background(p_sky)
            else
               write(*,*) "Command not defined for multi channel maps"
               cycle main_loop
            end if
         else
            write(*,*) 'No plot device open'
         end if

! 'Bayesian x-ray fit'
      else if (command.eq.25) then
         call bayesian_x_fit

! 'Graphics-operations'
      else if (command.eq.26) then
         if (plot_open) then
            call graphics_ops
         else
            write(*,*) 'No plot device open'            
         end if

! 'make-test-sky'
      else if (command.eq.27) then
         call select_map
         call make_test_sky

! 'Make-pointed-observation'
      else if (command.eq.28) then
         call make_pointed_observation(idum)

! 'Plot-likelihood'
      else if (command.eq.29) then
         call plot_likelihood

! 'display-sz-sky'
      else if (command.eq.30) then
         if (plot_open) then
            if (nchan.gt.1) then
               call io_geti('Which channel?','1',chan,status)
               write(*,*) 'Displaying image at',nu(chan)/1e9,'GHz'
               xsum = sum(szsky(:,:,chan))
            else
               chan = nchan
               xsum = sum(szsky(:,:,nchan))
            end if
            call display_map(szsky(:,:,chan),maxsize,graph_adj,&
                 cellsize,cellsize,0D0,0D0)

            write(*,*) 'Sum of brightness temperature pixels in map is ',&
                 xsum
            write(*,*) 'Total flux density in map is ',&
                 xsum*(cellsize*sec2rad)**2*2.0*k_b*nu(chan)**2/const_c2*1.0D26,' Jy'
         else
            write(*,*) 'No plot device open'
         end if

! 'Difference-map-profiles'
      else if (command.eq.31) then
         call select_map
         if (flag_nd.eq.2) then
            call make_map_profile(p_sky,maxsize,map_profile)
         else
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if
         call select_map
         if (flag_nd.eq.2) then
            call make_map_profile(p_sky,maxsize,map_profile2)
         else
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if
         if (autoscale) then
            call display_profile(map_profile-map_profile2,maxsize/2,cellsize)
         else
            call display_profile(map_profile-map_profile2,maxsize/2,cellsize,&
                 prof_x1, prof_x2, prof_y1, prof_y2)
         end if

! 'Plot-maximum-amplitude'
      else if (command.eq.32) then
         if (plot_open) then
            chan = 1
            if (nchan.gt.1) then
               call io_geti('Which channel','*',chan,status)
            end if
            call plot_max_amplitude(re(:,:,chan),im(:,:,chan),&
                 maxsize,cellsize)
         else
            write(*,*) 'No plot device open'
         end if

! 'Get-geometry'
      else if (command.eq.33) then
         call get_geometry

! 'Clear'
      else if (command.eq.34) then
         call pgpage

! 'Simple-clean'
      else if (command.eq.35) then
         call simple_clean

! 'Show-cluster-parameters'
      else if (command.eq.36) then
         call show_cluster_parms

! 'Read-uv-fits'
      else if (command.eq.37) then
         call read_uv_fits

! 'plot-cut'
      else if (command.eq.38) then
         call select_map
         if (flag_nd.eq.2) then
            call plot_cut(p_sky,maxsize)
         else
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if

! 'Mcg-xray-fit'
      else if (command.eq.39) then
         call io_getd('Background level:','*',back,status)
         call mcg_x_fit

! 'Combine-maps'
      else if (command.eq.40) then
         do_plot = io_yesno('Do plots','n',status)
         do_regrd = io_yesno('Regrid maps', 'y', status)
         if (io_yesno('Input map coordinates in galactic coordinates',&
                      'yes',status)) then
            call combine_maps_gal(do_regrd)
         else
            call combine_maps(do_regrd)
         end if

! 'Write-model'
      else if (command.eq.41) then
         select case(model_type)
         case('f','F')
            write(*,*) "Can not write for this type of model"
         case default
            call write_model
         end select

! 'Read-model'
      else if (command.eq.42) then
         call read_model

! 'Explain-profile'
      else if (command.eq.43) then
         do i = 1,15
            write(*,*) ExplainString(i)
         end do

! 'Explain-command'
      else if (command.eq.44) then

! only explain one command at a time....
         write (*,*) 'One one command can be explained at a time....'
         call io_getcmd('[Profile.Explain]>', comms, ncomm, Ecommand, status)
         write(*,*) 'Command: '//comms(Ecommand)
         write(*,*)
         do i = 1,ExplainLength(Ecommand)
            write(*,*) ExplainStrings(Ecommand,i)
         end do
         write(*,*)

! 'Stat-map'
      else if (command.eq.45) then
         call select_map
         if (flag_nd.ne.2) then
            write(*,*) 'Not available for multi-channel arrays'
         else
            call stat_map(p_sky,maxsize)
         end if

! 'Test-cmd'
      else if (command.eq.46) then
         prompt = 'Filename for ASCII beam pattern'
         extn = '.dat'
         which_dir = profile_data_dir
         rootname = 'testBeamPattern'
         call get_read_filename(filename)
         call read_ascii_beam(filename)

! 'Estimate-source-confusion'
      else if (command.eq.47) then
         call estimate_source_confusion(idum)

! 'Lensing'
      else if (command.eq.48) then
         call lensing(z)

! 'Adaptive-smooth'
      else if (command.eq.49) then

         call select_map 
! Allocate space to smooth map
         allocate(smooth_map(maxsize,maxsize))
         smooth_map = 0.0
         call io_getd('min conv size? (arcsec)','5.',mingauss,status)
         call io_getd('max conv size? (arcsec)','60.',maxgauss,status)
         mingauss = mingauss/cellsize
         maxgauss = maxgauss/cellsize
         call io_getc('Write original with smoothed map (y/n):','y',&
              & chr1,status)

         if (flag_nd.eq.2) then
            call adaptive_smooth(p_sky,smooth_map,maxsize)
            call display_map(p_sky,maxsize,graph_adj,&
                 cellsize,cellsize,0D0,0D0)
            if ((chr1=='y').or.(chr1=='Y')) then
               p_sky = smooth_map
            end if
         elseif (flag_nd.eq.3) then
            do chan = 1, nchan
               p_sky = p3_sky(:,:,chan)
               call adaptive_smooth(p_sky,smooth_map,maxsize)
               call display_map(p_sky,maxsize,graph_adj,&
                    cellsize,cellsize,0D0,0D0)
               if ((chr1=='y').or.(chr1=='Y')) then
                  p_sky = smooth_map
               end if
            end do
         end if
         deallocate(smooth_map)

! 'Fill-visibilities'
      else if (command.eq.50) then
         if (nchan.gt.1) then
            fill_single = io_yesno&
                 ('Draw visibilities all from one channel','no',status)
            if (fill_single) then
               call io_geti('Which channel ','1',chan,status)
            end if
            call fill_visibilities(idum,fill_single)
         else
            call fill_visibilities(idum,fill_single)
         end if

! 'Project-cluster'
      else if (command.eq.51) then
         call project_cluster

! 'Run'
      else if (command.eq.52) then

!  Get command file
         prompt = 'Filename for command file:'
         which_dir = profile_scripts_dir
         rootname = 'sparcs_sim'
         extn = '.cmd'
         call get_read_filename(filename)

         call run_cmd(filename,status)

! 'Read-sz-sky'
      else if (command.eq.53) then
         call read_sz_sky

! 'Plot-source-confusion'
      else if (command == 54) then 
         call plot_source_confusion

! Command 55 not currently in use
      else if (command == 55) then

! Command 56 not currently in use
      else if (command == 56) then

! Command 57 not currently in use
      else if (command.eq.57) then

! 'Multiple bayesian x-ray fit'
      else if (command.eq.59) then
         call multi_bayesian_x_fit

! Not currently in use
      else if (command.eq.60) then

! 'Monte-Carlo-X-fit'
      else if (command.eq.61) then

         true_back = back
         likely_sum = 0.0
         likely_sum2 = 0.0
         chisq_sum = 0.0
         chisq_sum2 = 0.0

         call io_getc('Show only fit region (y/n):','n',chr1,status)
         if (chr1.eq.'y') then               
            call read_flagpixel(ok_data,maxsize)
         else
            ok_data = .true.
         end if

         call io_geti('Number of iterations','*',n_iter,status)
         do i = 1, n_iter

! Make a mock rosat map
            back = true_back
            call make_model_xray_map(idum,.true.)
            rosat_map = xmap
            back = true_back 
            call make_model_xray_map(idum,.false.)
            call diff_map(rosat_map,xmap,chisq,likely,.false.)
            write(*,*) 'Iteration ',i,' chisq ',chisq,' likely ',likely
            likely_sum = likely_sum+likely
            likely_sum2 = likely_sum2+likely**2
            chisq_sum = chisq_sum+chisq
            chisq_sum2 = chisq_sum2+chisq**2
         end do
         chisq_sum = chisq_sum/n_iter
         chisq_sum2 = sqrt(chisq_sum2/n_iter-chisq_sum**2)
         likely_sum = likely_sum/n_iter
         likely_sum2 = sqrt(likely_sum2/n_iter-likely_sum**2)
         write(*,*) 'Mean Chisq ',chisq_sum,' +/- ',chisq_sum2
         write(*,*) 'Mean Likely ',likely_sum,' +/- ',likely_sum2

! Command 62 not in use
      else if (command.eq.62) then        

! 'Make-mosaic-observation'
      else if (command.eq.63) then
         u = 0.d0
         v = 0.d0
         w = 0.d0
         data_re = 0.d0
         data_im = 0.d0
         rms = 0.d0
         weight = 0.d0
         which_pointing = 0  
!         call io_getc('Write u,v AND w?','n',chr1,tatus)
!         if (chr1=='y') wis = .true.
!         if (chr1=='n') wis = .false.
         call make_mosaic_observation(idum)

! Command 64 not in use
      else if (command.eq.64) then

! 'Get-telescope'
      else if (command.eq.65) then
         call get_telescope

! Not currently in use
      else if (command.eq.66) then 

! 'weight-my-data':
      elseif (command.eq.67) then 
         offset = 20.27*60./cellsize
         sigma = 0.6*3600./(sqrt(2.*log(2.)))
         sigma = sigma/cellsize
         allocate(off_beam(maxsize,maxsize))
         do i = 1,maxsize
            do j = 1,maxsize
               offrad2 = (real(i-(maxsize/2+1))**2+real(j-(maxsize/2+1+offset))**2)
               off_beam(i,j) = dble(exp(-1.*offrad2/(2.*sigma**2)))
            end do
         end do
         mid_beam = off_beam(maxsize/2+1,maxsize/2+1)
         off_beam = off_beam/mid_beam
         call display_map(off_beam,maxsize,graph_adj,&
              cellsize,cellsize,0D0,0D0)
         do i = 1,maxsize
            do j = 1,maxsize
               if (off_beam(i,j).le.0.3) off_beam(i,j)=0.3
            end do
         end do
         szsky(:,:,nchan) = szsky(:,:,nchan)/off_beam
         call display_map(szsky(:,:,nchan),maxsize,graph_adj,&
              cellsize,cellsize,0D0,0D0)
         deallocate(off_beam)
         write(*,*) 'Use make-aperture to apply correct PB'
         pbcor = .true.


! Make-Spectral-Index-Map:
      elseif (command.eq.68) then
         call select_map
         if (flag_nd.eq.2) then
            temp_sky = p_sky
            if (associated(p_sky,target=szsky(:,:,chan))) then
               fr(1) = nu(chan)/1.d9
            else
               fr(1) = fr(1)/1.d9
            end if
            call io_getd('Frequency of map (GHz)','*',fr(1),status)
            fr(1) = fr(1)*1.d9
            call select_map
            if (associated(p_sky,target=szsky(:,:,chan))) then
               fr(2) = nu(chan)/1.d9
            else
               fr(2) = fr(2)/1.d9
            end if
            call io_getd('Frequency of map (GHz)','*',fr(2),status)
            fr(2) = fr(2)*1.d9
            call make_spectral_index_map(maxsize,temp_sky,p_sky,fr)
         else if (flag_nd.eq.3) then
            call make_spectral_index_map(maxsize,nchan,p3_sky,nu)
         end if

! Overlay mosaic:
      elseif (command.eq.69) then 
         call overlay_mosaic

! Pointing offset:
      elseif (command.eq.70) then
         call pointing_power

! evaluate SZ spectrum at a user defined frequency
      elseif (command.eq.71) then
         f_start = f_start/1.d9
         call io_getd('Frequency (GHz)','*',f_start,status)
         f_start = f_start*1.d9
         write(*,*)
         write(*,*) 'G_nu/x^2       = ',g_nu_x2(f_start)
         write(*,*) 'G_nu_thermo    = ',g_nu_thermo(f_start)
         write(*,*) 'Planck fn / RJ = ',&
              planck_fn(T0,f_start)/(T0*rayleigh_jeans_fn(f_start))
         write(*,*) 'DT RJ / thermo = ',&
              g_nu_x2(f_start)/g_nu_thermo(f_start)
         write(*,*) 'Rel_corr ratio = ',&
              rel_corr_sz(f_start)/&
              (g_nu_x2(f_start)*(2*k_b*T0)*(f_start/const_c)**2),&
              '         for T_e(K) = ',Te
         write(*,*)

! Remove SKA ptsrcs:

!---------------------------------------------------------------
! This is a project specific task... careful what you use it on.
!---------------------------------------------------------------

      elseif (command.eq.72) then
         write(*,*) ' '
         write(*,*) '------------------------------------------------------------------'
         write(*,*) ' This is a project specific task... be careful what you use it on.'
         write(*,*) '------------------------------------------------------------------'
         write(*,*) ' '
         call io_getc('Ptsrc list:','ps_list.dat',plist,status)
         call io_getr('Source subtraction limit (uJy):','86.30',slim,status)
         slim = slim*1.E-6 ! convert to Jy
         write(*,*) ' '
         write(*,*) 'Working...'
         write(*,*) ' '
         open(unit=72,file=plist)
         allocate(tempmap(maxsize,maxsize))
         adjust = ((maxsize/2)+1)-((mapsize/2)+1)
!         obsfreq = const_c/obslambda
         tempmap = 0.
         temp = 0.
777      read(72,*,end=401) xin,yin,s_30,al ! from michael peel's files

         yin=yin+real(adjust)+1. ! account for the fact that the map is now 2^n
         xin=xin+real(adjust)+1.
         s_30 = s_30*1.e-3 ! fluxes are in mJy at 30GHz -> convert to Jy
         flux = s_30*(obsfreq/(30.*1.e9))**(al) ! convert flux at 30 to flux at obsfreq
! account for the effect of multiple small sources:
         tempmap(int(yin),int(xin)) = tempmap(int(yin),int(xin)) + flux

         if (tempmap(int(yin),int(xin)).ge.(slim)) then
            src = tempmap(int(yin),int(xin))
! x,y -> y,x for F90/C convention change:
            szsky(int(yin),int(xin),nchan) = szsky(int(yin),int(xin),nchan) - src
            if (src.gt.temp) then
               maxx = yin
               maxy = xin
               maxf = src
               temp = src
            end if
            tempmap(int(yin),int(xin)) = 0.
         end if

         goto 777

401      continue
         close(unit = 72)
         write(*,*) ' '
         write(*,*) '          ...Done.'
         write(*,*) 'All pixels with flux greater than ',slim,' have been source subtracted.'
         write(*,*) ' '
         call io_geti('blc x-coord of submap:','600',blc_x,status)
         call io_geti('blc y-coord of submap:','700',blc_y,status)
         blc_x = blc_x+adjust
         blc_y = blc_y+adjust
         call io_geti('Dimension of submap:','1024',dim,status)
         deallocate(tempmap)
         allocate(tempmap(dim,dim))
         tempmap = 0.d0
         do i = 1,dim
            do j = 1,dim
               tempmap(i,j) = szsky(blc_x+(i-1),blc_y+(j-1),nchan)
            end do
         end do
! converted to nominal brightness temperature
! this gets converted back in make_aperture
         tempmap = tempmap/((cellsize*sec2rad)**2*2.0*k_b*obsfreq**2/const_c2*1.0D26)
         deallocate(szsky)
         maxsize = dim
         allocate(szsky(maxsize,maxsize,nchan))
         szsky = 0.d0
         szsky(:,:,nchan) = tempmap
         deallocate(tempmap)

! Not currently in use
      elseif (command.eq.73) then

! plot-scaling which plots z-r-T scaling relation
      elseif (command.eq.74) then
         call plot_scaling

! make-cmb-sky:
      elseif (command.eq.75) then
         exact_spec = io_yesno('Make CMB sky with exact Cl spectrum',&
              'no',status)
         overwrite = io_yesno('Overwrite previous sky model',&
              'no',status)
         do_plot = io_yesno('Do plots','n',status)
         call make_cmb_sky(idum,exact_spec)

! apply w-projection correction:
      elseif (command.eq.76) then

! grid visibilities into uv-space:
      elseif (command.eq.77) then
         write(*,*) '--------------------------------'
         write(*,*) 'Convolution kernel choices:     '
         write(*,*) '              SINC - like Mapper'
         write(*,*) '              PSWF - like AIPS  '
         write(*,*) '      W-PROJECTION - like CONRAD'
         write(*,*) '--------------------------------'
         call io_getc('Convolution kernel ((n)one, (s)inc, (p)swf, (w)-proj:','n',chr1,status)
         if (chr1=='n') kernel = 'none'
         if (chr1=='s') kernel = 'sinc'
         if (chr1=='p') kernel = 'pswf'
         if (chr1=='w') kernel = 'wpro'
         chan = 1
         call grid_vis(kernel)

! 'read-visibilities'
      elseif (command.eq.78) then
         u = 0.
         v = 0.
         w = 0.
         data_re = 0.
         data_im = 0.
         rms = 0.
         weight = 0.
         call read_vis_data

! write-uv-fits:
      elseif (command.eq.79) then
         prompt = 'FITS filename'
         which_dir = profile_data_dir
         rootname = 'uv'
         extn = '.fits'
         call get_write_filename(filename)

         ra = obs_raref*hr2rad
         call io_getra('RA','*',ra,status)
         obs_raref = ra/hr2rad
         obs_rarot = obs_raref

         dec = obs_decref*deg2rad
         call io_getdec('Declination','*',dec,status)
         obs_decref = dec/deg2rad
         obs_decrot = obs_decref

         call io_getc('Observation name:','*',obs_name,status)
         chr_epoch='J'
         call io_getc('Epoch [B/J]:','*',chr_epoch,status)
         if ((chr_epoch=='b').or.(chr_epoch=='B')) obs_epoch = 1
         if ((chr_epoch=='j').or.(chr_epoch=='J')) obs_epoch = 2
         call write_uv_fits(filename)
         write(*,*) 'FITS file written'

! map-visibilities:
      elseif (command.eq.80) then
         call map_vis(kernel)

! not currently in use
      elseif (command.eq.81) then

! subtract-point-sources
      elseif (command.eq.82) then
         call subtract_sources(bmaj,bmin)
         call io_getc('Subtract baselevel?','y',chr1,status)
         baselev = minval(szsky(:,:,chan))
         if (chr1=='y') szsky(:,:,chan) = szsky(:,:,chan) - baselev

! 'write-fits-map'
      elseif (command.eq.83) then
         call select_map
         if (flag_nd.eq.2) then

! Get filename for map
            prompt = 'FITS filename'
            extn = '.fits'
            which_dir = profile_data_dir
            rootname = 'map'
            call get_write_filename(filename)

            ra = obs_raref*hr2rad
            call io_getra('RA','*',ra,status)
            obs_raref = ra/hr2rad
            obs_rarot = obs_raref
            dec = obs_decref*deg2rad
            call io_getdec('Declination','*',dec,status)
            obs_decref = dec/deg2rad
            obs_decrot = obs_decref
            
            call write_fits_map(p_sky,maxsize,maxsize,cellsize,filename)
         else
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if

! 'correlate-maps'
      elseif (command.eq.84) then
         call select_map
         if (flag_nd.ne.2) then
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if
         temp_sky = p_sky

         call select_map
         if (flag_nd.ne.2) then
            write(*,*) "Command not defined for multi channel maps"
            cycle main_loop
         end if

         ok_data = .true.
         if (io_yesno('Define region for fit','no',status)) then
            call display_map(temp_sky,maxsize,graph_adj,&
                 cellsize,cellsize,0D0,0D0)                     
            call def_map_region
         end if

! Allocate working arrays
         npix = count(ok_data)
         allocate(array1(npix,2))
         allocate(xbar(2))
         allocate(std(2))
         allocate(ssp(2,npix))
         allocate(r(2,2))

         k = 0
         do i = 1,maxsize
            do j = 1,maxsize
               if (ok_data(i,j)) then
                  k = k+1
                  array1(k,1) = temp_sky(i,j)
                  array1(k,2) = p_sky(i,j)
               end if
            end do
         end do

         call pgenv(minval(real(array1(1:k,1))),maxval(real(array1(1:k,1))),&
              minval(real(array1(1:k,2))),maxval(real(array1(1:k,2))),0,1)
         call pgsci(2)
         call pgpt(k,real(array1(:,1)),real(array1(:,2)),17)
         call pgsci(1)

         ifail = 0
         call G02BAF(npix,2,array1,npix,xbar,std,ssp,2,r,2,ifail)
         write(*,*) 'Pearson co-efficient:     | Map 1 | Map 2|'
         write(*,*) '                     Map 1|',r(1,1), r(1,2)
         write(*,*) '                     Map 2|',r(2,1), r(2,2)

         deallocate(array1)
         deallocate(xbar)
         deallocate(std)
         deallocate(ssp)
         deallocate(r)

! command 85 currently unused
      elseif (command.eq.85) then

! command 86 currently unused
      elseif (command.eq.86) then

! test quadrature:
      elseif (command.eq.87) then
         allocate(ctempmap(maxsize,maxsize))
         open(12,file='quadtest_fourn.dat')
         open(13,file='quadtest_quad.dat')
         uvcell = 1./(cellsize*maxsize*sec2rad)
         do i = 1,maxsize
            do j = 1,maxsize
               upt = real(i-(maxsize/2+1))*uvcell
               vpt = real(j-(maxsize/2+1))*uvcell
!               call by_quadrature(upt,vpt,redat,imdat)
               write(12,*) re(i,j,chan),im(i,j,chan)
               write(13,*) redat,imdat
               ctempmap(i,j) = cmplx(redat,imdat)
            end do
         end do
         call display_map(real(ctempmap),maxsize,graph_adj,&
              cellsize,cellsize,0D0,0D0)
         call display_map(aimag(ctempmap),maxsize,graph_adj,&
              cellsize,cellsize,0D0,0D0)
         ctempmap = cmplx(real(ctempmap)-re(:,:,chan),aimag(ctempmap)-im(:,:,chan))
         call display_map(real(ctempmap),maxsize,graph_adj,&
              cellsize,cellsize,0D0,0D0)
         call display_map(aimag(ctempmap),maxsize,graph_adj,&
              cellsize,cellsize,0D0,0D0)
         deallocate(ctempmap)
         close(12)
         close(13)

! Get-xray-telescope
      elseif (command.eq.89) then
         call get_xray_telescope

! Get-cosmology
      elseif (command.eq.90) then
         call plot_cosmology

! Display-primary-beam
      elseif (command.eq.91) then
         call display_primary_beam

! End of commands case loop        
      end if

   end do main_loop

! Deallocate all common arrays
   call do_deallocation

! Close pgplot 
   call pgend

end program profile







