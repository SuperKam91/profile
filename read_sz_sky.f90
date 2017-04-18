   subroutine read_sz_sky
! Moved from main program
      
      use kind_def
      use sz_globals
      use physics
      use file_handling
      use fits_handling
      use map_handling

      implicit none

      character*1 :: which_form, which_spec, chr1
      character(len=fname_len) :: fits_file
      real :: bmaj, bmin
      real(kind=dp) :: scale, synth_beam, fits_freq, alpha, T1, xsum, incell
      real(kind=dp) :: norm 
!      real(kind=dp), dimension(:,:), allocatable :: fits_map 
      integer :: i
      logical :: stat_fl
      logical, external :: io_yesno

      overwrite = io_yesno('Overwrite existing SZ sky','no',status)

      write(*,*) 'Select type of input map'
      write(*,*) 
      write(*,*) 'FITS format map'
      write(*,*) '   (I)ntensity map in Jy/beam'
      write(*,*) '   Intensity map in Jy/(P)ixel'      
      write(*,*) '   Intensity map in (M)Jy/sr'
      write(*,*) '   (T)hermodynamic temperature'
      write(*,*) '   (B)rightness temperature'
      write(*,*) '   (Y) parameter'
      write(*,*) 
      write(*,*) '(D)aisuke format map'
      write(*,*) '(A)SCII y map'
      write(*,*) 
      call io_getc('Select one','B',which_form,status)

      scale = 1.0
      select case(which_form)
      case('i','I')
         write(*,*) 'Convolving beam size in arcsec'
         call io_getr('Major axis','15.0',bmaj,status)
         call io_getr('Minor axis','15.0',bmin,status)

! conversion T to Jy/bm:
         synth_beam = pi*bmaj*bmin*(sec2rad**2)/(4.*log(2.)) !!CHECK
         norm = 1.0D26*synth_beam*k_b*2.0*obsfreq**2/const_c2
         scale = 1.0/norm
         write(*,*) 'Appropriate Amplitude rescaling is ',scale
      case('p','P')
         norm = cellsize*sec2rad
         norm = norm*norm
         norm = norm*2.0*k_b*obsfreq**2/const_c2*1.0D26
         scale = 1./norm
         write(*,*) 'Appropriate Amplitude rescaling is ',scale
      case('m','M')
         norm = 2.0*k_b*obsfreq**2/const_c2*1.0D26
         scale = 1d6/norm
         write(*,*) 'Appropriate Amplitude rescaling is ',scale
      case('t','T')
         call io_getd('Map observing frequency (GHz)','15.0',&
              fits_freq,status)
         scale = dB_by_dT(T0,fits_freq)*const_c2/&
                 (fits_freq**2*2.0D0*k_b)
         call io_getc('Is your map in uK?','yes',chr1,status)
         if (chr1 == 'y') then
            scale = scale*1e-6 
         end if
         write(*,*) 'Appropriate Amplitude rescaling is ',scale
      case('y','Y','a','A')
         scale = g_nu_x2(obsfreq)*T0
         write(*,*) 'Appropriate Amplitude rescaling is ',scale
      case('d','D')
!Daisuke maps in y-parameter
         scale = g_nu_x2(obsfreq)*T0         
         write(*,*) 'Appropriate Amplitude rescaling is ',scale
      end select

! Select temp_map for loading in files
      p_sky => temp_sky
      
      select case(which_form)
      case('d','D')
         prompt = 'Daisuke format filename'
         which_dir = profile_data_dir
         rootname = 'tsz_a0'
         extn = '.626_x'
         call get_read_filename(fits_file)
         call read_daisuke_map(fits_file)
      case('i','I','p','P','m','M','t','T','b','B','y','Y')
         prompt = 'FITS filename'
         which_dir = profile_data_dir
         rootname = 'simclusterRaw'
         extn = '.fits'
         call get_read_filename(fits_file)
         call read_fits_map(fits_file)
         write(*,*) 'Sum of map =', sum(p_sky)*(cellsize/60d0)**2
      case('a','A')
         prompt = 'FITS filename'
         which_dir = profile_data_dir
         rootname = 'amiR10Y10_15arcsec'
         extn = '.txt'
         call get_read_filename(fits_file)
         call read_ascii_map(fits_file)
      end select

! Ask what spectrum to use for channel maps
      write(*,*) 'Select spectrum of input map'
      write(*,*) 
      write(*,*) '(S)Z'
      write(*,*) '(C)MB'
      write(*,*) '(T)hermal'
      write(*,*) '(P)ower law'
      write(*,*) 
      call io_getc('Select one','S',which_spec,status)       
      if ((which_spec.eq.'P').or.(which_spec.eq.'p'))then
         call io_getd('Spectral index','0.7',alpha,status)
      end if
      if ((which_spec.eq.'t').or.(which_spec.eq.'T'))then
         call io_getd('Temperature of emitter','70.',T1,status)
      end if

! Overwrite szsky if necessary
      if (overwrite) then
         szsky = 0.d0
      end if

! Write into szsky array
      do i = 1,nchan
         select case (which_spec)
         case('s','S')
            norm = scale*g_nu_x2(dble(nu(i)))/g_nu_x2(obsfreq)
         case('c','C')
            norm = scale*planck_fn(T0,dble(nu(i)))/planck_fn(T0,obsfreq)
         case('t','T')
            norm = scale*planck_fn(T1,dble(nu(i)))/planck_fn(T1,obsfreq)
         case('p','P')
            norm = scale*(nu(i)/obsfreq)**(-2.-alpha)
         end select
         szsky(:,:,i) = szsky(:,:,i)+temp_sky*norm
      end do

      if (verbose) then
         chan = nchan/2+1
         xsum = sum(szsky(:,:,chan))
         write(*,*) 'At frequency of ',nu(chan)
         write(*,*) '   MAXVAL in map:',maxval(szsky(:,:,chan)),' K (Tb)'
         write(*,*) '   Sum of brightness temperature pixels in map is ',xsum
         write(*,*) '   Total flux density in map is ',&
           xsum*(cellsize*sec2rad)**2*2.0*k_b*nu(chan)**2/const_c2*1.0D26,' Jy'
      end if

   end subroutine read_sz_sky




