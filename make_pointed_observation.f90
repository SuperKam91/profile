subroutine make_pointed_observation(idum)

   use kind_def
   use sz_globals
   use file_handling
   use fits_handling
   use maths

   implicit none

   character(len=fname_len) :: filename
   character*1 :: chr_epoch
   real(kind=dp) :: psi, inc, rad, sindec
   integer :: i, j, k, ndays, l_dir
   real(kind=dp) :: ha, ra, dec, ha_start, ha_stop, rms_s
   real(kind=dp), dimension(n_antennas,n_antennas) :: u_2d,v_2d,w_2d
   integer, dimension(n_antennas,n_antennas) :: shadow      
   integer :: idum
   logical :: add_noise, sep_vis_chan
   logical :: io_yesno
   external io_yesno
   logical :: ex
   character*1 :: chr1

   select case(which_tel)
   case('d','D')
      ha_start = -6.0*15.0*deg2rad
      ha_stop = 6.0*15.0*deg2rad
      ha_inc = 200.0d0/3600.0d0*15.0d0*deg2rad
      n_samp = 1+(ha_stop-ha_start)/ha_inc
      rms_s = 0.132d0

      ra = obs_raref*hr2rad
      call io_getra('RA','*',ra,status)
      obs_raref = ra/hr2rad
      obs_rarot = obs_raref

      dec = obs_decref*deg2rad
      call io_getdec('Declination','*',dec,status)
      obs_decref = dec/deg2rad
      obs_decrot = obs_decref

      ha_start = ha_start/(deg2rad*15.0)
      call io_getd('HA start','*',ha_start,status)
      ha_start = ha_start*deg2rad*15.0
      ha_stop = ha_stop/(deg2rad*15.0)
      call io_getd('HA stop','*',ha_stop,status)
      ha_stop = ha_stop*(deg2rad*15.0)
      ha_inc = ha_inc*3600.0d0/(15.0d0*deg2rad)
      call io_getd('Sample integration time (s) ','*',ha_inc,status)
      ha_inc = ha_inc/3600.0d0*15.0d0*deg2rad
      n_samp = 1+(ha_stop-ha_start)/ha_inc
      write(*,*) 'Taking ',n_samp,' samples'
      add_noise = io_yesno('Add Gaussian noise to the data','yes',status)
      if (add_noise) then
         write(*,*) 
         write(*,*) 'RMS noise on RT per baseline per second - 0.132 Jy'
         write(*,*) 'RMS noise on AMILA per baseline per second - 0.016 Jy'
         write(*,*) 'RMS noise on AMISA per baseline per second - 0.191 Jy'
         write(*,*)
         call io_getd('RMS noise (Jy) per channel per baseline per second',&
                      '*',rms_s,status)
      else
         rms_s = 1.0
      end if
      if (n_samp*n_antennas*(n_antennas-1)/2.gt.mvis) then
         write(*,*) 'Sample time too short/mvis too small'
         write(*,*) n_samp*n_antennas*(n_antennas-1)/2,' required; ',mvis,&
                   ' availaible'
         return
      end if
      call io_geti('Number of days of observations','1',ndays,status)

! Check total number of visibilities to be generated
      if (n_samp*nchan*n_antennas*(n_antennas-1)/2.gt.mvis) then
         write(*,*) 'Sample time too short/mvis too small'
         return
      end if

! Calculate rms on one sample
      rms = rms_s/sqrt(ndays*ha_inc*3600.0/(15.0*deg2rad))

! Calculate u, v and check for shadowing
      call calc_basel
      nvis = 1
      do k = 1, n_samp
         ha = (ha_start+(k-1)*ha_inc)
         call calc_uvw(ha,dec,u_2d,v_2d,w_2d,shadow)
         do i = 1, n_antennas-1
            do j = i+1, n_antennas
               if (shadow(i,j).eq.1) then
                  do chan = 1, nchan
                     u(nvis) = u_2d(i,j)*nu(chan)/obsfreq
                     v(nvis) = v_2d(i,j)*nu(chan)/obsfreq
                     w(nvis) = w_2d(i,j)*nu(chan)/obsfreq
                     call extract_visibility(re(:,:,chan),im(:,:,chan),&
                                          maxsize,-u(nvis),v(nvis),&
                                          cellsize,data_re(nvis),data_im(nvis))
                     weight(nvis) = 1.0
                     nvis = nvis+1
                  end do
               end if
            end do
         end do
      end do
      nvis = nvis-1

      write(*,*) nvis,' out of ',nchan*n_samp*n_antennas*(n_antennas-1)/2,&
                      ' unshadowed'

   case default
      write(*,*) 'Not currently supported'
      return
   end select

! Add noise if necessary
   if (add_noise) then
      do i = 1, nvis
         data_re(i) = data_re(i)+rms(i)/sqrt(2.)*gasdev(idum)
         data_im(i) = data_im(i)+rms(i)/sqrt(2.)*gasdev(idum)
         if (rms(i).ne.0.0) then
            weight(i) = 1/rms(i)**2
         else
            weight(i) = 1.0
         end if
      end do
   end if

   if (io_yesno('Write vis file','yes',status)) then
      if (nchan.gt.1) then
         sep_vis_chan = io_yesno('Write separate vis file for each channel',&
                        'yes',status)
      else
         sep_vis_chan = .false.
      end if

! Get filename
      if (sep_vis_chan) then
         prompt = 'File name (_n.vis added):'
      else
         prompt = 'File name (.vis added):'
      end if
      extn = ''
      which_dir = profile_data_dir
      rootname = 'default'
      call get_write_filename(filename)

      call write_vis_data(filename,sep_vis_chan)
   end if
      
   if (io_yesno('Write fits file','yes',status)) then

! Get filename
      prompt = 'FITS filename'
      extn = '.fits'
      which_dir = profile_data_dir
      rootname = 'uv'
      call get_write_filename(filename)

! Define other parameters for .fits file
      call io_getc('Observation name:','*',obs_name,status)
      chr_epoch='J'
      call io_getc('Epoch [B/J]:','*',chr_epoch,status)
      if ((chr_epoch=='b').or.(chr_epoch=='B')) obs_epoch = 1
      if ((chr_epoch=='j').or.(chr_epoch=='J')) obs_epoch = 2

      call write_uv_fits(filename)
      write(*,*) 'FITS file written'
   end if

end subroutine make_pointed_observation
