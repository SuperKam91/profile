subroutine make_mosaic_observation(idum)

! 08/08/02 - corrected bug in declination offset spacing
! 05/03/03 - filenames now always appended by four integer [KG]

   use kind_def
   use sz_globals
   use file_handling
   use fits_handling
   use kgraph
   use maths

   implicit none
   character*80 :: logfile, ps_file
   character*1 :: mos_type
   real(kind=dp) :: psi, inc, rad, sindec, ndays
   integer :: i, j, k, l, lunit, vis_len
   integer :: lp, ps_stat
   real(kind=dp) :: ha, dec, ha_start, ha_stop, rms_s, ra
   real(kind=dp), dimension(n_antennas,n_antennas) :: u_2d,v_2d,w_2d
   integer, dimension(n_antennas,n_antennas) :: shadow      
   integer, dimension(:), allocatable :: memb_in_row
   real(kind=dp) :: mos_spac
   real :: diam
   character :: chra*16, chdec*16
   integer :: prra, prdec, lr, ld
   character(len=fname_len) :: visname, basename, fitsname
   character*8 :: chr1
   character*1 :: chr_epoch
   character*80 :: text
   integer :: idum, n_flag, l_dir
   integer :: n_row, n_mosaic, n_order, n_obs, base_len, sl1, fl1
   integer :: n_RA, n_DEC
   logical :: add_noise, sing_phase_cent, exp_wei, sig2noi
   logical :: indiv_vis, indiv_fits, multi_fits, sep_vis_chan
   logical :: io_yesno
   external io_yesno
   integer :: chr_lenb
   external chr_lenb
   integer :: pgopen
   external pgopen

   if ((which_tel.ne.'d').or.(which_tel.ne.'d')) then
      write(*,*) 'Not supported for non-defined telescopes at this time'
      return
   end if

   prra = 1
   prdec = 0
   ha_start = -4.0*15.0*deg2rad
   ha_stop = 4.0*15.0*deg2rad
   ha_inc = 600.0/3600.0*15.0*deg2rad
   n_samp = 1+(ha_stop-ha_start)/ha_inc
   dec = obs_decref*deg2rad
   ra = obs_raref*hr2rad
   rms_s = 0.540
   n_flag = 0


   multi_fits = io_yesno('Write multi-fits file for all pointings',&
               'yes',status)
   if (.not. multi_fits) then
      indiv_vis = io_yesno('Write individual vis files for each pointing',&
               'yes',status)
      indiv_fits = io_yesno('Write individual fits files for each pointing',&
                'no',status)
      !write(*,*) indiv_fits
   else
      indiv_vis = .false. 
      indiv_fits = .false.
   end if

   if (indiv_fits.or.multi_fits) then
      call io_getc('Observation name:','*',obs_name,status)
      chr_epoch='J'
      call io_getc('Epoch [B/J]:','*',chr_epoch,status)
      if ((chr_epoch=='b').or.(chr_epoch=='B')) obs_epoch = 1
      if ((chr_epoch=='j').or.(chr_epoch=='J')) obs_epoch = 2
   end if

   if (indiv_vis) then
      if (nchan.gt.1) then
         sep_vis_chan = io_yesno('Write seperate vis file for each channel',&
                        'yes',status)
      else
         sep_vis_chan = .false.
      end if
   end if

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
   ha_inc = ha_inc*3600.0/(15.0*deg2rad)
   call io_getd('Time between uv samples (s) ','*',ha_inc,status)
   ha_inc = ha_inc/3600.0*15.0*deg2rad
   n_samp = 1+(ha_stop-ha_start)/ha_inc
   write(*,*) 'Taking ',n_samp,' samples around each uv track'
   add_noise = io_yesno('Add Gaussian noise to the data','yes',status)
   if (add_noise) then
      write(*,*) 
      write(*,*) 'RMS noise on RT per baseline per second - 0.132 Jy'
      write(*,*) 'RMS noise on upgraded RT per channel per baseline per second - 0.045 Jy'
      write(*,*) 'RMS noise on AMI per channel per baseline per second - 0.540 Jy'
      write(*,*)
      call io_getd('RMS noise (Jy) per channel per baseline per second','*',rms_s,status)
   else
      rms_s = 1.0
   end if

   if (n_samp*nchan*n_antennas*(n_antennas-1)/2.gt.mvis) then
      write(*,*) 'Sample time too short/mvis too small'
      return
   end if

   call io_getd('Number of days repeat in simulation',&
                '1.0',ndays,status)

   mos_spac = 1.0/(sqrt(2.0)*dish_diameter*obsfreq/const_c*sec2rad)
   write(*,*) 'Optimum sepn (for sensitivity) = FWHM/sqrt 2'
   write(*,*) 'Optimum sepn (for coherent mosaicing) = FWHM/2'
   call io_getd('Separation of mosaic points (arcsec)','*',mos_spac,status)

! Select mosaicing scheme
   write(*,*) 'Available mosaicing schemes'
   write(*,*) '(H)exangonal close packed'
   write(*,*) '(T)riangle of pointings'
   write(*,*) '(A)MI-type field'
   write(*,*)
   call io_getc('Select Type of mosaic','H',mos_type,status)

   select case(mos_type)

! A 3 field triangular mosaic 
   case('t','T')
      n_mosaic = 3
      n_obs = 3

! Allocate arrays
      allocate(ra_offset(n_mosaic))
      allocate(dec_offset(n_mosaic))
      if (io_yesno('Triangle pointing up','yes',status)) then
         ra_offset(1) = 0.d0
         dec_offset(1) = mos_spac/sqrt(3.)
         ra_offset(2) = -mos_spac/2.
         dec_offset(2) = -mos_spac/(2.*sqrt(3.))
         ra_offset(3) = mos_spac/2.
         dec_offset(3) = -mos_spac/(2.*sqrt(3.))
      else
         ra_offset(1) = -mos_spac/2.
         dec_offset(1) = mos_spac/(2.*sqrt(3.))
         ra_offset(2) = mos_spac/2.
         dec_offset(2) = mos_spac/(2.*sqrt(3.))
         ra_offset(3) = 0.d0
         dec_offset(3) = -mos_spac/sqrt(3.)
      end if

   case('A','a')
      call io_geti('Number of pointings per row (in RA)?', '3', n_RA, status)
      call io_geti('Number of pointings per column (in DEC)?', '2', n_DEC, status)
      ! Need statement if pointings are too large for current field
      n_mosaic = n_RA*n_DEC
      allocate(ra_offset(n_RA*n_DEC))
      allocate(dec_offset(n_RA*n_DEC))
      n_obs = 1
      
      do l= 1,n_DEC
         do k = 1, n_RA
            ra_offset(n_obs) = ((float(n_RA)-2.0*float(k))/2. + float(mod(l,2))/2.) * mos_spac
            dec_offset(n_obs) = float(l-n_DEC/2)*(mos_spac)*sqrt(3./4.)
            n_obs = n_obs+1
         enddo
      enddo

! A hexagonal close packed mosaic
   case default
      n_order = (cellsize*maxsize/mos_spac)/2+1
      call io_geti('Order of mosaic','*',n_order,status)
      n_mosaic = 3*n_order**2-3*n_order+1
      write(*,*) 'Doing ',n_mosaic,' observations'

! Allocate arrays
      allocate(memb_in_row(-n_order:n_order))
      allocate(ra_offset(n_mosaic))
      allocate(dec_offset(n_mosaic))

! loop over rows
      do n_row = -(n_order-1), 0
         memb_in_row(n_row) = 2*n_order-1+n_row
         memb_in_row(-n_row) = 2*n_order-1+n_row
      end do
      n_obs = 1

! loop over rows
      do n_row = -(n_order-1), n_order-1

! loop over observations in this row
         do i = 1,memb_in_row(n_row)
            dec_offset(n_obs) = n_row*mos_spac*sqrt(3.0)/2.0
            ra_offset(n_obs) = (i-((memb_in_row(n_row)+1.0)/2.0))*mos_spac
            n_obs = n_obs+1
         end do
      end do

! deallocate working array
      deallocate(memb_in_row)
   end select

! Get file names
   which_dir = profile_data_dir
   l_dir  = chr_lenb(which_dir)  
   basename = which_dir(1:l_dir)//'default'
   call io_getfil('File name (_n.vis / .fits added):','*',basename,status)
   sing_phase_cent = io_yesno('Phase centre for all maps coincident',&
                     'yes', status)

! Not a useful option any more?
   sig2noi = io_yesno('Set weight to s-to-n ratio **2','no',status)

! Open log file and write header
   call io_nxtlun(lunit,status)
   base_len = chr_lenb(basename)
   logfile = basename(1:base_len)//'.log'
   open (lunit,file=logfile,form='formatted')
   write(*,*) 'Writing log file to ',logfile
   write(lunit,*) 'Observation of ',n_mosaic,' fields centred on ',&
                  ra/hr2rad,'  ',  dec/deg2rad
   write(lunit,*) 'Hour angle ',ha_start/(deg2rad*15.0),' to ',&
                   ha_stop/(deg2rad*15.0)
   write(lunit,*) 'Each visibility represents ',ndays,' integration'
   write(lunit,*) 
   close(lunit)

! plot mosaicing strategy
   if (plot_open) then
      chan = 1
      call display_map(szsky(:,:,chan),maxsize,.false.,cellsize,cellsize,&
                       0.D0,0.D0)
      call pgsfs(2)
      diam = 1.0/(2.0*dish_diameter*obsfreq/const_c*sec2rad)
      if (diam/cellsize.lt.20) then
         call pgsch(0.5)
      else
         call pgsch(1.0)
      end if

      do i = 1, n_mosaic
         call pgsci(plot_col)
         call pgcirc(ra_offset(i),dec_offset(i),diam) 
         call pgnumb(i,0,1,text,status)
         call pgtext(ra_offset(i),dec_offset(i),text)
      end do
   end if

! plot postscript of mosaicing strategy
   ps_file = basename(1:base_len)//'.ps'
   lp = chr_lenb(ps_file)
   ps_file = ps_file(1:lp)//'/PS'
   ps_stat = pgopen(ps_file)
   call pgslct(ps_stat)
   call pgsfs(2)
   chan = 1
   call display_map(szsky(:,:,chan),maxsize,.false.,cellsize,cellsize,0D0,0D0)
   diam = 1.0/(2.0*dish_diameter*obsfreq/const_c*sec2rad)
   if (diam/cellsize.lt.20) then
      call pgsch(0.5)
   else
      call pgsch(1.0)
   end if

   do i = 1, n_mosaic
      call pgsci(0)
      call pgcirc(ra_offset(i),dec_offset(i),diam) 
      call pgnumb(i,0,1,text,status)
      call pgtext(ra_offset(i),dec_offset(i),text)
   end do
   call pgsci(plot_col)
   call pgclos

   if (plot_open) then
      call pgslct(pl_dev1)
      call pgsch(char_hgt)
   end if

! Calculate rms on one sample
   rms = rms_s/sqrt(ndays*ha_inc*3600.0/(15.0*deg2rad))

   if (rms_s.ne.0.0) then
      weight = 1/rms**2
   else
      weight = 1.0
   end if

! Calculate u, v and check for shadowing
   call calc_basel

   nvis = 1
   weight = 1.

! loop over all mosaic observations
   do n_obs = 1, n_mosaic

! restart nvis count if writing out individual vis or fits files
      if ((indiv_fits).or.(indiv_vis)) then
         nvis = 1
      end if

! prepare name of vis and fits files
      write(chr1,'(I4.4)') n_obs
      sl1 = chr_lenb(chr1)
      base_len = chr_lenb(basename)

! leave addition of '.vis' to writing routine
      visname = basename(1:base_len)//'_'//chr1(1:sl1)
      fitsname = basename(1:base_len)//'_'//chr1(1:sl1)//'.fits'

! setup offset appropriate to this observation
      pb_x = ra_offset(n_obs)
      pb_y = dec_offset(n_obs)
      if (sing_phase_cent) then
         ph_x = -ra_offset(n_obs)
         ph_y = -dec_offset(n_obs)
      else
         ph_x = 0.0
         ph_y = 0.0
      end if
      obs_decref = (dec+pb_y*sec2rad)/deg2rad
      obs_raref = (ra-(pb_x*sec2rad)/(cos(dec-ph_y*sec2rad)))/hr2rad
      call chr_chrtos(obs_raref,prra,chra,lr)
      call chr_chrtos(obs_decref,prdec,chdec,ld)
      if (sing_phase_cent) then
         obs_rarot = ra/hr2rad
         obs_decrot = dec/deg2rad
      else
         obs_rarot = obs_raref
         obs_decrot = obs_decref
      end if
      write(*,*)
      write(*,*) 'Pointing ',n_obs,' :'
      write(*,*) '          offset   ',pb_x, pb_y
      write(*,*) '          position ',chra(1:lr),'   ',chdec(1:ld)
      vis_len = chr_lenb(visname)
      open (lunit,file=logfile,form='formatted',position='append')
      write(lunit,*) visname(1:vis_len),' at ',chra(1:lr),'   ',chdec(1:ld)
      close(lunit)

! simulate observation of sky with interferometer
      verbose = .false.
      do_plot = .false.
      call make_aperture

! loop over all sample times
      do k = 1, n_samp
         ha = (ha_start+(k-1)*ha_inc)

! calculate uv coverage for this HA
         call calc_uvw(ha,dec,u_2d,v_2d,w_2d,shadow)

! loop over first antenna
         do i = 1, n_antennas-1

! loop over second antenna
            do j = i+1, n_antennas

! check whether telescope has more than one frequency channel
               if (nchan.gt.1) then

! loop over channels
                  do chan = 1, nchan
                     u(nvis) = u_2d(i,j)*nu(chan)/obsfreq
                     v(nvis) = v_2d(i,j)*nu(chan)/obsfreq
                     call extract_visibility(re(:,:,chan),im(:,:,chan),&
                          maxsize,-u(nvis),v(nvis),&
                          cellsize,data_re(nvis),data_im(nvis))
                     if (sig2noi) then
                        weight(nvis) = (data_re(nvis)**2+data_im(nvis)**2)/&
                             rms_s**2
                     end if
                     which_pointing(nvis) = n_obs

! flag visibility if shadowed baseline
                     if (shadow(i,j).ne.1) then
                        weight(nvis) = 0.0
                        data_re(nvis) = 0.0
                        data_im(nvis) = 0.0
                        n_flag = n_flag+1
                     end if
                     nvis = nvis+1
                  end do
               else
                  chan = 1
                  u(nvis) = u_2d(i,j)
                  v(nvis) = v_2d(i,j)
                  call extract_visibility(re(:,:,chan),im(:,:,chan),&
                       maxsize,-u(nvis),v(nvis),&
                       cellsize,data_re(nvis),data_im(nvis))
                  if (sig2noi) then
                     weight(nvis) = (data_re(nvis)**2+data_im(nvis)**2)/&
                          rms_s**2
                  end if
                  if (shadow(i,j).ne.1) then
                     weight(nvis) = 0.0
                     data_re(nvis) = 0.0
                     data_im(nvis) = 0.0
                  end if
                  which_pointing(nvis) = n_obs
                  nvis = nvis+1
               end if
               if (nvis.ge.mvis) then
                  write(*,*) 'mvis too small in globals'
                  return
               end if
            end do
         end do
      end do

! adjust nvis count if writing out individual vis or fits files
      if ((indiv_fits).or.(indiv_vis)) then
         nvis = nvis-1
      end if

! Add noise if necessary
      if (add_noise) then
         do i = 1, nvis
            data_re(i) = data_re(i)+rms(i)/sqrt(2.)*gasdev(idum)
            data_im(i) = data_im(i)+rms(i)/sqrt(2.)*gasdev(idum)
         end do
      end if

! Write out individual vis or fits fills is required
      if (indiv_vis) then
         call write_vis_data(visname,sep_vis_chan)
      end if
      if (indiv_fits) then 
         call write_uv_fits(fitsname)
         fl1 = chr_lenb(fitsname)
         write(*,*) nvis,' visibilities written to ',fitsname(1:fl1)
      end if

   end do

! Ajust nvis count 
   nvis = nvis-1

! Reset obs_raref and obs_decref
   obs_raref = ra/hr2rad
   obs_decref = dec/deg2rad

! Write multi-source fits file
   if (multi_fits) then
      obs_raref = ra/hr2rad
      obs_decref = dec/deg2rad
      fitsname = basename(1:base_len)//'.fits'
      call write_multi_fits(n_mosaic,fitsname,sing_phase_cent)
      fl1 = chr_lenb(fitsname)
      write(*,*) nvis,' visibilities written to ',fitsname(1:fl1)
      write(*,*) n_flag, ' visibilities flagged'
   end if

! Rewind log file
   rewind(lunit)
   status = 0

! Deallocate arrays
   deallocate(ra_offset)
   deallocate(dec_offset)

end subroutine make_mosaic_observation




