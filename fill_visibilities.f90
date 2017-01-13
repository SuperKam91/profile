subroutine fill_visibilities(idum,fill_single)

   use kind_def
   use sz_globals
   use file_handling
   use fits_handling
   use maths

   implicit none

   character(len=fname_len) :: filename
   character*1 :: chr1, chr_epoch
   integer :: idum, ivis
   real(kind=dp) :: visrms, scale
   real(kind=dp) :: tmp_re,tmp_im
   logical :: ex, fill_single, add_noise, userms, sep_vis_chan
   logical, external :: io_yesno

! Initialise Profile's sampled visibilities arrays
   u = 0.
   v = 0.
   w = 0.
   data_re = 0.
   data_im = 0.
   rms = 0.
   weight = 0.
   which_pointing = 0  
   nvis = 1

! get template
   call io_getc('Read template uv postions from a (V)is or (F)its file','v',&
                 chr1,status)
   select case(chr1)
   case('v','V')
      call read_vis_data
   case('f','F')
      call read_uv_fits
   case default
      write(*,*) 'No supported file type selected'
      return
   end select
 
   overwrite = io_yesno('Overwrite read in visibilities','yes',status)

! Add Gaussian random noise to the visibilities?
   add_noise = io_yesno('Add noise?:','y',status)
   if (add_noise) then
      userms = io_yesno('Use rms from file (scaled)?:','y',status)
      if (userms) then
         call io_getd('Noise scaling factor: ','1.0',scale,status)
      else
         call io_getd('rms (Jy) per visibility to use?:','1.0',scale,status)
      end if
   end if

   if (.not.fill_single) chan = 1

! loop over visibilies generating new sampled data at template positions
   do ivis = 1, nvis

! reset chan if required
      if (chan.gt.nchan) chan = 1

! find new real and imag NB extract visibility from -u instead of u to 
! be consistent with fit_visibilities
      call extract_visibility(re(:,:,chan),im(:,:,chan),maxsize,-u(ivis),&
                              v(ivis),cellsize,tmp_re,tmp_im)

      if (overwrite) then
         data_re(ivis) = tmp_re
         data_im(ivis) = tmp_im
      else
         data_re(ivis) = data_re(ivis)+tmp_re
         data_im(ivis) = data_im(ivis)+tmp_im
      end if

      if (add_noise) then
         if (userms) then
            visrms = rms(ivis)*scale
         else
            visrms = scale
         end if
         data_re(ivis) = data_re(ivis)+(visrms/sqrt(2.d0)*gasdev(idum))
         data_im(ivis) = data_im(ivis)+(visrms/sqrt(2.d0)*gasdev(idum))
      end if

      if (.not.fill_single) chan = chan+1 
   end do

! Set number of visibilities etc
   n_samp = nvis/(nchan*n_antennas*(n_antennas-1)/2)
   n_pol = 1

! Write vis file?
   if (io_yesno('Write vis file','yes',status)) then
      if (nchan.gt.1) then
         sep_vis_chan = io_yesno('Write seperate vis file for each channel',&
                        'yes',status)
      else
         sep_vis_chan = .false.
      end if

! Get filename
      prompt = 'File name (_n.vis added):'
      extn = ''
      which_dir = profile_data_dir
      rootname = 'default'
      call get_write_filename(filename)

      call write_vis_data(filename,sep_vis_chan)
   end if
      
! Write fits file?
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

end subroutine fill_visibilities








