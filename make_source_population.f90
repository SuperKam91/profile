subroutine make_source_population(idum)

   use kind_def
   use sz_globals
   use maths
   use kgraph
   use file_handling
   use astronomy

   implicit none
   integer :: idum
   character*1 :: chr1, chr2
   character(len=fname_len) :: filename
   real(kind=dp) :: k_parm, gamma, s1, s2, alpha, alpha_mean, alpha_rms, norm
   real(kind=dp) :: dec2, ra2, ra1, dec1
   real(kind=dp) :: base_flux, chan_flux, frac_ra, frac_dec, dx, dy
   integer :: centre, i, j, k, nn, src, temp1, temp2
   integer :: ra_pix_cent, dec_pix_cent, out_ra_pix, out_dec_pix
   integer :: iunit
   logical :: to_file, ex, src_file
   logical :: io_yesno
   external io_yesno

   centre = maxsize/2+1
   alpha_mean = 0.7d0
   alpha_rms = 0.3d0

   do_plot = io_yesno('Do plots','n',status)
 
   overwrite = io_yesno('Overwrite existing model sky','n',status)

   write(*,*) 'Available point source populations:'
   write(*,*) '(R)andom distribution of point sources from power law'
   write(*,*) '(L)ist of point sources from file'

   call io_getc('What type of test sky','r',chr1,status)

   if (overwrite) then
      szsky = 0.0D0
   end if

   what_sky : select case(chr1)
   case('r','R')
      write(*,*) '*****NOTE UNIT CHANGE*****'
      write(*,*)
      write(*,*) 'dN/dS = K S^(-gamma) Jy^-1 str^-1   S1<S<S2'
      write(*,*)
      write(*,*) &
           'Values from 9C RT counts:'
      write(*,*) 'K = 51.0'
      write(*,*) 'Gamma = 2.15'
      write(*,*) 
      write(*,*) &
           'Values from 10C AMI LA counts:'
      write(*,*) 'K = 340'
      write(*,*) 'Gamma = 1.81'
      write(*,*) 


!      k_parm = 51
!      gamma = 2.15
      k_parm = 340.0d0
      gamma = 1.81d0 
      s1 = 1.d-5
      s2 = 3.d-4
      verbose = .true.
      call io_getd('Value of K',&
           & '*',k_parm,status)
      call io_getd('Value of gamma',&
           & '*',gamma,status)
      call io_getd('Minimum flux, S1 (Jy) ',&
           & '*',s1,status)
      call io_getd('Maximum flux, S2 (Jy) ',&
           & '*',s2,status)
      call io_getd('Mean value of spectral index ',&
           & '*',alpha_mean,status)
      call io_getd('RMS of spectral index ',&
           & '*',alpha_rms,status)
      if (io_yesno('Write source list to file','no',status)) then
	 to_file = .true.
         src_file = io_yesno('Write in .src format','yes',status) 

! Get filename
         prompt = 'Source list filename'
         which_dir = profile_data_dir
         rootname = 'sources'
         status = 0
         if (src_file) then
            extn = '.src'
         else
            extn = '.list'
         end if
         call get_write_filename(filename)

! Open file
         call io_nxtlun(iunit,status)
         open (iunit,file=filename,form='formatted')
      else
         iunit = 60
         to_file = .false.
      end if
      call scatter_random_sources&
         (k_parm,gamma,s1,s2,alpha_mean,alpha_rms,to_file,src_file,iunit,idum)
      if (to_file) then
         close(iunit)
         rewind(iunit)
      end if

   case('l','L') 
      ra2 = obs_raref*hr2rad
      call io_getra('RA','*',ra2,status)
      obs_raref = ra2/hr2rad
      dec2 = obs_decref*deg2rad
      call io_getdec('Declination','*',dec2,status)
      obs_decref = dec2/deg2rad
      ra_pix_cent = maxsize/2+1
      dec_pix_cent = maxsize/2+1
      call read_ptsrc_list
      do src = 1,nsrc
         ra1 = rapos(src)*hr2rad
         dec1 = decpos(src)*deg2rad
         call sin_proj_forward(ra1,dec1,ra2,dec2,cellsize,dx,dy)
!         write(*,*) src,rapos(src),decpos(src),ra,dec,dx,dy
         dx = dx+ra_pix_cent
         dy = dy+dec_pix_cent
         if (dx.ge.0.0) then
            out_ra_pix = dx
         else
            out_ra_pix = dx-1
         end if
         if (dy.ge.0.0) then
            out_dec_pix = dy
         else
            out_dec_pix = dy-1
         end if
         frac_ra = 1.d0-(dx-out_ra_pix)
         frac_dec = 1.d0-(dy-out_dec_pix)

! Error checking
         if ((frac_ra.gt.1.0).or.(frac_ra.lt.0.0).or.&
              (frac_dec.gt.1.0).or.(frac_dec.lt.0.0)) then
            write(*,*) 'Error in make_source_population fraction calc'
            write(*,*) 'Fractions ',frac_ra,frac_dec
            write(*,*) 'Pixel offsets from projection ',dx, dy
            write(*,*) 'Actual pix offs ',out_ra_pix,out_dec_pix
            write(*,*) 'Centre of map ',ra_pix_cent,dec_pix_cent
         end if

! Add into output map if source within map
         if ((out_ra_pix.gt.1).and.&
              (out_ra_pix.lt.maxsize-1).and.&
              (out_dec_pix.gt.1).and.&
              (out_dec_pix.lt.maxsize-1)) then

! Change from Jy to brightness temperature
            base_flux = srcflux(src)/&
                   ((cellsize*sec2rad)**2*2.0*k_b*obsfreq**2/const_c2*1.0D26)

! Loop over channels
            do chan = 1, nchan
               chan_flux = base_flux*(nu(chan)/obsfreq)**(-2.-srcalpha(src))

! The source flux can appear in 4 possible pixels in the output map
               szsky(out_ra_pix,out_dec_pix,chan) = &
                    szsky(out_ra_pix,out_dec_pix,chan)+ &
                    frac_ra*frac_dec*chan_flux
               szsky(out_ra_pix+1,out_dec_pix,chan) = &
                    szsky(out_ra_pix+1,out_dec_pix,chan)+ &
                    (1.-frac_ra)*frac_dec*chan_flux
               szsky(out_ra_pix,out_dec_pix+1,chan) = &
                    szsky(out_ra_pix,out_dec_pix+1,chan)+ &
                    frac_ra*(1.-frac_dec)*chan_flux
               szsky(out_ra_pix+1,out_dec_pix+1,chan) = &
                    szsky(out_ra_pix+1,out_dec_pix+1,chan)+ &
                    (1.-frac_ra)*(1.-frac_dec)*chan_flux
            end do
         end if

!         radist(src)  = (ra/hr2rad)-rapos(src)
!         decdist(src) = decpos(src)-(dec/deg2rad) 
!         temp1 = 0
!         temp2 = 0
!         temp1 = int(radist(src)*15.*cos(dec)*3600./cellsize)+(maxsize/2)
!         temp2 = int(decdist(src)*3600./cellsize)+(maxsize/2)
!         src_x(src) = temp1
!         src_y(src) = temp2
!         radist(src) = radist(src)*hr2rad
!         decdist(src) = decdist(src)*deg2rad

! change from Jy to brightness temperature:
!         if ((temp1.gt.0.).and.(temp1.le.maxsize)) then
!            if ((temp2.gt.0.).and.(temp2.le.maxsize)) then
!               do chan = 1, nchan
!                  szsky(temp1,temp2,chan) = srcflux(src)/&
!                   ((cellsize*sec2rad)**2*2.0*k_b*obsfreq**2/const_c2*1.0D26)&
!                       *(nu(chan)/obsfreq)**(-2.-srcalpha(src))
!               end do
!            end if
!         end if
      end do      
   end select what_sky

   if (verbose) then
      write(*,*) 'Maximum radio brightness  ',maxval(szsky), ' K'
   end if

   if (do_plot) then
      write(*,*) 'Displaying radio sky'
      chan = nchan/2+1
      call display_map(szsky(:,:,chan),maxsize,graph_adj,&
                       cellsize,cellsize,0D0,0D0)
   end if

   status = 0

end subroutine make_source_population
