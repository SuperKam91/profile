! Reads in a list of .FITS maps and their associated rms coverages and 
! combines them together 
! Currently reads in FITS maps, puts them into square arrays and sums these 
! At some point in the future it may be necessary to work with rectagular
! maps all the way through.

subroutine combine_maps(do_regrd)

   use kind_def
   use sz_globals
   use file_handling
   use fits_handling
   use map_handling
   use kgraph
   use maths
   use astronomy

   implicit none
   
   real(kind=dp), allocatable ,dimension(:,:) :: comb_map, comb_weight
   real(kind=dp) :: delta_ra, delta_dec, mapsize_ra, mapsize_dec
   real(kind=dp) :: ra_cent, dec_cent, in_ra, in_dec, frac_ra, frac_dec
   real(kind=dp) :: cosdec, maxweight, cutlevel, scale, ra2, dec2
   real(kind=dp) :: sepn_ra, sepn_dec, mapsize, tmp_ra, tmp_dec, dx, dy
   real(kind=dp) :: min_fits_ra, max_fits_ra, min_fits_dec, max_fits_dec
   real(kind=dp) :: base_ra_cent, base_dec_cent, ra_spac, dec_spac
   real(kind=dp) :: sum_weight, sep, tol
   integer :: l_dir, l_root, l_ext, imap_ra, imap_dec, nmaps_ra, nmaps_dec
   integer :: npix_ra, npix_dec, ra_pix_cent, dec_pix_cent, l_file
   integer :: out_ra_pix, out_dec_pix, iunit, bksize
   character(len=fname_len) :: basename, filename, noise_ext
   character*80  :: ftcom
   integer :: fcount, flen, i, j, len_noise_ext
   logical :: debug, debug2, doall, do_regrd
   character :: chra*16, chdec*16
   integer :: prra, prdec, lr, ld

   integer, external :: chr_lenb
   logical, external :: io_yesno

   verbose = .true.
   debug = .false.
   debug2 = .true.
   cutlevel = 0.00001d0
   prra = 1
   prdec = 0

! Check whether imsize is correct
   if (maxsize.ne.128) then
      write(*,*) 'NB imsize not 128 - recommend that you run get-geometry'
   end if

! Read in list of FITS maps to be processed produced by
! ls *.????.fits | awk ' {print $1 "." $2 }' FS="." > test.flist
   call read_flist
   write(*,*) n_fits_files,' maps available to be combined'

   if (debug) then 
      do fcount = 1, n_fits_files
         filename = fits_list(fcount)
         flen = chr_lenb(filename)
         write(*,*) fcount, filename(1:flen)
      end do
   end if

! Read in headers or use previously prepared digest?
   if (io_yesno('Does a header digest exist','yes',status)) then
      call read_header_digest
   else

! Read in headers for all FITS maps
      write(*,*) 'Reading headers'
      do fcount = 1, n_fits_files
         filename = fits_list(fcount)
         flen = chr_lenb(filename)
         filename = filename(1:flen)//'.fits'
         flen = chr_lenb(filename)
         call read_fits_header(filename,fits_npix_ra(fcount),&
              fits_npix_dec(fcount),fits_ra(fcount),fits_dec(fcount),&
              fits_bmaj(fcount), fits_bmin(fcount), fits_bpa(fcount),&
              fits_incell(fcount))
         if (n_fits_files.gt.1000) then
            if (1000*(fcount/1000).eq.fcount) then
               write(*,*) 'Reading header for file ',fcount
            end if
         end if
         if (verbose) then
            call chr_chdtos(fits_ra(fcount)/hr2rad,prra,chra,lr)
            call chr_chdtos(fits_dec(fcount)/deg2rad,prdec,chdec,ld)
            write(*,*) fcount, filename(1:flen)
            write(*,*)' at ', chra(1:lr),' ',chdec(1:ld)
         end if
      end do

! Write a digest of the header information for the future
      call write_header_digest
   end if

! If not regridding, need the actual CRPIX and CRVAL values of the (pre-regridded) maps
   if (.not.do_regrd) then
     do fcount = 1, n_fits_files
         filename = fits_list(fcount)
         flen = chr_lenb(filename)
         filename = filename(1:flen)//'.fits'
         flen = chr_lenb(filename)
         fits_crpix_ra(fcount) = 0
         fits_crpix_dec(fcount) = 0
         status=0
         call io_nxtlun(iunit,status)
! Open the FITS file with readonly access
         call ftopen(iunit,filename,0,bksize,status)
! Check if file found correctly
         if (status.ne.0) then
            status = 0
            write(*,*) filename(1:flen),' file not found'
            write(*,*)
            goto 20
         end if
! Read the required keywords
         status = 0
         call ftgkye(iunit,'CRPIX1',fits_crpix_ra(fcount),ftcom,status)
         status = 0
         call ftgkye(iunit,'CRPIX2',fits_crpix_dec(fcount),ftcom,status)
         status = 0
         call ftgkye(iunit,'CRVAL1',fits_ra1(fcount),ftcom,status)
         fits_ra1(fcount) = fits_ra1(fcount)*deg2rad
         status = 0
         call ftgkye(iunit,'CRVAL2',fits_dec1(fcount),ftcom,status)
         fits_dec1(fcount) = fits_dec1(fcount)*deg2rad
! Close the file
         status = 0
         call ftclos(iunit,status)
 20      continue
     enddo
   endif

! Find extent of FITS files
   min_fits_ra = minval(fits_ra)
   max_fits_ra = maxval(fits_ra)
   min_fits_dec = minval(fits_dec)
   max_fits_dec = maxval(fits_dec)

! Extension for associated noise files
   noise_ext = '.noise.fits'
   call io_getc('Extension for noise files or q to quit','*',noise_ext,status)
   len_noise_ext = chr_lenb(noise_ext)
! Quit here after just writing the header digest
   if (trim(noise_ext).eq.'q') goto 50

! Get cutoff level
   call io_getd('Fractional weight cutoff level','*',cutlevel, status)

! At present output map is definitely square
   call io_getd('Output map size (deg)','6.0d0',mapsize,status)
   call io_getd('Spacing of maps in RA (mins)','25.0d0',ra_spac,status)
   ra_spac = ra_spac*hr2rad/60.d0
   call io_getd('Spacing of maps in Dec (deg)','6.0d0',dec_spac,status)
   dec_spac = dec_spac*deg2rad

! Get number of output maps desired
   nmaps_ra = (max_fits_ra-min_fits_ra)/ra_spac+1
   call io_geti('Number of maps along RA axis','*',nmaps_ra,status)
   nmaps_dec = (max_fits_dec-min_fits_dec)/dec_spac+1
   call io_geti('Number of maps along Dec axis','*',nmaps_dec,status)

! Calculate map sizes
   npix_ra = mapsize*3600.d0/cellsize+1
   npix_dec = npix_ra
   ra_pix_cent = npix_ra/2+1
   dec_pix_cent = npix_dec/2+1
   mapsize_ra = maxval(fits_npix_ra*fits_incell)
   mapsize_dec = maxval(fits_npix_dec*fits_incell)
   mapsize = mapsize*deg2rad

! Generate coordinates for first map, suggesting a rounded value
   base_ra_cent = int((min_fits_ra+mapsize/2.d0)/hr2rad*12.d0)/(12.d0/hr2rad)
   call io_getra('Central RA for first map','*',base_ra_cent,status)
   base_dec_cent = int((min_fits_dec-mapsize/2.d0)/deg2rad*6.d0)/(6.d0/deg2rad)
   call io_getdec('Central Dec for first map','*',base_dec_cent,status)

   if (verbose) then
      write(*,*) 'Creating ',npix_ra,' by ',npix_dec,' output maps'
   end if

! Allocate arrays
   allocate(comb_map(npix_ra,npix_dec))
   allocate(comb_weight(npix_ra,npix_dec))

   write(*,*) '-------------------'

! Loop over output maps
   do imap_ra = 1, nmaps_ra
      do imap_dec = 1, nmaps_dec

! Find centre of this output map
         ra_cent = base_ra_cent+(imap_ra-1)*ra_spac
         dec_cent = base_dec_cent+(imap_dec-1)*dec_spac

! Take care of wrapping round from 24hr -> 0hr
         if (ra_cent.gt.2.d0*pi) then
            ra_cent = ra_cent-2.d0*pi
         end if
            
! Construct rootname for this output map
         prra = 0
         prdec = 0
         call chr_chdtos(ra_cent/hr2rad,prra,chra,lr)
         call chr_chdtos(dec_cent/deg2rad,prdec,chdec,ld)
         if (dec_cent.gt.0.0) then
            rootname = 'J'//chra(3:4)//chra(6:7)//'+'//chdec(3:4)//chdec(6:7)
         else
            rootname = 'J'//chra(3:4)//chra(6:7)//'-'//chdec(3:4)//chdec(6:7)
         end if

         l_root = chr_lenb(rootname)
         write(*,*) 'Generating map centred at ',rootname(1:l_root)

! Determine which of the available input FITS files to use
         fits_use_file = .false.
         ra2 = ra_cent
         dec2 = dec_cent
         do fcount = 1, n_fits_files
! If regridding, check that the centre of the input map falls within the output raster
            if (do_regrd) then
              tmp_ra = fits_ra(fcount)
              tmp_dec = fits_dec(fcount)
              call sin_proj_forward(tmp_ra,tmp_dec,ra2,dec2,cellsize,dx,dy)
              if ((abs(dx).lt.(mapsize+mapsize_ra)/(2.*cellsize*sec2rad)).and.&
                 (abs(dy).lt.(mapsize+mapsize_dec)/(2.*cellsize*sec2rad))) then
                 fits_use_file(fcount) = .true.
              end if
! If not regridding, check that the file has been regridded to the correct centre
            else
              tol = cellsize/10.d0*sec2rad
              tmp_ra = fits_ra1(fcount)
              tmp_dec = fits_dec1(fcount)
              call calc_sepn(tmp_ra, tmp_dec, ra2, dec2, sep)
              if (sep.lt.tol) then
                fits_use_file(fcount) = .true.
              endif
            endif
         end do
         write(*,*) 'Using ',count(fits_use_file),' out of ',&
                    n_fits_files,' files'
         if (verbose) then
            write(*,*)
            call chr_chdtos(ra_cent/hr2rad,prra,chra,lr)
            call chr_chdtos(dec_cent/deg2rad,prdec,chdec,ld)
            write(*,*) 'Map centre at ',chra(1:lr),chdec(1:ld),&
                 ' at pixel ',ra_pix_cent,dec_pix_cent
         end if

! Initialise output maps
         comb_map = 0.0
         comb_weight = 0.0

! Loop over all input maps 
         verbose = .false.
         do fcount = 1, n_fits_files
            if (fits_use_file(fcount)) then
               basename = fits_list(fcount)
               flen = chr_lenb(basename)
               sum_weight = 0.0d0

! Read in input rms - note need to read rms *first* because otherwise temp_sky is overwritten in read_fits_map...
               filename = basename(1:flen)//noise_ext(1:len_noise_ext)
               p_sky => sky_weight
               call read_fits_map(filename)
               if (debug) then
                  write(*,*) 'rms ',fcount
                  call display_map(p_sky,maxsize,graph_adj,&
                       cellsize,cellsize,0D0,0D0)         
               end if

! Read in input map
               filename = basename(1:flen)//'.fits'
               if (debug2) then
                  write(*,*) 'Processing file ',fcount
                  write(*,*) '  ',basename(1:flen)
               end if
               p_sky => temp_sky
               call read_fits_map(filename)
               if (debug) then
                  write(*,*) 'Map ',fcount
                  call display_map(p_sky,maxsize,graph_adj,&
                       cellsize,cellsize,0D0,0D0)         
               end if

! Fill fits_noise entry
               fits_noise(fcount) = sky_weight(maxsize/2+1,maxsize/2+1)

! Convert into weight
! A rather horrid way to convert to do check that rms is not zero
!               where (sky_weight.gt.sky_weight(maxsize/2,maxsize/2)/2.)
!                  sky_weight = 1.d0/sky_weight**2
!               elsewhere
!                  sky_weight = 0.0
!               end where
! For drift scan maps the lowest noise pixel is *not necessarily* the centre!  What's wrong with just checking that it's not 0?
               where (sky_weight.ne.0.0)
                  sky_weight = 1.d0/sky_weight**2
               end where

! Plot input maps for debugging
               if (debug) then
                  write(*,*) 'weight ',fcount
                  call display_map(p_sky,maxsize,graph_adj,&
                       cellsize,cellsize,0D0,0D0)      
                  write(*,*) 'flux map'
                  call display_map(temp_sky*sky_weight,maxsize,graph_adj,&
                       cellsize,cellsize,0D0,0D0)                      
               end if

! Loop over all pixels in input map
               do i = 1, maxsize
                  do j = 1, maxsize

                     if (sky_weight(i,j).ne.0.0) then
                        if (do_regrd) then
! Find RA and Dec of this pixel in the input map
                          dx = (i-(maxsize/2+1))*1.d0
                          dy = (j-(maxsize/2+1))*1.d0
                          ra2 = fits_ra(fcount)
                          dec2 = fits_dec(fcount)
                          call sin_proj_back(in_ra,in_dec,ra2,dec2,&
                                           cellsize,dx,dy)

! This pixel can appear in 4 possible pixels in the output map
! Find bottom left pixel in output map to which this corresponds and
! find what fraction of input pixel to place in each of the 4 output pixels
                          ra2 = ra_cent
                          dec2 = dec_cent
                          call sin_proj_forward(in_ra,in_dec,ra2,dec2,&
                                              cellsize,dx,dy)

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
                             write(*,*) 'Error in combine map fraction calc'
                             write(*,*) 'Fractions ',frac_ra,frac_dec
                             write(*,*) 'Pixel offsets from projection ',dx, dy
                             write(*,*) 'Actual pix offs ',out_ra_pix,out_dec_pix
                             write(*,*) 'Centre of map ',ra_pix_cent,dec_pix_cent
                             write(*,*) 'Pixel position ',in_ra, in_dec
                             write(*,*) 'Map centre ',ra_cent, dec_cent
                          end if

! Add into output map
                          if ((out_ra_pix.gt.1).and.&
                            (out_ra_pix.lt.npix_ra-1).and.&
                            (out_dec_pix.gt.1).and.&
                            (out_dec_pix.lt.npix_dec-1)) then
                             comb_map(out_ra_pix,out_dec_pix) = &
                                comb_map(out_ra_pix,out_dec_pix)+ &
                                frac_ra*frac_dec* &
                                temp_sky(i,j)*sky_weight(i,j) 
                             comb_weight(out_ra_pix,out_dec_pix) = &
                                comb_weight(out_ra_pix,out_dec_pix)+ &
                                frac_ra*frac_dec* &
                                sky_weight(i,j)
                             comb_map(out_ra_pix+1,out_dec_pix) = &
                                comb_map(out_ra_pix+1,out_dec_pix)+ &
                                (1.-frac_ra)*frac_dec* &
                                temp_sky(i,j)*sky_weight(i,j)
                             comb_weight(out_ra_pix+1,out_dec_pix) = &
                                comb_weight(out_ra_pix+1,out_dec_pix)+ &
                                (1.-frac_ra)*frac_dec* &
                                sky_weight(i,j)
                             comb_map(out_ra_pix,out_dec_pix+1) = &
                                comb_map(out_ra_pix,out_dec_pix+1)+ &
                                frac_ra*(1.-frac_dec)* &
                                temp_sky(i,j)*sky_weight(i,j)
                             comb_weight(out_ra_pix,out_dec_pix+1) = &
                                comb_weight(out_ra_pix,out_dec_pix+1)+ &
                                frac_ra*(1-frac_dec)* &
                                sky_weight(i,j)
                             comb_map(out_ra_pix+1,out_dec_pix+1) = &
                                comb_map(out_ra_pix+1,out_dec_pix+1)+ &
                                (1.-frac_ra)*(1.-frac_dec)* &
                                temp_sky(i,j)*sky_weight(i,j)
                             comb_weight(out_ra_pix+1,out_dec_pix+1) = &
                                comb_weight(out_ra_pix+1,out_dec_pix+1)+ &
                                (1.-frac_ra)*(1.-frac_dec)* &
                                sky_weight(i,j)
                             sum_weight = sum_weight+sky_weight(i,j)
                          end if
                        else
                          ! No regridding required, just add into the output grid
                          out_ra_pix = npix_ra/2+1 - (fits_crpix_ra(fcount)-i)
                          out_dec_pix = npix_dec/2+1 - (fits_crpix_dec(fcount)-j)
                          if ((out_ra_pix.ge.1).and.(out_ra_pix.le.npix_ra).and. &
                            (out_dec_pix.ge.1).and.(out_dec_pix.le.npix_dec)) then
                            comb_map(out_ra_pix,out_dec_pix) = &
                              comb_map(out_ra_pix,out_dec_pix)+temp_sky(i,j)*sky_weight(i,j)
                            comb_weight(out_ra_pix,out_dec_pix) = &
                              comb_weight(out_ra_pix,out_dec_pix) + sky_weight(i,j)
                            sum_weight = sum_weight + sky_weight(i,j)
                          endif
                        end if
                     end if
                  end do
               end do
               write(*,*) 'Weight contributed by this map ',sum_weight

! If this map is not actually contributing at all then reset its fits_use_file
! which will be used in the output FITS table
               if (sum_weight.eq.0.0d0) then
                  fits_use_file(fcount) = .false.
               end if

            end if
         end do

! Do plotting if requested
         if (do_plot) then
            write(*,*) 'Combined map'
            p_sky => temp_sky
            scale = 1.d0
            call subim(comb_map,p_sky,cellsize,cellsize,&
                       npix_ra,npix_dec,maxsize,scale,.true.)
            call display_map(p_sky,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
            write(*,*) 'Summed Weights'
            call subim(comb_weight,p_sky,cellsize,cellsize,&
                       npix_ra,npix_dec,maxsize,scale,.true.)
            call display_map(p_sky,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
         end if

! Convert comb_map into Jy and comb_weight into an rms map
         maxweight = maxval(comb_weight)
         where (comb_weight.gt.maxweight*cutlevel)
            comb_map = comb_map/comb_weight
            comb_weight = 1./sqrt(comb_weight)
         elsewhere
            comb_map = MAGICBLANK
            comb_weight = MAGICBLANK
         end where 

! Do plotting if requested
         if (do_plot) then
            write(*,*) 'Signal map'
            p_sky => temp_sky
            scale = 1.d0
            call subim(comb_map,p_sky,cellsize,cellsize,&
                       npix_ra,npix_dec,maxsize,&
                 scale,.true.)
            call display_map(p_sky,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
            write(*,*) 'rms'
            call subim(comb_weight,p_sky,cellsize,cellsize,&
                       npix_ra,npix_dec,maxsize,&
                 scale,.true.)
            call display_map(p_sky,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
         end if

! Write out FITS file
         obs_raref = ra_cent/hr2rad
         obs_decref = dec_cent/deg2rad

! Get filename for map
         extn = '.fits'
         l_dir  = chr_lenb(profile_data_dir)
         l_root = chr_lenb(rootname)
         l_ext  = chr_lenb(extn)
         filename = profile_data_dir(1:l_dir)//rootname(1:l_root)//&
                    extn(1:l_ext)
         l_file = chr_lenb(filename)
         if (count(fits_use_file).gt.0) then
             write(*,*) 'Writing ',filename(1:l_file)
             write_fits_beam_table = .true.
             call write_fits_map(comb_map,npix_ra,npix_dec,cellsize,filename)

! Get filename for rms
             extn = noise_ext(1:len_noise_ext)
             l_ext  = chr_lenb(extn)
             filename = profile_data_dir(1:l_dir)//rootname(1:l_root)//&
                    extn(1:l_ext)
             l_file = chr_lenb(filename)
             write(*,*) 'Writing ',filename(1:l_file)
             call write_fits_map(comb_weight,npix_ra,npix_dec,cellsize,filename)
             write_fits_beam_table = .false.
         else
             write(*,*) 'No data for ',rootname(1:l_root),', not writing files'
         endif
         write(*,*) '-------------------'
      end do
   end do

! Do deallocation
   deallocate(comb_map)
   deallocate(comb_weight)
   call deallocate_combine_fits
   
50 continue

end subroutine combine_maps


