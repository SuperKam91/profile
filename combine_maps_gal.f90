! Reads in a list of .FITS maps and their associated rms coverages and 
! combines them together 
! Currently reads in FITS maps, puts them into square arrays and sums these 
! At some point in the future it may be necessary to work with rectagular
! maps all the way through.

subroutine combine_maps_gal(do_regrd)

   use kind_def
   use sz_globals
   use file_handling
   use fits_handling
   use map_handling
   use kgraph
   use maths
   use astronomy

   implicit none
   
   real(kind=dp), allocatable, dimension(:) :: fits_l, fits_b
   real(kind=dp), allocatable ,dimension(:,:) :: comb_map, comb_weight
   real(kind=dp) :: delta_ra, delta_dec, mapsize_ra, mapsize_dec
   real(kind=dp) :: ra_cent, dec_cent, in_ra, in_dec, frac_ra, frac_dec
   real(kind=dp) :: l_cent, b_cent, in_l, in_b
   real(kind=dp) :: cosdec, maxweight, cutlevel, scale, ra2, dec2
   real(kind=dp) :: sepn_ra, sepn_dec, mapsize, tmp_l, tmp_b, dx, dy
   real(kind=dp) :: min_fits_ra, max_fits_ra, min_fits_dec, max_fits_dec
   real(kind=dp) :: min_fits_l, max_fits_l, min_fits_b, max_fits_b
   real(kind=dp) :: min_l, max_l, min_b, max_b, mapsize_true
   real(kind=dp) :: min_ra, max_ra, min_dec, max_dec, max_sepn1, max_sepn2
   real(kind=dp) :: base_l_cent, base_b_cent, l_spac, b_spac
   real(kind=dp) :: sum_weight, tmp_ra, tmp_dec, sep, tol
   real(kind=dp), dimension(4) :: map_ra, map_dec
   integer :: l_dir, l_root, l_ext, imap_l, imap_b, nmaps_l, nmaps_b
   integer :: npix_ra, npix_dec, ra_pix_cent, dec_pix_cent, l_file
   integer :: out_ra_pix, out_dec_pix, iunit, bksize
   character(len=fname_len) :: basename, filename, noise_ext
   character*80  :: ftcom
   integer :: fcount, flen, i, j, len_noise_ext
   logical :: debug, debug2, doall, do_regrd
   character :: chra*16, chdec*16, chl*16, chb*16
   integer :: prra, prdec, lr, ld, ll, lb

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

   if (allocated(fits_l)) deallocate(fits_l)
   if (allocated(fits_b)) deallocate(fits_b)
   allocate(fits_l(n_fits_files))
   allocate(fits_b(n_fits_files))
   do fcount = 1, n_fits_files
     call sla_eqgal(dble(fits_ra(fcount)),dble(fits_dec(fcount)),fits_l(fcount),fits_b(fcount))
   end do
! NB sla_eqgal MUST have double precision inputs/outputs otherwise it doesn't work!

! Find extent of FITS files
   min_fits_ra = minval(fits_ra)
   max_fits_ra = maxval(fits_ra)
   min_fits_dec = minval(fits_dec)
   max_fits_dec = maxval(fits_dec)
   min_fits_l = minval(fits_l)
   max_fits_l = maxval(fits_l)
   min_fits_b = minval(fits_b)
   max_fits_b = maxval(fits_b)

! Extension for associated noise files
   noise_ext = '.noise.fits'
   call io_getc('Extension for noise files, or q to quit here','*',noise_ext,status)
   len_noise_ext = chr_lenb(noise_ext)
! Quit here after just writing the header digest
   if (trim(noise_ext).eq.'q') goto 50

! Get cutoff level
   call io_getd('Fractional weight cutoff level','*',cutlevel, status)

! At present output map is definitely square
   call io_getd('Output map size (deg)','6.0d0',mapsize,status)
   call io_getd('Spacing of maps in l (deg)', '6.0', l_spac, status)
   l_spac = l_spac*deg2rad
   call io_getd('Spacing of maps in b (deg)', '6.0', b_spac, status)
   b_spac = b_spac*deg2rad
   mapsize=mapsize*deg2rad

! Get number of output maps desired
   nmaps_l = (max_fits_l-min_fits_l)/l_spac+1
   call io_geti('Number of maps along l axis','*',nmaps_l,status)
   nmaps_b = (max_fits_b-min_fits_b)/b_spac+1
   call io_geti('Number of maps along b axis','*',nmaps_b,status)

! Calculate map sizes - now need to calculate true map sizes in RA, dec.
! Need to get central coords for first map first

! Generate coordinates for first map, suggesting a rounded value
   base_l_cent = (min_fits_l+mapsize/2.d0)/deg2rad
   call io_getd('Central l for first map', '*', base_l_cent, status)
   base_l_cent = base_l_cent*deg2rad
   base_b_cent = (min_fits_b+mapsize/2.d0)/deg2rad
   call io_getd('Central b for first map', '*', base_b_cent, status)
   base_b_cent = base_b_cent*deg2rad
   
   mapsize_ra = maxval(fits_npix_ra*fits_incell)
   mapsize_dec = maxval(fits_npix_dec*fits_incell)

! Mapsize will in general be different for each map
   do imap_l = 1, nmaps_l
     do imap_b = 1, nmaps_b
        l_cent = base_l_cent + l_spac*(imap_l-1)
        b_cent = base_b_cent + b_spac*(imap_b-1)
        call sla_galeq(l_cent, b_cent, ra_cent, dec_cent)
        write(*,*) l_cent, b_cent, ra_cent, dec_cent
        min_l = l_cent - mapsize/2.d0
        max_l = l_cent + mapsize/2.d0
        min_b = b_cent - mapsize/2.d0
        max_b = b_cent + mapsize/2.d0
        call sla_galeq(min_l, min_b, map_ra(1), map_dec(1))
        call sla_galeq(min_l, max_b, map_ra(2), map_dec(2))
        call sla_galeq(max_l, min_b, map_ra(3), map_dec(3))
        call sla_galeq(max_l, max_b, map_ra(4), map_dec(4))
        min_ra=minval(map_ra)
        max_ra=maxval(map_ra)
        do while (abs(min_ra-max_ra+2.d0*pi).lt.abs(max_ra-min_ra))
          map_ra(maxloc(map_ra)) = map_ra(maxloc(map_ra))-2.d0*pi
          min_ra=minval(map_ra)
          max_ra=maxval(map_ra)
        end do
        min_dec=minval(map_dec)
        max_dec=maxval(map_dec)
        call calc_sepn(min_ra, min_dec, ra_cent, min_dec, max_sepn1)
        call calc_sepn(max_ra, min_dec, ra_cent, min_dec, max_sepn2)
        mapsize_true = max(max_sepn1, max_sepn2, max_dec-dec_cent, dec_cent-min_dec)
        mapsize_true = 2.d0*mapsize_true
        write(*,*) 'Mapsize in RA, dec is ', mapsize_true/deg2rad
        npix_ra = mapsize_true/deg2rad*3600.d0/cellsize+1
        npix_dec = npix_ra
        ra_pix_cent = npix_ra/2+1
        dec_pix_cent = npix_dec/2+1

        if (verbose) then
           write(*,*) 'Creating ',npix_ra,' by ',npix_dec,' output map'
        end if

! Allocate arrays
        allocate(comb_map(npix_ra,npix_dec))
        allocate(comb_weight(npix_ra,npix_dec))

        write(*,*) '-------------------'

! Construct rootname for this output map
        write(chl, "(F16.2)"), l_cent/deg2rad
        write(chb, "(F16.2)"), b_cent/deg2rad

        if (b_cent.ge.0.0) then
          rootname = 'G'//trim(adjustl(chl))//'+'//trim(adjustl(chb))
        else
          rootname = 'G'//trim(adjustl(chl))//trim(adjustl(chb))
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
              call sla_eqgal(tmp_ra-mapsize_ra/2, tmp_dec-mapsize_dec/2, tmp_l, tmp_b)
              if ((tmp_l.gt.min_l).and.(tmp_l.lt.max_l).and.&
                 (tmp_b.gt.min_b).and.(tmp_b.lt.max_b)) then
                 fits_use_file(fcount) = .true.
              end if
              call sla_eqgal(tmp_ra+mapsize_ra/2, tmp_dec+mapsize_dec/2, tmp_l, tmp_b)
              if ((tmp_l.gt.min_l).and.(tmp_l.lt.max_l).and.&
                 (tmp_b.gt.min_b).and.(tmp_b.lt.max_b)) then
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
           prra = 0
           prdec = 0
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

! Read in input rms
               filename = basename(1:flen)//noise_ext(1:len_noise_ext)
               p_sky => sky_weight
               call read_fits_map(filename)
               if (debug) then
                  write(*,*) 'rms ',fcount
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
                                           
! Work out if the pixel falls inside the l, b ranges
                          call sla_eqgal(in_ra, in_dec, in_l, in_b)
                          if ((in_l.lt.min_l).or.(in_l.gt.max_l).or.&
                             (in_b.lt.min_b).or.(in_b.gt.max_b)) then
                              cycle
                          end if
                        
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
                          out_ra_pix = ra_pix_cent - (fits_crpix_ra(fcount)-i)
                          out_dec_pix = dec_pix_cent - (fits_crpix_dec(fcount)-j)
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

! Make sure edges of the output map are properly blanked
        do i = 1, npix_ra
          do j = 1, npix_dec
            dx = (i-(npix_ra/2+1))*1.d0
            dy = (j-(npix_dec/2+1))*1.d0
            ra2 = ra_cent
            dec2 = dec_cent
            call sin_proj_back(in_ra,in_dec,ra2,dec2,&
                                           cellsize,dx,dy)
            call sla_eqgal(in_ra, in_dec, in_l, in_b)
            if ((in_l.lt.min_l).or.(in_l.gt.max_l).or.&
                 (in_b.lt.min_b).or.(in_b.gt.max_b)) then
                comb_map(i,j) = 0.d0
                comb_weight(i,j) = 0.d0
            end if
          enddo
        enddo

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

         deallocate(comb_map)
         deallocate(comb_weight)
      end do
   end do

! Do deallocation
   deallocate(fits_l)
   deallocate(fits_b)
   call deallocate_combine_fits

50 continue

end subroutine combine_maps_gal


