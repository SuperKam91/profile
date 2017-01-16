module fits_handling

! Module which contains routines for reading / writing FITS format images 
! and uv data sets. Currently depends on a whole set of F77 fits routines.

   use kind_def
   use sz_globals
   use file_handling

   implicit none

   public :: read_fits_map
   public :: write_fits_map
   public :: read_uv_fits
   public :: write_uv_fits
   public :: write_multi_fits
   public :: read_daisuke_map
   public :: read_fits_header
   public :: printerror
   public :: write_beam_table

contains

!***************************************************************************

   subroutine read_fits_map(file)

      use kind_def
      use sz_globals
      use map_handling

      implicit none

      character(len=fname_len), intent(in) :: file
      character*1 :: chr1
      real(kind=dp), dimension(:,:), allocatable :: fits_map
      real, dimension(:,:), allocatable :: ivalue
      character*80 :: ftcom
      real(kind=dp) :: scale, incell
      integer:: iunit, nx, ny
      integer, parameter :: maxdim=99
      integer, dimension(maxdim) :: naxes
      integer :: group,fpixel,i,j,rwstat,bksize
      real    :: obsra, obsdec, delta_ra, delta_dec
      logical :: simple, extend, anyflg, surface_brightness
      
      status=0
      call io_nxtlun(iunit,status)

! Open the FITS file with readonly access
      rwstat=0

      call ftopen(iunit,file,rwstat,bksize,status)

! Check if file found correctly
      if (status.ne.0) then
         status = 0
         write(*,*) file,' file not found'
         write(*,*)
         return
      end if
! Read the required primary array keywords
! ftghpr is redundant - use ftghpr instead
!      call ftgprh(iunit,simple,bitpix,naxis,naxes,pcount,gcount,extend,status)

      call ftghpr(iunit,maxdim,simple,bitpix,naxis,naxes,pcount,gcount,&
                  extend,status)


! Get size of map to be read
      nx = naxes(1)
      ny = naxes(2)

      allocate (fits_map(nx,ny))
      allocate (ivalue(nx,ny))

      status = 0
      call ftgkye(iunit,'CRVAL1',obsra,ftcom,status)
      status = 0
      call ftgkye(iunit,'CRVAL2',obsdec,ftcom,status)
      status = 0
      call ftgkye(iunit,'CDELT1',delta_ra,ftcom,status)
      status = 0
      call ftgkye(iunit,'CDELT2',delta_dec,ftcom,status)
      status = 0
      !!call ftgkys(iunit,'TELESCOP',tel_tmp,ftcom,status)
      

      if (verbose) then
         write(*,*) 'Map at ',obsra,obsdec
         write(*,*) nx,' by ',ny,' pixel map'
         write(*,*) 'deg/pixel ',delta_ra, delta_dec
      end if
      if (delta_dec.eq.0.0) then
         write(*,*) 'Unable to read cellsize for FITS image - using ',cellsize
         incell = cellsize
      else
         incell = abs(delta_dec)*3600.d0
      end if

! Find whether incell and cellsize are the same to within 0.01% and if so
! set them to be exactly equal
      if (abs(cellsize-incell)/cellsize.lt.0.0001) then
         incell = cellsize
      end if

      if (abs(delta_ra).ne.(abs(delta_dec))) then
         write(*,*) 'NB coordinate increments not the same in RA and Dec'
      end if

      if (verbose) write(*,*) 'FITS map cellsize ',real(incell)

      surface_brightness = .true.
      if (real(incell).ne.real(cellsize)) then
         write(*,*) real(incell), real(cellsize)
         write(*,*) 'Converting to different pixel size'
         write(*,*) 'Is this a map of surface brightness (i.e. keep the & 
                     & mean the same'
         write(*,*) 'during conversion) or of flux (i.e. keep total number &
                     & of counts the same'
         call io_getc('Select (F)lux or (S)urface Brightness','s',chr1,status)
         select case(chr1)
         case('f','F')
            surface_brightness = .false.
         case default
            surface_brightness = .true.
         end select
      end if

      group = 0
      fpixel = 1

! Read the data
      status = 0
      call ftg2de(iunit,group,0.0,nx,naxes(1),naxes(2),ivalue,&
           & anyflg,status)

! Close the file
      status = 0
      call ftclos(iunit,status)

! Write the data into the array to be returned
      fits_map = ivalue

! Set BLANKed regions to 0
      where (fits_map.ne.fits_map) fits_map=0.0

      scale = 1.0
      call subim(fits_map, temp_sky, incell, cellsize, nx, ny, &
                 maxsize,scale,surface_brightness)
!      write(*,*) 'Exiting subim'
      deallocate (fits_map)
      deallocate (ivalue)

! Set BLANKed regions to 0
      where (temp_sky.ne.temp_sky) temp_sky=0.0

! Put map into all channels if required
      if (flag_nd.eq.3) then
        do i = 1, nchan
          p3_sky(:,:,i) = temp_sky
        end do
      else
        p_sky = temp_sky
      end if

   end subroutine read_fits_map

!**************************************************************************

   subroutine write_fits_map(array,nx,ny,cell,filename)!, keep_unit, iunit)

      use kind_def

      implicit none
	      
      integer  :: nx, ny
      real(kind=dp) :: array(nx,ny)
      real(kind=dp) :: cell
      real(kind=dp) :: nullval
      real(kind=dp) :: xrefpix, yrefpix, bscale, bzero, datamax, datamin
      real(kind=dp) :: ra, dec, zero, one, skycell, equinox, rf_freq, dval1
      integer :: unit, blocksize, group, fpixel, nelements, firstpix
      integer :: bitpix, naxis, decimals
      integer, allocatable :: naxes(:)
      logical :: simple, extend, ex, lopen
      character(len=fname_len) :: filename
      character*32 :: comment
      character*16 :: kword
      character*16 :: stringa
      !logical, optional, intent(in) :: keep_unit
      !integer, optional, intent(out) :: iunit
      character*100 :: name

      status = 0
      blocksize = 1

! Initialize variables
      group =  1
      fpixel = 1
      nelements = nx*ny
      simple = .true.
      bitpix = -32
      naxis = 4
      allocate(naxes(naxis))

      naxes(1) = nx
      naxes(2) = ny
      naxes(3) = 1
      naxes(4) = 1

      extend = .true.
      firstpix = 1
      nullval = MAGICBLANK
      rf_freq = obsfreq

! Header (world coordinate system)
      ra = obs_raref*15.
      dec = obs_decref      
      skycell = real(cell/3600.0)
      xrefpix = nx/2+1
      yrefpix = ny/2+1
      zero = 0.0
      one = 1.0
      decimals = 9
      equinox = 2000.0
      datamax = maxval(array)
      datamin = minval(array)

      call ftgiou (unit, status)
      call ftinit (unit, filename, blocksize, status)
      if (status.ne.0) call printerror(status)

      call FTFLNM(unit, name, status)
      !call FTFLMD(unit, iunit, status)
      !write(*,*) name, iunit

! Write compulsory part of header:
      call ftphpr (unit, simple, bitpix, naxis, naxes, 0, 1, extend, status)
      if (status.ne.0) call printerror(status)

! Define data scaling.
      bscale = 1.0d0                                                          
      bzero = 0.0d0
      call ftpscl (unit, bscale, bzero, status)
      if (status.ne.0) call printerror(status)
	
! Write the optional keywords:
      if (tel_name.eq.'LA') kword='AMI-LA'
      if (tel_name.eq.'SA') kword='AMI-SA'
      call ftpkys(unit,'OBJECT', 'Simul','Source name', status)
!      call ftpkys(unit,'TELESCOP', 'Profile','Radio Telescope', status)
      
      call ftpkys(unit,'TELESCOP', kword,'Radio Telescope', status)
      call ftpkys(unit,'INSTRUME', kword,'Receiver', status)
      call ftpkys(unit,'OBSERVER', '','', status)

! Write optional keywords to the header
      
      comment = 'REAL = TAPE * BSCALE + BZERO'
      kword = 'BSCALE'
      call ftpkyd (unit, kword, bscale, decimals, comment, status)

      comment = ''
      kword = 'BZERO'
      call ftpkyd (unit, kword, bzero, decimals, comment, status)

      comment = 'Units of flux density'
      kword = 'BUNIT'
      stringa = 'JY/BEAM '
      call ftpkys (unit, kword, stringa, comment, status)

      comment = 'Epoch of RA DEC'
      kword = 'EQUINOX'
      dval1 = 2.0d3
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      comment = 'Antenna pointing RA'
      kword = 'OBSRA'
      dval1 = ra
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      comment = 'Antenna pointing Dec'
      kword = 'OBSDEC'
      dval1 = dec
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      comment = 'Maximum pixel value'
      kword = 'DATAMAX'
      call ftpkyd (unit, kword, datamax, decimals, comment, status)

      comment = 'Minimum pixel value'
      kword = 'DATAMIN'
      call ftpkyd (unit, kword, datamin, decimals, comment, status)

      comment = ''
      kword = 'CTYPE1'
      stringa = 'RA---SIN'
      call ftpkys (unit, kword, stringa, comment, status)

      kword = 'CRVAL1'
      dval1 = ra
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CDELT1'
      dval1 = -skycell
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CRPIX1'
      dval1 = xrefpix
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CROTA1'
      dval1 = 0.0d0
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CTYPE2'
      stringa  = 'DEC--SIN'
      call ftpkys (unit, kword, stringa, comment, status)

      kword = 'CRVAL2'
      dval1 = dec
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CDELT2' 
      dval1 = skycell
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CRPIX2'
      dval1 = yrefpix
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CROTA2'
      dval1 = 0.0d0
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CTYPE3'
      stringa = 'FREQ    '
      call ftpkys (unit, kword, stringa, comment, status)

      kword = 'CRVAL3'
      dval1 = obsfreq
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CDELT3'
      dval1 = 1.0
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CRPIX3'
      dval1 = 1.0d0
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CROTA3'
      dval1 = 0.0d0
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CTYPE4'
      stringa = 'STOKES  '
      call ftpkys (unit, kword, stringa, comment, status)

      kword = 'CRVAL4'
      dval1 = 1.0d0
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CDELT4'
      dval1 = 1.0d0
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CRPIX4'
      dval1 = 1.0d0
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

      kword = 'CROTA4'
      dval1 = 0.0d0
      call ftpkyd (unit, kword, dval1, decimals, comment, status)

! Write data to file
      call ftppnd (unit, group, fpixel, nelements, array, nullval, status)
      
! If necessary, write beam table
      if (write_fits_beam_table) then
         call write_beam_table(unit)
      end if
      
! Close FITS file      
!!$      if (present(keep_unit)) then
!!$         if (.not.keep_unit) then
!!$            write(*,*) 'Not keeping da unit open', unit
!!$            call ftclos(unit, status)
!!$            call ftfiou(unit, status)
!!$         endif
!!$      else
         call ftclos (unit, status)
         call ftfiou(unit, status)
!!$      end if
      
      if (status.ne.0) call printerror(status)
      !iunit = unit
      

! Deallocate array
       deallocate(naxes)

   end subroutine write_fits_map

!****************************************************************************

   subroutine read_uv_fits

      use kind_def
      use sz_globals

      implicit none

      character(len=fname_len) :: filename
      character*80 :: name, ftcom
      integer :: iunit
      integer, parameter :: maxdim=7      
      integer :: naxes(maxdim)
      integer :: group, i, j, rwstat, bksize, a, ifreq
      integer :: fpixels(maxdim), lpixels(maxdim), inc(maxdim)
      integer :: axes(6), npoints, count, b, igroup, l_tel
      real    :: ivalue(105), parms(6), dfreq, freq, u1, v1, w1
      real    :: obsra, obsdec, pb_ra, pb_dec
      logical :: simple, extend, anyflg, debug
      integer, external :: chr_lenb

      debug = .false.

      prompt = 'FITS filename'
      which_dir = profile_data_dir
      rootname = 'uv'
      extn = '.fits'
      call get_read_filename(filename)
      
      status = 0
      call io_nxtlun(iunit,status)

! open the existing FITS file with readonly access
      rwstat=0

      call ftopen(iunit,filename,rwstat,bksize,status)
      if (debug) write(*,*) 'status', status,' ;bksize ',bksize

! read the required primary array keywords
!      call ftgprh(iunit,simple,bitpix,naxis,naxes,pcount,gcount,&
!             extend,status)
      call ftghpr(iunit,maxdim,simple,bitpix,naxis,naxes,pcount,gcount,&
                  extend,status)
      if (debug) write(*,*)iunit,simple,bitpix,naxis,naxes,pcount,gcount,&
             extend, status

      call ftgkys(iunit,'telescop',name,ftcom,status)
      call ftgkye(iunit,'obsra',obsra,ftcom,status)
      status = 0
      call ftgkye(iunit,'obsdec',obsdec,ftcom,status)
      status = 0
      call ftgkye(iunit,'crval5',pb_ra,ftcom,status)
      status = 0
      call ftgkye(iunit,'crval6',pb_dec,ftcom,status)
      status = 0
      call ftgkye(iunit,'crval4',freq,ftcom,status)
      status = 0
      call ftgkye(iunit,'cdelt4',dfreq,ftcom,status)
      status = 0
      call ftgkye(iunit,'pscal1',uscale,ftcom,status)
      status = 0
      call ftgkye(iunit,'pscal2',vscale,ftcom,status)

! Set the pointing and phase centres, RAs stored in hrs
      obs_raref = obsra/15.
      obs_decref = obsdec
      obs_rarot = pb_ra/15.
      obs_decrot = pb_dec

      l_tel = chr_lenb(name)
      write(*,*) 'Telescope          ',name(1:l_tel)
      write(*,*) 'Chan 1 frequency  ',freq,' ;channel width ',dfreq
      write(*,*) 'Pointing centre   ',pb_ra,pb_dec
      write(*,*) 'Phase centre      ',obsra,obsdec
      write(*,*)

      do i = 1,naxis-1
         fpixels(i) = 1
         axes(i) = naxes(i+1)
         lpixels(i) = axes(i)
         inc(i) = 1
      end do

! Find number of data points per FITS group
      npoints = 1
      do i = 2,naxis
         if (debug) write(*,*) 'naxes ',i,' is ',naxes(i)
         npoints = npoints * naxes(i)
      end do

      write (*,*) npoints,' data points per group'

! Some checks
      if (naxis.eq.7) then 
         write(*,*) 'This appears to be a multi-pointing FITS file - &
                     &not supported'
      end if
      if (naxes(4).ne.nchan) then 
         write(*,*) &
       'WARNING: this FITS file has a different number of frequency channels',&
       'to the currently selected telescope; continuing but results might be',&
       'odd....'
      end if

      write(*,*) 'Reading uv data'
      b = 1
      do igroup = 1, gcount
         call ftggpe(iunit,igroup,1,naxis,parms,status)
         call ftgsve(iunit,igroup,naxis-1,axes,fpixels,lpixels,inc,0,&
               ivalue,anyflg,status)
         u1 = parms(1)
         v1 = parms(2)
         w1 = 0.
         do j = 1,npoints, 3
            ifreq = j/3+1
            u(b) = u1*nu(ifreq)
            v(b) = v1*nu(ifreq)
            w(b) = 0.0
            data_re(b) = ivalue(j)
            data_im(b) = ivalue(j+1)
            weight(b) = ivalue(j+2)
            rms(b) = 1./sqrt(weight(b))
            if (debug) then
               write(*,*) b,ifreq,u(b),v(b),data_re(b),data_im(b),weight(b)
            end if
            b = b+1
         end do
      end do

      b = b-1

      nvis = b
      write(*,*) 'Read ',nvis,' visibilities'
      write(*,*) sum(weight),' sum of weights'

      call ftclos(iunit,status)

   end subroutine read_uv_fits


!**************************************************************************

   subroutine read_daisuke_map(infile)

      use kind_def
      use sz_globals
      use cosmology
      use map_handling

      implicit none

      character(len=fname_len) :: infile
      real(kind=dp), dimension(:,:), allocatable :: daisuke_map
      real(kind=dp),dimension(:), allocatable :: temp
      real(kind=dp) :: scale, incell
      integer  :: i,j,k,nx(2)
      real     :: aexpn,dx(2)

      open(unit=15,file=infile)

      read(15,*) aexpn,nx(1),nx(2)
      read(15,*) dx(1),dx(2)

      if (nx(1).ne.nx(2)) then
         write(*,*) 'read_daisuke_map only for square input maps'
         return
      end if

      allocate(temp(nx(1)*nx(2)))
      allocate(daisuke_map(nx(1),nx(2)))

      call io_getd('Hubble parameter:','*',H0,status)
      call io_getd('Redshift:','*',z,status)
      call angdist
      write(*,*) 'Angular distance:',D_theta,'Mpc.'
      scale = 1.0

      incell = dx(1)/(D_theta*(H0/100.))
      incell = incell/sec2rad
      write(*,*) dx(1),'/h Mpc = ',incell,' arcsec'
!      write(*,*) 'Resizing map cells to:',cellsize,'arcsec.'

      if (aexpn.gt.1.0) then
         write(*,*) 'Expansion factor, a, reset from',aexpn,'to 1.0'
         aexpn=1.0
      end if

      write(*,*) 'Reading Daisuke file format...'
      j = 1
      k = -4
      do i = 1,nx(1)
         do j = 1,nx(2),5
            k = k+5
            if (k.lt.((nx(1)*nx(2))-4)) then
               read(15,*) temp(k),temp(k+1),temp(k+2),temp(k+3),temp(k+4)
            else
               read(15,*) temp(k),temp(k+1),temp(k+2),temp(k+3)
               goto 100
            end if
         end do
      end do

100   daisuke_map = 0.
      k = 0
      do i = 1,nx(1)
         do j = 1,nx(2)
            k = k+1
            daisuke_map(i,j) = temp(k)
         end do
      end do

      close(unit=15)

      call subim(daisuke_map, p_sky, incell, cellsize, nx(1), nx(2),&
                 maxsize,scale,.true.)

      deallocate(daisuke_map)
      deallocate(temp)

   end subroutine read_daisuke_map

!***************************************************************************

   subroutine read_fits_header(file,npix_ra,npix_dec,obsra,obsdec,&
                               bmaj, bmin, bpa, incell)

      use kind_def
      use sz_globals
      use map_handling
      use astronomy

      implicit none

      character(len=fname_len), intent(in) :: file
      character*80 :: ftcom
      real    :: obsra, obsdec, delta_ra, delta_dec, cent_ra, cent_dec
      real    :: bmaj, bmin, bpa 
      real(kind=dp) :: incell, obsra1, obsdec1, dx, dy
      integer:: iunit, npix_ra, npix_dec
      integer, parameter :: maxdim=99
      integer, dimension(maxdim) :: naxes
      integer :: group, fpixel, i, j, rwstat, bksize
      integer :: nkeys, nkeys2, nspace
      logical :: simple, extend, anyflg, beam_found
      integer, parameter :: nkeymax = 24
      character*24 :: keys(nkeymax),bmaj_string,bmin_string,bpa_string
      character*24 :: match      
      character*128 :: record   

      status=0
      call io_nxtlun(iunit,status)

! Open the FITS file with readonly access
      rwstat=0

      call ftopen(iunit,file,rwstat,bksize,status)

! Check if file found correctly
      if (status.ne.0) then
         status = 0
         write(*,*) file,' file not found'
         write(*,*)
         return
      end if

! Read the required primary array keywords
!      call ftgprh(iunit,simple,bitpix,naxis,naxes,pcount,gcount,extend,status)
      call ftghpr(iunit,maxdim,simple,bitpix,naxis,naxes,pcount,gcount,&
                  extend,status)

      npix_ra = naxes(1)
      npix_dec = naxes(2)

      status = 0
      call ftgkye(iunit,'CRVAL1',obsra,ftcom,status)
      status = 0
      call ftgkye(iunit,'CRVAL2',obsdec,ftcom,status)
      status = 0
      call ftgkye(iunit,'CRPIX1',cent_ra,ftcom,status)
      status = 0
      call ftgkye(iunit,'CRPIX2',cent_dec,ftcom,status)
      status = 0
      call ftgkye(iunit,'CDELT1',delta_ra,ftcom,status)
      status = 0
      call ftgkye(iunit,'CDELT2',delta_dec,ftcom,status)

! Find pixel size in radians
      if (delta_dec.eq.0.0) then
         write(*,*) 'Unable to read cellsize for FITS image - using ',cellsize
         incell = cellsize*sec2rad
      else
         incell = abs(delta_dec)*deg2rad
      end if

      if (abs(delta_ra).ne.(abs(delta_dec))) then
         write(*,*) 'NB coordinate increments not the same in RA and Dec'
      end if

! Find RA and Dec of central pixel 
! Convert RA and Dec into radians
!      obsra = (obsra-(cent_ra-npix_ra/2-1)*delta_ra)*deg2rad
!      obsdec = (obsdec-(cent_dec-npix_ra/2-1)*delta_dec)*deg2rad

! Do a proper sin projection for cases where the reference pixel is far away from the central pixel
      dx = npix_ra/2+1-cent_ra
      dy = npix_dec/2+1-cent_dec
      obsra1 = obsra*deg2rad
      obsdec1 = obsdec*deg2rad
      call sin_proj_back(obsra1, obsdec1, obsra1, obsdec1, incell/sec2rad, dx, dy)
      obsra = obsra1
      obsdec = obsdec1

! Attempt to find beam parameters

! Determine the number of keywords in the header
      call ftghsp(iunit,nkeys2,nspace,status)     
      
! Search for 80-character keyword record containing BMAJ
      i = nkeys2      
      beam_found = .false.
      do while((.not.beam_found).and.(i.gt.0))         
        call ftgrec(iunit,i,record,status)

! Parse RECORD for list of KEYS separated by spaces, return number NKEYS
	call get_keys(record,nkeymax,keys,nkeys)
	read(keys(4),*,err=9) match
 9      continue
	if (match.eq.'BMAJ=') then 
           beam_found = .true.
        end if
	i = i-1
      enddo

      if (beam_found) then

!  Read values of bmaj and bmin     
         read(keys(4),*,err=10) bmaj_string
         read(keys(6),*,err=10) bmin_string
         read(keys(8),*,err=10) bpa_string
         if((bmaj_string.eq.'BMAJ=').and.(bmin_string.eq.'BMIN=').and.&
            (bpa_string.eq.'BPA=')) then
            read(keys(5),*,err=10) bmaj
            read(keys(7),*,err=10) bmin
            read(keys(9),*,err=10) bpa	  
         endif
      else
         bmaj = 0.0
         bmin = 0.0
         bpa = 0.0
      end if

! Close the file
      status = 0
      call ftclos(iunit,status)

      return 
      
 10   write(*,*) 'ERROR in reading beam parameters'
      return

   end subroutine read_fits_header

!***************************************************************************

   subroutine printerror(status)

      implicit none

      integer                        :: status
      character*30                   :: errtext
      character*80                   :: errmessage

! Return if status OK
      if (status.le.0) return

! Otherwise, print error messages.

      call ftgerr(status,errtext)
      write(*,*) '*** FITSIO Error Status = ', status, ': ', errtext

      call ftgmsg(errmessage)
      do while (errmessage.ne.' ')
         write(*,*) errmessage
         call ftgmsg(errmessage)
      end do

   end subroutine printerror

!***************************************************************************

   subroutine get_keys(record,nkeymax,keys,nkeys)
!
!  Parse RECORD for list of KEYS separated by spaces, return number NKEYS
!   - end line on # or ! for comments
!
      implicit none

      integer nkeymax,nkeys,i,j,lr,ibeg,iend
      character*(*) record
      character keys(nkeymax)*24      

      i = 0
      nkeys = 0
      lr = len(record)
      do while (i.lt.lr)
         j = 0
         i = i+1
         do while (record(i:i).eq.' ') 
            if(i.eq.lr) return
            i = i+1
         enddo
         if(record(i:i).eq.'#' .or. record(i:i).eq.'!') return

         do while (i+j.le.lr .and. record(i+j:i+j).ne.' ')
            j = j+1
         end do
         nkeys = nkeys+1
         ibeg = i
         iend = i+j-1
         i = iend

         keys(nkeys) = record(ibeg:iend)
         ! Ignore single quotes, can cause problems in read_fits_header
         if (keys(nkeys).eq.'''') nkeys = nkeys - 1
	 
      enddo    

   end subroutine get_keys
  
!**************************************************************************

   subroutine write_beam_table(unit)

      use kind_def
      use maths

      implicit none

      integer :: unit, status
      integer :: nspace, colnum, frow, felem, nrows, tabcount, fcount
      integer :: rowlen
      integer, parameter :: tfields=11
      integer :: tbcol(tfields)
      character*16 :: ttype(tfields),tform(tfields),tunit(tfields), extname
      real ::  tmp
      integer, allocatable, dimension (:) :: pnt,ra_hrs,ra_min
      integer, allocatable, dimension (:) :: dec_deg,dec_min
      real, allocatable, dimension (:) :: ra_sec,dec_sec
      real, allocatable, dimension (:) :: bmaj,bmin,bpa,noise 
      real(kind=dp) :: tmp_ra, tmp_dec
      character :: chra*16, chdec*16
      integer :: prra, prdec, lr, ld
      logical :: debug

      data ttype/'POINTING','RAH','RAM','RAS','DECD','DECM','DECS',&
     'BMAJ','BMIN','BPA','NOISE'/
      data tform/'I10','I5','I5','F10.3','I5','I5','F10.3','F10.2',&
     'F10.2','F10.2','F12.7'/
      data tunit/' ','hrs','min','sec','deg','min','sec','arcsec',&
     'arcsec','deg','Jy'/

      status = 0
      debug = .false.
      
! Do allocation
      nrows=count(mask=fits_use_file)
      allocate(pnt(nrows))
      allocate(ra_hrs(nrows))
      allocate(ra_min(nrows))
      allocate(ra_sec(nrows))
      allocate(dec_deg(nrows))
      allocate(dec_min(nrows))
      allocate(dec_sec(nrows))
      allocate(bmaj(nrows))
      allocate(bmin(nrows))
      allocate(bpa(nrows))
      allocate(noise(nrows))

! Prepare tables
      tabcount = 0
      prra = 3
      prdec = 3
      do fcount = 1, n_fits_files
         if (fits_use_file(fcount)) then
            tabcount = tabcount+1
            pnt(tabcount) = fcount
            tmp = fits_ra(fcount)/hr2rad
            ra_hrs(tabcount) = int(tmp)
            tmp = 60.*(tmp-ra_hrs(tabcount))
            ra_min(tabcount) = int(tmp)
            tmp = 60.*(tmp-ra_min(tabcount))
            ra_sec(tabcount) = tmp
            tmp = fits_dec(fcount)/deg2rad
            dec_deg(tabcount) = int(tmp)
            tmp = 60.*(tmp-dec_deg(tabcount))
            dec_min(tabcount) = int(tmp)
            tmp = 60.*(tmp-dec_min(tabcount))
            dec_sec(tabcount) = tmp
            bmaj(tabcount) = fits_bmaj(fcount)*3600.d0
            bmin(tabcount) = fits_bmin(fcount)*3600.d0
            bpa(tabcount) = fits_bpa(fcount)
            noise(tabcount) = fits_noise(fcount)
            tmp_ra = fits_ra(fcount)/hr2rad
            tmp_dec = fits_dec(fcount)/deg2rad
            call chr_chdtos(tmp_ra,prra,chra,lr)
            call chr_chdtos(tmp_dec,prdec,chdec,ld)
            if (debug) then
               write(*,*) tabcount,&
                       ra_hrs(tabcount),ra_min(tabcount),ra_sec(tabcount),&
                       dec_deg(tabcount),dec_min(tabcount),dec_sec(tabcount),&
                       bmaj(tabcount),bmin(tabcount),bpa(tabcount),&
                       noise(tabcount)
               write(*,*) chra(1:lr),chdec(1:ld)
            end if
         end if
      end do
      
! Error checking
      if (nrows.ne.tabcount) write(*,*) 'ERROR ', nrows,tabcount

! FTCRHD creates a new empty FITS extension following the current
! extension and moves to it.  In this case, FITSIO was initially
! positioned on the primary array when the FITS file was first opened, so
! FTCRHD appends an empty extension and moves to it.  All future FITSIO
! calls then operate on the new extension (which will be an ASCII
! table).
      call ftcrhd(unit,status)
      
! define parameters for the ASCII table (see the above data statements)
      extname='BEAM_ASCII'
      
! FTGABC is a convenient subroutine for calculating the total width of
! the table and the starting position of each column in an ASCII table.
! Any number of blank spaces (including zero)  may be inserted between
! each column of the table, as specified by the NSPACE parameter.
      nspace=1
      call ftgabc(tfields,tform,nspace,rowlen,tbcol,status)
      
! FTPHTB writes all the required header keywords which define the
! structure of the ASCII table. NROWS and TFIELDS give the number of
! rows and columns in the table, and the TTYPE, TBCOL, TFORM, and TUNIT
! arrays give the column name, starting position, format, and units,
! respectively of each column. The values of the ROWLEN and TBCOL parameters
! were previously calculated by the FTGABC routine.
      call ftphtb(unit,rowlen,nrows,tfields,ttype,tbcol,tform,tunit,&
                 extname,status)
      
! Write names to the first column, diameters to 2nd col., and density to 3rd
! FTPCLS writes the string values to the NAME column (column 1) of the
! table.  The FTPCLJ and FTPCLE routines write the diameter (integer) and
! density (real) value to the 2nd and 3rd columns.  The FITSIO routines
! are column oriented, so it is usually easier to read or write data in a
! table in a column by column order rather than row by row.
      frow=1
      felem=1
      colnum=1
      call ftpclj(unit,colnum,frow,felem,nrows,pnt,status)
      colnum=2
      call ftpclj(unit,colnum,frow,felem,nrows,ra_hrs,status)  
      colnum=3
      call ftpclj(unit,colnum,frow,felem,nrows,ra_min,status) 
      colnum=4
      call ftpcle(unit,colnum,frow,felem,nrows,ra_sec,status)      
      colnum=5
      call ftpclj(unit,colnum,frow,felem,nrows,dec_deg,status)	 
      colnum=6
      call ftpclj(unit,colnum,frow,felem,nrows,dec_min,status)	 
      colnum=7
      call ftpcle(unit,colnum,frow,felem,nrows,dec_sec,status)	 
      colnum=8
      call ftpcle(unit,colnum,frow,felem,nrows,bmaj,status)     
      colnum=9
      call ftpcle(unit,colnum,frow,felem,nrows,bmin,status)      
      colnum=10	   
      call ftpcle(unit,colnum,frow,felem,nrows,bpa,status)
      colnum=11	   
      call ftpcle(unit,colnum,frow,felem,nrows,noise,status)  
      
! Do deallocation
      deallocate(pnt)
      deallocate(ra_hrs)
      deallocate(ra_min)
      deallocate(ra_sec)
      deallocate(dec_deg)
      deallocate(dec_min)
      deallocate(dec_sec)
      deallocate(bmaj)
      deallocate(bmin)
      deallocate(bpa)
      deallocate(noise)

   end subroutine write_beam_table

!**********************************************************************

   subroutine write_multi_fits(n_src,filename,sing_phase_cent)

      use kind_def
      use sz_globals

      implicit none

      character*80, intent(in) :: filename
      integer, intent(in) :: n_src
      logical, intent(in) :: sing_phase_cent
      logical :: ps_cent
      integer  :: isrc,i,b1,b2,iunit
      !!integer   :: iy,id,im
      real(dp)   :: utsec,jd,sum_weight
      integer    :: base_id, date1, date2, src
      integer :: idate(3)
      real*8  :: epoch
      integer    :: basel, nbytes, ngroup, n_basel
      integer    :: ig,a,b,p,iq,ant1,ant2,samp
      integer, parameter :: ng = 7+MAX_FREQ*3
      real*4     :: group(ng)
      real*4     :: raref, decref
      real(dp)   :: obs_ut_start
      real(dp)   :: obs_mjd_start
      integer    :: smooth, n_smooth
      real(dp)   :: sm_inc
      character*1  :: chr1
      integer l, ls
      real*8  :: freq, dfreq, fd
      
      logical, parameter ::  simple = .true.
      logical, parameter ::  extend = .true.  
      integer, parameter :: blocksize = 1

      integer   chr_lenb
      external  chr_lenb

      logical   chr_match
      external  chr_match

      character*16 :: kword
      character*16 :: value
      character*72 :: comment
      character*32 :: chstr

      real*8 nullval
      logical anyf

      obs_ut_start = 0.
      obs_mjd_start = 0.
      jdscal1 = 1.d0/128.d0
      jdscal2 = 1.d0/6.d6
      jdzero = int(obs_mjd_start+2400000.5d0)


      basel = 1
      n_basel = n_antennas*(n_antennas-1)/2

      ngroup = n_samp*n_basel*n_src

      write(*,*) 'Writing ',ngroup,' visibilities to FITS'

      sum_weight = 0.
      status = 0
      bscale = 1.
      isrc = 1

      n_basel = n_antennas*(n_antennas-1)/2
      ngroup = n_samp*n_basel*n_src

      bitpix = -32
      naxis = 6 
      pcount = 7
      gcount = ngroup
      naxisn(1) =  0
      naxisn(2) =  3
      naxisn(3) =  1
      naxisn(4) =  nchan
      naxisn(5) =  1
      naxisn(6) =  1



!     Open and initialise new fits file
      call ftgiou(iunit, status)
      !write(*,*) 'ftgiou', status
      call ftinit(iunit,filename, blocksize, status)
      !write(*,*) 'ftinit', status
!      call ftnopn(iunit,filename,1,stat)
      raref = obs_raref
      decref = obs_decref


      call ftpkyl(iunit,'SIMPLE', .true.,'Standard FITS format',status)      
      call ftpkyj(iunit,'BITPIX', bitpix,'Bits per pixel',status)     
      call ftpkyj(iunit,'NAXIS',  naxis,'Number of axes',status)    
      call ftpkyj(iunit,'NAXIS1', naxisn(1),'No image: uv data', status)      
      call ftpkyj(iunit,'NAXIS2', naxisn(2),'Complex: cos, sin, weight',status)
      call ftpkyj(iunit,'NAXIS3',  naxisn(3),'Number of polarisations', status)
      call ftpkyj(iunit,'NAXIS4',  naxisn(4),'Number of frequencies', status) 
      call ftpkyj(iunit,'NAXIS5',  naxisn(5),'RA', status)
      call ftpkyj(iunit,'NAXIS6',  naxisn(6),'Dec', status)

      call ftpkyl(iunit,'BLOCKED', .true., 'Tape may be blocked',status)
      call ftpkyl(iunit,'GROUPS ', .true., 'Groups data structure',status)
      call ftpkyj(iunit,'PCOUNT ', pcount, 'Parameters per group',status)
      call ftpkyj(iunit,'GCOUNT ', gcount, 'Number of groups',status)
      call ftpkyl(iunit,'EXTEND ', .true., 'Extension is antenna table',status)



!  Object, telescope, date of observation
      ! Sort date/timestamps and epoch
      call util_enqdat(idate)
      idate(3)=idate(3)-2000
      call sla_CALDJ(idate(3), idate(2), idate(1),obs_mjd_start, i)
      call sla_djcl(obs_mjd_start,idate(3), idate(2), idate(1),fd,i)
      if (obs_epoch.eq.1) epoch=1950.0d0
      if (obs_epoch.eq.2) epoch=2000.0d0
      !Get telescope name
      if (chr_match(tel_name,'LA')) then
         write (value, '(A)') 'AMI-LA'
      elseif (chr_match(tel_name,'LA')) then
         value='AMI-SA'
      else
         value='AMI-SA'
      endif

      
      call ftpkys(iunit,'OBJECT', obs_name(1:16),'Source/field name',status)     
      call ftpkys(iunit,'TELESCOP', value,' Profile',status)      
      call ftpkys(iunit,'INSTRUME', value,' Sim',status)      
      write(value, '(I4 ,("-",I2.2),("-",I2.2))')  idate(3),idate(2),idate(1)  
      call ftpkys(iunit,'DATE-OBS', value,'Start date of observation',status)
      call ftpkys(iunit,'BUNIT', ' ','Units of data',status)     
      call ftpkyd(iunit,'BZERO', 0.0d0, -10,'Data offset',status)
      call ftpkyd(iunit,'BSCALE', 1.0d0, -10,'Data = FITS*BSCALE + BZERO',status)
      call ftpkyd(iunit,'EPOCH', epoch, -10,'Equinox of RA, Dec',status)
      call ftpkye(iunit,'OBSRA', raref, -10,'Antenna pointing RA',status)
      call ftpkye(iunit,'OBSDEC', decref, -10,'Antenna pointing DEC',status)



      !  Frequency channels

      freq = nu(1)
      dfreq = 0.0d0
      if (nchan.gt.1) dfreq = nu(2)-nu(1)

!  Details of group data
      

      call ftpkys(iunit,'CTYPE2', 'COMPLEX ','',status)
      call ftpkyd(iunit,'CRVAL2', 1.0d0, -10,'',status)            
      call ftpkyd(iunit,'CDELT2', 1.0d0, -10,'',status)      
      call ftpkyd(iunit,'CRPIX2', 1.0d0, -10,'',status)     
      call ftpkyd(iunit,'CROTA2', 0.0d0, -10,'',status)      

      call ftpkys(iunit,'CTYPE3', 'STOKES  ','Stokes parameter : I+Q' ,status)
      call ftpkyd(iunit,'CRVAL3', -6.0d0, -10,'',status)
      call ftpkyd(iunit,'CDELT3', 1.0d0, -10,'',status)
      call ftpkyd(iunit,'CRPIX3', 1.0d0, -10,'',status)
      call ftpkyd(iunit,'CROTA3', 0.0d0, -10,'',status)      

      call ftpkys(iunit,'CTYPE4', 'FREQ','Frequency, Hertz',status)
      call ftpkyd(iunit,'CRVAL4', freq, -10,'',status)
      call ftpkyd(iunit,'CDELT4', dfreq, -10,'',status)
      call ftpkyd(iunit,'CRPIX4', 1.0d0, -10,'',status)
      call ftpkyd(iunit,'CROTA4', 0.0d0, -10,'',status)      

      call ftpkys(iunit,'CTYPE5', 'RA','Right Ascension, degrees',status)
      call ftpkye(iunit,'CRVAL5', raref, -10,'',status)      
      call ftpkyd(iunit,'CDELT5', 1.0d0, -10,'',status)      
      call ftpkyd(iunit,'CRPIX5', 1.0d0, -10,'',status)
      call ftpkyd(iunit,'CROTA5', 0.0d0, -10,'',status)       

      call ftpkys(iunit,'CTYPE6', 'DEC','Declination, degrees',status)
      call ftpkye(iunit,'CRVAL6', decref, -10,'',status)
      call ftpkyd(iunit,'CDELT6', 1.0d0, -10,'',status)
      call ftpkyd(iunit,'CRPIX6', 1.0d0, -10,'',status)
      call ftpkyd(iunit,'CROTA6', 0.0d0, -10,'',status)

!  Details of group parameters
      uscale = 1.
      uzero  = 0.0d0
      vscale = uscale
      vzero  = 0.0d0
      wscale = uscale
      wzero  = 0.0d0


      call ftpkys(iunit,'PTYPE1', 'UU','U coordinate, seconds',status)
      call ftpkyd(iunit,'PSCAL1', uscale, -10,'Real = FITS*PSCAL + PZERO',status)
      call ftpkyd(iunit,'PZERO1', uzero,  -10, '',status)

      call ftpkys(iunit,'PTYPE2', 'VV','V coordinate, seconds',status)
      call ftpkyd(iunit,'PSCAL2', vscale, -10,'Real = FITS*PSCAL + PZERO',status)
      call ftpkyd(iunit,'PZERO2', vzero,  -10,'',status)      

      call ftpkys(iunit,'PTYPE3', 'WW','W coordinate, seconds',status)
      call ftpkyd(iunit,'PSCAL3', wscale, -10,'Real = FITS*PSCAL + PZERO',status)
      call ftpkyd(iunit,'PZERO3', wzero,  -10,'',status)     

      call ftpkys(iunit,'PTYPE4', 'DATE','Julian Date',status)
      call ftpkyd(iunit,'PSCAL4', jdscal1 , -10,'',status)
      call ftpkyd(iunit,'PZERO4', jdzero, -10,'',status)           
      
      call ftpkys(iunit,'PTYPE5', 'DATE','Julian Date',status)
      call ftpkyd(iunit,'PSCAL5', jdscal2, -10,'',status)
      call ftpkyd(iunit,'PZERO5', 0.0d0, -10, '',status)      
      
      call ftpkys(iunit,'PTYPE6', 'BASELINE','Ant1*256 + Ant2',status)
      call ftpkyd(iunit,'PSCAL6', 1.0d0, -10,'',status)
      call ftpkyd(iunit,'PZERO6', 0.0d0, -10,'',status)

      call ftpkys(iunit,'PTYPE7', 'SOURCE  ','Source ID. NO. (in SU table)',status)
      call ftpkyd(iunit,'PSCAL7', 1.0d0, -10,'',status)
      call ftpkyd(iunit,'PZERO7', 0.0d0, -10,'',status)

! Write comments and history
      comment='Profile Simulation'
      call ftpcom(iunit, comment, status)
      chstr = ' observation'
      ls = chr_lenb(chstr)
      l = chr_lenb(obs_name)
      write(comment,'(A X A)') obs_name(1:l),chstr(1:ls)
      call ftpcom(iunit, comment, status)

      do i = 1,lhist
         call ftphis(iunit,history(i), status)
      enddo

      history(1) = "HISTORY AIPS  SORT ORDER = ''TB''"  
      history(2) = '                   / Where T is time (IAT), B is baseline num'
      call ftphis(iunit,history(1), status)
      call ftphis(iunit,history(2), status)
      write(value,'(E16.9)'),1.0 
      history(3) = "HISTORY AIPS WTSCAL = " // value // "//Complex wts = WTSCAL*(FITS*BSCALE+BZERO)"
      call ftphis(iunit,history(3), status)
!Header done

     
! loop over number of samples in each track:
      b=0
      do samp = 1,n_samp*n_src
         utsec = (real(samp)*ha_step)-obs_ut_start
         jd = obs_mjd_start + utsec*const_st2day + 2400000.5d0
         date1=nint((jd-jdzero)/jdscal1)
         date2=nint((jd-jdzero-date1*jdscal1)/jdscal2)

! loop over baselines:
! write visibilities for each baseline as a single group of data:
         do ant1 = 1,n_antennas-1
            do ant2 = ant1+1,n_antennas               
               b = b+1
               base_id = ant1*256 + ant2
               i = ((samp-1)*n_basel)+basel
! Initialise group
               !do ig = 1,ng
               group(:) = blank
               !end do

! uvw in seconds of time:
               group(1) = real(u(b)/(nu(1))/uscale, 4)               
               group(2) = real(v(b)/(nu(1))/vscale, 4)
               group(3) = w(b)/(nu(1))/wscale
               group(4) = date1
               group(5) = date2
               group(6) = base_id
               group(7) = which_pointing(b)
               ig = 7! loop over frequency:               
               do chan =1,nchan
                  group(ig+1) = data_re(b)/bscale
                  group(ig+2) = data_im(b)/bscale
                  group(ig+3) = weight(b)/bscale
                  sum_weight = sum_weight+weight(b)
                  b=b+1
                  ig=ig+3
               end do
               
               b=b-1
! Convert to standard representation
               !call util_r4cvt(group,ng)
               !write(*,*) group

!
               call ftpgpe(iunit, i, 1, ng, group, status)
               if (status.ne.0)then
                  call FTGERR(status, comment)
                  write(*,*) 'ftpgpe error',status, comment
                  stop
               endif
               basel = basel + 1
            end do            
         end do
         basel = 1 ! Reset baseline count
      end do
     
      call fits_wrant(iunit, n_antennas,x_pos,y_pos,z_pos,obsfreq,status)

      call fits_wrsu(iunit, obs_name,obs_epoch,raref,decref,n_src, &
           ra_offset,dec_offset,ps_cent,status)


      call ftclos(iunit, status) 
      call ftfiou(iunit, status)

      write(*,*) b,' total visibilities; ',sum_weight,&
           &' sum of weights'

    end subroutine write_multi_fits

!****************************************************************************
   
   subroutine write_uv_fits(filename)

      use kind_def
      use sz_globals

      implicit none

      character*80, intent(in)   :: filename
      real(dp)   :: utsec,jd,sum_weight
      integer    :: base_id, date1, date2
      integer    :: basel, nbytes, ngroup, n_basel
      integer    :: ig,a,b,p,iq,ant1,ant2,samp, i, l ,ls
      integer, parameter :: ng = 6+MAX_FREQ*3
      real*4     :: group(ng)
      real*4     :: raref, decref
      real(dp)   :: obs_ut_start
      real(dp)   :: obs_mjd_start
      character*1  :: chr1
      integer :: idate(3)
      real*8  :: epoch
      real*8  :: freq, dfreq

      integer :: iunit
      integer :: naxes(ng)
         
      integer   chr_lenb
      external  chr_lenb

      logical   chr_match
      external  chr_match


!c Header items

      !integer, dimension(naxis) :: naxisn
      
      logical, parameter ::  simple = .true.
      logical, parameter ::  extend = .true.  
      integer, parameter :: blocksize = 1

      character*72 :: comment
      character*16 :: kword
      character*16 :: value
      real*8 :: fvalue
      character*32 :: chstr

      real*8 nullval
      logical anyf

      basel = 1
      n_basel = n_antennas*(n_antennas-1)/2
      ngroup = n_samp*n_basel

      bitpix = -32
      naxis = 6 
      pcount = 6
      gcount = ngroup
      naxisn(1) =  0
      naxisn(2) =  3
      naxisn(3) =  1
      naxisn(4) =  nchan
      naxisn(5) =  1
      naxisn(6) =  1


      obs_ut_start = 0.
      obs_mjd_start = 0.
      jdscal1 = 1.d0/128.d0
      jdscal2 = 1.d0/6.d6
      jdzero = int(obs_mjd_start+2400000.5d0)
      uscale = 1.0d0
      vscale = uscale
      wscale = uscale      

      sum_weight = 0.
      status = 0
      bscale = 1.

     
!     Open and initialise new fits file
      call ftgiou(iunit, status)
      call ftinit(iunit,filename, blocksize, status)

!     Convert RA to degrees
      raref = obs_raref*15d0 
      decref = obs_decref

      call ftpkyl(iunit,'SIMPLE', .true.,'Standard FITS format',status)      
      call ftpkyj(iunit,'BITPIX', bitpix,'Bits per pixel',status)     
      call ftpkyj(iunit,'NAXIS',  naxis,'Number of axes',status)    
      call ftpkyj(iunit,'NAXIS1', naxisn(1),'No image: uv data', status)      
      call ftpkyj(iunit,'NAXIS2', naxisn(2),'Complex: cos, sin, weight',status)
      call ftpkyj(iunit,'NAXIS3',  naxisn(3),'Number of polarisations', status)
      call ftpkyj(iunit,'NAXIS4',  naxisn(4),'Number of frequencies', status) 
      call ftpkyj(iunit,'NAXIS5',  naxisn(5),'RA', status)
      call ftpkyj(iunit,'NAXIS6',  naxisn(6),'Dec', status)
      call ftpkyl(iunit,'BLOCKED', .true., 'Tape may be blocked',status)
      call ftpkyl(iunit,'GROUPS ', .true., 'Groups data structure',status)
      call ftpkyj(iunit,'PCOUNT ', pcount, 'Parameters per group',status)
      call ftpkyj(iunit,'GCOUNT ', gcount, 'Number of groups',status)
      call ftpkyl(iunit,'EXTEND ', .true., 'Extension is antenna table',status)

      !     

      !call ftphpr(iunit,simple,bitpix,naxisn,naxis,pcount,gcount,extend,status)
      !write(*,*) 'ftphpr', status, naxis, naxisn

      


!  Object, telescope, date of observation

      call util_enqdat(idate)
!     call sla_cldj(idate(3),idate(2),idate(1),mjd1,j) 
!     obs_mjd_start = mjd1
!     source = obs_name
      if (obs_epoch.eq.1) epoch=1950.0d0
      if (obs_epoch.eq.2) epoch=2000.0d0

      if (chr_match(tel_name,'LA')) then
         write (value, '(A)') 'AMI-LA'
      elseif (chr_match(tel_name,'LA')) then
         value='AMI-SA'
      else
         value='AMI-SA'
      endif

      
      call ftpkys(iunit,'OBJECT', obs_name(1:16),'Source/field name',status)
      
      call ftpkys(iunit,'TELESCOP', value,' Profile',status)
      
      call ftpkys(iunit,'INSTRUME', value,' Sim',status)
      
      write(value, '(I4 ,("-",I2.2),("-",I2.2))')  idate(3),idate(2),idate(1)  
      call ftpkys(iunit,'DATE-OBS', value,'Start date of observation',status)
      call ftpkys(iunit,'BUNIT', ' ','Units of data',status)
      
!      write(value,'(E16.9)'), 0.0d0
      call ftpkyd(iunit,'BZERO', 0.0d0, -10,'Data offset',status)
      
 !     write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'BSCALE', 1.0d0, -10,'Data = FITS*BSCALE + BZERO',status)
      
!      write(value,'(F10.4)'), epoch
      call ftpkyd(iunit,'EPOCH', epoch, -10,'Equinox of RA, Dec',status)
 !     write(value,'(E16.9)'), raref
      call ftpkye(iunit,'OBSRA', raref, -10,'Antenna pointing RA',status)
!      write(value,'(F10.4)'), decref
      call ftpkye(iunit,'OBSDEC', decref, -10,'Antenna pointing DEC',status)



!  Phase centre [degrees]
      raref = obs_rarot*15d0
      decref = obs_decrot


!  Frequency channels

      freq = nu(1)
      dfreq = 0.0d0
      if (nchan.gt.1) dfreq = nu(2)-nu(1)

!  Details of group data
      

      call ftpkys(iunit,'CTYPE2', 'COMPLEX ','',status)
      write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'CRVAL2', 1.0d0, -10,'',status)
      
      
      write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'CDELT2', 1.0d0, -10,'',status)
      
      write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'CRPIX2', 1.0d0, -10,'',status)
      
      write(value,'(E16.9)'), 0.0d0
      call ftpkyd(iunit,'CROTA2', 0.0d0, -10,'',status)
      


      call ftpkys(iunit,'CTYPE3', 'STOKES  ','Stokes parameter :',status)

      write(value,'(E16.9)'), -6.0d0
      call ftpkyd(iunit,'CRVAL3', -6.0d0, -10,'',status)
      
      write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'CDELT3', 1.0d0, -10,'',status)
      
      write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'CRPIX3', 1.0d0, -10,'',status)
      
      write(value,'(E16.9)'), 0.0d0
      call ftpkyd(iunit,'CROTA3', 0.0d0, -10,'',status)
      


      call ftpkys(iunit,'CTYPE4', 'FREQ','Frequency, Hertz',status)
      write(value,'(E16.9)'), freq
      call ftpkyd(iunit,'CRVAL4', freq, -10,'',status)
      
      write(value,'(E16.9)'), dfreq
      call ftpkyd(iunit,'CDELT4', dfreq, -10,'',status)
      
      write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'CRPIX4', 1.0d0, -10,'',status)
      
      write(value,'(E16.9)'), 0.0d0
      call ftpkyd(iunit,'CROTA4', 0.0d0, -10,'',status)
      


      call ftpkys(iunit,'CTYPE5', 'RA','Right Ascension, degrees',status)
      write(value,'(E16.9)'), raref
      call ftpkye(iunit,'CRVAL5', raref, -10,'',status)
      
      write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'CDELT5', 1.0d0, -10,'',status)
      
      write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'CRPIX5', 1.0d0, -10,'',status)
      
      write(value,'(E16.9)'), 0.0d0
      call ftpkyd(iunit,'CROTA5', 0.0d0, -10,'',status)
       


      call ftpkys(iunit,'CTYPE6', 'DEC','Declination, degrees',status)
      write(value,'(E16.9)'), decref
      call ftpkye(iunit,'CRVAL6', decref, -10,'',status)
      
      write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'CDELT6', 1.0d0, -10,'',status)
      
      write(value,'(E16.9)'), 1.0d0
      call ftpkyd(iunit,'CRPIX6', 1.0d0, -10,'',status)
      
      write(value,'(E16.9)'), 0.0d0
      call ftpkyd(iunit,'CROTA6', 0.0d0, -10,'',status)
       


!  Details of group parameters

!     uscale = 1.0d-10
      uscale = 1.0d0
      uzero  = 0.0d0

      call ftpkys(iunit,'PTYPE1', 'UU','U coordinate, seconds',status)
      call ftpkyd(iunit,'PSCAL1',uscale, -10,'Real = FITS*PSCAL + PZERO',status)
      call ftpkyd(iunit,'PZERO1',  uzero, -10, '',status)


      vscale = uscale
      vzero = 0.0d0

      call ftpkys(iunit,'PTYPE2', 'VV','V coordinate, seconds',status)
      call ftpkyd(iunit,'PSCAL2',vscale, -10,'Real = FITS*PSCAL + PZERO',status)
      call ftpkyd(iunit,'PZERO2', vzero, -10,'',status)
      

      wscale = uscale
      wzero = 0.0d0

      call ftpkys(iunit,'PTYPE3', 'WW','W coordinate, seconds',status)
      call ftpkyd(iunit,'PSCAL3',wscale, -10,'Real = FITS*PSCAL + PZERO',status)
      call ftpkyd(iunit,'PZERO3', wzero, -10,'',status)
      

      jdzero = int(obs_mjd_start+2400000.5d0)
      jdscal1 = 1.d0/128.d0
      jdscal2 = 1.d0/6.d6

      call ftpkys(iunit,'PTYPE4', 'DATE','Julian Date',status)
      call ftpkyd(iunit,'PSCAL4', jdscal1 , -10,'',status)
      call ftpkyd(iunit,'PZERO4', jdzero, -10,'',status)
            
      
      call ftpkys(iunit,'PTYPE5', 'DATE','Julian Date',status)
      call ftpkyd(iunit,'PSCAL5', jdscal2, -10,'',status)
      call ftpkyd(iunit,'PZERO5', 0.0d0, -10, '',status)
      

      
      call ftpkys(iunit,'PTYPE6', 'BASELINE','Ant1*256 + Ant2',status)
      call ftpkyd(iunit,'PSCAL6', 1.0d0, -10,'',status)
      call ftpkyd(iunit,'PZERO6', 0.0d0, -10,'',status)


! Write comments and history
      comment='Profile Simulation'
      call ftpcom(iunit, comment, status)
      chstr = ' observation'
      ls = chr_lenb(chstr)
      l = chr_lenb(obs_name)
      write(comment,'(A X A)') obs_name(1:l),chstr(1:ls)
      call ftpcom(iunit, comment, status)
      

      !!write(history(1), '(A)') 
      history(1) = "HISTORY AIPS  SORT ORDER = ''TB''"  
      history(2) = '                   / Where T is time (IAT), B is baseline num'
      call ftphis(iunit,history(1), status)
      call ftphis(iunit,history(2), status)
      write(value,'(E16.9)'),1.0 
      history(1) = "HISTORY AIPS WTSCAL = " // value // "//Complex wts = WTSCAL*(FITS*BSCALE+BZERO)"
      call ftphis(iunit,history(1), status)
      !call ftpkys(iunit,'END', '','',status)
      


! loop over number of samples in each track:
      b = 0
      do samp = 1,n_samp
         utsec = (real(samp)*ha_step)-obs_ut_start
         jd = obs_mjd_start + utsec*const_st2day + 2400000.5d0
         date1=nint((jd-jdzero)/jdscal1)
         date2=nint((jd-jdzero-date1*jdscal1)/jdscal2)


!! Check whether AMI type telescope with + and - spacings
!         if ((tel_name.eq.'LA').or.(tel_name.eq.'SA')) n_pol = 1


! loop over baselines:
! write visibilities for each baseline as a single group of data:
        ! write(*,*) n_samp, n_pol, n_antennas, nchan
         do p =1,n_pol
            do ant1 = 1,n_antennas-1
               do ant2 = ant1+1,n_antennas
                  b = b+1
                  base_id = ant1*256 + ant2
                  !do ig = 1,ng
                     group(:) = blank
                  !end do
                  i = ((samp-1)*n_basel)+basel
! uvw in seconds of time:
                  group(1) = u(b)/(nu(1))/uscale
                  group(2) = v(b)/(nu(1))/vscale
                  group(3) = w(b)/(nu(1))/wscale
                  group(4) = date1
                  group(5) = date2
                  group(6) = base_id
                  ig=6

! loop over frequency:
                  do chan =1,nchan
                     group(ig+1) = data_re(b)/bscale
                     group(ig+2) = data_im(b)/bscale
                     group(ig+3) = weight(b)/bscale
                     b=b+1
                     ig=ig+3
                     !write(*,*) b,u(b),v(b),data_re(b),data_im(b),weight(b)
                  end do
                  b=b-1
!    Convert to standard representation                 
                  call ftpgpe(iunit, i, 1, ng, group, status)
                  if (status.ne.0)then
                     call FTGERR(status, comment)
                     write(*,*) 'ftpgpe error',status, comment
                     stop
                  endif                  
                  basel = basel + 1
               end do
            end do
         end do
         basel = 1 ! Reset baseline count
      end do


      call fits_wrant(iunit, n_antennas,x_pos,y_pos,z_pos,obsfreq,status)

      sum_weight = sum(weight)

      call ftclos(iunit, status) 
      call ftfiou(iunit, status)
      

    end subroutine write_uv_fits




end module fits_handling
