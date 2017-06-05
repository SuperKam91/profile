  module file_handling

  use kind_def
   use sz_globals

   interface read_model
      module procedure read_model
   end interface

   interface read_vis_data
      module procedure read_vis_data
   end interface

   interface write_model
      module procedure write_model
   end interface

   public write_weights

   interface write_vis_data
      module procedure write_vis_data_nof
      module procedure write_vis_data_givenf
   end interface

   interface write_flagpixel
      module procedure write_flagpixel
   end interface

   interface read_flagpixel
      module procedure read_flagpixel
   end interface

   interface read_ptsrc_list
      module procedure read_ptsrc_list
   end interface

   public read_nt
   
   interface read_telescope
      module procedure read_telescope
      module procedure read_specified_telescope
   end interface

   public write_telescope

   interface read_gridded_pb
      module procedure read_gridded_pb
   end interface

   public init_directory
   public get_directory
   public get_read_filename
   public get_write_filename
   public write_header_digest
   public read_header_digest
   public read_ascii_beam
   public read_ascii_map

contains

!***************************************************************************

!  Reads in a particular model
   subroutine read_model

      use kind_def
      use sz_globals
      use cosmology

      implicit none

      real(kind=dp) :: d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d6a
      real(kind=dp) :: d11,d12,d13,d14,d15,d16,d17,d18,d19
      real(kind=dp) :: d20, d21, d22, d23, d24, d25, d26
      real(kind=dp), dimension(3) :: da1, da2
      logical :: l1, l2, ex
      character(len=fname_len) :: filename, file_type, dch1
      integer :: iunit

      prompt = 'Model filename'
      which_dir = model_dir
      rootname = 'model'
      extn = '.mod'
      call get_read_filename(filename)

! Open the file
      call io_nxtlun(iunit,status)
      open(iunit,file = filename,form = 'formatted')

! Read the saved data
      read(iunit,'(A)') file_type
      if (file_type(1:5).ne.'model') then
         write(*,*) 'Not a model file'
         return
      else
         read(iunit,*) dch1
         read(iunit,*) d1,d2,d3,d4,d5,d6,d6a
         read(iunit,*) d7,d8,d9,d10,d11,d12,d13,d14
         read(iunit,*) d15,d16,d17,d18
         read(iunit,*) da1(1), da1(2), da1(3), l2
         read(iunit,*) da2(1), da2(2), da2(3)
         read(iunit,*) d19, d20, d21, d22
         read(iunit,*) d23, d24, d25, d26
         read(iunit,*) a_GNFW, b_GNFW, c_GNFW, c500_GNFW, Ytot, theta_s

! Write the data to the correct variables
         beta = d1
         theta_cx = d2
         n0 = d3
         z = d4
         H0 = d5
         T_central = d6
         back = d6a
         T_halo = d7
         Tcorerad = d8
         cellsize = d9
         pb_x = d10
         pb_y = d11
         sz_x = d12
         sz_y = d13
         sz_z = d14
         x_exp = d15
         const = d16
         Te = d17
         el_alpha = d18
         ang(1) = el_alpha
         theta_c = da1
!         temp_halo = l1
         ellipsoid = l2
         ang = da2
         nx_2 = d19
         beta2 = d20
         cutoff1 = d21
         cutoff2 = d22
         poly_a = d23
         poly_b = d24
         poly_c = d25
         poly_d = d26

! Rewind and close file
         rewind (iunit)
         close (iunit)
      end if

! Calculate the angular size distance and luminosity distance in metres
      call lumdist
      D_lum = D_lum*Mpc         ! D_lum in metres now....
      call angdist
      D_theta = D_theta*Mpc
      write(*,*) 'Luminosity distance is : ',D_lum/Mpc,' Mpc'
      write(*,*) 'Angular distance is    : ',D_theta/Mpc,' Mpc'

   end subroutine read_model

!***************************************************************************

   subroutine read_vis_data

      use kind_def
      use sz_globals

      implicit none

      real, dimension(mvis) :: su,sv,sw,srms,sdata_re,sdata_im
      character(len=fname_len) :: filename
      integer :: dummy, iunit, i
      logical :: ex
      character*1 chr1

      su = 0.0
      sv = 0.0
      sw = 0.0
      sdata_re = 0.0
      sdata_im = 0.0
      srms = 0.0

! Get filename
      prompt = 'vis filename'
      which_dir = profile_data_dir
      rootname = 'default'
      extn = '.vis'
      call get_read_filename(filename)

! Open file
      call io_nxtlun(iunit,status)
      open (iunit,file = filename)

! Get visibility data
      if (wis) then
         do i = 1,mvis
            read (iunit,*,end=101) dummy,su(i),sv(i),sw(i),sdata_re(i),&
                 & sdata_im(i),srms(i)
         end do
      else
         do i = 1,mvis
            read (iunit,*,end=101) dummy,su(i),sv(i),sdata_re(i),&
                 & sdata_im(i),srms(i)
         end do
      end if

! Come here at end of file and find how many visibilities have been read in
101   nvis = i-1
      write(*,*) 'Read ',nvis,' visibilities.'

! Rewind and close file
      close (unit=iunit)
      rewind (iunit)

! Write data to correct variables
      u = su
      v = sv
      rms = srms
      call io_getc('Use (N)atural or (W)hitening weighting?','*',chr1,status)
      if ((chr1=='n').or.(chr1=='N')) then
         do i = 1,nvis
            data_re(i) = sdata_re(i)
            data_im(i) = sdata_im(i)
            weight(i) = 1./rms(i)**2
         end do
      else
         do i = 1,nvis
            data_re(i) = sdata_re(i)
            data_im(i) = sdata_im(i)
            weight(i) = 1./rms(i)
         end do
      end if

      which_pointing = 1
      n_samp = nvis/(nchan*n_antennas*(n_antennas-1)/2)
      n_pol = 1

   end subroutine read_vis_data

!*************************************************************************

   subroutine write_weights(wtfile)
   
      use kind_def
      use sz_globals
      
      implicit none

      character(len=fname_len) :: wtfile, comment
      integer, parameter :: maxgroup = 8
      integer :: igroup(maxgroup)
      real :: group(maxgroup)
      doubleprecision :: group8(maxgroup/2)
      integer :: iwfile, l, iptr, nbytes, i, chr_lenb
      logical :: ex
      character*1 :: chr1
      integer :: vis_per_samp, j
      real*4 :: sum_wt, ra
      external chr_lenb
      equivalence (igroup,group,group8)
      
      call io_nxtlun(iwfile,status)
      open(iwfile,file=wtfile,status='UNKNOWN',iostat=status)
      if (status.ne.0) then
        l = chr_lenb(wtfile)
        write(*,*) '*** problem opening '//wtfile(1:l)
      else
        iptr = 0
        nbytes = 32
        vis_per_samp = nvis/n_samp
        do i = 1, n_samp
          igroup(1) = i
          igroup(2) = which_pointing((i-1)*vis_per_samp + 1)
          ra = obs_raref*hr2rad+(i-1)*ha_inc
          if (ra.gt.2.0*pi) then
            ra = ra - 2.0*pi
          endif
          group8(2) = ra
          group8(3) = obs_decref*deg2rad
          sum_wt = 0.0
          do j = 1, vis_per_samp
            sum_wt = sum_wt + weight((i-1)*vis_per_samp + j)
          end do
          group8(4) = sum_wt
          call util_i4cvt(igroup,2)
          call util_r8cvt(group8(2),3)
          call io_fwrite(iwfile,iptr,group,nbytes,status)
          iptr = iptr+nbytes
        enddo
        call io_clrerr(status)
        call io_close(iwfile,status)
      endif
      
   end subroutine write_weights

!*************************************************************************

!  writes out a particular model
   subroutine write_model

      use kind_def
      use sz_globals

      implicit none

      character(len=fname_len) :: filename, comment
      integer :: iunit
      logical :: ex
      character*1 :: chr1
      integer :: write_time(3), write_date(3),j

! Get filename
      prompt = 'Model filename'
      extn = '.mod'
      which_dir = model_dir
      rootname = 'model'
      call get_write_filename(filename)

! Open file
      call io_nxtlun(iunit,status)
      open (iunit,file=filename,form='formatted')

! Write data to file
      write(iunit,'(A)') 'model'
      write(iunit,*) model_type
      write(iunit,*) beta, theta_cx, n0, z, H0, T_central, back
      write(iunit,*) T_halo, Tcorerad, cellsize, pb_x, pb_y, sz_x, sz_y, sz_z
      write(iunit,*) x_exp, const, Te, el_alpha
      write(iunit,*) theta_c(1), theta_c(2), theta_c(3), ellipsoid
      write(iunit,*) ang(1),ang(2),ang(3)
      write(iunit,*) nx_2, beta2, cutoff1, cutoff2
      write(iunit,*) poly_a, poly_b, poly_c, poly_d
      write(iunit,*) a_GNFW, b_GNFW, c_GNFW, c500_GNFW, Ytot, theta_s

      call itime(write_time)
      call idate(write_date)
      write(iunit,'(A,2(I2,A),I2,4X,2(I2,A),I4)') &
            '* File written ', &
            (write_time(j),'.',j=1,2),write_time(3),&
            (write_date(j),'/',j=1,2),write_date(3)

      call io_getstr('comment>', ' ', comment, status)
      write(iunit,'(2A)') '* ',comment

! Close the file
      close (iunit)
      rewind (iunit)

   end subroutine write_model


! ***********************************************************************

   subroutine write_flagpixel(flag,mapsize)

      use kind_def

      implicit none
      logical, dimension(mapsize,mapsize), intent(in) :: flag
      integer, intent(in) :: mapsize
      character(len=fname_len) :: filename
      integer :: iunit
      logical :: ex
      character*1 :: chr1

! Get filename
      prompt = 'Region filename'
      extn = '.reg'
      which_dir = profile_data_dir
      rootname = 'region'
      call get_write_filename(filename)

! Open file
      call io_nxtlun(iunit,status)
      open (iunit,file=filename,form='formatted')
      write(iunit,'(A)') 'region'
      write(iunit,*) mapsize
      write(iunit,*) flag
      close(iunit)
      rewind(iunit)

   end subroutine write_flagpixel

! **************************************************************************

   subroutine read_flagpixel(flag,mapsize)

      use kind_def
      use sz_globals

      implicit none
      logical, dimension(mapsize,mapsize), intent(out) :: flag
      integer, intent(in) :: mapsize
      character(len=fname_len) :: filename, file_type
      integer :: filesize, iunit
      logical :: ex

      prompt = 'Region filename'
      which_dir = profile_data_dir
      rootname = 'region'
      extn = '.reg'
      call get_read_filename(filename)

! Open file
      call io_nxtlun(iunit,status)
      open (iunit,file=filename,form='formatted')
      read(iunit,'(A)') file_type
      if (file_type(1:6).ne.'region') then
         write(*,*) 'Not a region file'
         return
      else
         read(iunit,*) filesize
         if (filesize.ne.mapsize) then
            write(*,*) 'Error reading file'
            write(*,*) filename,' has array of size ',filesize
            flag = .true.
         else
            read(iunit,*) flag
         end if
         close(iunit)
         rewind(iunit)
      end if

   end subroutine read_flagpixel

! ***********************************************************************

   subroutine read_nt(ndat,Tdat,num_nt,dat_read)

      implicit none
      integer, intent(in) :: num_nt
      integer, intent(out) :: dat_read
      real(kind=dp), dimension(num_nt), intent(out) :: ndat, Tdat
      character(len=fname_len) :: filename
      integer :: i,iunit
      logical :: ex

! Get filename
      prompt = 'File name'
      which_dir = profile_data_dir
      rootname = 'nt'
      extn = '.dat'
      call get_read_filename(filename)

      call io_nxtlun(iunit,status)
      open(iunit,file=filename, action='read')
      do i=1,num_nt
         read(iunit,*,end=113) ndat(i),Tdat(i)
      end do

113   dat_read=i-1

      rewind (iunit)
      close (unit=iunit)

   end subroutine read_nt

!**************************************************************************

   subroutine read_telescope

      use sz_globals

      implicit none

      integer :: iunit, i, old_n_antennas, old_nchan, io_stat
      real(kind=dp) :: dummy1
      character(len=fname_len) :: filename, file_type, temp_line
      logical :: ex

! Ask for file name
      prompt = 'Telescope filename'
      which_dir = telescope_dir
      rootname = 'amidc_sa'
      extn = '.tel'
      call get_read_filename(filename)

! Open the file
      call io_nxtlun(iunit,status)
      open(iunit,file = filename,form = 'formatted')

! Read the saved data
      read(iunit,'(A)') file_type
      if (file_type(1:4).ne.'tele') then
         write(*,*) 'Not a telescope file'
         close(unit=iunit)
         return
      else
         old_n_antennas = n_antennas
         old_nchan = nchan
         read(iunit,*) n_antennas
         read(iunit,*) nchan

         if ((n_antennas.ne.old_n_antennas).or.&
             (nchan.ne.old_nchan)) then
            call deallocate_telescope
            call allocate_telescope
         end if

         if (nchan.ne.old_nchan) then
            call do_deallocation
            call do_allocation
         end if

         read(iunit,*) obsfreq, dish_diameter, tel_lat
         
         do i = 1, n_antennas
            read(iunit,*) x_pos(i), y_pos(i), z_pos(i)
         end do
         do i = 1, nchan
            read(iunit,*) nu(i)
         end do

         ! Check for a telescope name - not present in older files
         do while (io_stat==0)
            read(iunit,'(A)',iostat=io_stat) temp_line
            if (temp_line(1:8).eq.'tel_name') then
              tel_name=trim(temp_line(9:fname_len))
              !write(*,*) 'Found telescope name', trim(tel_name)
            endif
         enddo
      end if
     
! Rewind and close file
      close (unit=iunit)
      rewind (iunit)

   end subroutine read_telescope

!**************************************************************************

   subroutine read_specified_telescope(filename,errorcode)

      use sz_globals

      implicit none

      integer :: errorcode
      integer :: iunit, i, old_n_antennas, old_nchan
      real(kind=dp) :: dummy1
      character*80 :: filename, file_type

! Initialise errorcode
      errorcode = 0

! Open the file
      call io_nxtlun(iunit,status)
      open(iunit,file = filename,form = 'formatted',status = 'OLD',iostat=status)

! Check whether file exists
      if (status.ne.0) then
         write(*,*) 'Can not find file ',filename
         close(unit=iunit)
         errorcode = 1
         return
      end if

! Read the saved data
      read(iunit,'(A)') file_type
      if (file_type(1:4).ne.'tele') then
         write(*,*) 'Not a telescope file'
         close(unit=iunit)
         errorcode = 2
         return
      else
         old_n_antennas = n_antennas
         old_nchan = nchan
         read(iunit,*) n_antennas
         read(iunit,*) nchan
         if ((n_antennas.ne.old_n_antennas).or.&
             (nchan.ne.old_nchan)) then
            call deallocate_telescope
            call allocate_telescope
         end if
         read(iunit,*) obsfreq, dish_diameter, tel_lat
         
         do i = 1, n_antennas
            read(iunit,*) x_pos(i), y_pos(i), z_pos(i)
         end do
         do i = 1, nchan
            read(iunit,*) nu(i)
         end do
      end if
      close (unit=iunit)
      rewind (iunit)

   end subroutine read_specified_telescope

!**************************************************************************

   subroutine write_telescope

      implicit none

      integer :: iunit, i
      character(len=fname_len) :: filename, comment
      logical :: ex
      character*1 :: chr1
      integer :: write_time(3), write_date(3),j

! Get filename
      prompt = 'Telescope filename'
      extn = '.tel'
      which_dir = telescope_dir
      rootname = 'default'
      call get_write_filename(filename)

! Open file
      call io_nxtlun(iunit,status)
      open (iunit,file=filename,form='formatted')

! Write data to file
      write(iunit,'(A)') 'tele'
      write(iunit,*) n_antennas
      write(iunit,*) nchan
      write(iunit,*) obsfreq, dish_diameter, tel_lat
      do i = 1, n_antennas
         write(iunit,*) x_pos(i), y_pos(i), z_pos(i)
      end do
      do i = 1, nchan
         write(iunit,*) nu(i)
      end do

      call itime(write_time)
      call idate(write_date)
      write(iunit,'(A,2(I2,A),I2,4X,2(I2,A),I4)') &
            '* File written ', &
            (write_time(j),'.',j=1,2),write_time(3),&
            (write_date(j),'/',j=1,2),write_date(3)

      call io_getstr('comment>', ' ', comment, status)
      write(iunit,'(2A)') '* ',comment

! Close the file
      close (iunit)
      rewind (iunit)

   end subroutine write_telescope

!**************************************************************************

   subroutine write_vis_data_nof

      use kind_def
      use sz_globals

      implicit none
      real, dimension(mvis) :: su,sv,sw,srms,sdata_re,sdata_im
      character(len=fname_len) :: filename
      integer :: iunit, i
      logical :: ex
      character*1 :: chr1

! Get filename
      prompt = 'vis filename'
      if (wis) then
         extn = '.wis'
      else
         extn = '.vis'
      end if
      which_dir = profile_data_dir
      rootname = 'default'
      call get_write_filename(filename)

! Open file
      call io_nxtlun(iunit,status)
      open (iunit,file = filename,form = 'formatted')

! Write variables to data arrays to be outputed
      su = u
      sv = v
      sw = w
      sdata_re = data_re
      sdata_im = data_im
      srms = rms

! Write visibility data
      if (wis) then
         do i = 1,nvis
            write (iunit,*) i,su(i),sv(i),sw(i),sdata_re(i),&
                 & sdata_im(i),srms(i)
         end do
      else
         do i = 1,nvis
            write (iunit,*) i,su(i),sv(i),sdata_re(i),&
                 & sdata_im(i),srms(i)
         end do
      end if

! Rewind and close file
      close (unit=iunit)
      rewind (iunit)

   end subroutine write_vis_data_nof

!**************************************************************************

   subroutine write_vis_data_givenf(visname,sep_vis_chan)

      use kind_def
      use sz_globals

      implicit none
      character(len=fname_len), intent(in) :: visname
      logical, intent(in) :: sep_vis_chan
      real, dimension(mvis) :: su,sv,sw,srms,sdata_re,sdata_im
      character(len=fname_len), dimension(nchan) :: filename
      character(len=fname_len) :: chr1
      character*2 :: chr2
      integer, dimension(nchan) :: iunit, vis_count
      integer :: i, chn, sl1, base_len, fl1, extn_len, dumi
      logical :: ex
      integer :: chr_lenb
      external chr_lenb

! Define filename and open file(s) 
      call io_nxtlun(dumi,status)

      if (dumi.lt.10) dumi = 10

      do chn = 1, nchan
         iunit(chn) = dumi+chn-1
      end do

! Get extension
      if (wis) then 
         extn = '.wis'
      else
         extn = '.vis'
      end if

      if (sep_vis_chan) then
         do chn = 1, nchan
            write(chr2,'(I2.2)') chn
            sl1 = chr_lenb(chr2)
            base_len = chr_lenb(visname)
            extn_len = chr_lenb(extn)
            filename(chn) = visname(1:base_len)//'_'//chr2(1:sl1)//&
                            extn(1:extn_len)
            open (iunit(chn),file = filename(chn))
         end do
      else
         base_len = chr_lenb(visname)
         extn_len = chr_lenb(extn)
         filename(1) = visname(1:base_len)//extn(1:extn_len)
         open (iunit(1),file = filename(1))
      end if

      do chn = 1, nchan
         vis_count(chn) = 0
      end do

! Write variables to data arrays to be outputed
      su = u
      sv = v
      sw = w
      sdata_re = data_re
      sdata_im = data_im
      srms = rms

! Write visibility data
      chn = 1
      if (wis) then
         do i = 1,nvis
            write (iunit(chn),*) i,su(i),sv(i),sw(i),sdata_re(i),&
                 & sdata_im(i),srms(i)
            vis_count(chn) = vis_count(chn)+1
            if (sep_vis_chan) then 
               chn = chn+1
               if (chn.gt.nchan) chn = 1
            end if
         end do
      else
         do i = 1,nvis
            write (iunit(chn),*) i,su(i),sv(i),sdata_re(i),&
                 & sdata_im(i),srms(i)
            vis_count(chn) = vis_count(chn)+1
            if (sep_vis_chan) then 
               chn = chn+1
               if (chn.gt.nchan) chn = 1
            end if
         end do
      end if

      write(*,*) 'Data written'

! Rewind and close file
      if (sep_vis_chan) then 
         do chn = 1,nchan
            chr1 = (filename(chn))
            fl1 = chr_lenb(chr1)
            write(*,*) vis_count(chn),' visibilities written to ',chr1(1:fl1)
            close (unit=iunit(chn))
            rewind (iunit(chn))
         end do
      else
         chr1 = (filename(1))
         fl1 = chr_lenb(chr1)
         write(*,*) vis_count(1),' visibilities written to ',chr1(1:fl1)
         close (unit=iunit(1))
         rewind (iunit(1))
      end if
         
   end subroutine write_vis_data_givenf

!**************************************************************************

   subroutine read_ptsrc_list

      use kind_def
      use sz_globals

      implicit none

      integer  :: src, dummy, iunit
      character(len=fname_len) :: filename, file_type
      real(dp) :: hh, mm, ss, dd, dm, ds, flux, alpha
      logical :: ex

      prompt = 'Source list filename'
      extn = '.src'
      which_dir = profile_data_dir
      rootname = 'list'
      call get_read_filename(filename)

! Open the file
      call io_nxtlun(iunit,status)
      open(iunit,file = filename,form = 'formatted')

! Count the number of sources
      src = 0
100   read (iunit,*,end=101) dummy,hh,mm,ss,dd,dm,ds,flux,alpha
      src = src+1
      goto 100
      
101   rewind(iunit)
      nsrc = src
      
! Allocate arrays
      call deallocate_sources
      call allocate_sources

! Loop over sources and read data
      do src = 1,nsrc
         read (iunit,*) dummy,hh,mm,ss,dd,dm,ds,flux,alpha
         rapos(src) = (hh+mm/60.+ss/3600.)
         decpos(src) = (dd+dm/60.+ds/3600.)
         srcflux(src) = flux
         srcalpha(src) = alpha
      end do

      rewind(iunit)
      close(iunit)

      write(*,*) '>>>Read ',nsrc,' sources from file<<<' 

   end subroutine read_ptsrc_list

!**************************************************************************

! Reads in a primary beam model from a set of gridded measurements
   subroutine read_gridded_pb

      use kind_def

      implicit none

      integer :: iunit, i, j
      character(len=fname_len) :: filename
      logical :: ex

! Ask for file name

      prompt = 'Primary beam filename'
      extn = '.pb'
      which_dir = profile_data_dir
      rootname = 'ryle'
      call get_read_filename(filename)

! Open the file
      call io_nxtlun(iunit,status)
      open(iunit,file = filename,form = 'formatted')

      pb_points = 23
      pb_spacing = 60.0
      call io_geti('Size of pb array','*',pb_points,status)
      call io_getd('Spacing of pb array (arcsec)','*',pb_spacing,status)

! Allocate array
      if (allocated(pb_data)) deallocate(pb_data)
      allocate(pb_data(pb_points,pb_points))

      do i = 1,pb_points
         do j = 1,pb_points
            read(iunit,*) pb_data(i,j)

! Need to interpolate this primary beam data, so make sure it goes negative
! beyond first null - will take absolute value later
!         if (sqrt(1.0*((i-12)**2+(j-12)**2)).gt.6.35) then
!            pb_data(i,j)=-pb_data(i,j)
!         end if
         end do
      end do

   end subroutine read_gridded_pb

! -------------------------------------------------------------------------

   subroutine init_directory

      use kind_def
      use sz_globals

      implicit none
      
      integer :: l
      integer, external :: chr_lenb
      character*100 :: exec_dir

      ! Set default directory as the one where the executable is called from
      call getarg(0, exec_dir)
      l = scan(exec_dir, '/', .true.)
      if (l.gt.0) exec_dir=exec_dir(1:l)

      call getenv('PROFILE_DATA', profile_data_dir)
      l = chr_lenb(profile_data_dir)
      if (l.eq.0) then
         profile_data_dir = trim(exec_dir)//'data/'
         l = chr_lenb(profile_data_dir)
      end if
      if (profile_data_dir(l:l).ne.'/') profile_data_dir(l+1:)='/'

      call getenv('PROFILE_MODEL', model_dir)
      l = chr_lenb(model_dir)
      if (l.eq.0) then
         model_dir = trim(exec_dir)//'models/'
         l = chr_lenb(model_dir)
      end if
      if (model_dir(l:l).ne.'/') model_dir(l+1:)='/'

      call getenv('PROFILE_TELESCOPE', telescope_dir)
      l = chr_lenb(telescope_dir)
      if (l.eq.0) then
         telescope_dir = trim(exec_dir)//'tel/'
         l = chr_lenb(telescope_dir)
      end if
      if (telescope_dir(l:l).ne.'/') telescope_dir(l+1:)='/'

      call getenv('PROFILE_POSTSCRIPT', postscript_dir)
      l = chr_lenb(postscript_dir)
      if (l.eq.0) then
         postscript_dir = trim(exec_dir)//'plots/'
         l = chr_lenb(postscript_dir)
      end if
      if (postscript_dir(l:l).ne.'/') postscript_dir(l+1:)='/'

      call getenv('PROFILE_SCRIPTS', profile_scripts_dir)
      l = chr_lenb(profile_scripts_dir)
      if (l.eq.0) then
         profile_scripts_dir = trim(exec_dir)//'scripts/'
         l = chr_lenb(profile_scripts_dir)
      end if
      if (profile_scripts_dir(l:l).ne.'/') profile_scripts_dir(l+1:)='/'

   end subroutine init_directory

! -------------------------------------------------------------------------

   subroutine get_directory

      use kind_def
      use sz_globals

      implicit none

      integer :: l
      integer, external :: chr_lenb

      call io_getc('Profile data directory:', '*', profile_data_dir, status)
      l = chr_lenb(profile_data_dir)
      if (profile_data_dir(l:l).ne.'/') profile_data_dir(l+1:)='/'

      call io_getc('Model directory:', '*', model_dir, status)
      l = chr_lenb(model_dir)
      if (model_dir(l:l).ne.'/') model_dir(l+1:)='/'

      call io_getc('Telescope directory:', '*', telescope_dir, status)
      l = chr_lenb(telescope_dir)
      if (telescope_dir(l:l).ne.'/') telescope_dir(l+1:)='/'

      call io_getc('Postscript directory:', '*', postscript_dir, status)
      l = chr_lenb(postscript_dir)
      if (postscript_dir(l:l).ne.'/') postscript_dir(l+1:)='/'

      call io_getc('Scripts directory:', '*', profile_scripts_dir, status)
      l = chr_lenb(profile_scripts_dir)
      if (profile_scripts_dir(l:l).ne.'/') profile_scripts_dir(l+1:)='/'


   end subroutine get_directory

! -------------------------------------------------------------------------

   subroutine get_read_filename(filename)

      use kind_def
      use sz_globals

      implicit none

      character(len=fname_len) :: filename
      character*130 :: s_command, message
      integer :: l_dir, l_root, l_ext
      logical :: ex
      integer, external :: chr_lenb
      logical, external :: io_yesno

! Ask for file name
      get_file : do
         status = 0
         l_dir  = chr_lenb(which_dir)
         l_root = chr_lenb(rootname)
         l_ext  = chr_lenb(extn)
         if (l_root.eq.0) then
            filename = which_dir
         else
            filename = which_dir(1:l_dir)//rootname(1:l_root)//extn(1:l_ext)
         end if
         call io_getfil(prompt,'*',filename,status)

! Check if this file actually exists
         inquire(file=filename,exist=ex)

         if (.not. ex) then
            write(*,*) 'File does not exist'
            if (io_yesno('Try again','yes',status)) then 
               write(*,*)
               message  = 'Available '//extn(1:l_ext)//' files in '//&
                           which_dir(1:l_dir)//' are:-'
               write(*,*) message
               s_command = 'ls '//which_dir(1:l_dir)//'*'//extn(1:l_ext)
               call io_system(s_command,status)
               write(*,*)
            else
               exit get_file
            end if
         else
            exit get_file
         end if

      end do get_file

   end subroutine get_read_filename

! -------------------------------------------------------------------------

   subroutine get_write_filename(filename)

      use kind_def
      use sz_globals

      implicit none

      character(len=fname_len) :: filename
      integer :: l_dir, l_root, l_ext
      logical :: ex
      integer, external :: chr_lenb
      logical, external :: io_yesno

      get_file : do
         status = 0
         l_dir  = chr_lenb(which_dir)
         l_root = chr_lenb(rootname)
         l_ext  = chr_lenb(extn)
         if (l_root.eq.0) then
            filename = which_dir
         else
            filename = which_dir(1:l_dir)//rootname(1:l_root)//extn(1:l_ext)
         end if
         call io_getfil(prompt,'*',filename,status)

! Check if this file already exists
         inquire(file=filename,exist=ex)

         if (ex) then
            if (io_yesno('File already exists. Overwrite (y/n):',&
                         'n',status)) then
               exit get_file
            end if
         else
            exit get_file
         end if

      end do get_file
      
   end subroutine get_write_filename

!***************************************************************************

   subroutine read_flist

      use kind_def
      use sz_globals

      implicit none

      character(len=fname_len) :: filename, fitsname
      integer :: fcount, iunit

      prompt = 'Filename for list of FITS files'
      extn = '.flist'
      which_dir = profile_data_dir
      rootname = 'default'
      call get_read_filename(filename)

 ! Open the file
      call io_nxtlun(iunit,status) 
      open(iunit,file = filename,form = 'formatted')

! Count the number of files
      fcount = 0
100   read (iunit,'(A)',end=101) fitsname
      fcount = fcount+1
      goto 100
      
101   rewind(iunit)
      n_fits_files = fcount

      call deallocate_combine_fits
      call allocate_combine_fits

! Read in filenames
      do fcount = 1, n_fits_files
         read (iunit,'(A)') fitsname
         fits_list(fcount) = fitsname
      end do

      rewind(iunit)
      close(iunit)

   end subroutine read_flist

!***************************************************************************
 
   subroutine write_header_digest

      use kind_def
      use sz_globals

      implicit none

      character(len=fname_len) :: filename
      integer :: fcount, iunit

      prompt = 'Filename for digest of FITS file headers'
      extn = '.digest'
      which_dir = profile_data_dir
      rootname = 'default'
      call get_write_filename(filename)

! Open the file
      call io_nxtlun(iunit,status) 
      open(iunit,file = filename,form = 'formatted')

! Write the number of files
      write(iunit,*) n_fits_files

! Write list of fits files
      do fcount = 1, n_fits_files
         write(iunit,*) fits_list(fcount)
         write(iunit,*) fits_ra(fcount)
         write(iunit,*) fits_dec(fcount)
         write(iunit,*) fits_npix_ra(fcount)
         write(iunit,*) fits_npix_dec(fcount)
         write(iunit,*) fits_bmaj(fcount)
         write(iunit,*) fits_bmin(fcount)
         write(iunit,*) fits_bpa(fcount)
         write(iunit,*) fits_incell(fcount)
      end do

      close(iunit)
      rewind(iunit)

   end subroutine write_header_digest

!***************************************************************************
 
   subroutine read_header_digest

      use kind_def
      use sz_globals

      implicit none

      character(len=fname_len) :: filename
      integer :: fcount, iunit

      prompt = 'Filename for digest of FITS file headers'
      extn = '.digest'
      which_dir = profile_data_dir
      rootname = 'default'
      call get_read_filename(filename)

! Open the file
      call io_nxtlun(iunit,status) 
      open(iunit,file = filename,form = 'formatted')

! Read the number of files
      read(iunit,*) n_fits_files

! Do array allocation
      call deallocate_combine_fits
      call allocate_combine_fits

! Read in data
      do fcount = 1, n_fits_files
         read(iunit,'(a)') fits_list(fcount)
         read(iunit,*) fits_ra(fcount)
         read(iunit,*) fits_dec(fcount)
         read(iunit,*) fits_npix_ra(fcount)
         read(iunit,*) fits_npix_dec(fcount)
         read(iunit,*) fits_bmaj(fcount)
         read(iunit,*) fits_bmin(fcount)
         read(iunit,*) fits_bpa(fcount)
         read(iunit,*) fits_incell(fcount)
      end do

      close(iunit)
      rewind(iunit)

   end subroutine read_header_digest     

!***************************************************************************
 
   subroutine read_ascii_beam(filename)

      use kind_def
      use sz_globals
      use astronomy
!      use kgraph

      implicit none

      character(len=fname_len) :: filename, dummy
      integer :: iunit, fcount, x_p, y_p
      integer :: max_x, max_y, min_x, min_y, npix, cent_x, cent_y
      real(kind = dp) :: ra2, dec2, az, el, tmp_re, tmp_im, insize, dx, dy
      real(kind = dp), allocatable, dimension(:,:) :: beam_re, beam_im

      ra2 = 0.0
      dec2 = pi/2.0
      insize = 720.
      max_x = -1e9
      max_y = -1e9
      min_x = 1e9
      min_y = 1e9

! Open the file
      call io_nxtlun(iunit,status) 
      open(iunit,file = filename,form = 'formatted')

!  Find required mapsize

!  Read in first line and ignore
      read (iunit,'(A)') dummy

!  Read in data points for beam
100   read (iunit,*,end=101) az, el, tmp_re, tmp_im
      az = az*deg2rad
      el = el*deg2rad
      call sin_proj_forward(az,el,ra2,dec2,insize,dx,dy)
      x_p = nint(dx)
      y_p = nint(dy)
      if (x_p.gt.max_x) max_x = x_p
      if (x_p.lt.min_x) min_x = x_p
      if (y_p.gt.max_y) max_y = y_p
      if (y_p.lt.min_y) min_y = y_p
      fcount = fcount+1
      goto 100
     
101   rewind(iunit)
 
      write(*,*) min_x,min_y,max_x,max_y,fcount

      npix = 1
      do while ((npix.lt.(max_x-min_x+1)).or.(npix.lt.(max_y-min_x+1)))
         npix = npix*2
      end do

      cent_x = npix/2+1
      cent_y = npix/2+1

      write(*,*) npix

      if (allocated(beam_re)) deallocate(beam_re)
      if (allocated(beam_im)) deallocate(beam_im)
      allocate(beam_re(npix,npix))
      allocate(beam_im(npix,npix))
      beam_re = 0.0d0
      beam_im = 0.0d0

!  Read in first line and ignore
      read (iunit,'(A)') dummy

!  Read in data points for beam
200   read (iunit,*,end=201) az, el, tmp_re, tmp_im
      az = az*deg2rad
      el = el*deg2rad
      call sin_proj_forward(az,el,ra2,dec2,insize,dx,dy)
      x_p = nint(dx)+cent_x
      y_p = nint(dy)+cent_y
      beam_re(x_p,y_p) = tmp_re
      beam_im(x_p,y_p) = tmp_im

      goto 200      

201   close(iunit)
      rewind(iunit)

!      call display_map(beam_re,npix,graph_adj,insize,insize,0D0,0D0)

   end subroutine read_ascii_beam

!***************************************************************************
 
   subroutine read_ascii_map(filename)

      use kind_def
      use sz_globals
      use astronomy
      implicit none

      character(len=fname_len) :: filename, dummy
      integer :: iunit, nx, ny
      real(kind = dp) :: val

! Initialise array
      p_sky = 0.d0

! Open the file
      call io_nxtlun(iunit,status) 
      open(iunit,file = filename,form = 'formatted')

!  Read in data points
100   read (iunit,*,end=101) nx, ny, val
      if ((nx.ge.1).and.(nx.le.maxsize).and.(ny.ge.1).and.(ny.le.maxsize))then
         p_sky(nx,ny) = val
      end if

      goto 100

101   rewind(iunit)   

   end subroutine read_ascii_map

 !***************************************************************************

end module file_handling



