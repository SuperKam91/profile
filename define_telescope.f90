module define_telescope
   
   use kind_def
   use sz_globals
   use kgraph
   use file_handling

   implicit none

   public :: def_telescope
   public :: edit_telescope
   public :: show_telescope
   public :: plot_telescope
   public :: get_telescope
   public :: def_pb

   logical :: write_num

contains

!***************************************************************************

   subroutine def_telescope

      use file_handling

      integer :: i, old_n_antennas, old_nchan
      real(kind=dp) :: bandwidth, f_step, f_start
      character*1 :: chr1
      logical :: io_yesno
      external io_yesno

      write_num = .true.

      call io_getc('Enter telescope positions by hand, from file or edit','f',&
                   chr1,status)     

      if (chr1.eq.'h') then

         old_n_antennas = n_antennas
         old_nchan = nchan
         call io_geti('Number of antennas for telescope','*',n_antennas,status)
         call io_geti('Number of frequency channels','*',nchan,status)

         if ((n_antennas.ne.old_n_antennas).or.&
             (nchan.ne.old_nchan)) then
            call deallocate_telescope
            call allocate_telescope
         end if

         if (nchan.ne.old_nchan) then
            call do_deallocation
            call do_allocation
         end if

         obsfreq = obsfreq/1.0D9
         call io_getd('What frequency to observe at (GHz):','*',obsfreq,status)
         obsfreq = obsfreq*1.0D9
         if (nchan.eq.1) then
            nu(1) = obsfreq
         else
            call io_getd('Total bandwidth (GHz)','4.5',bandwidth,status)
            bandwidth = bandwidth*1.0D9
            f_step = bandwidth/(2.0*nchan)
            f_start = obsfreq-bandwidth/2.0
            do i = 1, nchan
               nu(i) = f_start+(2*i-1)*f_step
               write(*,*) 'Channel ',i,' centre frequency ',nu(i)/1.0D9
            end do
         end if
         call io_getd('What dish diameter','*',dish_diameter,status)       
         tel_lat = tel_lat/deg2rad
         call io_getd('What telescope latitude','*',tel_lat,status)
         tel_lat = tel_lat*deg2rad

         do i = 1, n_antennas
            write(*,*) 
            write(*,*) 'Antenna ',i
            call io_getr('x position ','*',x_pos(i),status)
            call io_getr('y position ','*',y_pos(i),status)
            call io_getr('z position ','*',z_pos(i),status)
         end do
         if (io_yesno('write telescope file (y/n)','y',status)) then
            call write_telescope
         end if
      else if (chr1.eq.'f') then
         call read_telescope
      else
         call edit_telescope
         if (io_yesno('write telescope file (y/n)','y',status)) then
            call write_telescope
         end if
      end if

   end subroutine def_telescope

!***************************************************************************

   subroutine edit_telescope

      use file_handling

      implicit none

      integer :: indev, pgopen, i, j, l1, istat, k, n_uv_samp
      external pgopen
      real(kind=dp) :: bandwidth, f_step, f_start, mean_x, mean_y
      character*1 :: chr1
      character*80 :: plot_device, text
      character(len=fname_len) :: ps_file
      real :: minpos, maxpos, max_basel, x, y, diam, datmin, datmax
      real(kind=dp) :: ha, dec, ha_start, ha_stop, inc_ha
      logical :: uv_sampling
      integer :: nbin

      integer, external :: chr_lenb
      logical, external :: io_yesno

      datmin = 0.0
      datmax = 1000.0
      dec = obs_decref*deg2rad
      nbin = 10
      uv_sampling = .true.

      obsfreq = obsfreq/1.0D9
      call io_getd('What frequency to observe at (GHz):','*',obsfreq,status)
      obsfreq = obsfreq*1.0D9
      if (nchan.eq.1) then
         nu(1) = obsfreq
      else
         call io_getd('Total bandwidth (GHz)','6.0',bandwidth,status)
         bandwidth = bandwidth*1.0D9
         f_step = bandwidth/(2.0*nchan)
         f_start = obsfreq-bandwidth/2.0
         do i = 1, nchan
            nu(i) = f_start+(2*i-1)*f_step
            write(*,*) 'Channel ',i,' centre frequency ',nu(i)/1.0D9
         end do
      end if
      call io_getd('What dish diameter','*',dish_diameter,status)       
      tel_lat = tel_lat/deg2rad
      call io_getd('What telescope latitude','*',tel_lat,status)
      tel_lat = tel_lat*deg2rad
      call io_getc('Input plot device:','/xwindow',plot_device,status)
      indev = pgopen(plot_device)
      call io_getr('Enter max baseline (m): ','60',max_basel,status)

      call pgsci(plot_col)
      call pgslct(indev)
      call pgscf(2)
      call pgsfs(2)
      call pgask(.false.)
      call pgpap(8.0,1.0)

      maxpos = max_basel/2+diam
      minpos = -max_basel/2-diam
      ha_start = -3.5*15.0*deg2rad
      ha_stop = 3.5*15.0*deg2rad
      n_uv_samp = 20
      write(*,*) 'For plotting uv coverage....'
      dec = dec/deg2rad
      call io_getd('What obs. declination?','*',dec,status)
      dec = dec*deg2rad
      ha_start = ha_start/(deg2rad*15.0)
      call io_getd('HA start','*',ha_start,status)
      ha_start = ha_start*deg2rad*15.0
      ha_stop = ha_stop/(deg2rad*15.0)
      call io_getd('HA stop','*',ha_stop,status)
      ha_stop = ha_stop*(deg2rad*15.0)
      call io_geti('Resolution of uv tracks','*',n_uv_samp,status)
      if (n_uv_samp.ne.1) then
         inc_ha = (ha_stop-ha_start)/(n_uv_samp-1)
      else
         inc_ha = 0.
      end if

      call calc_basel
      call pgsci(plot_col)
      call pgslct(pl_dev1)
      call plot_uv(max_basel,n_uv_samp,ha_start,inc_ha,dec,.false.)
      call pgslct(indev)

      write(*,*) 'Press H for Help'

      edit_loop : do
         status = 0
         call plot_telescope(max_basel)
         call pgcurs(x,y,chr1)
         j = iachar(chr1)-48
         if ((j.ge.49).and.(j.le.74)) then
            j = j-39
         end if
         if ((j.gt.0).and.(j.le.9).and.(j.le.n_antennas)) then
            x_pos(j) = x
            y_pos(j) = y

            call calc_basel
            call pgslct(pl_dev1)
            if (uv_sampling) then
               call plot_uv(max_basel,n_uv_samp,ha_start,inc_ha,dec,.false.)
            else
               call  plot_uvhist(n_uv_samp,ha_start,inc_ha,dec,datmin,datmax,nbin)
            end if   
            call pgslct(indev)

         end if

         if ((chr1=='Q').or.(chr1=='q')) then
            exit edit_loop 

         else if ((chr1=='A').or.(chr1=='a')) then
            call io_geti('Which antenna','1',j,status)
            if ((j.gt.0).and.(j.le.n_antennas)) then
               x_pos(j) = x
               y_pos(j) = y

               call calc_basel
               call pgslct(pl_dev1)
               if (uv_sampling) then
                  call plot_uv(max_basel,n_uv_samp,ha_start,inc_ha,dec,.false.)
               else
                call plot_uvhist(n_uv_samp,ha_start,inc_ha,dec,datmin,datmax,nbin)
               end if
               call pgslct(indev)
            end if
            
         else if ((chr1=='H').or.(chr1=='h')) then
            write(*,*) 'Help'
            write(*,*)
            write(*,*) 'Q        quit'
            write(*,*) 'P        write postscript file'
            write(*,*) 'O        change observing parameters'
            write(*,*) 'B        change maximum baseline'
            write(*,*) 'C        centre telescope'
            write(*,*) 'N        toggle writing of antenna numbers'
            write(*,*) 'A        place an antenna at this point'
            write(*,*) 'V        write a template .vis file'
            write(*,*) 'S        plot sampling in uv plane'
            write(*,*) 'G        plot histogram of uv distances'
            write(*,*) '(1 - 9)  shortcut to place relevant antenna'

         else if ((chr1=='O').or.(chr1=='o')) then
            dec = dec/deg2rad
            call io_getd('Declination','*',dec,status)
            dec = dec*deg2rad
            ha_start = ha_start/(deg2rad*15.0)
            call io_getd('HA start','*',ha_start,status)
            ha_start = ha_start*deg2rad*15.0
            ha_stop = ha_stop/(deg2rad*15.0)
            call io_getd('HA stop','*',ha_stop,status)
            ha_stop = ha_stop*(deg2rad*15.0)
            call io_geti('Number of points','*',n_uv_samp,status)
            inc_ha = (ha_stop-ha_start)/(n_uv_samp-1)
            call calc_basel
            call pgsci(plot_col)
            call pgslct(pl_dev1)
            if (uv_sampling) then
               call plot_uv(max_basel,n_uv_samp,ha_start,inc_ha,dec,.false.)
            else
               call plot_uvhist(n_uv_samp,ha_start,inc_ha,dec,datmin,datmax,nbin)
            end if
            call pgslct(indev)

         else if ((chr1=='B').or.(chr1=='b')) then
            call io_getr('Enter max baseline (m): ','*',max_basel,status)
            call pgsci(plot_col)
            call pgslct(pl_dev1)
            if (uv_sampling) then
               call plot_uv(max_basel,n_uv_samp,ha_start,inc_ha,dec,.false.)
            else
               call plot_uvhist(n_uv_samp,ha_start,inc_ha,dec,datmin,datmax,nbin)
            end if
            call pgslct(indev)

         else if ((chr1=='C').or.(chr1=='c')) then
            mean_x = 0.0
            mean_y = 0.0
            do j = 1,n_antennas
               mean_x = mean_x+x_pos(j)
               mean_y = mean_y+y_pos(j)
            end do
            mean_x = mean_x/n_antennas
            mean_y = mean_y/n_antennas
            do j = 1,n_antennas
               x_pos(j) = x_pos(j)-mean_x
               y_pos(j) = y_pos(j)-mean_y
            end do

         else if ((chr1=='N').or.(chr1=='n')) then

            write_num = .not. write_num

         else if ((chr1=='P').or.(chr1=='p')) then
            
            prompt = 'Filename for postscript plot'
            which_dir = postscript_dir
            rootname = 'plot'
            extn = '.ps'
            call get_write_filename(ps_file)

            l1 = chr_lenb(ps_file)
            ps_file = ps_file(1:l1)//'/CPS'
            istat = pgopen(ps_file)
            if (istat.le.0) then
               write(*,*) 'Error opening postscript file'
               call pgslct(indev)
            else
               call plot_telescope(max_basel)
               if (uv_sampling) then
                  call plot_uv(max_basel,n_uv_samp,ha_start,inc_ha,dec,.false.)
               else
                  call plot_uvhist(n_uv_samp,ha_start,inc_ha,dec,datmin,datmax,&
                                   nbin)
               end if
               call pgclos
               call pgslct(indev)
            end if
         else if ((chr1=='V').or.(chr1=='v')) then

            wis = io_yesno('Write out w coordinate as well','no',status)
            call plot_uv(max_basel,n_uv_samp,ha_start,inc_ha,dec,.true.)

         else if ((chr1=='S').or.(chr1=='s')) then
            uv_sampling = .true.
            call pgslct(pl_dev1)
            call plot_uv(max_basel,n_uv_samp,ha_start,inc_ha,dec,.false.)
            call pgslct(indev)

         else if ((chr1=='G').or.(chr1=='g')) then
            uv_sampling = .false.
            call io_getr('Minimum uv distance','*',datmin,status)
            call io_getr('Maximum uv distance','*',datmax,status)
            call io_geti('Number of bins','*',nbin,status)
            call pgslct(pl_dev1)
            call plot_uvhist(n_uv_samp,ha_start,inc_ha,dec,datmin,datmax,nbin)
            call pgslct(indev)

         end if

      end do edit_loop

      call pgslct(indev)
      call pgclos
      call pgslct(pl_dev1)
      call pgsci(plot_col)
     
   end subroutine edit_telescope

!***************************************************************************

   subroutine show_telescope
      
      integer :: i, n_uv_samp
      real :: minpos, maxpos, max_basel, x, y, diam
      real(kind=dp) :: ha, dec, ha_start, ha_stop, inc_ha

      write(*,*) 'Number of antennas           ',n_antennas
      write(*,*) 'Number of frequency channels ',nchan
      write(*,*) 'Frequency (GHz)              ',obsfreq/1.0D9
      write(*,*) 'Dish diameter (m)            ',dish_diameter
      write(*,*) 'Latitude                     ',tel_lat/deg2rad
      write(*,*)
      if (nchan.ne.1) then
         do i = 1, nchan
            write(*,*) 'Channel ',i,' centre frequency ',nu(i)/1.0D9
         end do
      end if

   end subroutine show_telescope

!***************************************************************************

   subroutine plot_telescope(max_basel)
      
      integer :: i, k
      real :: minpos, maxpos, diam, max_basel
      character*80 :: text
      
      diam = dish_diameter
      maxpos = max_basel/2+diam
      minpos = -max_basel/2-diam
      call pgsci(plot_col)
      call pgslw(line_width)
      call pgsch(char_hgt)
      call pgvsize(2., 6.5, 1., 5.5)
      call pgscf(2)
      call pgsfs(2)
      call pgenv(minpos,maxpos,minpos,maxpos,1,1)
      call pgscf(3)
      call pglabel('East - West','North - South','Antenna Positions')
      do i = 1,n_antennas
         k = i+1
         colour_alloc: do
            if (k.gt.12) then
               k=k-10
            else
               exit colour_alloc
            end if
         end do colour_alloc
         call pgsci(k)
         call pgcirc(x_pos(i),y_pos(i),diam/2.0)
         if (write_num) then
            call pgnumb(i,0,1,text,status)
            status = 0
            call pgtext(x_pos(i),y_pos(i),text)
         end if
      end do
      call pgsci(plot_col)
      
   end subroutine plot_telescope

!***************************************************************************

!Changes the S-Z and X-ray telescopes whose observations are being modelled
   subroutine get_telescope

      use kind_def
      use sz_globals

      implicit none

      character*1 :: chr1
      logical, external :: io_yesno

      status = 0

      write(*,*) 'Possible S-Z telescopes :-'
      write(*,*) 'User (D)efined telescope'
      write(*,*) '(R)yle'
      write(*,*) '(V)LA'
      write(*,*) '(C)AT'
      write(*,*) 'VSA:- (M)ain array of VSA (compact)'
      write(*,*) '    - (E)xtended Main array of VSA'
      write(*,*) '    - SuperEx(T)ended VSA'
      write(*,*) '    - (S)ource subtractor for VSA'
      write(*,*) 'Cyg(X) ext array mosaic'
      write(*,*) 'S(Z)I'
      call io_getc('Select one (D/R/V/C/M/E/T/S/X/Z):','d',which_tel,status)

! Set dish diameter and frequency for each telescope
      select case(which_tel)
      case('v','V')
         obsfreq = obsfreq/1.0D9
         call io_getd('What frequency to observe at (GHz):','*',obsfreq,status)
         obsfreq = obsfreq*1.0D9
         dish_diameter = vla_dish_diameter
         nchan = 1
         call deallocate_frequency
         call allocate_frequency
         nu(1) = obsfreq   
         fwhm(1) = const_c/(obsfreq*dish_diameter)
         which_pb = 'v'
      case('r','R')
         obsfreq = rt_freq
         dish_diameter = rt_dish_diameter
         nchan = 1
         call deallocate_frequency
         call allocate_frequency
         nu(1) = obsfreq      
         fwhm(1) = const_c/(obsfreq*dish_diameter)
         which_pb = 'r'
      case('c','C')
         obsfreq = cat_freq
         dish_diameter = cat_dish_diameter
         nchan = 1
         call deallocate_frequency
         call allocate_frequency
         nu(1) = obsfreq
         fwhm(1) = const_c/(obsfreq*dish_diameter)
         which_pb = 'c'
      case('m','M')
         obsfreq = vsa_freq/1.0D9
         call io_getd('What frequency to observe at (GHz):','*',obsfreq,status)
         obsfreq = obsfreq*1.0D9
         dish_diameter = vsa_ma_dish_diameter
         nchan = 1
         call deallocate_frequency
         call allocate_frequency
         nu(1) = obsfreq
         which_pb = 'c'
         fwhm(1) = vsa_ma_fwhm
      case('e','E')
         obsfreq = vsa_freq/1.0D9
         call io_getd('What frequency to observe at (GHz):','*',obsfreq,status)
         obsfreq = obsfreq*1.0D9
         dish_diameter = vsa_ea_dish_diameter
         nchan = 1
         call deallocate_frequency
         call allocate_frequency
         nu(1) = obsfreq
         which_pb = 'c'
         fwhm(1) = vsa_ea_fwhm
      case('t','T')
         obsfreq = vsa_freq/1.0D9
         call io_getd('What frequency to observe at (GHz):','*',obsfreq,status)
         obsfreq = obsfreq*1.0D9
         dish_diameter = vsa_sa_dish_diameter
         nchan = 1
         call deallocate_frequency
         call allocate_frequency
         nu(1) = obsfreq
         which_pb = 'c'
         fwhm(1) = vsa_ma_fwhm
      case('x','X')
         obsfreq = vsa_freq/1.0D9
         call io_getd('What frequency to observe at (GHz):','*',obsfreq,status)
         obsfreq = obsfreq*1.0D9
         dish_diameter = vsa_ea_dish_diameter
         nchan = 1
         call deallocate_frequency
         call allocate_frequency
         nu(1) = obsfreq
         fwhm(1) = const_c/(obsfreq*dish_diameter)
         which_pb = 'g'
         fwhm(1) = vsa_ea_fwhm
      case('s','S')
         obsfreq = vsa_freq/1.0D9
         call io_getd('What frequency to observe at (GHz):','*',obsfreq,status)
         obsfreq = obsfreq*1.0D9
         dish_diameter = vsa_ss_dish_diameter
         nchan = 1
         call deallocate_frequency
         call allocate_frequency
         nu(1) = obsfreq
         which_pb = 'r'
         fwhm(1) = vsa_ss_fwhm
      case('z','Z')
         call io_geti('Number of receivers in fp array:','19',fparray,status)
         call io_getd('Spacing for primary beams [fwhm]:','1.',pbspace,status)
         call io_getc('Are you planning to do a jiggle map?','y',chr1,status)
         if (chr1=='y') jpoints=4
         jpoint = 1
         which_pb = 'g'
         call def_telescope
      case default
         which_pb = 'g'
         call def_telescope
      end select

! Warning if field of view to small
      if (maxsize*cellsize*sec2rad/2.0.lt.&
          1.0/(dish_diameter*obsfreq/const_c)) then
         write(*,*)
         write(*,*) &
'** WARNING ** Field of view to small - increase cellsize or number of cells'
         write(*,*)
      end if

      if (nchan==1) then
         nu(1) = obsfreq
      end if

      if ((which_tel.eq.'d').or.(which_tel.eq.'D')) then
         call def_pb
      else
         if (io_yesno('Change from default primary beam','no',status)) then
            call def_pb
         end if
      end if

   end subroutine get_telescope

!***************************************************************************

   subroutine def_pb

      implicit none

      integer :: i
      character*1 :: chr1
      real :: cent_fwhm
      real(dp) ::  conv_fact
      logical, external :: io_yesno
 
      write(*,*) 'Possible Primary Beam models :-'
      write(*,*) '(P)olynomial'
      write(*,*) '(G)aussian'
      write(*,*) '(E)lliptical'
      write(*,*) 'Spherical (H)armonic'
      write(*,*) '(C)AT scaled'
      write(*,*) '(R)yle scaled'
      write(*,*) '(V)LA scaled'
      write(*,*) 'Gridded from (F)ile'
      write(*,*) 'Approximate (S)urvey beam'
      write(*,*) 'S(Z)I'
      write(*,*) 'Cyg(X) ext array mosaic'
      write(*,*) '(U)nity beam'
      write(*,*) '(A)MI-SA beam'
      write(*,*) 'AMI-(L)A beam'
      call io_getc('Select one:','*',which_pb,status)

      select case(which_pb)
      case('p','P')
         call io_geti('What order polynomial','*',n_poly_pb,status)
         if (allocated(poly_pb)) deallocate(poly_pb)
         allocate(poly_pb(n_poly_pb))
         do i = 1, n_poly_pb
            call io_getr('Polynomial coeff','*',poly_pb(i),status)
         end do
      case('g','G') 
         if ((nchan.eq.1).or.&
              (io_yesno('Scale FWHM by 1/frequency','yes',status))) then
            cent_fwhm = const_c/(obsfreq*dish_diameter*sec2rad)
            write(*,*) 'Central frequency is', obsfreq, 'Hz'
            call io_getr('FWHM (arcsec)','*',cent_fwhm,status)
            do chan = 1,nchan
               fwhm(chan) = cent_fwhm*obsfreq/nu(chan)
            end do
         else
            do chan = 1,nchan
               fwhm(chan) = const_c/(nu(chan)*dish_diameter*sec2rad)
               write(*,*) 'Channel', chan, nu(chan), 'Hz'
               call io_getr('FWHM (arcsec)', '*',fwhm(chan),status)
            end do
         end if
      case('h','H')
         write(*,*) 'Available spherical harmonic beams:'
         write(*,*) '1 - SPH10'
         write(*,*) '2 - SPH11'
         write(*,*) '3 - SPH22'
         write(*,*) '4 - SPH31'
         write(*,*) '5 - SPH52'
         call io_geti('Select one ','5',which_sph_harm,status)
         if (which_sph_harm.eq.5) then
            call io_getd('Rotation angle ','45.0',sph_theta,status)
         end if
      case('e','E')
         call io_getd('sigma X [arcmin]:','*',el_sigma1,status)
         call io_getd('sigma Y [arcmin]:','*',el_sigma2,status)
         call io_getd('theta [degrees]:','*',el_theta,status)   
      case('f','F') 
         call read_gridded_pb
      case('a', 'A')
         conv_fact = 2.*sqrt(2.*log(2.))*60.
         do chan = 1, nchan
            fwhm(chan) = (101.1/(nu(chan)*1d-9) + 1.89)*conv_fact
         end do
      case('l', 'L')
         conv_fact = 2.*sqrt(2.*log(2.))*60.
         do chan = 1, nchan
            fwhm(chan) = (24.905/(nu(chan)*1d-9) + 0.79)*conv_fact
         end do
      case('s','S')
         call io_getd('survey beam size','3600.0',pb_sbeam,status)
         call io_getd('rolloff scale','600.0',pb_rolloff,status)
      case('u','U')
         write(*,*) '*** Setting primary beam to unity - for testing only'
      end select
      

   end subroutine def_pb

!***************************************************************************

end module define_telescope
            







