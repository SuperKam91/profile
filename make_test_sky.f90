! makes a test sky
subroutine make_test_sky

   use kind_def
   use sz_globals
   use maths
   use kgraph

   implicit none

   character*1 :: chr1
   real(kind=dp) :: alpha, sig_lev, d_ra, d_dec, x1, y1, gau_val
   real(kind=dp) :: sigma, rad, gau_fwhm, c_rad, ring_rad, ring_wid
   real(kind=dp) :: dec, ra, norm
   real(kind=dp), dimension(maxsize) :: gau_tab
   integer :: dis_ra, dis_dec, centre, i, j, k
   logical :: do_jansky
   logical, external :: io_yesno

   centre = maxsize/2+1

   do_plot = io_yesno('Do plots','n',status)
 
   overwrite = io_yesno('Overwrite existing model sky','y',status)

   write(*,*) 'Available test skies:'
   write(*,*) '(P)oint source'
   write(*,*) '(G)aussian source'
   write(*,*) '(D)C sky'
   write(*,*) '(R)ing source'
   call io_getc('What type of test sky','p',chr1,status)

   if (overwrite) then
      write(*,*) 'Overwriting old model sky'
      if (flag_nd.eq.2) then
         p_sky = 0.d0
      else if (flag_nd.eq.3) then
         p3_sky = 0.d0
      end if
   end if

   what_sky : select case(chr1)
   case('p','P')
      call io_geti('Displacement of point in x (pixels)','0',dis_ra,status)
      call io_geti('Displacement of point in y (pixels)','0',dis_dec,status)

! Set a flag which will allow definition of source flux in Janskys if szsky 
! chosen
      do_jansky = .false.
      if (flag_nd.eq.3) then
         if (associated(p3_sky, target=szsky)) then
            do_jansky = .true.
         end if
      end if
      if (flag_nd.eq.2) then 
         if (associated(p_sky, target=szsky(:,:,chan))) then
            do_jansky = .true.
         end if
      end if

      if (do_jansky) then

! Deal with single channel and multi channel differently
         if (flag_nd.eq.2) then
            call io_getd('What flux density for point (Jy)','1.0',&
                       sig_lev,status)

! Convert from flux density to brightness temperature
            sig_lev = sig_lev/((cellsize*sec2rad)**2*2.0*k_b*nu(chan)**2&
                  /const_c2*1.0D26)

            p_sky(centre+dis_ra,centre+dis_dec) = &
                  p_sky(centre+dis_ra,centre+dis_dec) &
                 +sig_lev
         else
            write(*,*) 'Flux density defined at ',obsfreq/1.D9
            call io_getd('What flux density for point (Jy)','1.0',&
                       sig_lev,status)
            call io_getd('Spectral index:','0.7',alpha,status)
            if (verbose) then
               do chan=1,nchan
                  write(*,*) 'Flux in channel ',chan,':',sig_lev*&
                                       (nu(chan)/obsfreq)**(-alpha)
               end do
            end if

! Convert from flux density to brightness temperature
            sig_lev = sig_lev/((cellsize*sec2rad)**2*2.0*k_b*obsfreq**2&
                  /const_c2*1.0D26)
            do chan=1,nchan
               p3_sky(centre+dis_ra,centre+dis_dec,chan) = &
                     p3_sky(centre+dis_ra,centre+dis_dec,chan) &
                    +sig_lev*(nu(chan)/obsfreq)**(-2.-alpha)
            end do
         end if
      else if (flag_nd.eq.2) then
         call io_getd('What level for point source','1.0',sig_lev,status)
         p_sky(centre+dis_ra,centre+dis_dec) =&
             p_sky(centre+dis_ra,centre+dis_dec)+sig_lev
      end if

   case('g','G')
      call io_getd('Displacement of gaussian in RA (arcsec)',&
           & '0.0',d_ra,status)
      call io_getd('Displacement of gaussian in Dec (arcsec)',&
           & '0.0',d_dec,status)
      call io_getd('What central level for gaussian',&
           & '1.0',sig_lev,status)
      if (flag_nd.eq.3) then
         call io_getd('Radio spectral index:','0.7',alpha,status)
      end if
      call io_getd('What FWHM for Gaussian (arcsec)','1.0',gau_fwhm,status)

! width sigma of a Gaussian is FWHM/(2sqrt(2ln2))
      sigma = gau_fwhm/(2.0*sqrt(2.0*log(2.0)))

! tabulate Gaussian - exp((x/sigma)^2)/2
      do i = 1, maxsize
         rad = (i-1)*cellsize
         rad = ((rad/sigma)**2)/2.0
         gau_tab(i) = exp(-rad)
      end do

! put on sky
      do i = 1, maxsize
         do j = 1, maxsize

! Calculate offset in ra and dec of this point in the map
            x1 = d_ra/cellsize+(i-centre)
            y1 = d_dec/cellsize+(j-centre)

! Add 1 to radius because gau_tab(1) is value of Gaussian at zero radius
            rad = sqrt(x1**2+y1**2)+1.0
            call interp(gau_tab,maxsize,rad,gau_val)
            if (flag_nd.eq.3) then
               do chan = 1, nchan
                 p3_sky(i,j,chan) = p3_sky(i,j,chan)+sig_lev*gau_val*&
                                    (nu(chan)/obsfreq)**(-2.-alpha) 
              end do
            else if (flag_nd.eq.2) then
               p_sky(i,j) = p_sky(i,j)+sig_lev*gau_val
            end if
         end do
      end do

   case('c','C')
      call io_getd('Displacement of circle in RA (arcsec)',&
           & '0.0',d_ra,status)
      call io_getd('Displacement of circle in Dec (arcsec)',&
           & '0.0',d_dec,status)
      call io_getd('What level for circle',&
           & '1.0',sig_lev,status)
      if (flag_nd.eq.3) then
         call io_getd('Radio spectral index:','0.7',alpha,status)
      end if
      call io_getd('What radius for circle (arcsec)','1.0',c_rad,status)
      
! put on sky
      do i = 1, maxsize
         do j = 1, maxsize

! Calculate offset in ra and dec of this point in the map
            x1 = d_ra/cellsize+(i-centre)
            y1 = d_dec/cellsize+(j-centre)
            rad = sqrt(x1**2+y1**2)
            if (rad.le.c_rad) then
               if (flag_nd.eq.3) then
                  do chan = 1, nchan
                    p3_sky(i,j,chan) = p3_sky(i,j,chan)+sig_lev*&
                                       (nu(chan)/obsfreq)**(-2.-alpha) 
                  end do
               else if (flag_nd.eq.2) then
                  p_sky(i,j) = p_sky(i,j)+sig_lev
               end if
            end if
         end do
      end do

   case('r','R')
      call io_getd('Inner radius of ring (arcsec)',&
           & '0.0',ring_rad,status)
      call io_getd('Width of ring (arcsec)',&
           & '1.0', ring_wid, status)
      call io_getd('Displacement of ring in RA (arcsec)',&
           & '0.0',d_ra,status)
      call io_getd('Displacement of ring in Dec (arcsec)',&
           & '0.0',d_dec,status)
      call io_getd('Total flux density in ring',&
           & '1.0',sig_lev,status)
      if (flag_nd.eq.3) then
         call io_getd('Radio spectral index:','0.7',alpha,status)
      end if

! Calculate flux per pixel
      sig_lev = sig_lev/pi/(2*ring_rad*ring_wid+ring_wid**2)*cellsize**2
! put on sky
      do i = 1, maxsize
        do j = 1, maxsize
! Calculate offset in ra and dec of this point in the map
          x1 = d_ra/cellsize+(i-centre)
          y1 = d_dec/cellsize+(j-centre)
          rad = sqrt(x1**2+y1**2)*cellsize
          if ((rad.gt.ring_rad).and.(rad.lt.ring_rad+ring_wid)) then
            if (flag_nd.eq.3) then
              do chan = 1, nchan
! Convert to brightness temp
                norm = (cellsize*sec2rad)**2*2.0*k_b*nu(chan)**2/const_c2*1.0D26
                p3_sky(i,j,chan) = p3_sky(i,j,chan)+sig_lev*&
                                  (nu(chan)/obsfreq)**(-alpha)/norm
              end do
            else if (flag_nd.eq.2) then
              p_sky(i,j) = p_sky(i,j)+sig_lev
            end if
          endif
        enddo
      enddo

   case default
      call io_getd('What DC level','1.0',sig_lev,status)
      if (flag_nd.eq.2) then
         p_sky = p_sky+sig_lev
      else if (flag_nd.eq.3) then
         p3_sky = p3_sky+sig_lev
      end if
   end select what_sky

   if (do_plot) then
      write(*,*) 'Displaying sky'
      if (flag_nd.eq.2) then
         call display_map(p_sky,maxsize,graph_adj,&
                          cellsize,cellsize,0D0,0D0)
      else if (flag_nd.eq.3) then
         chan = nchan/2+1
         call display_map(p3_sky(:,:,chan),maxsize,graph_adj,&
                          cellsize,cellsize,0D0,0D0)
      end if
   end if

end subroutine make_test_sky




