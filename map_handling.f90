module map_handling

   use kind_def
   use sz_globals

   implicit none

   public :: adaptive_smooth 
   public :: diff_map
   public :: def_map_region
   public :: estimate_map_background
   public :: remsource
   public :: subim
   public :: mod_map_region
   public :: find_centroid
   public :: smoo_map
   public :: simple_clean
   public :: select_map
   public :: stat_map
   interface make_spectral_index_map
      module procedure make_spectral_index_map_2inp
      module procedure make_spectral_index_map_n_inp
   end interface
   public :: make_map_profile

contains

! **************************************************************************

! Smooth a map with a gaussian filter whose sigma varies as 1/(pixel value) 
! from mingauss pixels at the  brightest pixel to a maximum of maxgauss

   subroutine adaptive_smooth(inmap,outmap,mapsize)

      use kind_def
      use sz_globals

      implicit none

      integer, parameter :: max_conv = 30
      integer, intent(in) :: mapsize
      real(kind=dp), dimension(mapsize, mapsize), intent(in) :: inmap
      real(kind=dp), dimension(mapsize, mapsize), intent(out) :: outmap
      real(kind=dp) :: max, sigma, norm
      real(kind=dp), dimension(-max_conv:max_conv,-max_conv:max_conv) :: conv
      integer :: ncell
      integer :: npix 
      integer :: i, j, ii, jj, iii, jjj, box

      npix = mapsize**2
      ncell = mapsize

      max = maxval(inmap)
      outmap = 0.0

      do i = 1, mapsize
         do j = 1, mapsize
            if (inmap(j,i).gt.0) then
               sigma = (max/inmap(j,i)) * mingauss
               if (sigma.gt.maxgauss) sigma = maxgauss
               box = 3*sigma
               if (box.gt.max_conv) then
                  box = max_conv
               end if
               sigma = 1/(2*sigma**2)
               norm = 0
               do ii = -box,box
                  do jj = -box,box
                     iii = i+ii
                     jjj = j+jj
                     if ((iii.gt.0).and.(iii.le.mapsize).and.&
                         (jjj.gt.0).and.(jjj.le.mapsize)) then
                        conv(jj,ii) = exp(-(ii**2+jj**2)*sigma)
                        norm = norm+conv(jj,ii)
                     end if
                  end do
               end do
               norm = inmap(j,i)/norm
               do ii = -box, box
                  do jj = -box, box
                     iii = i+ii
                     jjj = j+jj
                     if ((iii.gt.0).and.(iii.le.mapsize).and.&
                         (jjj.gt.0).and.(jjj.le.mapsize)) then
                        outmap(jjj,iii)=outmap(jjj,iii)+conv(jj,ii)*norm
                     end if
                  end do
               end do
            end if
         end do
      end do
      norm = sum(outmap)
      write(*,*) 'total flux = ', norm

   end subroutine adaptive_smooth

! ***************************************************************************

   subroutine diff_map(map1,map2,chisq,likely,report_stat)
!     puts map1 - map2 in diff
!     calculates chisq = (map1-map2)**2/map)
!     therefore map1 should be the observed results and map2 the model results

      use kind_def
      use maths
      use kgraph

      implicit none
      real(kind=dp), dimension(maxsize,maxsize), intent(in) :: map1, map2
      logical, intent(in) :: report_stat
      real(kind=dp), dimension(maxsize,maxsize) :: ln_prob
      real(kind=dp), intent(inout) :: chisq, likely
      real(kind=dp), dimension(maxsize,maxsize) :: chi_work
      logical, dimension(maxsize,maxsize) :: use_region
      integer :: i,j
      real(kind=dp) :: size, dummy1,dummy2

! Check that map2 is non-zero; remove from region if not
      use_region = ((map2 > 0.0).and.ok_data)

! Calculate chi-squared
      where(use_region)
         diff = map1-map2
         chi_work = diff**2/map2
      elsewhere
         chi_work = 0.0
      end where

      chisq = sum(chi_work)
      size = count(use_region)

      if (size.ne.0) then
         chisq = chisq/size
      else
         chisq = 0.0
      end if

! Calculate poisson likelihood
! map1 should be the observed results and map2 the model results 
      do i = 1,maxsize
         do j = 1,maxsize
            if (use_region(i,j)) then 
               dummy1 = (1.0*nint(map1(i,j)))
               ln_prob(i,j) = map1(i,j)*log(map2(i,j))-1.0*map2(i,j)-&
                              log_factorial(dummy1)

            else
               ln_prob(i,j) = 0.0
            endif
         end do
      end do

      if (size.ne.0) then
         likely = -sum(ln_prob)/size
      else
         likely = 0.0
      end if

      if (plot_open) then
         write(*,*) 'Showing ',count(ok_data),' pixels'
         call display_map(diff,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
      end if

      if (report_stat) then
         write(*,*) 'Max in difference map is    :',maxval(diff)
         write(*,*) 'Minimum in difference map is:',minval(diff)
         write(*,*) 'Total flux in difference map:',sum(diff)
         write(*,*) 'Chi-squared:',chisq
         write(*,*) 'NB Chi-squared is not the right statistic to use&
                       & - handle with care'
         write(*,*) 'Normalised likelihood', likely
      end if

   end subroutine diff_map

! *************************************************************************

   subroutine def_map_region

      use kind_def
      use sz_globals
      use kgraph

      implicit none
      integer :: imin, imax, jmin, jmax, num_pix, ximin, ximax, xjmin, xjmax
      real, dimension(4) :: map_rectangle
      real, dimension(2) :: circle_cent
      real :: circle_radius
      real :: new_rad_squared, true_rad_squared
      integer :: i, j
      character*1 :: chr1

      real,dimension(maxsize) :: x,y

      ok_data = .false.
      status = 0
      inc_region : do
         call io_getc('Define (c)ircular or (r)ectangular region?',&
              'c',chr1,status)
         if (chr1 == 'c') then
            call get_map_circle(circle_radius,circle_cent,maxsize,&
                 cellsize,cellsize,0D0,0D0)
            do i = 1, maxsize 
               do j = 1, maxsize
                  imin = i   ! just because they're already
                  imax = j   ! declared above
                  jmin = int((i - circle_cent(1))**2)
                  jmax = int((j - circle_cent(2))**2)
                  new_rad_squared = jmin + jmax
                  true_rad_squared = circle_radius**2
                  if (new_rad_squared .le. true_rad_squared) then
                     ok_data(imin,imax) = .true.
                  end if
               end do
            end do
         else
            write(*,*) 'Draw region to include'
            call get_map_rectangle(map_rectangle,maxsize,&
                                cellsize,cellsize,0D0,0D0)
            imin = int(map_rectangle(1))
            jmin = int(map_rectangle(2))
            imax = int(map_rectangle(3))
            jmax = int(map_rectangle(4))
            ok_data(imin:imax,jmin:jmax) = .true.        
         end if

         call io_getc('Another region to include (y/n):','n',chr1,status)
         if (chr1=='n') then
            exit inc_region
         end if
      end do inc_region

      exc_region : do
         call io_getc('Draw region to exclude (y/n):','n',chr1,status)
         if (chr1=='n') then
            exit exc_region
         end if
         write(*,*) 'Draw region to exclude'
         call get_map_rectangle(map_rectangle,maxsize,&
                                cellsize,cellsize,0D0,0D0)
         ximin = int(map_rectangle(1))
         xjmin = int(map_rectangle(2))
         ximax = int(map_rectangle(3))
         xjmax = int(map_rectangle(4))
         ok_data(ximin:ximax,xjmin:xjmax) = .false.
      end do exc_region

   end subroutine def_map_region

! ***************************************************************************

   subroutine estimate_map_background(map1)

      use kind_def
      use sz_globals
      use kgraph

      implicit none

      real(kind=dp), dimension(maxsize,maxsize), intent(in) :: map1
      integer :: num_pix
      real(kind=dp) :: tot_flux
      logical :: old_graph_adj

! Display rosat map and allow user to adjust the greyscale level
      old_graph_adj = graph_adj
      graph_adj = .true.
      write(*,*) 'Adjust grey-scale level and then press q'
      call display_map(map1,maxsize,graph_adj,cellsize,cellsize,0.D0,0.D0)

! Determine which area of the map to use to estimate background
      call def_map_region

      tot_flux = sum(map1,mask=ok_data)
      num_pix = count(ok_data)
      back = tot_flux/num_pix

      write(*,*) 'Using a total of ',num_pix, ' pixels'

      write(*,*) 'Background level in this region is ', back
      
      graph_adj = old_graph_adj

   end subroutine estimate_map_background

! **************************************************************************

!   this subroutine replaces the point sources with blank patches of 
!   sky of level back

   subroutine remsource(inmap,flag,insize)

      use kind_def
      use sz_globals
      use file_handling
      use kgraph

      implicit none

      real(kind=dp), dimension(insize,insize), intent(inout) :: inmap
      logical, dimension(insize,insize), intent(inout) :: flag
      integer, intent(in) :: insize
      real, dimension(4) :: map_rectangle
      integer :: i, j, field, ximin, ximax, xjmin, xjmax
      character*1 :: chr1

      call display_map(inmap,insize,graph_adj,cellsize,cellsize,0D0,0D0)
      exc_region : do
         write(*,*) 'Draw source to exclude'
         call get_map_rectangle(map_rectangle,maxsize,&
              cellsize,cellsize,0D0,0D0)
         ximin = int(map_rectangle(1))
         xjmin = int(map_rectangle(2))
         ximax = int(map_rectangle(3))
         xjmax = int(map_rectangle(4))
         flag(ximin:ximax,xjmin:xjmax) = .false.
         call io_getc('Another region to exclude (y/n):','n',chr1,status)
         if (chr1=='n') then
            exit exc_region
         end if
      end do exc_region
      where(.not.(flag))
         flag = .false.
      end where

   end subroutine remsource

! *************************************************************************

   subroutine subim(inmap, outmap, incell, outcell, nx, ny, outsize,&
                    scale, surface_brightness)

! cuts a nx x ny map down to outsize x outsize, interpolates onto
! a outcell grid

! 12/03/02 change to condition for FITS map being too small
! 18/06/02 changes to allow processing of surface brightness map (no change 
!          with change of pixel size, eg szsky) as well as flux maps (change
!          proportional to ratio of new to old pixel area) [KG]
! 20/06/02 changed to allow any combination of insize, outsize, incell, 
!          outcell [KG]
!          Centroid fitting commented out [KG]
! 02/04/07 changed to same coordinate convention as AIPS [AMS]
     use kind_def
     use sz_globals
     use kgraph

     implicit none

     integer, intent(in) :: nx, ny, outsize
     real(kind=dp), dimension(nx,ny), intent(inout) :: inmap
     real(kind=dp), dimension(outsize,outsize), intent(out) :: outmap
     real(kind=dp), dimension(outsize,outsize) :: temp_map
     real(kind=dp), intent(inout) :: scale
     logical, intent(in) :: surface_brightness
     real(kind=dp) :: tot
     real(kind=dp), intent(in) :: outcell, incell

     integer :: i,j,k,l, p, q, outcentre, box, p2, q2
     integer :: imin, imax, jmin, jmax, size, ncontrib
     integer :: ximin, ximax, xjmin, xjmax, adjust_x, adjust_y
     real :: pcentre,qcentre,pc,qc
     real, dimension(4) :: map_region
     real(kind=dp) :: pr, qr, t, u0, t1, u1, x, y, pr2, qr2, p_com, q_com
     real(kind=dp) :: temp, psum, qsum, flux_contrib, common_area
     real(kind=dp) :: tot_area, x2, y2

     character*1 :: chr1, chr2


     outmap = 0.d0

!     min = (insize-outsize)/2
!     max = min+outsize

     if (.not.surface_brightness) then
        scale = scale*(outcell/incell)**2
     end if

! copy maps to new array
     
! check whether we have to do integration or interpolation 
     if (outcell==incell) then
        outmap = 0.
        adjust_x = ((outsize/2)+1)-((nx/2)+1)
        adjust_y = ((outsize/2)+1)-((ny/2)+1)
        do i = 1,outsize
           do j = 1,outsize
              if (((i-adjust_x).ge.1).and.((i-adjust_x).le.nx).and.&
                  ((j-adjust_y).ge.1).and.((j-adjust_y).le.ny)) then
                 outmap(i,j) = inmap(i-adjust_x,j-adjust_y)     
              end if
           end do
        end do
        
        outmap = outmap*scale
        
     elseif (outcell.lt.incell) then

! interpolate onto finer grid
        outcentre = (outsize/2)+1
        pcentre = (nx/2)+1
        qcentre = (ny/2)+1

! loop over all pixels in output map
        do i = 1,outsize
           do j = 1,outsize

! calculate where this corresponds to on input map
              x = (i-outcentre)*outcell
              y = (j-outcentre)*outcell
              pr = x/incell+pcentre
              qr = y/incell+qcentre
              p = int(pr)
              q = int(qr)
              t = pr-p
              u0 = qr-q
              t1 = 1-t
              u1 = 1-u0

! Check whether the desired region exists on input map
              if ((p.ge.1).and.(q.ge.1).and.&
                  (p.le.nx+1).and.(q.le.ny+1)) then
                 outmap(i,j) = t1*u1*inmap(p,q)&
                              +t*u1*inmap(p+1,q)&
                              +t*u0*inmap(p+1,q+1)&
                              +t1*u0*inmap(p,q+1)

! Apply user defined rescaling
                 outmap(i,j) = outmap(i,j)*scale

              else

! If region doesn't exits on input map
                 outmap(i,j) = 0.0
              end if
           end do
        end do
     else

! integrate onto courser grid
        outcentre = (outsize/2)+1
        pcentre = (nx/2)+1
        qcentre = (ny/2)+1
        ncontrib = int(outcell/incell)+2

! loop over all pixels in output map
        do i = 1,outsize
           do j = 1,outsize

! calculate where this corresponds to on input map
              x = (i - outcentre)*outcell
              y = (j - outcentre)*outcell
              x2 = (i+1 - outcentre)*outcell
              y2 = (j+1 - outcentre)*outcell
              pr = x/incell + pcentre
              qr = y/incell + qcentre
              p = int(pr)
              q = int(qr)
              pr2 = (x2)/incell + pcentre
              qr2 = (y2)/incell + qcentre
              p2 = int(pr2)
              q2 = int(qr2)

! calculate area in arcsec^2 of this input pixel which is in common with
! the output pixel
              flux_contrib = 0.0
              tot_area = 0.0

! Check whether this point on the output map has any points contributing to it
              if ((p2.ge.1).and.(p.le.nx).and.&
                  (q2.ge.1).and.(q.le.ny)) then
                 
! loop over all pixels which could contribute to this output pixel
                 do k = 0, ncontrib-1

! determine common p length

! the first pixel will contribute a fractional amount
                    if (k.eq.0) then
                       p_com = (p+1-pr)*incell

! if the end of the pixel is beyond the end of the output pixel it won't 
! contribute a full amount
                    else if ((p+k+1).gt.pr2) then

! need to check whether this pixel contributes at all
                       if (p+k.lt.pr2) then
                          p_com = (pr2-(p+k))*incell
                       else
                          p_com = 0.0
                       end if

! otherwise this pixel contributes fully
                    else
                       p_com = incell
                    end if

                    do l = 0, ncontrib-1

! determine common q length

! the first pixel will contribute a fractional amount
                       if (l.eq.0) then
                          q_com = (q+1-qr)*incell

! if the end of the pixel is beyond the end of the output pixel it won't 
! contribute a full amount
                       else if ((q+l+1).gt.qr2) then

! need to check whether this pixel contributes at all
                          if (q+l.lt.qr2) then
                             q_com = (qr2-(q+l))*incell
                          else
                             q_com = 0.0
                          end if

! otherwise this pixel contributes fully
                       else
                          q_com = incell
                       end if

                       common_area = p_com*q_com

! check that this pixel exists in inmap
                       if ((p.ge.1).and.(q.ge.1).and.&
                            (p.le.nx).and.(q.le.ny)) then
                          flux_contrib = flux_contrib+common_area*inmap(p,q)
                          tot_area = tot_area+common_area
                       end if
                    end do
                 end do
                 outmap(i,j) = flux_contrib/(outcell**2)*scale
              end if
           end do
        end do
     end if
   end subroutine subim

! ***************************************************************************

   subroutine mod_map_region

      use kind_def
      use sz_globals
      use kgraph

      implicit none
      integer :: imin, imax, jmin, jmax, num_pix, ximin, ximax, xjmin, xjmax
      real, dimension(4) :: map_rectangle
      character*1 :: chr1

      status = 0
      inc_region : do
         write(*,*) 'Draw region to include'
         call get_map_rectangle(map_rectangle,maxsize,&
                                cellsize,cellsize,0D0,0D0)
         imin = int(map_rectangle(1))
         jmin = int(map_rectangle(2))
         imax = int(map_rectangle(3))
         jmax = int(map_rectangle(4))
         ok_data(imin:imax,jmin:jmax) = .true.
         call io_getc('Another region to include (y/n):','n',chr1,status)
         if (chr1=='n') then
            exit inc_region
         end if
      end do inc_region

      exc_region : do
         write(*,*) 'Draw region to exclude'
         call get_map_rectangle(map_rectangle,maxsize,&
                                cellsize,cellsize,0D0,0D0)
         ximin = int(map_rectangle(1))
         xjmin = int(map_rectangle(2))
         ximax = int(map_rectangle(3))
         xjmax = int(map_rectangle(4))
         ok_data(ximin:ximax,xjmin:xjmax) = .false.
         call io_getc('Another region to exclude (y/n):','n',chr1,status)
         if (chr1=='n') then
            exit exc_region
         end if
      end do exc_region

   end subroutine mod_map_region

! ***************************************************************************

   subroutine find_centroid(inmap,incell,insize,pcentre,qcentre)
      use kind_def
      use sz_globals
      use kgraph

      implicit none

      integer, intent(in) :: insize
      real(kind=dp), dimension(insize,insize), intent(inout) :: inmap
      real(kind=dp), intent(in) :: incell
      integer :: i, j, incentre
      real(kind=dp) :: xsum, psum, qsum, pcentre, qcentre, temp

      
      write(*,*) 'Determining map centroid'
      write(*,*) 'Adjust grey scale of map (press q when finished)'
      if (graph_adj) then 
         call display_map(inmap,insize,graph_adj,incell,incell,0.D0,0.D0)
      end if

! Determine which area of the map to use to estimate centroid
      write(*,*) ' '
      write(*,*) 'Draw region to use for centroid fitting'
      write(*,*) ' '
      call def_map_region

      xsum = 0.0
      psum = 0.0
      qsum = 0.0
      incentre = (insize/2)+1

! Centroid fitting.
      do j = 1,insize
         do i = 1,insize
            if (ok_data(i,j)) then 
               temp = inmap(i,j)
               xsum = xsum + temp
               psum = psum+(i*temp)
               qsum = qsum+(j*temp)
            end if
         end do
      end do

      write(*,*) 'centroid = ',psum/xsum,qsum/xsum 
      pcentre = ((psum/xsum) - incentre)*incell 
      qcentre = ((qsum/xsum) - incentre)*incell 
      write(*,*) 'Centroid = ',pcentre,qcentre, ' arcsec from centre' 
      write(*,*) 'Total flux used for fitting ',xsum


   end subroutine find_centroid

! ***************************************************************************
   subroutine smoo_map(inmap,mapsize,smoo_scale)

      use kind_def
      use sz_globals
      use kgraph
      use maths

      implicit none
 
      integer, intent(in) :: mapsize
      real(kind=dp), dimension(mapsize,mapsize), intent(inout) :: inmap
      real(kind=dp), intent(in) :: smoo_scale
      real(kind=dp), dimension(mapsize,mapsize) :: gaumap
      real(kind=dp), dimension(mapsize,mapsize) :: tempre, tempim
      real(kind=dp) :: norm2, dist2, synth_beam, norm, xsum
      integer :: i, j, cen_pix

      tempre = inmap
      tempim = 0.0
      
      if (do_plot) then
         write(*,*) 'Input sky'
         call display_map(tempre, mapsize, graph_adj,cellsize,cellsize,0D0,0D0)
      end if
      
!     FT it to get aperture
      call do_2d_fft(tempre,tempim,mapsize,mapsize,1)

      if (do_plot) then
         write(*,*) 'Displaying real part of aperture'
         call display_map(tempre,mapsize,graph_adj,cellsize,cellsize,0D0,0D0)
      end if

      if (do_plot) then
         write(*,*) 'Displaying imaginary part of aperture'
         call display_map(tempim,mapsize,graph_adj,cellsize,cellsize,0D0,0D0)
      end if

! Multiply aperture plane by Gaussian
      cen_pix = mapsize/2+1
      norm2 = (mapsize*cellsize/smoo_scale)**2
      do i = 1, mapsize
         do j = 1, mapsize
            dist2 = (i-cen_pix)**2+(j-cen_pix)**2
            gaumap(i,j) = exp(-dist2/(2*norm2))
         end do
      end do

      tempre = tempre*gaumap
      tempim = tempim*gaumap

      if (do_plot) then
         write(*,*) 'Displaying real part of filtered aperture'
         call display_map(tempre,mapsize,graph_adj,cellsize,cellsize,0D0,0D0)
      end if

      if (do_plot) then
         write(*,*) 'Displaying imaginary part filtered of aperture'
         call display_map(tempim,mapsize,graph_adj,cellsize,cellsize,0D0,0D0)
      end if

! FT back to map plane
      call do_2d_fft(tempre,tempim,mapsize,mapsize,-1)

      if (do_plot) then
         write(*,*) 'Displaying real part of smoothed map'
         call display_map(tempre,mapsize,graph_adj,cellsize,cellsize,0D0,0D0)
      end if

      if (do_plot) then
         write(*,*) 'Displaying imaginary part of smoothed map'
         call display_map(tempim,mapsize,graph_adj,cellsize,cellsize,0D0,0D0)
      end if

! Overwrite input map with smoothed map
      inmap = tempre

! Rescale map to preserve flux density
      inmap = inmap/mapsize**2

   end subroutine smoo_map

! ***************************************************************************
   
   subroutine simple_clean
      
      use kind_def
      use sz_globals
      use kgraph

      implicit none

      integer  :: i,j,k,l,m,n,niter,xmax
      real(dp) :: rad,fwhm_maj,fwhm_min,theta_maj, theta_min
      real(dp) :: sigma_maj,sigma_min,temp,radmax,gain,half_power
      real(dp) :: cl_beam(maxsize,maxsize),temparray(maxsize,maxsize)
      integer,allocatable  :: source(:,:)
      real, allocatable :: amp(:)
      integer  :: disp_rate
      real     :: x(maxsize),bmin,bmax,cc_x,cc_y
      logical  :: show_progress
      logical, external :: io_yesno
      
      if (io_yesno('Restart CLEANing from scratch','no',status)) then
         resid_map = dirty_map
         clean_map = 0.d0
         comp_map = 0.d0
      end if

      gain = 1.d-1
      call io_geti('Number of iterations:','100',niter,status)
      call io_getd('Gain','*',gain,status)

      ok_data = .true.
      if (io_yesno('Define CLEAN region','yes',status)) then
         call display_map(resid_map,maxsize,.false.,&
                          cellsize,cellsize,0D0,0D0)      
         call def_map_region
      end if

      show_progress = io_yesno('Show progress during CLEAN','no',status)
      if (show_progress) then
         disp_rate = 20
         call io_geti('Number of iterations between plots','*',&
                       disp_rate,status)
      end if

! determine Gaussian approximation to beam
      write(*,*) 'Estimating restoring beam'
      radmax = 0.
      half_power = maxval(real(beam))/2.
      do i = 1,maxsize
         do j = 1,maxsize
            if (beam(i,j).ge.half_power) then
               rad = sqrt(real(i-(maxsize/2)-1)**2+real(j-(maxsize/2)-1)**2)
               if (rad.ge.radmax) radmax = rad
            end if
         end do
      end do
      write(*,*) 'CLEAN BEAM: MAJ = ',radmax,' MIN = ',radmax,' PA = ',0.
      sigma_maj = radmax/(sqrt(2.*log(2.)))

! make clean beam:
      cl_beam = 0.
      do i = 1,maxsize
         do j = 1,maxsize
            temp = sqrt(real(i-maxsize/2)**2+real(j-maxsize/2)**2)
            cl_beam(i,j) = exp(-1.*real(temp)**2/(2.*(real(sigma_maj))**2))  
         end do
      end do

      write(*,*) 'Clean beam determined'

! Allocate and initialise arrays
      allocate(source(niter,2))
      allocate(amp(niter))
      source = 0
      amp = 0.d0

      write(*,*) 'Starting Simple CLEAN'

! subtract PSF iteratively:
      do k = 1,niter
         source(k,1:2) = maxloc(resid_map,MASK=ok_data)
         i = source(k,1)
         j = source(k,2)
         amp(k) = resid_map(i,j)*gain
         do l = 1,maxsize
            do m = 1,maxsize
               if (((maxsize/2+1+(l-i)).gt.0)&
                    .and.((maxsize/2+1+(l-i)).le.maxsize)) then
                  if (((maxsize/2+1+(m-j)).gt.0)&
                       .and.((maxsize/2+1+(m-j)).le.maxsize)) then
                     resid_map(l,m) = resid_map(l,m) &
                          - amp(k)*beam(maxsize/2+1+(l-i),maxsize/2+1+(m-j))
                  end if
               end if
            end do
         end do

! Display progress
         if (show_progress.and.(1.*(k/disp_rate).eq.(1.*k)/disp_rate)) then
            write(*,*) 'Iteration: ',k            
            call display_map(resid_map,maxsize,.false.,&
                             cellsize,cellsize,0D0,0D0)
            call pgsci(2)
            call pgsch(0.5)            
            do i = 1,k
               cc_x = (source(i,1)-maxsize/2)*cellsize
               cc_y = (source(i,2)-maxsize/2)*cellsize
               call pgpt1(cc_x,cc_y,2)
            end do
            write(*,*) sum(amp),' Jy CLEANed'
            call pgsci(plot_col)
            call pgsch(char_hgt)
         end if
      end do
      write(*,*)
      write(*,*) '**** CLEAN completed; restoring'
      write(*,*)

! restore sources using CLEAN beam
! this convolution done in a very inefficient way
      do k = 1,niter
         do i = 1,maxsize
            do j = 1,maxsize
               temp = dble(exp(-1.*real(i - source(k,1))**2&
                    /(2.*(sigma_maj)**2)))
               temp = temp*dble(exp(-1.*real(j-source(k,2))**2&
                    /(2.*(sigma_maj)**2)))
               comp_map(i,j) = comp_map(i,j)+temp*amp(k)
            end do
         end do
      end do

      clean_map = comp_map+resid_map
 
      write(*,*)
      write(*,*) 'Restored map calculated'
      write(*,*)

      deallocate(source)
      deallocate(amp)

      write(*,*) 'CLEAN residuals map'
      call display_map(resid_map,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
      write(*,*) 'CLEAN components map'
      call display_map(comp_map,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
      write(*,*) 'CLEAN map'
      call display_map(clean_map,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)

   end subroutine simple_clean

! ***************************************************************************

   subroutine select_map

      use kind_def
      use sz_globals

      implicit none

     write(*,*) 'Available maps:'
      write(*,*)
      write(*,*) '(S)Z sky'
      write(*,*) 'Compton (Y) sky'
      write(*,*) '(E)lectron density sky'
      write(*,*) '(R)eal part of aperture plane'
      write(*,*) '(I)maginary part of aperture plane'
      write(*,*) '(G)ridded real part of aperture plane'
      write(*,*) 'Gridded imaginary part of (A)perture plane'
      write(*,*) '(U)V coverage'
      write(*,*) '(D)irty map'
      write(*,*) 'Synthesised (B)eam'
      write(*,*) '(C)LEANed map'
      write(*,*) 'CLEAN com(P)onents map'
      write(*,*) 'CLEAN residua(L)s map'
      write(*,*) '(X)-ray sky'
      write(*,*) '(M)odel X-ray map'
      write(*,*) '(O)bserved X-ray map'
      write(*,*) 'Di(F)ference map'
      write(*,*)

      call io_getc('Select a map','*',map_type,status)

      flag_nd = 2
      select case (map_type)
      case('s','S')
         call io_geti('Which channel number (0 for all)','1',chan,status)
         if(chan.ne.0) then
            p_sky => szsky(:,:,chan)
         else
            flag_nd = 3
            p3_sky => szsky
         end if
      case('y','Y')
         p_sky => compton_y
      case('e','E')
         p_sky => nsky
      case('r','R')
         call io_geti('Which channel number (0 for all)','1',chan,status)
         if(chan.ne.0) then
            p_sky => re(:,:,chan)
         else
            flag_nd = 3
            p3_sky => re
         end if
      case('i','I')
         call io_geti('Which channel number (0 for all)','1',chan,status)
         if(chan.ne.0) then
            p_sky => im(:,:,chan)
         else
            flag_nd = 3
            p3_sky => im
         end if
      case('g','G')
         p_sky => gr_re        
      case('a','A')
         p_sky => gr_im
      case('u','U')
         p_sky => uvcov
      case('d','D')
         p_sky => dirty_map
      case('b','B')
         p_sky => beam         
      case('c','C')
         p_sky => clean_map
      case('p','P')
         p_sky => comp_map
      case('l','L')
         p_sky => resid_map
      case('x','X')
         p_sky => xsky
      case('m','M')
         p_sky => xmap
      case('o','O')
         p_sky => rosat_map
      case('f','F')
         p_sky => diff
      end select

   end subroutine select_map

! ***************************************************************************

! NB routine uses ok_data so mapsize had better be equal to maxsiz(!)
   subroutine stat_map(inmap,mapsize)

      use kind_def
      use sz_globals
      use kgraph

      implicit none

      integer, intent(in) :: mapsize
      real(kind=dp), dimension(mapsize,mapsize), intent(in) :: inmap
      real(kind=dp), dimension(mapsize,mapsize) ::  in2
      real(kind=dp) :: mean, var
      integer, dimension(2) :: locn

      logical, external :: io_yesno

      ok_data = .true.
      if (io_yesno('Define region for statistics','yes',status)) then
         call display_map(inmap,mapsize,graph_adj,cellsize,cellsize,0D0,0D0)
         call def_map_region
      endif
      write(*,*) 'Number of pixels in region is ',count(mask=ok_data)
      write(*,*) 'Sum over this region is ', sum(inmap,mask=ok_data)
      locn = maxloc(inmap)
      write(*,*) 'Max of the region is ',maxval(inmap),' at ',locn
      locn = minloc(inmap)
      write(*,*) 'Min of the region is ',minval(inmap),' at ',locn
      mean = sum(inmap,mask=ok_data)/count(mask=ok_data)
      in2 = inmap**2
      var = sum(in2,mask=ok_data)/count(mask=ok_data)-mean**2
      write(*,*) 'Mean value in region is ',mean
      write(*,*) 'RMS in region is ',sqrt(var)

   end subroutine stat_map

! ***************************************************************************

   subroutine make_spectral_index_map_2inp(mapsize,map1,map2,fr)

      use kind_def
      use sz_globals
      use kgraph

      implicit none

      integer, intent(in)   :: mapsize
      real(dp), dimension(mapsize,mapsize), intent(in) :: map1, map2
      real(dp), dimension(2), intent(in) :: fr

      real(dp), dimension(mapsize,mapsize) :: index_map
      real                  :: min_ind, max_ind
      integer               :: i, j
      real(dp)              :: mincut, maxcut, blanked, tmp1, tmp2
      logical, external     :: io_yesno

      blanked = 3.0
      mincut = -2.5
      maxcut = 2.5

! Define limits for spectral index map
      call io_getd('Minimum allowed spectral index','*',mincut,status)
      call io_getd('Maximum allowed spectral index','*',maxcut,status)
      call io_getd('Blanking level','*',blanked,status)

      ok_data = .true.
      if (io_yesno("Define region for calculation","no",status)) then
         call display_map(map1,mapsize,graph_adj,&
           cellsize,cellsize,0D0,0D0)
         call def_map_region
      end if

      index_map = blanked

! Only allow regions which are positive
!      ok_data = ok_data.and.(map1.gt.0.0).and.(map2.gt.0.0)

! Calculate spectral index
      do i = 1, mapsize
         do j = 1, mapsize
            if (ok_data(i,j)) then
               tmp1 = -1.d0*(log(abs(map1(i,j)))-log(abs(map2(i,j))))
               tmp2 = (log(fr(1))-log(fr(2)))
               if (tmp2.ne.0.0) then
                  index_map(i,j) = tmp1/tmp2-2.d0
                  if (index_map(i,j).gt.maxcut) index_map(i,j) = maxcut
                  if (index_map(i,j).lt.mincut) index_map(i,j) = mincut
               else
                  index_map(i,j) = blanked
               end if
            else
               index_map(i,j) = blanked               
            end if
         end do
      end do

      call display_map(index_map,mapsize,graph_adj,&
           cellsize,cellsize,0D0,0D0)
      min_ind = minval(index_map)
      max_ind = maxval(index_map)
      call pgwedg('rg',0.5,4.,max_ind,min_ind,&
                  'spectral index')
      call pgsci(2)
      call superpose_cntr(map1,mapsize,cellsize,cellsize,&
           0.D0,0.D0)
      call pgsci(1)

   end subroutine make_spectral_index_map_2inp

! ***************************************************************************

   subroutine make_spectral_index_map_n_inp(mapsize,mchan,freqsky,fr)

      use kind_def
      use sz_globals
      use maths
      use kgraph

      implicit none

      integer, intent(in)   :: mapsize, mchan
      real(dp), dimension(mapsize,mapsize,mchan), intent(in) :: freqsky     
      real, dimension(nchan), intent(in) :: fr
      real(dp), dimension(mapsize,mapsize) :: index_map
      real(dp), dimension(mchan) :: xpts, ypts
      real(dp)              :: m_fit, c_fit
      real                  :: min_ind, max_ind
      integer               :: i, j, k
      real(dp)              :: mincut, maxcut, blanked, tmp1, tmp2
      logical, external     :: io_yesno

      blanked = 3.0
      mincut = -2.5
      maxcut = 2.5

! Define limits for spectral index map
      call io_getd('Minimum allowed spectral index','*',mincut,status)
      call io_getd('Maximum allowed spectral index','*',maxcut,status)
      call io_getd('Blanking level','*',blanked,status)

      ok_data = .true.
      if (io_yesno("Define region for calculation","no",status)) then
         call display_map(freqsky(:,:,1),mapsize,graph_adj,&
           cellsize,cellsize,0D0,0D0)
         call def_map_region
      end if

      index_map = blanked

! Only allow regions for which all channels are positive
!      do i = 1, mchan
!         ok_data = ok_data.and.(freqsky(:,:,i).gt.0.0)
!      end do

! Do one pixel at a time
      do i = 1, mapsize
         do j = 1, mapsize
            if (ok_data(i,j)) then

! Loop over channels finding log flux and log frquency
               do k = 1, mchan
                  xpts(k) = log(abs(fr(k)))
                  ypts(k) = log(abs(freqsky(i,j,k)))
               end do

! Calculate line of best fit
               call best_fit_line(xpts,ypts,mchan,m_fit,c_fit)

               index_map(i,j) = -m_fit-2.d0
               if (index_map(i,j).gt.maxcut) index_map(i,j) = maxcut
               if (index_map(i,j).lt.mincut) index_map(i,j) = mincut
            else
               index_map(i,j) = blanked
            end if
         end do
      end do

      call display_map(index_map,mapsize,graph_adj,&
           cellsize,cellsize,0D0,0D0)
      min_ind = minval(index_map)
      max_ind = maxval(index_map)
      call pgwedg('rg',0.5,4.,max_ind,min_ind,&
                  'spectral index')
      call pgsci(2)
      call superpose_cntr(freqsky(:,:,1),mapsize,cellsize,cellsize,&
           0.D0,0.D0)
      call pgsci(1)

   end subroutine make_spectral_index_map_n_inp

! ***************************************************************************
   
! smooths a 2-d map azimuthally into a 1-d profile
   subroutine make_map_profile(inmap,mapsize,out)

      use kind_def
      use sz_globals
      use maths

      implicit none
      integer, intent(in) :: mapsize
      real(kind=dp), dimension(mapsize,mapsize), intent(in) :: inmap
      real(kind=dp), dimension(mapsize/2), intent(out) :: out
      real(kind=dp) :: psi, inc, x, y, val
      integer :: nstep, centre, i, j, centre1, centre2
      integer, parameter :: res = 200

      centre = mapsize/2+1

      nstep = mapsize/2
      inc = 2*pi/float(res)
      psi = 0
      out = 0.0

! Loop over all radii
      do i = 0, nstep-1

! Loop over all angles at this radius
         do j = 0, res-1
            psi = j*inc
            x = i*sin(psi)+centre    
            y = i*cos(psi)+centre
            call interp(inmap,mapsize,mapsize,x,y,val)
            out(i+1) = out(i+1)+val
         end do
      end do
      out = out/res

   end subroutine make_map_profile

! ***************************************************************************

end module map_handling





