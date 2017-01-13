subroutine monte_carlo_primordial_effect

! relies on user having inputted geometry and cluster parameters
! (g-g, g-c-p)

  use kind_def
  use sz_globals
  use maths
  use kgraph
  use cosmology
  use file_handling
  ! use map_handling
  use physics

  ! ma_p_size = size of fits file that is read in
  ! ma_x_size = size of the FoV of the telescope (ie sky arrays)

  implicit none
  integer, parameter :: resol=220
  integer, parameter :: npoints = 50
  integer :: i,j, mapsize
! real (kind=dp), dimension(maxsize,maxsize) :: re_map, im_map
  real (kind=dp) :: cell
  real (kind=dp) :: incell
  real (kind=dp) :: upt,vpt,psi,inc,rad_inc
  real (kind=dp), dimension(npoints) :: re_pt, im_pt, amp
  real (kind=dp), dimension(resol) :: max_amp, mean_amp
  real (kind=dp), dimension(resol) :: mean_cluster_amp 
  real (kind=dp), dimension(resol) :: mean_cluster_cmb_amp
  real (kind=dp), dimension(resol) :: max_cluster_amp, rad, plotable
  real (kind=dp), dimension(resol) :: noise, decr
  real (kind=dp), dimension(:,:), allocatable :: fits_map
  real (kind=dp) :: rzerofl,izerofl, norm
  character :: char1
  character(len=fname_len)  :: outfile, fits_file
  logical :: stat_fl, sky_made, aperture_made, map_exists

  telescope : select case(which_tel)
  case('c','C','m','M')
     rad_inc = 150./resol
  case('e','E')
     rad_inc = 300./resol
  case('d','D')
     rad_inc = 5.*dish_diameter/resol
  case default
     rad_inc = 2750./resol
  end select telescope

  call angdist
  D_theta = D_theta * Mpc
  from_file=.false.
  from_poly=.false.
  do_plot = .false.
  overwrite = .true.

  write(*,*) 'Assumes cluster parameters have been entered'
  write(*,*)
  write(*,*) 'Enter CMB FITS filename, and details'
  write(*,*)

  ! Prompt user for details of CMB map

  status = 0
  fits_file = ''
  getfileloop:  do 
     
     call io_getfil('FITS filename>','/home/kl247/data/true.fits'&
          &,fits_file,status)
     write(*,*) fits_file
     inquire(file=fits_file,exist=map_exists)
     if (map_exists) then 
        exit getfileloop
     else
        write(*,*) 
        write(*,*) 'File doesn''t exist'
     end if
  end do getfileloop
  call io_getd('Input pixel size (arcsec):','225.0',incell,status)
  call io_geti('Input map size (pixels):','512',mapsize,status)
  allocate (fits_map(mapsize,mapsize))

  ! Make cluster skies

  call make_skies
  write(*,*) 'Skies for cluster observation made'
  write(*,*)

  ! Make Aperture

  do_plot = .false.           
  call make_aperture
  write(*,*) 'Aperture for cluster observation made'
  write(*,*)

  ! Calculate Max/Mean Amplitudes

  inc = pi/npoints
  do j = 1, resol
     rad(j) = j*rad_inc
     psi = 0.
     do i = 1, npoints
        upt = rad(j)*sin(psi)
        vpt = rad(j)*cos(psi)
        call extract_visibility&
             & (re,im,maxsize,upt,vpt,cellsize,re_pt(i),im_pt(i))
        psi = psi+inc
     end do
     amp = sqrt(re_pt**2+im_pt**2)
     max_cluster_amp(j) = maxval(amp)
     mean_cluster_amp(j) = sum(amp)/npoints
  end do
  write(*,*) 'Max/Mean visibility amplitudes calculated'
  write(*,*)
  
  ! Read in CMB fits map

! HACK !!!!
  call read_fits_map(fits_file)
  write(*,*) 'Assumes that FITS map will be in brightness &
       & temperature (K)'
  norm = cellsize*sec2rad
  norm = norm*norm
  norm = norm*2.0*k_b*obsfreq**2/const_c2*1.0D26
  write(*,*) 'Amplitude scaling factor if in flux density (Jy) ', &
       1.0/norm
  norm = -g_nu_x2(obsfreq)*T0
  write(*,*) 'Amplitude scaling factor if in y-parmeter', norm
  write(*,*) 'incell=',incell, 'cellsize=',cellsize, 'mapsize=',mapsize

! loop starts here!

  call katy_subim&
  &(fits_map, szsky, incell, cellsize, mapsize, maxsize)
  sky_made = .true.
  from_file = .false.
  
  call angdist
  D_theta = D_theta * Mpc
  from_file=.false.
  from_poly=.false.
  do_plot = .false.
  overwrite = .false.

  ! 'Make-skies' (don't overwrite previous sky)

  call make_skies
  write(*,*) 'Cluster + CMB Skies made'
  write(*,*)

  !  Make Aperture

  do_plot = .false.           
  call make_aperture
  aperture_made = .true.
  write(*,*) 'Cluster + CMB Aperture Made'
  write(*,*)

  !  Calculate Max/Mean Amplitude

  inc = pi/npoints
  do j = 1, resol
     rad(j) = j*rad_inc
     psi = 0.
     do i = 1, npoints
        upt = rad(j)*sin(psi)
        vpt = rad(j)*cos(psi)
        call extract_visibility&
             & (re,im,maxsize,upt,vpt,cellsize,re_pt(i),im_pt(i))
        psi = psi+inc
     end do
     amp = sqrt(re_pt**2+im_pt**2)
     max_amp(j) = maxval(amp)
     mean_amp(j) = sum(amp)/npoints
  end do

! loop ends here!

! now workout the average and st.dev. before plotting...

!!$  write(*,*) 'Displaying cluster'
!!$  call pgpage
!!$  call display_profile(rad,mean_cluster_amp,resol,.true.)
!!$  call pglabel('Observing baseline ( \gl )','Predicted flux (&
!!$       & Jy)','Mean amplitude: cluster only')

  ! Plot signal to noise

  do j = 1, resol
     plotable(j) = 1./ (-mean_cluster_amp(j)/(mean_amp(j)-mean_cluster_amp(j)))
     noise(j) = (mean_amp(j))
     decr(j) = mean_cluster_amp(j)
  end do
  write(*,*) 'Displaying cluster + noise'
  call pgpage
  call display_profile(rad,mean_amp,resol,.true.)
  call pglabel('Observing baseline ( \gl )','Predicted flux (&
        &Jy)','Mean amplitude: cluster + CMB')

  call pgpage
  call display_profile(rad,plotable,resol,.true.)
  call pglabel('Observing baseline ( \gl )','ratio of fluxes'&
       &,'Mean amplitude: CMB / cluster')


!!$
!!$  write(*,*) 'Displaying signal to noise'
!!$  call pgpage
!!$  call display_profile(rad,plotable,resol,.true.)
!!$
!!$  write(*,*) 'Displaying cluster decrement'
!!$  call pgpage
!!$  call display_profile(rad,mean_cluster_amp,resol,.true.)
!!$  call pglabel('Observing baseline ( \gl )','Predicted flux (&
!!$        &Jy)','Mean amplitude')

  deallocate (fits_map)


end subroutine monte_carlo_primordial_effect









