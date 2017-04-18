module sz_globals

   use kind_def

! Constants
   real(kind=dp), parameter :: pi = 3.1415926535D0
   real(kind=dp), parameter :: deg2rad = pi/180.0D0
   real(kind=dp), parameter :: min2rad = deg2rad/60.0D0
   real(kind=dp), parameter :: sec2rad = deg2rad/3600.0D0
   real(kind=dp), parameter :: hr2rad = pi/12.0D0
   real(kind=dp), parameter :: Mpc = 3.086D22
   real(kind=dp), parameter :: m2Mpc = 1.d0/Mpc
   real(kind=dp), parameter :: kpc = 3.086D19
   real(kind=dp), parameter :: T0 = 2.725D0
   real(kind=dp), parameter :: k_b = 1.380658D-23
   real(kind=dp), parameter :: G_c = 6.67259D-11
   real(kind=dp), parameter :: m_p = 1.6726231D-27
   real(kind=dp), parameter :: mu_e = 1.14d0*m_p !average gas mass per electron
   real(kind=dp), parameter :: mu_n = 0.6d0*m_p  !average gas mass per particle
   real(kind=dp), parameter :: sigma_T = 6.652D-29
   real(kind=dp), parameter :: m_e = 9.1093897D-31
   real(kind=dp), parameter :: const_h = 6.6260755D-34
   real(kind=dp), parameter :: const_c = 2.99792458D8
   real(kind=dp), parameter :: const_c2 = const_c**2
   real(kind=dp), parameter :: const_e = 1.60217733D-19
   real(kind=dp), parameter :: const_st2day = 0.1157407407407407D-04 
   real(kind=dp), parameter :: m_0 = 1.989D30 ! solar mass
   real(kind=dp), parameter :: keVtoK = 1.D3*(const_e/k_b)
   real(kind=dp), parameter :: omega_e = 7.292D-5  ! ang vel of the earth
   real(kind=dp), parameter :: I_0 = 2.D0*(k_b*T0)**3/(const_h*const_c)**2

! Telescope specific parameters
   real(kind=dp), parameter :: rt_dish_diameter = 12.8D0
   real(kind=dp), parameter :: rt_freq = 1.55D10
   real(kind=dp), parameter :: cat_dish_diameter = 0.7D0
   real(kind=dp), parameter :: cat_freq = 1.545D10
   real(kind=dp), parameter :: vla_dish_diameter = 25.0D0
   real(kind=dp), parameter :: vsa_freq = 3.3D10
   real(kind=dp), parameter :: vsa_ss_dish_diameter = 3.7D0
   real(kind=dp), parameter :: vsa_ss_fwhm = 6.D2
   real(kind=dp), parameter :: vsa_ma_dish_diameter = 0.143D0
   real(kind=dp), parameter :: vsa_ma_fwhm = 1.485D4
   real(kind=dp), parameter :: vsa_ea_dish_diameter = 0.32D0
   real(kind=dp), parameter :: vsa_ea_fwhm = 6.636D3
   real(kind=dp), parameter :: vsa_sa_dish_diameter = 0.55D0
   real(kind=dp), parameter :: vsa_sa_fwhm = 3.981D3

   integer, parameter       :: convsize = 3
   integer, parameter       :: MAX_FREQ=10

! Global status
   integer :: status

! Global arrays
!
!  re and im are the simulated visibilities in the aperture plane
!  xmap is the simulated xray sky seen by an X-ray telescope
!  fits_map is the original X-ray fits image
!  rosat_map is the corrected rosat image
!  diff is the diffence of rosat_map and xmap
!  szsky and xsky are the modelled skies before observation
!  nsky is the number density map (used for lens models)

! variables for observed visibilities
   integer, parameter :: mvis = 5000000
   real(kind=dp), dimension(mvis) :: u,v,w,jd1,jd2,baseline
   real(kind=dp), dimension(mvis) :: weight,rms,data_re,data_im
   integer, dimension(mvis) :: which_pointing

! variable for mcg fitting
   integer, parameter :: shape(2) = 1.

! 2d arrays for maps, aperture planes etc
   integer, save :: maxsize
   character*1, save :: map_type
   real(kind=dp), dimension(:,:), pointer :: p_sky
   real(kind=dp), dimension(:,:), allocatable, target :: xmap, temp_sky,&
        rosat_map, diff, xsky, nsky, azelbeam, gr_re, gr_im, uvcov
   real(kind=dp), dimension(:,:), allocatable, target :: compton_y, &
        dirty_map, clean_map, beam, comp_map, resid_map, sky_weight

! 3d arrays for maps, aperture planes etc with frequency channels
   real(kind=dp), dimension(:,:,:), pointer :: p3_sky
   real(kind=dp), dimension(:,:,:), allocatable, target :: szsky,re,im,szibeam

! Flag for whether addressing 2d or 3d map   
   integer :: flag_nd

! Logical array for indicating flagged pixels
   logical, dimension(:,:), allocatable, save :: ok_data

! 1d arrays for profiles
   real(kind=dp), dimension(:), allocatable, save :: map_profile, map_profile2

! other arrays...
   real(kind=dp), dimension(:,:), allocatable :: szioverlay
   complex, dimension(:), allocatable, save :: aperture

! variables for defining observation
   real(kind=dp), save :: pb_x, pb_y, back, mingauss, maxgauss
   real(kind=dp), save :: ph_x, ph_y, ha_inc, ha_step, pbspace
   integer, save :: nvis, jpoints, jpoint
   real  :: joff                               

! variables for X-ray telescope
   character*2, save :: x_tel
   logical, save :: x_beam
   real(kind=dp), dimension(:,:), allocatable, save :: ros_beam_re, ros_beam_im
   real(kind=dp), save :: rosat_beam_sum

! flags
   logical, save :: do_plot, overwrite, ellipsoid, verbose, pbcor
   logical, save :: smooth_x_model, graph_adj, autoscale, plot_open, post
   logical, save :: vary_beta, from_file, from_poly, vary_posn
   logical, save :: match_ami, ami, wis, ptlist, wlight, dorot, overlay
   real(kind=dp), save :: prof_x1, prof_x2, prof_y1, prof_y2
   integer,save :: new_centre(2),fparray
   
! globals for the cluster model
   character*1, save :: model_type
   character*1, save :: DM_type !adapted for Einasto DM profile kj 26/02/17

! general cluster parameters
   real(kind=dp), save :: beta, theta_cx, z, n0, cellsize, Te, const
   real(kind=dp), save :: beta2, nx_2, theta_cx2a, gsig, el_alpha
   real(kind=dp), save :: T_central, T_halo, Tcorerad, x_exp, sz_x, sz_y, sz_z
   real(kind=dp), dimension(3), save :: theta_c, ang
   real(kind=dp), save :: a_GNFW, b_GNFW, c_GNFW, c500_GNFW, Ytot, theta_s, f_GNFW
   real(kind=dp), save :: aEin_DM, rm2_DM, rhom2_DM !adapted for Einasto DM profile, needed globally kj 26/02/17
   real(kind=dp), save :: thetai
   real(kind=dp), parameter :: thetalimit = 20.0
   real(kind=dp), parameter :: thetamin = 0.2
   real(dp), allocatable, dimension(:) :: logyintegrand,logtheta,logyarray
   real(kind=dp), save :: MT200, fg200, c200, r200, rp_GNFW, rs_DM
   real(kind=dp), parameter :: rmin=0.001, rmax=100.0
   real(kind=dp), save :: cutoff1, cutoff2
   logical, save :: hanning

! globals for model type File
   integer, save :: num_nt
   real(kind=dp) :: ntstep
   real(kind=dp), dimension(:), allocatable :: ndat, Tdat   

! globals for model type Poly
   real(kind=dp) :: poly_a, poly_b, poly_c, poly_d

   logical, save :: loud_fitting

! variables for cosmology
   real(kind=dp), save :: omegal, omegam, omegab, w_de
   real(kind=dp), save :: D_lum, D_theta, H0, q0, rhocrit0 

! variables for defining a telescope
   integer, save :: nchan, chan
   real(kind=dp), save :: tel_lat, obsfreq, dish_diameter
   real, dimension(:), allocatable, save :: x_pos, y_pos, z_pos, nu, fwhm
   real, dimension(:,:), allocatable, save :: basel_x, basel_y, basel_z
   integer, save :: n_antennas, n_samp, n_pol, n_poly_pb
   real, dimension(:), allocatable, save :: poly_pb
   real, save :: poly_cut
   character*1, save :: which_tel, which_pb
   character*128 :: tel_name
   integer :: pb_points, which_sph_harm
   real(kind=dp), allocatable, dimension(:,:), save :: pb_data
   real(kind=dp) :: pb_spacing, pb_sbeam, pb_rolloff
   real(kind=dp) :: el_sigma1, el_sigma2, el_theta, sph_theta

! arrays for defining point sources
   real(dp), allocatable, dimension(:)  :: rapos, decpos, srcflux,srcalpha
   real(dp), allocatable, dimension(:)  :: radist, decdist, jradist, jdecdist
   integer, allocatable, dimension(:)  :: src_x, src_y
   integer :: nsrc

! parameters associated with FITS output format.
!  FITS header items.
   integer :: bitpix
   integer :: naxis, naxisn(7)
   integer :: gcount, pcount
   real(kind=dp) :: blank
   real(kind=dp) :: bscale, bzero
   real(kind=dp) :: uscale, uzero
   real(kind=dp) :: vscale, vzero
   real(kind=dp) :: wscale, wzero
   real(kind=dp) :: jdscal1, jdscal2, jdzero

!  Output buffer for FITS records
   integer, parameter :: blocks=1
   integer, parameter :: bsize=2880
   integer, parameter :: mbuff=blocks*bsize/4
   real*4 :: buffer(mbuff)
   integer :: ibuff, nblock

!  Output buffer for FITS records
   integer, parameter :: MAX_SAMP = 43200
   integer*4 :: utss(MAX_SAMP), stss(MAX_SAMP)
   integer :: obs_epoch
   character*32 :: obs_name
   real*4 :: obs_raref,obs_decref
   real*4 :: obs_rarot,obs_decrot
   integer :: ifile, filep
   integer, parameter :: MAX_HIST = 32
   character*72 :: history(MAX_HIST)
   integer*4 :: histkey
   integer :: lhist


!  Magic blanking value for FITS files
   real(kind=dp), parameter :: MAGICBLANK=1.23456d-9

!  String length for file names
   integer, parameter :: fname_len=120

!  Arrays for mosaiced observations
   real*4, dimension(:), allocatable :: ra_offset, dec_offset

!  Graphics definitions
   integer, save :: pl_dev1
   integer ::  plot_col, line_sty, line_width
   real :: char_hgt  

!  Directories and other variables for file_handling
   character*120 :: profile_data_dir, model_dir, telescope_dir, postscript_dir
   character*120 :: profile_scripts_dir
   character*120 :: prompt, which_dir, rootname, extn

!  Variables for stitching together maps in combine_maps
   integer :: n_fits_files
   character*120, dimension(:), allocatable :: fits_list
   real, allocatable, dimension(:) :: fits_ra, fits_dec, fits_ra1, fits_dec1
   integer, allocatable, dimension(:) :: fits_npix_ra, fits_npix_dec 
   real, allocatable, dimension(:) :: fits_bmaj, fits_bmin, fits_bpa
   real, allocatable, dimension(:) :: fits_crpix_ra, fits_crpix_dec
   real, allocatable, dimension(:) :: fits_noise
   real(dp), allocatable, dimension(:) :: fits_incell
   logical, allocatable, dimension(:) :: fits_use_file
   logical :: write_fits_beam_table

!  Top-level commands
   integer, parameter :: ncomm = 92
   character*40, dimension(ncomm) :: comms
   data comms(1) /'quit'/
   data comms(2) /'read-fits-map'/
   data comms(3) /'display-map'/
   data comms(4) /'make-aperture'/
   data comms(5) /'difference-maps'/
   data comms(6) /'copy-map'/
   data comms(7) /'make-cluster'/
   data comms(8) /'help'/
   data comms(9) /'get-cluster-model'/
   data comms(10) /'plot-map-profile'/
   data comms(11) /'histogram-map'/
   data comms(12) /'make-model-xray-map'/
   data comms(13) /'flag-region'/
   data comms(14) /'get-directory'/
   data comms(15) /'sum-maps'/
   data comms(16) /'open-plot-device'/
   data comms(17) /'close-plot-device'/
   data comms(18) /'rescale-map'/
   data comms(19) /'smooth-map'/
   data comms(20) /'plot-model-visibility'/
   data comms(21) /'make-source-population'/
   data comms(22) /'plot-spectra'/
   data comms(23) /'get-cosmology'/
   data comms(24) /'estimate-map-background'/
   data comms(25) /'bayesian-xray-fit'/
   data comms(26) /'graphics-operations'/
   data comms(27) /'make-test-sky'/
   data comms(28) /'make-pointed-observation'/
   data comms(29) /'plot-likelihood'/
   data comms(30) /'display-sz-sky'/
   data comms(31) /'difference-map-profiles'/
   data comms(32) /'plot-maximum-amplitude'/
   data comms(33) /'get-geometry'/
   data comms(34) /'clear'/
   data comms(35) /'simple-clean'/
   data comms(36) /'show-cluster-parameters'/
   data comms(37) /'read-uv-fits'/
   data comms(38) /'plot-cut'/
   data comms(39) /'mcg-xray-fit'/
   data comms(40) /'combine-maps'/
   data comms(41) /'write-model'/
   data comms(42) /'read-model'/
   data comms(43) /'explain-profile'/
   data comms(44) /'explain-command'/
   data comms(45) /'stat-map'/
   data comms(46) /'test-cmd'/
   data comms(47) /'estimate-source-confusion'/
   data comms(48) /'lensing'/
   data comms(49) /'adaptive-smooth'/
   data comms(50) /'fill-visibilities'/
   data comms(51) /'project-cluster'/
   data comms(52) /'run'/
   data comms(53) /'read-sz-sky'/
   data comms(54) /'plot-source-confusion'/
!   data comms(55) /'make-3d-skies'/
!   data comms(56) /'wills-make-skies'/  
!   data comms(57) /'load-sim-skies'/    ! a quick hacky version...
   data comms(58) /'exit'/
   data comms(59) /'multiple-bayessian-xray-fit'/
!   data comms(60) /'read-helen-sky'/
   data comms(61) /'monte-carlo-xray-fit'/
!   data comms(62) /'read-sky'/
   data comms(63) /'make-mosaic-observation'/
!   data comms(64) /'monte-carlo-primordial-effect'/
   data comms(65) /'get-telescope'/
!   data comms(66) /'measure-fit-quality-bulk'/ ! experimental command
   data comms(67) /'weight-my-data'/ ! uv sampling
   data comms(68) /'make-spectral-index-map'/
   data comms(69) /'overlay_mosaic'/ !plot fwhm of mosaic fields
   data comms(70) /'pointing_power'/ !calc change in power for offset
   data comms(71) /'evaluate'/
   data comms(72) /'remove-ska-ptsrcs'/ !specific to ska simulations
!   data comms(73) /'write-sim-skies'/ ! writes out sky as txt file
   data comms(74) /'plot-scaling'/ ! plots z-r-T scaling relation
   data comms(75) /'make-cmb-sky'/ ! creates cmb visibilities
   data comms(76) /'w-projection'/ ! applies w correction
   data comms(77) /'grid-visibilities'/ ! grids visibilities
   data comms(78) /'read-visibilities'/ ! read data from .vis file
   data comms(79) /'write-uv-fits'/ ! write uv fits file
   data comms(80) /'map-visibility-data'/
!   data comms(81) /'set-wis'/ ! set visibility files to include w
   data comms(82) /'subtract-point-sources'/ 
   data comms(83) /'write-fits-map'/
   data comms(84) /'correlate-maps'/
!   data comms(85) /'vla-true'/
!   data comms(86) /'vla-false'/
   data comms(87) /'test-quadrature'/
   data comms(88) /'correct-single-dish-beam'/
   data comms(89) /'get-xray-telescope'/
   data comms(90) /'plot-cosmology'/
   data comms(91) /'display-primary-beam'/
   data comms(92) /'make-drift-observation'/

! Modified 'help'
! prints commands in logical order

   integer :: jparm
   integer, parameter :: nsep=6, n_unused=16, n_index = ncomm+nsep-n_unused
   integer index(ncomm+nsep)

   data (index(jparm), jparm=1, n_index)    /&
       & 16,17,34,26,36,43,44,8,1,58    ,0,  & ! General control commands
       & 2,3,5,6,13,15,18,19,24,31,38,45  ,  & ! map commands
       & 49,68,84,11                    ,0,  & ! map commands
       & 4,7,21,27,28,50,63,92,75,77,80,12 ,0,  & ! make command
       & 9,23,33,65,89                  ,0,  & ! get commands
       & 41,42,53,78,79,83              ,0,  & ! file handling commands
       & 10,20,22,29,90,91              ,0,  & ! plot commands
       & 25,30,32,35,39,47,48,51,59,61,71 ,  & ! misc commands
       & 67,69,70,72,74,76,81,82,87,88       & ! misc commands
       & /

   data (index(jparm), jparm=n_index+1, n_index+n_unused)   / &
       & 99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99    /

   interface do_allocation
      module procedure do_allocation
   end interface

   interface do_deallocation
      module procedure do_deallocation
   end interface

   public :: allocate_telescope
   public :: deallocate_telescope
   public :: allocate_frequency
   public :: deallocate_frequency
   public :: allocate_sources
   public :: deallocate_sources
   public :: allocate_combine_fits
   public :: deallcate_combine_fits

contains

! *************************************************************************

! NB try and keep order of allocation, initialisation and deallocation of
! arrays the same
   subroutine do_allocation

      implicit none

      logical :: reduced

! A flag to only use szsky rather than all the available arrays - useful for
! simulating huge sky areas requiring huge numbers of pixels
      reduced = .false.

! Allocate global arrays
      if (reduced) then
         allocate (szsky(maxsize,maxsize,nchan))
      else
         allocate (re(maxsize,maxsize,nchan))
         allocate (im(maxsize,maxsize,nchan))
         allocate (xmap(maxsize,maxsize))
         allocate (rosat_map(maxsize,maxsize))
         allocate (diff(maxsize,maxsize))
         allocate (szsky(maxsize,maxsize,nchan))
         allocate (compton_y(maxsize,maxsize))
         allocate (temp_sky(maxsize,maxsize))
         allocate (xsky(maxsize,maxsize))
         allocate (ros_beam_re(maxsize,maxsize))
         allocate (ros_beam_im(maxsize,maxsize))
         allocate (nsky(maxsize,maxsize))
         allocate (ok_data(maxsize,maxsize))
         allocate (gr_re(maxsize,maxsize))
         allocate (gr_im(maxsize,maxsize))
         allocate (uvcov(maxsize,maxsize))
         allocate (dirty_map(maxsize,maxsize))
         allocate (clean_map(maxsize,maxsize))
         allocate (beam(maxsize,maxsize))
         allocate (comp_map(maxsize,maxsize))
         allocate (resid_map(maxsize,maxsize))
         allocate (sky_weight(maxsize,maxsize))
      end if
     
! Allocate global arrays with dimension maxsize/2
      allocate (map_profile(maxsize/2))
      allocate (map_profile2(maxsize/2))

! Allocate global arrays with dimension (maxsize*(maxsize/2+1))
      allocate (aperture(maxsize*(maxsize/2+1)))

! Initialise arrays
      if (reduced) then
         szsky = 0.0
      else
         re = 0.d0
         im = 0.d0
         xmap = 0.d0
         rosat_map = 0.d0
         diff = 0.d0
         szsky = 0.d0
         compton_y = 0.d0
         xsky = 0.d0
         ros_beam_re = 0.d0
         ros_beam_im = 0.d0
         nsky = 0.d0
         ok_data = .false.
         gr_re = 0.d0
         gr_im = 0.d0
         uvcov = 0.d0
         dirty_map = 0.d0
         clean_map = 0.d0
         beam = 0.d0
         comp_map = 0.d0
         resid_map = 0.d0
         sky_weight = 0.d0
      end if
      map_profile = 0.d0
      map_profile2 = 0.d0
      aperture = 0.d0

      flag_nd = 0

   end subroutine do_allocation

! ***************************************************************************

! NB try and keep order of allocation, initialisation and deallocation of
! arrays the same
   subroutine do_deallocation

      implicit none

      if(associated(p_sky)) nullify(p_sky)
      if(associated(p3_sky)) nullify(p3_sky)
      if(allocated(re)) deallocate (re)
      if(allocated(im)) deallocate (im)
      if(allocated(xmap)) deallocate (xmap)
      if(allocated(rosat_map)) deallocate (rosat_map)
      if(allocated(diff)) deallocate (diff)
      if(allocated(szsky)) deallocate (szsky)
      if(allocated(compton_y)) deallocate (compton_y)
      if(allocated(temp_sky)) deallocate (temp_sky)
      if(allocated(xsky)) deallocate (xsky)
      if(allocated(ros_beam_re)) deallocate (ros_beam_re)
      if(allocated(ros_beam_im)) deallocate (ros_beam_im)
      if(allocated(nsky)) deallocate (nsky)
      if(allocated(ok_data)) deallocate (ok_data)
      if(allocated(gr_re)) deallocate (gr_re)
      if(allocated(gr_im)) deallocate (gr_im)
      if(allocated(uvcov)) deallocate (uvcov)
      if(allocated(dirty_map)) deallocate (dirty_map)
      if(allocated(clean_map)) deallocate (clean_map)
      if(allocated(beam)) deallocate (beam)
      if(allocated(comp_map)) deallocate (comp_map)
      if(allocated(resid_map)) deallocate (resid_map)
      if(allocated(sky_weight)) deallocate (sky_weight)
      if(allocated(map_profile)) deallocate (map_profile)
      if(allocated(map_profile2)) deallocate (map_profile2)
      if(allocated(aperture)) deallocate (aperture)

   end subroutine do_deallocation

!***************************************************************************

   subroutine allocate_telescope

      allocate(x_pos(n_antennas))
      allocate(y_pos(n_antennas))
      allocate(z_pos(n_antennas))
      allocate(basel_x(n_antennas,n_antennas))
      allocate(basel_y(n_antennas,n_antennas))
      allocate(basel_z(n_antennas,n_antennas))
      call deallocate_frequency
      call allocate_frequency

   end subroutine allocate_telescope

!***************************************************************************

   subroutine deallocate_telescope

      if (allocated(x_pos)) deallocate(x_pos)
      if (allocated(y_pos)) deallocate(y_pos)
      if (allocated(z_pos)) deallocate(z_pos)
      if (allocated(basel_x)) deallocate(basel_x)
      if (allocated(basel_y)) deallocate(basel_y)
      if (allocated(basel_z)) deallocate(basel_z)
      call deallocate_frequency

   end subroutine deallocate_telescope

!***************************************************************************

   subroutine allocate_frequency

      allocate(nu(nchan))
      allocate(fwhm(nchan))

   end subroutine allocate_frequency

!***************************************************************************

   subroutine deallocate_frequency

      if (allocated(nu)) deallocate(nu)
      if (allocated(fwhm)) deallocate(fwhm)

   end subroutine deallocate_frequency

!***************************************************************************

   subroutine allocate_sources

      allocate(rapos(nsrc))
      allocate(decpos(nsrc))
      allocate(srcflux(nsrc))
      allocate(srcalpha(nsrc))
      allocate(radist(nsrc))
      allocate(decdist(nsrc))
      allocate(jradist(nsrc))
      allocate(jdecdist(nsrc))
      allocate(src_x(nsrc))
      allocate(src_y(nsrc))

   end subroutine allocate_sources

!***************************************************************************

   subroutine deallocate_sources

      if (allocated(rapos)) deallocate(rapos)
      if (allocated(decpos)) deallocate(decpos)
      if (allocated(srcflux)) deallocate(srcflux)
      if (allocated(srcalpha)) deallocate(srcalpha)
      if (allocated(radist)) deallocate(radist)
      if (allocated(decdist)) deallocate(decdist)
      if (allocated(jradist)) deallocate(jradist)
      if (allocated(jdecdist)) deallocate(jdecdist)
      if (allocated(src_x)) deallocate(src_x)
      if (allocated(src_y)) deallocate(src_y)

   end subroutine deallocate_sources

!***************************************************************************

   subroutine allocate_combine_fits

      allocate(fits_list(n_fits_files))
      allocate(fits_ra(n_fits_files))
      allocate(fits_dec(n_fits_files))
      allocate(fits_ra1(n_fits_files))
      allocate(fits_dec1(n_fits_files))
      allocate(fits_npix_ra(n_fits_files))
      allocate(fits_npix_dec(n_fits_files))
      allocate(fits_crpix_ra(n_fits_files))
      allocate(fits_crpix_dec(n_fits_files))
      allocate(fits_bmaj(n_fits_files))
      allocate(fits_bmin(n_fits_files))
      allocate(fits_bpa(n_fits_files))
      allocate(fits_noise(n_fits_files))
      allocate(fits_incell(n_fits_files))
      allocate(fits_use_file(n_fits_files))

   end subroutine allocate_combine_fits

!***************************************************************************

   subroutine deallocate_combine_fits

      if (allocated(fits_list)) deallocate(fits_list)
      if (allocated(fits_ra)) deallocate(fits_ra)
      if (allocated(fits_dec)) deallocate(fits_dec)
      if (allocated(fits_ra1)) deallocate(fits_ra1)
      if (allocated(fits_dec1)) deallocate(fits_dec1)
      if (allocated(fits_npix_ra)) deallocate(fits_npix_ra)
      if (allocated(fits_npix_dec)) deallocate(fits_npix_dec)
      if (allocated(fits_crpix_ra)) deallocate(fits_crpix_ra)
      if (allocated(fits_crpix_dec)) deallocate(fits_crpix_dec)
      if (allocated(fits_bmaj)) deallocate(fits_bmaj)
      if (allocated(fits_bmin)) deallocate(fits_bmin)
      if (allocated(fits_bpa)) deallocate(fits_bpa)
      if (allocated(fits_noise)) deallocate(fits_noise)
      if (allocated(fits_incell)) deallocate(fits_incell)
      if (allocated(fits_use_file)) deallocate(fits_use_file)

   end subroutine deallocate_combine_fits

!***************************************************************************

end module sz_globals


























