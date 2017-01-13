! subroutine that projects a given cluster model back in redshift
subroutine project_cluster

   use kind_def
   use sz_globals
   use cosmology
   use kgraph

   implicit none

   integer :: i, j, k, n_iter
   real(kind=dp) :: z_max, z_min, r_core
   integer, parameter :: resol=2
   integer, parameter :: npoints = 10
   real (kind=dp) :: upt,vpt,psi,inc,rad_inc
   real (kind=dp), dimension(:), allocatable :: max_amp1, max_amp2,z1
   real (kind=dp), dimension(npoints) :: re_pt, im_pt, amp
   real (kind=dp), dimension(resol) :: max_amp, rad
   character*80 :: title, x_axis, y_axis

   z_max = 0.4
   r_core = 0.25*Mpc
   ellipsoid = .false.
   n_iter = 100

   call io_getd('Maximum redshift ','*',z_max,status)
   call io_geti('Number of redshifts ','*',n_iter,status)

! Allocate arrays
   allocate(max_amp1(n_iter))
   allocate(max_amp2(n_iter))
   allocate(z1(n_iter))

! Loops over all redshifts
   do k = 1, n_iter

! Logarithmic increments
      z_min = 0.002
      z = exp((k-1.0)*(log(z_max)-log(z_min))/(n_iter-1.0)+log(z_min))

! Linear increments
!      z = (zmax*k)/n_iter
      z1(k) = z
      call lumdist
      D_lum = D_lum*Mpc         ! D_lum in metres now....
      call angdist
      D_theta = D_theta*Mpc

! Calculate core radius in radians
      theta_cx  = r_core/ D_theta

! Convert to arcsec
      theta_cx = theta_cx/sec2rad

      verbose = .false.
      from_file = .false.
      do_plot = .false.
      from_poly = .false.
      overwrite = .true.

      call make_cluster

      call make_aperture

      telescope : select case(which_tel)
      case('c','C','m','M')
         rad_inc = 100./resol
      case('e','E')
         rad_inc = 200./resol
      case default
         rad_inc = 2000./resol
      end select telescope


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
      end do

      max_amp1(k) = max_amp(1)
      max_amp2(k) = max_amp(2)

      write(*,*) 'At redshift',z,'core radius is ',theta_cx
      write(*,*) 'Angular diameter distance is ',D_theta
      write (*,*)'Max amp at ',rad(1),' lambda is ',max_amp1(k)
      write (*,*)'Max amp at ',rad(2),' lambda is ',max_amp2(k)
      write(*,*)

   end do

   x_axis = 'Redshift'
   y_axis = 'Flux Density'
   title = 'Cluster projected in redshift'
   call pgpage 
   call display_profile(z1,max_amp1,n_iter,.true.) 
   call pglabel(x_axis,y_axis,title)
   call pgpage
   call display_profile(z1,max_amp2,n_iter,.true.)
   call pglabel(x_axis,y_axis,title)

! Deallocate arrays
   deallocate(z1)
   deallocate(max_amp1)
   deallocate(max_amp2)

end subroutine project_cluster
