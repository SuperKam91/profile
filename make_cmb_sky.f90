! New version of this code. The original was *very* strange and possibly
! incorrect.... [KG]

subroutine make_cmb_sky(idum,exact_spec)

   use kind_def
   use sz_globals
   use map_handling
   use kgraph
   use physics
   use maths

   integer :: idum
   logical :: exact_spec
   integer :: i, j, k, uvmapsize, centre
   integer, parameter :: nl=10000
   integer :: ell_filter, ell_width
   real(dp) :: uvcell, ell(nl), cl(nl), power, uvdist, clval, f_start
   real(dp),allocatable,dimension(:,:) :: uvreal, uvimag, sigmacl
   logical :: apply_hpf
   character*80 :: inspec
   logical, external :: io_yesno
   character*100 :: exec_dir
   
   ell_filter = 1000
   ell_width = 100

! Read input power spectrum
   call getarg(0, exec_dir)
   l = scan(exec_dir, '/', .true.)
   if (l.gt.0) exec_dir=exec_dir(1:l)
   call io_getc('Input power spectrum:',trim(exec_dir)//'cl8000.dat',inspec,status)
   open(15,file = inspec)
   
   i = 1
777 read(15,*,end=401) ell(i),cl(i)
   maxi = i
   i = i+1
   
   goto 777
         
401 close(unit=15)
         
   write(*,*) 'Cl spectrum read: ', maxi,' lines.'
   
   if (io_yesno('Apply high pass filter','no',status)) then
      call io_geti('Filter out ell less than','*',ell_filter,status)
      do i = 1, ell_filter
         cl(i) = 0.d0
      end do
      call io_geti('Apply Hanning window of width','*',ell_width,status)
      do i = 1, ell_width
         j = i+ell_filter
         cl(j) = cl(j)*sin(pi/2.*i/ell_filter)
      end do

   end if

! Calculate cl spectrum in mu K^2
   do i = 1,maxi
      cl(i) = cl(i)*(T0*1.d6)**2*2.*pi/(ell(i)*(ell(i)+1))
   end do

   write(*,*) 'Total power in spectrum (muK) = ',&
              sqrt(sum((2.*ell+1.)*cl/4./pi))

   if (verbose) then
      uvcell = 1./(maxsize*cellsize*sec2rad)
      write(*,*) 'uv cell = ', uvcell, 'lambda'
      uvmapsize = (1./(cellsize*sec2rad))/(uvcell)
      write(*,*) 'uv mapsize = ', uvmapsize,' pixels square.'
   end if

   centre = uvmapsize/2+1

! Allocate arrays and initialise
   if (allocated(uvreal)) deallocate(uvreal)
   allocate(uvreal(uvmapsize,uvmapsize))
   if (allocated(uvimag)) deallocate(uvimag)
   allocate(uvimag(uvmapsize,uvmapsize))
   if (allocated(sigmacl)) deallocate(sigmacl)
   allocate(sigmacl(uvmapsize,uvmapsize))

   sigmacl = 0.0
   uvreal = 0.0
   uvimag = 0.0

! Generate Cl amplitudes for half of the uv-plane
   do i = 1,centre
      do j = 1,uvmapsize

! Calculate distance in uv plane and multiply by 2 pi to get in l distance
         uvdist = 2.*pi*sqrt(real(i-uvmapsize/2-1)**2+&
                             real(j-uvmapsize/2-1)**2)*uvcell

! Check with Mike that it makes sense to do this
         call interp(cl,nl,uvdist,clval)
         sigmacl(i,j) = sqrt(clval)      
      end do
   end do

! Generate a realisation of the Cl spectrum if required
   if (.not.exact_spec) then
      do i = 1,centre
         do j = 1,uvmapsize
            sigmacl(i,j) = sigmacl(i,j)*gasdev(idum) 
         end do
      end do
   end if

! Generate real and imaginary parts with randomised phases
   do i = 1,centre
      do j = 1,uvmapsize
         phase=2.*pi*ran2(idum)
         uvreal(i,j) = cos(phase)*sigmacl(i,j)
         uvimag(i,j) = sin(phase)*sigmacl(i,j)
      end do
   end do

! Symmetrise Fourier plane
   do i=centre,uvmapsize
      do j=1,uvmapsize
         uvreal(i,j) = uvreal(uvmapsize+2-i,uvmapsize+2-j)
         uvimag(i,j) = -uvimag(uvmapsize+2-i,uvmapsize+2-j)
      end do
   end do
   do i=centre+1,uvmapsize
      uvreal(i,1) = uvreal(uvmapsize+2-i,1)
      uvimag(i,1) = -uvimag(uvmapsize+2-i,1)
      uvreal(1,i) = uvreal(1,uvmapsize+2-i)
      uvimag(1,i) = -uvimag(1,uvmapsize+2-i)
   end do

   if (do_plot) then
      write(*,*) 'Power level in (half) uv plane'
      call display_map(sigmacl,uvmapsize,graph_adj,uvcell,uvcell,0D0,0D0)
   end if

! blank out DC term
   !uvreal(centre,centre) = 0.
   !uvimag(centre,centre) = 0.
 
! normalization
   uvreal = uvreal/uvcell
   uvimag = uvimag/uvcell
   if (verbose) then
      power = 0.
      do i = 1,uvmapsize
         do j = 1,uvmapsize
            power = power+uvreal(i,j)**2+uvimag(i,j)**2
         end do
      end do
      write(*,*) 'UV power (muK) = ',sqrt(real(power))*uvcell**2
   end if

   if (do_plot) then
      write(*,*) 'Real part of aperture plane'
      call display_map(uvreal,uvmapsize,graph_adj,uvcell,uvcell,0D0,0D0)
      write(*,*) 'Imaginary part of aperture plane'
      call display_map(uvimag,uvmapsize,graph_adj,uvcell,uvcell,0D0,0D0)
   end if

! FFT to get sky
   call do_2d_fft(uvreal,uvimag,uvmapsize,uvmapsize,-1)

   uvreal = uvreal*uvcell**2
   if (do_plot) then
      call display_map(uvreal,uvmapsize,graph_adj,cellsize,cellsize,0D0,0D0)   
      write(*,*) 'Real map, RMS(muK) = ',&
              sqrt(sum(uvreal**2)/(real(uvmapsize)**2))
   end if

! Write to SZ sky, converting from mu K of thermodynamic temperature to
! K brightness temperature
! Convert to intensity by multiplying by dB_by_dT; convert to brightness
! temperature by dviding by 2k/lambda^2
   do chan = 1, nchan
      f_start = nu(chan)
      write(*,*) (dB_by_dT(T0,f_start))/rayleigh_jeans_fn(f_start)
      if (overwrite) then
         szsky(:,:,chan) = uvreal/1.d6*(dB_by_dT(T0,f_start))/&
               rayleigh_jeans_fn(f_start)
      else
         szsky(:,:,chan) = szsky(:,:,chan)+uvreal/1.d6*(dB_by_dT(T0,f_start))/&
              rayleigh_jeans_fn(f_start)
      end if
   end do
   
! Deallocate arrays
   if (allocated(uvreal)) deallocate(uvreal)
   if (allocated(uvimag)) deallocate(uvimag)
   if (allocated(sigmacl)) deallocate(sigmacl)

end subroutine make_cmb_sky



