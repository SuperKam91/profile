! makes fake Rosat map from model, convolving with Rosat beam and adding
! background
subroutine make_model_xray_map(idum,do_poisson)

   use kind_def
   use sz_globals
   use maths

   implicit none
   integer :: idum
   logical, intent(in) :: do_poisson
   real(kind=dp), dimension(maxsize,maxsize) :: xmap_im, temp_re, temp_im
   integer :: i, j, k, kmax
   real(kind=dp) :: temp, tot, k_x, dl, pixsize

! Make model X-ray map by convolving with beam if observing with PSPC
   if (x_beam) then

      xmap = xsky*rosat_beam_sum

      xmap_im = 0.0

! To convolve rosat map, fft it, multiply by the fft of the beam, and then fft 
! back

! FFT raw map
      call do_2d_fft(xmap,xmap_im,maxsize,maxsize,1)

      temp_re = xmap
      temp_im = xmap_im
      xmap = temp_re*ros_beam_re-temp_im*ros_beam_im
      xmap_im = temp_re*ros_beam_im+temp_im*ros_beam_re

! FFT back to get xmap
      call do_2d_fft(xmap,xmap_im,maxsize,maxsize,-1)

! make a *real* map, ie with poisson noise
      if (do_poisson) then
         do j = 1,maxsize
            do i = 1,maxsize
               xmap(i,j) = poidev(xmap(i,j),idum)
               xmap(i,j) = xmap(i,j)+poidev(back,idum)
            end do
         end do

! otherwise just make smooth map
      else
         xmap = xmap+back
      end if

      tot = sum(xmap)

      if (verbose) then
         write(*,*) 'Total X-ray flux = ',tot,' counts'
!         pixsize in 1000 km squares: convert to cm^2 later
!         pixsize = 1.0D-6*cellsize*sec2rad*D_theta

!         length integral in cm
!         dl = pixsize*1.0D8
!         pixsize = pixsize**2
         k_x = (Mpc/D_lum)
         k_x = k_x**2
         k_x = k_x*const
!         k_x = k_x*pixsize*const
         k_x = k_x*x_exp 
!         write(*,*) 'This corresponds (?) to a luminosity of ', tot/k_x
      end if

   else

! Make HRI model x-ray map equal to the observed sky plus background

      xmap = xsky

! if back < 0, make a *real* map, ie with poisson noise
      if (back.lt.0.0) then
         write(*,*) 'back = ', back
         back = -back
         do j = 1,maxsize
            do i = 1,maxsize
               xmap(i,j) = poidev(xmap(i,j)+back,idum)
!               xmap(i,j) = xmap(i,j)+poidev(back,idum)
            end do
         end do

! otherwise just make smooth map
      else
         xmap = xmap+back
      end if

! find total number x-ray counts
      tot = sum(xmap)

      if (verbose) then
         write(*,*) 'Total X-ray flux = ',tot,' counts'
      end if
    end if

end subroutine make_model_xray_map

! **************************************************************************


