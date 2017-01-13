!     Subroutine to make aperture plane from sz sky

subroutine make_aperture

   use kind_def
   use sz_globals
   use maths
   use kgraph
   use map_handling
   use primary_beam

   implicit none

   real(kind=dp) :: temp1, temp2, taper, x, y
   real(kind=dp) :: r
   real(kind=dp) :: max, norm, rvalue
   real(kind=dp) :: rolloff, sbeam, sigma_pb
   integer :: i, j, k, maxi, maxj, ncell, nr, ni
   integer, dimension(2) :: pos_max
   real(kind=dp), dimension(maxsize,maxsize) :: retmp, ampl
   ncell = maxsize

! loop over frequency channels
   do chan = 1,nchan

      if (verbose) then
         write(*,*) 'Processing channel ',chan
      end if

! szsky is brightness temperatures therefore use R-J approx
! and multiply by solid angle of pixel size to get flux density in Jy
      write(*,*) cellsize, sec2rad, k_b, nu(chan), const_c2
      norm = (cellsize*sec2rad)**2*2.0*k_b*nu(chan)**2/const_c2*1.0D26
      if (verbose) then
         write(*,*) &
         'Converting from brightness temperature to Jy - factor = ', norm
      end if

! set up model sky distribution
      do j = 1,maxsize
         do i = 1, maxsize            
            temp1 = (float(i-maxsize/2)-1.0)*cellsize-pb_x
            temp2 = (float(j-maxsize/2)-1.0)*cellsize-pb_y
            call get_taper(temp1,temp2,taper)
            re(i,j,chan) = taper*szsky(i,j,chan)
            re(i,j,chan) = re(i,j,chan)*norm
            im(i,j,chan) = 0.0
         end do
      end do

      max = maxval(re(:,:,chan))
      pos_max = maxloc(re(:,:,chan))
      if (verbose) write(*,*)'Maximum ',max,' at ',pos_max
      if (do_plot) then
         write(*,*) 'Displaying sz sky'
         call display_map(re(:,:,chan),ncell,graph_adj,&
                          cellsize,cellsize,0D0,0D0)
      end if

! apply shift for phase center not coincident with array centre
      if ((pb_x.ne.0.).or.(ph_x.ne.0.).or.(pb_y.ne.0.).or.(ph_y.ne.0.))then
         if (verbose) write(*,*) 'Shifting phase and pointing centers'
         retmp = re(:,:,chan)
         do j = 1,maxsize
            do i = 1,maxsize
               x = float(i)+(pb_x+ph_x)/cellsize
               y = float(j)+(pb_y+ph_y)/cellsize
               call interp(retmp,maxsize,maxsize,x,y,rvalue)
               re(i,j,chan) = rvalue
            end do
         end do
         if (verbose) write(*,*) 'Shift done'
         if (do_plot) then
            write(*,*) 'Displaying displaced sz sky'
            call display_map(re(:,:,chan), ncell, graph_adj,&
                             cellsize,cellsize,0D0,0D0) 
         end if
      end if

! FT it to get aperture
      call do_2d_fft(re(:,:,chan),im(:,:,chan),maxsize,maxsize,1)

      max = maxval(re(:,:,chan))
      pos_max = maxloc(re(:,:,chan))
      if (verbose) write(*,*)'Maximum in real part of aperture ',max, &
                          ' at ',pos_max   

      if (do_plot) then
         write(*,*) 'Displaying real part of aperture'
         call display_map(re(:,:,chan),maxsize,graph_adj,&
                          cellsize,cellsize,0D0,0D0)
         write(*,*) 'Displaying imaginary part of aperture'
         call display_map(im(:,:,chan),maxsize,graph_adj,&
                          cellsize,cellsize,0D0,0D0)
      end if

   end do

end subroutine make_aperture



