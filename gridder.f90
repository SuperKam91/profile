
! this module contains code for gridding lists of visibilities onto a uv grid,
! gridding the w-projection and gridding the prolate spheroidal wavefunction.
! --------------------------------------------------------------------------

subroutine grid_vis(kernel)

   ! adopted from MAPPER

   use kind_def
   use sz_globals
   use pswf_globals
   use w_proj
   use map_handling
   use file_handling
   use kgraph
   use mikes_fft

   implicit none

   integer  :: i,j,k,centre,plane,fracu,fracv,n,m,count
   integer  :: job,iform,iw,nwplanes,pcen,file_len,lunit
   real     :: func,l,udist,vdist,temp
   real(dp) :: wei(maxsize,maxsize),uvcell,cutoff,wmax
   character*80, intent(in)  :: kernel
   character(len=fname_len) :: filename
   real(dp) :: spconvfunc(nplanes,csize,csize)
   integer, parameter :: maxsupport = 10
   logical  :: doplot, dofile
   character*1 :: chr1
   integer :: chr_lenb
   external chr_lenb
      
   uvcell = 1./(maxsize*cellsize*sec2rad)
   centre = maxsize/2+1

! Invert u coordinates of visibilities
   u = -u

   wei = 0.d0
   gr_re = 0.d0
   uvcov = 0.d0
   gr_im = 0.d0

! Convolve uv point with a sinc
   if (kernel=='sinc') then
      do k = 1,nvis
         call decide(i,u(k))
         call decide(j,v(k))
         if (((i.le.0).or.(i.gt.maxsize)).or.&
             ((j.le.0).or.(j.gt.maxsize))) then
            weight(k) = 0.d0
         else
            if (weight(k).gt.0.) then
               do n = -convsize,convsize
                  do m = -convsize,convsize
                     udist = real((i+n-centre)-(u(k)/uvcell))
                     vdist = real((j+m-centre)-(v(k)/uvcell))
                     l = sqrt(real(udist*udist+vdist*vdist))
                     call sinckernel(l,func)
                     if ((((j+m).ge.0).and.((j+m).lt.maxsize)).and.&
                          (((i+n).ge.0).and.((i+n).lt.maxsize))) then
                        gr_re(i+n,j+m) = gr_re(i+n,j+m)+func*data_re(k)*weight(k)
                        gr_re(maxsize+2-(i+n),maxsize+2-(j+m)) &
                             = gr_re(maxsize+2-(i+n),maxsize+2-(j+m)) &
                             + func*data_re(k)*weight(k)
                        gr_im(i+n,j+m) = gr_im(i+n,j+m)+func*data_im(k)&
                             *weight(k)
                        gr_im(maxsize+2-(i+n),maxsize+2-(j+m)) &
                             = gr_im(maxsize+2-(i+n),maxsize+2-(j+m)) &
                             - func*data_im(k)*weight(k)
                        uvcov(i+n,j+m) = uvcov(i+n,j+m)+func*weight(k)
                        uvcov(maxsize+2-(i+n),maxsize+2-(j+m)) &
                             = uvcov(maxsize+2-(i+n),maxsize+2-(j+m)) &
                             + func*weight(k)
                        wei(i+n,j+m) = wei(i+n,j+m)+weight(k)
                        wei(maxsize+2-(i+n),maxsize+2-(j+m)) &
                             = wei(maxsize+2-(i+n),maxsize+2-(j+m))+weight(k)
                     end if
                  end do
               end do
            end if
         end if
      end do
   elseif (kernel=='none') then
      do k = 1,nvis
         call decide(i,u(k))
         call decide(j,v(k))
         func = 1. 
         if (((i.le.0).or.(i.gt.maxsize)).or.&
             ((j.le.0).or.(j.gt.maxsize))) weight(k) = 0.
         if (weight(k).gt.0.) then
            gr_re(i,j) = gr_re(i,j)+func*data_re(k)*weight(k)
            gr_re(maxsize-i+2,maxsize-j+2) = gr_re(maxsize-i+2,maxsize-j+2) &
                 +func*data_re(k)*weight(k)
            gr_im(i,j) = gr_im(i,j)+func*data_im(k)*weight(k)
            gr_im(maxsize-i+2,maxsize-j+2) = gr_im(maxsize-i+2,maxsize-j+2) &
                 -func*data_im(k)*weight(k)
            uvcov(i,j) = uvcov(i,j)+func*weight(k)
            uvcov(maxsize-i+2,maxsize-j+2) = uvcov(maxsize-i+2,maxsize-j+2) &
                 + func*weight(k)
	    wei(i,j) = wei(i,j)+weight(k)
            wei(maxsize-i+2,maxsize-j+2) = wei(maxsize-i+2,maxsize-j+2) &
                 +weight(k)        
         end if
      end do
   elseif (kernel=='pswf') then
      call pswfgridder(spconvfunc)
      do k = 1,nvis
         call decide(i,u(k))
         call decide(j,v(k))
         if (((i.le.0).or.(i.gt.maxsize)).or.&
             ((j.le.0).or.(j.gt.maxsize))) weight(k) = 0.
         if (weight(k).gt.0.) then
            temp = int(u(k)+sign(0.5d0,u(k)))-u(k)-0.5
            temp = real(oversample)*temp
            fracu = int(temp+sign(0.5,temp))
            fracu = fracu+oversample
            temp = int(v(k)+sign(0.5d0,v(k)))-v(k)-0.5
            temp = temp*real(oversample)
            fracv = int(temp+sign(0.5,temp))
            fracv = fracv+oversample
            plane = fracu+oversample*fracv+1
            do n = -support,support
               do m = -support,support
                  if ((((j+m).ge.1).and.((j+m).le.maxsize)).and.&
                       (((i+n).ge.1).and.((i+n).le.maxsize))) then
                     func = spconvfunc(plane,n+ccentre,m+ccentre)
                     gr_re(i+n,j+m) = gr_re(i+n,j+m)+func*data_re(k)*weight(k)
                     gr_re(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          = gr_re(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          + func*data_re(k)*weight(k)
                     gr_im(i+n,j+m) = gr_im(i+n,j+m)+func*data_im(k)&
                          *weight(k)
                     gr_im(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          = gr_im(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          - func*data_im(k)*weight(k)
                     uvcov(i+n,j+m) = uvcov(i+n,j+m)+func*weight(k)
                     uvcov(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          = uvcov(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          + func*weight(k)
                     wei(i+n,j+m) = wei(i+n,j+m)+weight(k)
                     wei(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                           = wei(maxsize-(i+n)+2,maxsize-(j+m)+2)+weight(k)
                  end if
               end do
            end do
         end if
      end do
   elseif( kernel=='wpro') then
      wmax = dble(maxval(abs(real(w))))
100   call io_getd('wmax:','*',wmax,status)
      nwplanes = 21
      call io_geti('nwplanes:','*',nwplanes,status)
      cutoff = 1d-8
      call io_getd('cutoff:','*',cutoff,status)
      wstatus = 0
      if (allocated(convfunc)) deallocate(convfunc)
      allocate(convfunc(nwplanes*oversample*oversample,maxsupport,maxsupport))
      call wprojgridder(wmax,nwplanes,cutoff)
      write(*,*) 'Using wsupport = ',wsupport
      if (wstatus==1) then
         call io_getc('Try again?','y',chr1,status)
         if (chr1=='y') then
            wstatus = 0
            goto 100
         else
            return
         end if
      end if
      chan = 1
      wcentre = (nwplanes+1)/2 ! centre for w planes
      pcen = 2*(wsupport+1)/2  ! centre for pswf kernel
      wscale = wmax/(real(nwplanes-1)/2.)
      do k = 1,nvis
         call decide(i,u(k))
         call decide(j,v(k))
         if (((i.le.0).or.(i.gt.maxsize)).or.&
             ((j.le.0).or.(j.gt.maxsize))) weight(k) = 0.
         if (weight(k).gt.0.) then
            temp = int(u(k)+sign(0.5d0,u(k)))-u(k)-0.5
            temp = real(oversample)*temp
            fracu = int(temp+sign(0.5,temp))
            fracu = fracu+oversample
            temp = int(v(k)+sign(0.5d0,v(k)))-v(k)-0.5
            temp = temp*real(oversample)
            fracv = int(temp+sign(0.5,temp))
            fracv = fracv+oversample
! Why????
            w(k) = w(k)*nu(chan)/obsfreq
            iw = int((w(k)/wscale)+sign(0.5d0,w(k)))+wcentre
            plane = fracu+oversample*(fracv+oversample*iw)+1
            do n = -wsupport,wsupport
               do m = -wsupport,wsupport
                  if ((((j+m).ge.1).and.((j+m).le.maxsize)).and.&
                       (((i+n).ge.1).and.((i+n).le.maxsize))) then
                     func = convfunc(plane,n+pcen,m+pcen)
                     gr_re(i+n,j+m) = gr_re(i+n,j+m)+func*data_re(k)*weight(k)
                     gr_re(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          = gr_re(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          + func*data_re(k)*weight(k)
                     gr_im(i+n,j+m) = gr_im(i+n,j+m)+func*data_im(k)&
                          *weight(k)
                     gr_im(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          = gr_im(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          - func*data_im(k)*weight(k)
                     uvcov(i+n,j+m) = uvcov(i+n,j+m)+func*weight(k)
                     uvcov(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          = uvcov(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                          + func*weight(k)
                     wei(i+n,j+m) = wei(i+n,j+m)+weight(k)
                     wei(maxsize-(i+n)+2,maxsize-(j+m)+2) &
                           = wei(maxsize-(i+n)+2,maxsize-(j+m)+2)+weight(k)
                  end if
               end do
            end do
         end if
      end do
   end if

   dofile = .false.
   call io_getc('Write to file?','n',chr1,status)
   if (chr1=='y') dofile = .true.
   if (dofile) then
      prompt = 'File name'
      extn = '.grid'
      which_dir = profile_data_dir
      rootname = 'default'
      call get_write_filename(filename)
      call io_nxtlun(lunit,status)
      open (lunit,file=filename)
      write(*,*) 'Writing gridded data to ',filename
      write(lunit,*) 'DIMENSION:',maxsize
      write(lunit,*) 'CELLSIZE [lambda]:',uvcell
      write(lunit,*) 'CENTRE:', centre
      write(lunit,*) 'CONVOLUTION KERNEL:',kernel
      do i = 1,maxsize
         do j = 1,maxsize
            if (wei(i,j).ne.0.) then
               write(lunit,*) i,j,gr_re(i,j),gr_im(i,j),uvcov(i,j)
            else
               write(lunit,*) i,j,0.,0.,0.
            end if
         end do
      end do
   end if

   doplot = .false.
   call io_getc('Do plots?','n',chr1,status)
   if (chr1=='y') doplot = .true.

   if (doplot) then
      write(*,*) 'Displaying real part of visibilities'
      call display_map(gr_re,maxsize, graph_adj,cellsize,&
           cellsize,0D0,0D0)
      write(*,*) 'Displaying imaginary part of visibilities'
      call display_map(gr_im,maxsize, graph_adj,cellsize,&
           cellsize,0D0,0D0)
   end if
   
   write(*,*) nvis,' visibilities gridded'

! Re-invert u coords back again
   u = -u

end subroutine grid_vis

! -------------------------------------------------------------------------

! subroutine grid_correct(kernel,map)

!    use kind_def
!    use sz_globals
!    use pswf_globals
!    use w_proj
!    use kgraph

!    implicit none

!    integer  :: i,j,centre,job,iform,fracu,fracv,plane,m,n,pcen
!    complex  :: temparray(maxsize,maxsize)
!    real,allocatable  :: work(:)
!    real(dp)     :: map(maxsize,maxsize)
!    real     :: ker,temp,ucen,vcen,func
! !   real(dp) :: nux,nuy
! !   real(dp), allocatable  :: ccfx(:),ccfy(:)
!    character*80, intent(in)  :: kernel
!    real(dp) :: spconvfunc(nplanes,csize,csize)


   
! ! allocate dynamic arrays:
!    allocate(work(2*maxsize))
! !   allocate(ccfx(maxsize))
! !   allocate(ccfy(maxsize))

!    centre = (maxsize/2)
!    map = 0.
!    map = gr_re
!    gr_re = 0.

!    if (kernel=='sinc') then  
!       do i = 1,maxsize
!          do j = 1,maxsize
!             call sinckernel(sqrt(real((i-centre)**2+(j-centre)**2)),ker)
!             gr_re(i,j) = ker
!          end do
!       end do
!    elseif (kernel=='pswf') then
!       call pswfgridder(spconvfunc)
!       fracu = 0
!       fracv = 0
!       plane = fracu+oversample*fracv+1
!       do n = -support,support
!          do m = -support,support            
!             gr_re(maxsize/2+n+1,maxsize/2+m+1) = &
!                      spconvfunc(plane,n+ccentre,m+ccentre)
!          end do
!       end do
!    elseif (kernel=='wpro') then
!       fracu = 0
!       fracv = 0
!       plane = fracu+oversample*(fracv+oversample*wcentre)+1
!       pcen = 2*(wsupport+1)/2  ! centre for pswf kernel
!       do n = -wsupport,wsupport
!          do m = -wsupport,wsupport            
!             gr_re(maxsize/2+n+1,maxsize/2+m+1) = convfunc(plane,n+pcen,m+pcen)
!          end do
!       end do      
!    end if
   
!    temparray = cmplx(gr_re,0.)
   
! ! ... inverse fft of product
!    job=-1                     ! backward transform (-ve exponential)
!    iform=1                    ! data are complex
!    call makefft(maxsize,maxsize,temparray,job,iform,work)
!    temparray = temparray*real(maxsize)

!    gr_re = real(temparray)

!    ok_data = .true.
! 100 do i = 1,maxsize
!       do j = 1,maxsize
!          if (gr_re(i,j).ne.0.) then
!             map(i,j) = map(i,j)/gr_re(i,j)
!          end if
!          if (gr_re(i,j).lt.0.0001) then
!             map(i,j) = 0.0
!             ok_data(i,j) = .false.
!          end if
!       end do
!    end do
   
!    deallocate(work)
!  !  deallocate(ccfx)
!  !  deallocate(ccfy)

! end subroutine grid_correct

! -------------------------------------------------------------------------

! subroutine beam_correct(kernel,map)

!    use kind_def
!    use sz_globals
!    use pswf_globals
!    use kgraph
!    use w_proj

!    implicit none

!    integer  :: i,j,centre,job,iform,m,n,fracu,fracv,plane,pcen
!    complex  :: temparray(maxsize,maxsize)
!    real(kind=dp), dimension(:,:), allocatable :: fmap
!    real,allocatable     :: work(:)
!    real(dp)     :: map(maxsize,maxsize)
!    real     :: ker,ucen,vcen,func,temp
! !   real(dp) :: nux,nuy
! !   real(dp), allocatable  :: ccfx(:),ccfy(:)
!    character*80, intent(in)  :: kernel
!    real(dp) :: spconvfunc(nplanes,csize,csize)

! ! allocate dynamic arrays:
!    allocate(work(2*maxsize))
!    allocate(fmap(maxsize,maxsize))
! !   allocate(ccfx(maxsize))
! !   allocate(ccfy(maxsize))

!    centre = (maxsize/2)
!    map = 0.
!    map = uvcov      
!    fmap = 0.d0

!    if (kernel=='sinc') then
!       do i = 1,maxsize
!          do j = 1,maxsize
!             call sinckernel(sqrt(real((i-centre)**2+(j-centre)**2)),ker)
!             fmap(i,j) = ker
!          end do
!       end do
!    elseif (kernel=='pswf') then
!       call pswfgridder(spconvfunc)
!       fracu = 0
!       fracv = 0
!       plane = fracu+oversample*fracv+1
!       do n = -support,support
!          do m = -support,support
!             fmap(maxsize/2+n+1,maxsize/2+m+1) = spconvfunc(plane,n+ccentre,m+ccentre)
!          end do
!       end do
!    elseif (kernel=='wpro') then
!       fracu = 0
!       fracv = 0
!       plane = fracu+oversample*(fracv+oversample*wcentre)+1
!       pcen = 2*(wsupport+1)/2  ! centre for pswf kernel
!       do n = -wsupport,wsupport
!          do m = -wsupport,wsupport            
!             fmap(maxsize/2+n+1,maxsize/2+m+1) = convfunc(plane,n+pcen,m+pcen)
!          end do
         
!       end do
      
!    end if

!    temparray = cmplx(fmap(:,:),0.)
   
! ! ... inverse fft of product
!    job=-1                     ! backward transform (-ve exponential)
!    iform=1                    ! data are complex
!    call makefft(maxsize,maxsize,temparray,job,iform,work)
!    temparray = temparray*real(maxsize)

!    fmap(:,:) = real(temparray(:,:))

! 100 do i = 1,maxsize
!       do j = 1,maxsize
!          if (fmap(i,j).ne.0.) then
!             map(i,j) = map(i,j)/fmap(i,j)
!          end if
!          if (fmap(i,j).lt.0.0001) map(i,j) = 0.0
!       end do
!    end do

!    deallocate(work)

! end subroutine beam_correct

! -------------------------------------------------------------------------

subroutine kernel_correct(kernel)
! Replaces beam_correct and grid_correct

   use kind_def
   use sz_globals
   use maths
   use kgraph
   use pswf_globals
   use w_proj

   implicit none

   integer  :: i,j,centre,job,iform,m,n,fracu,fracv,plane,pcen
   real(kind=dp), dimension(:,:), allocatable :: fmap
   real     :: ker,ucen,vcen,func,temp
   character*80, intent(in)  :: kernel
   real(dp) :: spconvfunc(nplanes,csize,csize)

! allocate dynamic arrays:
   allocate(fmap(maxsize,maxsize))

   centre = (maxsize/2)+1
   fmap = 0.d0

   if (kernel=='sinc') then
      do i = 1, maxsize
         do j = 1, maxsize
            call sinckernel(sqrt(real((i-centre)**2+(j-centre)**2)),ker)
            fmap(i,j) = ker
         end do
      end do
   elseif (kernel=='pswf') then
      call pswfgridder(spconvfunc)
      fracu = 0
      fracv = 0
      plane = fracu+oversample*fracv+1
      do n = -support,support
         do m = -support,support
            fmap(maxsize/2+n+1,maxsize/2+m+1) = spconvfunc(plane,n+ccentre,m+ccentre)
         end do
      end do
   elseif (kernel=='wpro') then
      fracu = 0
      fracv = 0
      plane = fracu+oversample*(fracv+oversample*wcentre)+1
      pcen = 2*(wsupport+1)/2  ! centre for pswf kernel
      do n = -wsupport,wsupport
         do m = -wsupport,wsupport            
            fmap(maxsize/2+n+1,maxsize/2+m+1) = convfunc(plane,n+pcen,m+pcen)
         end do         
      end do
   end if

! Plot fft kernel 
   if (do_plot) then
      write(*,*) 'Displaying kernel'
      call display_map(fmap,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)  
   end if

! Inverse FFT
   temp_sky = 0.d0
   call do_fft(fmap,temp_sky,maxsize,maxsize,-1)
   ok_data = .true.

! Plot fft kernel 
   if (do_plot) then
      write(*,*) 'Displaying FFT of kernel'
      call display_map(fmap,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)  
   end if

   do i = 1,maxsize
      do j = 1,maxsize
         if (fmap(i,j).ne.0.) then
            dirty_map(i,j) = dirty_map(i,j)/fmap(i,j)
            beam(i,j) = beam(i,j)/fmap(i,j)
         end if
         if (fmap(i,j).lt.0.0001) then
            dirty_map(i,j) = 0.d0
            beam(i,j) = 0.d0
            ok_data(i,j) = .false.
         end if
      end do
   end do

   deallocate(fmap)

end subroutine kernel_correct

! -------------------------------------------------------------------------

subroutine map_vis(kernel)

   use sz_globals
   use maths
   use kgraph
   use map_handling
   use mikes_fft

   implicit none

   integer  :: job,iform,i,j,centre, nbin
   complex, allocatable, dimension(:,:)  :: temp
   real     :: r, x, y, width, datmin, datmax
   real(dp)  :: bmax
   character*80, intent(in)  :: kernel
   real,allocatable  :: work(:)
   logical  :: doplot, use_mike_fft
   logical, external :: io_yesno

   use_mike_fft = .false.

   centre = (maxsize/2)+1

   if (use_mike_fft) then
! FFT routine has some issues....
      allocate(temp(maxsize,maxsize))
      temp = 0.
      temp = cmplx(gr_re,gr_im)

! ... inverse fft of product
      allocate(work(2*maxsize))
      job=-1                     ! backward transform (-ve exponential)
      iform=1                    ! data are complex
      call makefft(maxsize,maxsize,temp,job,iform,work)

! like in norm in mapper:
      temp = temp*(real(maxsize)**2)
      temp = temp/real(maxsize)
      dirty_map = real(temp)
  
! do beam:
      temp = cmplx(uvcov,0.)

! ... inverse fft of product
      job=-1                     ! backward transform (-ve exponential)
      iform=1                    ! data are complex
      call makefft(maxsize,maxsize,temp,job,iform,work)
      temp = temp*(real(maxsize)**2)
      temp = temp/real(maxsize)
   
      uvcov = real(temp)

      deallocate(temp)
      deallocate(work)
   else

! Preferred FFT routine
      dirty_map = gr_re
      temp_sky = gr_im
      call do_fft(dirty_map,temp_sky,maxsize,maxsize,-1)
      beam = uvcov
      temp_sky = 0.d0
      call do_fft(beam,temp_sky,maxsize,maxsize,-1)
   end if

! using no kernel is like convolving with a DiracDelta so the 
! FFT is just one everywhere...
   if (kernel.ne.'none') then

! like in grid_correct:
! these routines caused things to go wrong in a spectacular fashion - again
! a problem with the fft routine? Use new version instead
!      call grid_correct(kernel,dirty_map)

! ...and for the beam:
!      call beam_correct(kernel,beam)
      do_plot = .true.
      call kernel_correct(kernel)
   end if

! Normalisation
   bmax = beam(centre,centre)
   beam = beam/bmax
   dirty_map = dirty_map/bmax
     
   doplot = io_yesno('Do plots?','n',status)
   if (doplot) then
      write(*,*) 'Displaying beam'
      call display_map(beam,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
      write(*,*) 'Displaying map'
      call display_map(dirty_map,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
      write(*,*) 'Max at:',maxloc(dirty_map)

      if (io_yesno('Display only area within PB FWHM?','y',status)) then
         temp_sky = dirty_map
         call io_getr('FWHM (amin):','60.',width,status)
         width = width/2.
         do i = 1,maxsize
            do j = 1,maxsize
               r = sqrt(real(i-centre)**2+real(j-centre)**2)
               r = r*cellsize/60.
               if (r.gt.width) then
                  temp_sky(i,j) = 0.
                end if
            end do
         end do
         write(*,*) 'Displaying map'
         call display_map(temp_sky,maxsize,graph_adj,&
                          cellsize,cellsize,0D0,0D0) 
      end if
      if (nsrc.gt.0) then
         if (io_yesno('Overlay source list?','y',status)) then
            call pgsci(2)
            write(*,*) 'Overlaying ',nsrc,' sources.'
            do i = 1,nsrc
               x = real(src_x(i))*cellsize - real(maxsize/2-1)*cellsize
               y = real(src_y(i))*cellsize - real(maxsize/2-1)*cellsize
               call pgpt1(x,y,2)
            end do
            call pgsci(1)
         end if
      end if
   end if

   if (io_yesno('Histogram fluxes?','n',status)) then
      ok_data = .true.
      datmin = minval(dirty_map)
      datmax = maxval(dirty_map)
      call io_geti('Number of bins','1000',nbin,status)
      call map_histogram(dirty_map,maxsize,datmin,datmax,nbin)
   end if

! Initialise CLEANing maps
   resid_map = dirty_map
   clean_map = 0.d0
   comp_map = 0.d0

end subroutine map_vis

! -------------------------------------------------------------------------
! Subroutine which determines the pixel number corresponding to a diplacement
! in the u or v directions
subroutine decide(i,upt)

   use kind_def
   use sz_globals

   integer, intent(out)  ::  i
   real(dp), intent(in)  ::  upt
   real(dp)              ::  uu,uvcell
   integer               ::  centre

   centre = maxsize/2+1
   uvcell = 1./(maxsize*cellsize*sec2rad)
   uu  = upt/uvcell
  
   i = int(uu+sign(0.5d0,uu))
   i = i+centre

   ! pad the edges for convolution:
   if (i.le.convsize) i = -2
   if (i.ge.(maxsize-convsize)) i = maxsize+1

end subroutine decide

! --------------------------------------------------------------------------

subroutine wprojgridder(wmax,nwplanes,cutoff)

   use kind_def
   use sz_globals
   use w_proj
   use pswf_globals  
   use maths
   use kgraph
   use mikes_fft

   implicit none

   real(dp), intent(in)   :: wmax
   integer, intent(in)    :: nwplanes
   real(dp), intent(in)   :: cutoff
   integer    :: i
   integer,parameter  :: maxsupport = 10
   real(dp) :: spconvfunc(nplanes,csize,csize)
   character*1 :: chr1

   integer, parameter    :: npol = 1

   integer    :: cmap(n_samp,npol,nchan),wsupport1,wsupport2,wsupport3
   integer    :: cenw,pol,plane,wsize
   real(dp)   :: wdash,wvar,wt,nux,nuy
   real(dp)   :: cellx,celly,ccellx,ccelly
   real(dp)   :: x2,y2,r2,wre,wim,freq
   integer    :: ix,iy,iw,nx,ny,fracu,fracv
   complex,allocatable :: thisplane(:,:)
   real       :: work(2*maxsize)
   integer    :: job,iform
   real       :: temp(maxsize,maxsize)
   real(dp), allocatable :: ccfx(:),ccfy(:)

  
   cenw = (nwplanes+1)/2
   
! Run some checks...

   if (wmax.le.0.) then
      write(*,*) 'Baseline length must be greater than zero'
      wstatus = 1
      return
   end if
   if (nwplanes.le.0.) then
      write(*,*) 'Number of w planes must be greater than zero'
      wstatus = 1
      return
   end if
   if (((real(nwplanes)/2.) - int(real(nwplanes)/2.)).eq.0.) then 
      write(*,*) 'Number of w planes must be odd.'
      wstatus = 1
      return
   end if
   if (oversample.le.0) then
      write(*,*) 'Over sampling must be greater than zero.'
      wstatus = 1
      return
   end if
   if (cutoff.lt.0.) then
      write(*,*) 'Cutoff must be positive.'
      wstatus = 1
      return
   end if
   if (cutoff.gt.1.) then
      write(*,*) 'Cutoff must be less than 1.'
      wstatus = 1
      return
   end if
   if (maxsupport.le.0) then
      write(*,*) 'Maximum support must be greater than zero.'
      wstatus = 1
      return
   end if

   wsupport = 0
   wscale = wmax/(real(nwplanes-1)/2.)

! Decide where visibility samples lie in cube.(Basically a check on number of w planes being sufficient.)
   
   do i = 1,n_samp
      wdash = w(i)
      do chan = 1,nchan
         do pol = 1,npol
            freq = nu(chan)/obsfreq         ! scale w with channel
                                            ! need to check this doesn't conflict with m-m-o
            if (nwplanes.gt.1) then
               cmap(i,pol,chan) = cenw+nint((wdash*freq/wscale)+0.5)
            else
               cmap(i,pol,chan) = 0
            end if

            if (cmap(i,pol,chan).lt.0) then
               write(*,*) 'w scaling error 1: recommend allowing larger range of w.'
               call io_getc('Do you want to start over again?','y',chr1,status)
               if (chr1=='y') then
                  wstatus = 1
                  return
               end if
            end if

            if (cmap(i,pol,chan).ge.nwplanes) then
               write(*,*) 'w scaling error 2: recommend allowing larger range of w.'
               call io_getc('Do you want to start over again?','y',chr1,status)
               if (chr1=='y') then
                  wstatus = 1
                  return
               end if
            end if
            if (cmap(i,pol,chan).le.-1) then
               write(*,*) 'w scaling error 3: recommend allowing larger range of w.'
               call io_getc('Do you want to start over again?','y',chr1,status)
               if (chr1=='y') then
                  wstatus = 1
                  return
               end if
            end if

         end do
      end do
   end do

! Initialize convolution function in MAP plane:

   wsupport = 0

! cellx = 1./maxsize*uvcell = cellsize in rads
   cellx = cellsize*sec2rad
   celly = cellsize*sec2rad

! could be: nx = min(maxsupport,maxsize)
   nx = maxsize
   ny = maxsize

   qnx = nx/oversample   ! ie. maxsize=1024 -> 8; 7?
   qny = ny/oversample

! make sure support for PSWF is ODD otherwise centre will 
! be in wrong place and FFT will go screwy:
   if ((real(qnx)/2.- real(int(real(qnx)/2.))).eq.0.) then
      qnx = qnx-1
   end if
   if ((real(qny)/2.- real(int(real(qny)/2.))).eq.0.) then
      qny = qny-1
   end if

   if (qnx.lt.3) qnx = 3
   if (qny.lt.3) qny = 3

   allocate(ccfx(qnx))
   allocate(ccfy(qny))
   
! ccellx = maxsize*cellx/qnx
   ccellx = real(cellx*maxsize)/real(qnx)
   ccelly = real(celly*maxsize)/real(qny)
   
!   call pswfgridder(spconvfunc)
!   plane = 1  

   do ix = 1,qnx
      nux = abs(real(ix-((qnx+1)/2)))/real(qnx-1)/2.
      ccfx(ix) = 0.
      call pswf(nux,ccfx(ix))
      ccfx(ix) = ccfx(ix)/real(qnx)
   end do
   do iy = 1,qny
      nuy = abs(real(iy-((qny+1)/2)))/real(qny-1)/2.
      ccfy(iy) = 0.
      call pswf(nuy,ccfy(iy))
      ccfy(iy) = ccfy(iy)/real(qny)
   end do
   
! Start stepping through w planes, starting with the furthest out. We calculate the support for that plane and use it for all the others.
   do iw = 1,nwplanes

      allocate(thisplane(nx,ny))
      thisplane = 0.
      wvar = real(iw-cenw)*wscale
      do iy = 1,qny
         y2 = (real(iy-((qny+1)/2))*ccelly)**2
         do ix = 1,qnx
            x2 = (real(ix-((qnx+1)/2))*ccellx)**2
            r2 = x2+y2
            if (r2.lt.1.0) then
!wt = spconvfunc(plane,ix,iy)
               wt = ccfx(ix)*ccfy(iy)
               call wproj((ix-((qnx+1)/2))*ccellx,(iy-((qny+1)/2))*ccelly,&
                    wvar,wre,wim)
               thisplane((ix-((qnx+1)/2))+(nx/2)+1,(iy-((qny+1)/2))+(ny/2)+1) &
                    = cmplx(wt*wre,-1.*wt*wim)
               
            end if
         end do
      end do

! We now have the phase screen multiplied by the PSWF, sampled on larger cellsize, in IMAGE SPACE. Only the central qny,qny pixels are non-zero. We calculate the FFT to get the convolution in UV SPACE:
      
      
!      call display_map(dble(real(thisplane)),maxsize,graph_adj,cellsize,cellsize,0D0,0D0)

      job=1                      ! forward transform (+ve exponential)
      iform=1                    ! data are complex
      call makefft(nx,ny,thisplane,job,iform,work)

!      call display_map(dble(real(thisplane)),maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
!      call display_map(dble(abs(thisplane)),maxsize,graph_adj,cellsize,cellsize,0D0,0D0)


      if (wsupport==0) then
        
        do ix = 1,(nx/2)+1
           
           if (abs(thisplane(ix,ny/2+1)).gt.cutoff) then
              wsupport1 = iabs(ix-(nx/2+1))/oversample
              exit
           end if

        end do
        do ix = 1,(nx/2)+1
          
           if (abs(thisplane(ix,ix)).gt.cutoff) then
              wsupport2 = int(abs(real(int(1.414*real(ix))-(nx/2+1))))/oversample
              exit
           end if

        end do

        do ix = 1,(nx/2)+1
           
           if (nx==ny) then
              if (abs(thisplane(nx/2+1,ix)).gt.cutoff) then
                 wsupport3 = int(abs(real(ix-(ny/2+1)))/oversample)
                 exit
              end if
           end if
           
        end do
     end if
     
     wsupport = max(wsupport1,wsupport2)
     wsupport = max(wsupport,wsupport3)

     ! Some checks:
     if (wsupport.le.0) then
        write(*,*) 'Unable to determine support of convolution function.'
        wstatus = 1
        return
     elseif (wsupport*oversample.ge.(nx/2+1)) then
        write(*,*) 'Overflowing convolution function - increase maxsize or decrease oversample.'
        wstatus  = 1
        return
     end if
     
     wsize = 2*wsupport+1
     wcentre = (wsize+1)/2
     
     do fracu = 0,oversample
        do fracv = 0,oversample
           plane = fracu+oversample*(fracv+oversample*iw)+1
           convfunc(plane,:,:) = 0.
           
           do iy = -wsupport,wsupport
              do ix = -wsupport,wsupport
                 convfunc(plane,(ix+wcentre),(iy+wcentre)) = thisplane((ix*oversample+fracu+nx/2),(iy*oversample+fracv+ny/2))
              end do
           end do
           
        end do
     end do
     ! write(*,*) 'Done w-plane:',iw
     deallocate(thisplane)


  end do

  write(*,*) 'W-projection convolved with PSWF is now gridded into uvw-space.'

end subroutine wprojgridder

! ---------------------------------------------------------------------------

subroutine pswfgridder(spconvfunc)

! this is the Fourier domain function - the PSWF is its own FT but has the multiplying term (1-t^2) in front in Fourier space.

   use kind_def
   use sz_globals
   use pswf_globals
   use w_proj

   implicit none

   integer  :: ix,iy,i
   integer  :: fracv, fracu, plane
   real(dp) :: nux,nuy,fx,fy
   real(dp) :: spconvfunc(nplanes,csize,csize)
   real(dp) :: volume

   real(dp)     :: nn(1000),fun(1000)
     
  
   do fracv = 0,(oversample-1)
      do fracu = 0,(oversample-1)

         plane = fracu+(oversample*fracv)+1 
         spconvfunc(plane,:,:) = 0.
        
         do ix = 1,csize
            nux = abs(real(oversample*(ix-ccentre)+fracu))/real(support*oversample)
            fx = 0.
            call pswf(nux,fx)
            fx = fx*(1.-nux**2)

            do iy = 1,csize
               nuy = abs(real(oversample*(iy-ccentre)+fracv))/real(support*oversample)
               fy = 0.
               call pswf(nuy,fy)
               fy = fy*(1.-nuy**2)

               spconvfunc(plane,ix,iy) = fx*fy
            end do
         end do

      end do
   end do


   volume = sum(real(spconvfunc(1,:,:)))
   if (volume.eq.0.) write(*,*) 'Integral of convolution function is ZERO.'
!   if (volume.gt.0.) write(*,*) 'Volume of kernel:',volume
   spconvfunc = spconvfunc/volume ! Normalization
!   write(*,*) 'Normalized volume:',sum(spconvfunc(1,:,:))


end subroutine pswfgridder

! -------------------------------------------------------------------------
