
subroutine overlay_mosaic

   use sz_globals
   use kgraph

   implicit none


   integer           :: i,j,np,count
   real              :: sfwhm,spacing,count_deg,sigma
   integer           :: n_beam
   real              :: sensi(maxsize,maxsize)
   real              :: centre(2,100), position(2,100)
   real              :: distance, gaussian, radius
   integer           :: h,m,d,dm, new_h, new_d
   integer           :: new_m,new_dm
   real              :: s,ds,new_s,new_ds
   real              :: dec,ra,new_dec,new_ra
   real              :: distance_x, distance_y

   real              :: maxlev,minlev,tr(6),x1,x2,y1,y2
   integer           :: just ,ppost
   real              :: array(maxsize,maxsize)
   character*80      :: chr1
   integer  pgopen
   external pgopen


   graph_adj = .true.
   ppost = 0

   call io_getr('FWHM of primary beam(amin):','85',sfwhm,status)
   call io_getr('Centre spacing (amin):','50',spacing,status)

! convert to pixels:
   sfwhm = sfwhm*60/cellsize
! convert to pixels:
   spacing = spacing*60./cellsize
! sigma in pixels:
   sigma = sfwhm/(2*sqrt(2*log(2.0)))

   do i = 1,maxsize
      do j = 1,maxsize
         sensi(i,j) = 0.0
      end do
   end do

   call io_geti('Number of beams (1/3/5/7):','3',n_beam,status)


   select case(n_beam)      
   case(1)
      centre(1,1) = maxsize/2
      centre(2,1) = maxsize/2      
   case(3)
      centre(1,1) =  maxsize/2-spacing/2.0
      centre(2,1) =  maxsize/2-((spacing**2-(spacing/2.0)**2)**0.5)/2.0
      centre(1,2) =  maxsize/2+spacing/2.0
      centre(2,2) =  maxsize/2-((spacing**2-(spacing/2.0)**2)**0.5)/2.0
      centre(1,3) =  maxsize/2
      centre(2,3) =  maxsize/2+((spacing**2-(spacing/2.0)**2)**0.5)/2.0
   case(7)
      centre(1,1) =  maxsize/2
      centre(2,1) =  maxsize/2
      centre(1,2) =  maxsize/2-spacing
      centre(2,2) =  maxsize/2
      centre(1,3) =  maxsize/2+spacing
      centre(2,3) =  maxsize/2
      centre(1,4) =  maxsize/2-spacing/2.0
      centre(2,4) =  maxsize/2-((spacing**2-(spacing/2.0)**2)**0.5)
      centre(1,5) =  maxsize/2+spacing/2.0
      centre(2,5) =  maxsize/2-((spacing**2-(spacing/2.0)**2)**0.5)
      centre(1,6) =  maxsize/2-spacing/2.0
      centre(2,6) =  maxsize/2+((spacing**2-(spacing/2.0)**2)**0.5)
      centre(1,7) =  maxsize/2+spacing/2.0
      centre(2,7) =  maxsize/2+((spacing**2-(spacing/2.0)**2)**0.5)
   end select

   do np = 1, n_beam
      do i = 1,maxsize
         do j = 1,maxsize
            position(1,np) = i
            position(2,np) = j
! find the distance from the centre
            distance = ((position(1,np)-centre(1,np))**2+&
                        (position(2,np)-centre(2,np))**2)**0.5
! calculate value of beam at this position
            gaussian = exp((-(distance**2))/(2*sigma**2))
! add in contribution from this pointing to sensitivity map
            sensi(i,j) = (sensi(i,j)**2 + gaussian**2)**0.5 
         end do
      end do
   end do

! count how many pixels have sensitivities above 0.5
   count = 0

   do i = 1,maxsize
      do j = 1,maxsize
         if (sensi(i,j) > 0.5) then
            count = count + 1
         else
            count = count
         endif

      end do
   end do

   count_deg = count/3600.0

   write (*,*) 'coverage in square degrees = ', count_deg

! display original map:
   maxlev = maxval(rosat_map)
   minlev = 0.
   just = 1

   tr(1) = 0
   tr(4) = 0
   tr(2) = 1.0
   tr(6) = 1.0

   x1 = 1.0
   x2 = maxsize
   y1 = 1.0
   y2 = maxsize

   call pgenv(x1,x2,y1,y2,just,0)

   array = real(rosat_map)
   call pggray(array,maxsize,maxsize,1,maxsize,1,maxsize,&
        maxlev,minlev,tr)

   radius = sfwhm/2.
   call pgsci(1)
   call pgsfs(2)
   call pgslw(2)
   do np = 1, n_beam
      call pgpt1(centre(1,np),centre(2,np), 1)
      call pgsls(4)
      call pgcirc(centre(1,np),centre(2,np),radius)
   end do

   call io_getc('Do you want to write this to ppostscript?','n',chr1,status)
   if (chr1 == 'y') then
      i = pgopen('mosaic.ps/PS')
      call pgenv(x1,x2,y1,y2,just,0)
      call pggray(array,maxsize,maxsize,1,maxsize,1,maxsize,&
           maxlev,minlev,tr)
      call pgsci(0)
      call pgsfs(2)
      call pgslw(2)
      do np = 1, n_beam
         call pgpt1(centre(1,np),centre(2,np), 1)
         call pgsls(4)
         call pgcirc(centre(1,np),centre(2,np),radius)
      end do
      call pgclos()
      i = pgopen('/xwindow')
   end if

   call io_getc('Would you like to output your mosaic centres?','y',chr1,status)
   if (chr1 == 'y') then
      write(*,*) 'Central co-ordinate:'
      read(*,*) h,m,s,d,dm,ds
      write(*,*) 'Centre at:',h,m,s,'',d,dm,ds
      do i = 1, n_beam         
         dec = (ds/3600.)+(dm/60.)+d
         ra = (s/3600.)+(m/60.)+h
         ra = ra*15
! distance in pixels:         
         distance_x = maxsize/2 - centre(1,i)
         distance_y = maxsize/2 - centre(2,i)
! distance in degrees:
         distance_x = distance_x*15./3600.
         distance_y = distance_y*15./3600.
! centre's are in pixels
         new_ra = (distance_x + ra)/15.
         new_h = int(new_ra)
         new_m = int((new_ra - new_h)*60.)
         new_s = ((new_ra - new_h) - (real(new_m)/60.))*3600.

         new_dec = distance_y + dec
         new_d  = int(new_dec)
         new_dm = int((new_dec - new_d)*60.)
         new_ds = ((new_dec - new_d) - (real(new_dm)/60.))*3600. 

         if (new_s .le. 0.01) then
            new_s = 0.
         end if
         
         if (new_ds .le. 0.01) then
            new_ds = 0.
         end if

         write(*,*) new_h,new_m,new_s,'',new_d,new_dm,new_ds
      end do
   end if
   call pgsci(1)
   call pgsls(1)

end subroutine overlay_mosaic


