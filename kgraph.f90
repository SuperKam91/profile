module kgraph

   use kind_def
   use sz_globals
   use file_handling

   implicit none

   integer, save :: BLCX,BLCY   ! The position of the graph and the size of 
   integer, save :: PLTSZ       ! the map after zooming........

   interface display_map
      module procedure display_map
   end interface

   interface superpose_cntr
      module procedure superpose_cntr
   end interface


! profile plotting routines. Allows:
! 1) plotting an equal spaced array, with autoscaled or defined axes
! 2) plotting 1 array against another, with autoscaled or defined axes
   interface display_profile
      module procedure display_profile_lin
      module procedure display_profile_2d
      module procedure display_def_profile_lin
      module procedure display_def_profile_2d
   end interface

   public :: get_map_rectangle
   public :: display_points
   public :: plot_uv
   public :: plot_uvhist
   public :: map_histogram

   private :: plot_map

contains

!****************************************************************************

   subroutine display_map(map,mapsize,interact,pix_xd,pix_yd,cent_xd,cent_yd)

      use kind_def

      implicit none

      real(kind=dp), dimension(mapsize,mapsize), intent(in) :: map
      real, allocatable, dimension(:) :: cont
      integer, intent(in) :: mapsize
      logical, intent(in) :: interact
      integer :: blc_x, blc_y, plotsize, ncont, status, i,j
      character(len=fname_len) :: ps_file
      character*3 :: switch
      real(kind=dp) :: y_parm, x_parm
      real :: x, y, maxlev, minlev, delta, centrelev, startlev, inclev
      real(kind=dp), intent(in) :: pix_xd, pix_yd, cent_xd, cent_yd
      real :: pix_x, pix_y, cent_x, cent_y
      real(kind=dp) :: max, min
      real :: x1, x2, y1, y2, x_world, y_world
      character*1 :: ch, char1
      integer :: istat, l1, just
      logical :: circle
      integer :: pgopen
      external pgopen
      integer :: chr_lenb
      external chr_lenb

      
! Set up initial plot
      plotsize = mapsize
!      blc_x = 0
!      blc_y = 0

      circle = .false.
      max = maxval(map)
      min = minval(map)

      maxlev = max
      minlev = min
      delta = max-min
      switch = 'gre'
      pix_x = pix_xd
      pix_y = pix_yd
      cent_x = cent_xd
      cent_y = cent_yd
      blc_x = 1
      blc_y = 1
      BLCX = blc_x
      BLCY = blc_y
      if (pix_x.eq.pix_y) then
         just = 1
      else
         just = 0
      end if

! Do initial plot
      call pgvsize(2., 6.5, 1., 5.5)
      x1 = (blc_x-mapsize/2-1.5)*pix_x+cent_x
      x2 = (blc_x-mapsize/2+plotsize-1.5)*pix_x+cent_x
      y1 = (blc_y-mapsize/2-1.5)*pix_y+cent_y
      y2 = (blc_y-mapsize/2+plotsize-1.5)*pix_y+cent_y
      call pgenv(x1,x2,y1,y2,just,0)
      call pgask(.false.)

      if (.not.interact) then 
! A subtle change in default behaviour when plotting to a postscript file
!         switch = 'cnt'
!         if (allocated(cont)) then
!            deallocate(cont)
!         end if
!         ncont = 10
!         allocate (cont(ncont))
!         cont(1) = minlev
!         inclev = (maxlev-minlev)/(ncont+1)
!         do i = 2,ncont
!            cont(i) = cont(i-1)+inclev
!            if (cont(i)==0.0) then
!               cont(i) = cont(i)+inclev
!            end if
!         end do
      end if

      
      call plot_map(map, mapsize, switch, blc_x, blc_y, plotsize, &
           ncont,cont,maxlev,minlev,pix_x,pix_y,cent_x,cent_y)
      
! Allow interactive change of grey scale, contours, zooming?
      if (interact) then
         write(*,*) 'Type h for help'

! Allow user to change plot
         plot_loop : do

! Get input from cursor
            call pgcurse(x_world,y_world,ch)
            x = (x_world-cent_x)/pix_x+mapsize/2+1.5-blc_x
            y = (y_world-cent_y)/pix_y+mapsize/2+1.5-blc_y

! Quit
            if ((ch=='q').or.(ch=='Q')) then
               exit plot_loop

! Reset plot
            else if ((ch=='r').or.(ch=='R')) then
               blc_x = 1
               blc_y = 1
               maxlev = max
               minlev = min         
               plotsize = mapsize
               switch = 'gre'
               x1 = (blc_x-mapsize/2-1.5)*pix_x+cent_x
               x2 = (blc_x-mapsize/2+plotsize-1.5)*pix_x+cent_x
               y1 = (blc_y-mapsize/2-1.5)*pix_y+cent_y
               y2 = (blc_y-mapsize/2+plotsize-1.5)*pix_y+cent_y
               call pgenv(x1,x2,y1,y2,just,0)

! Fiddle with grey scale        
            else if (ch=='A') then
               x_parm = x/plotsize
               y_parm = y/plotsize
               y_parm = (y_parm-0.5)*2.0
               x_parm = exp((x_parm-0.5)*2.0)
               centrelev = (max+min)/2.0+y_parm*delta
               maxlev = centrelev+x_parm*delta/2.0
               minlev = centrelev-x_parm*delta/2.0

! Zoom in
            else if ((ch=='z').or.(ch=='Z').or.(ch=='D')) then
               blc_x = blc_x+int(x)-plotsize/4
               blc_y = blc_y+int(y)-plotsize/4
               BLCX = blc_x
               BLCY = blc_y
               plotsize = plotsize/2
               x1 = (blc_x-mapsize/2-1.5)*pix_x+cent_x
               x2 = (blc_x-mapsize/2+plotsize-1.5)*pix_x+cent_x
               y1 = (blc_y-mapsize/2-1.5)*pix_y+cent_y
               y2 = (blc_y-mapsize/2+plotsize-1.5)*pix_y+cent_y
               call pgenv(x1,x2,y1,y2,just,0)

! Centre plot 
            else if (ch==' ') then
!               write(*,*) blc_x,int(x)
               blc_x = blc_x+int(x)-plotsize/2
               blc_y = blc_y+int(y)-plotsize/2
!               write(*,*) blc_x
               BLCX = blc_x
               BLCY = blc_y
               x1 = (blc_x-mapsize/2-1.5)*pix_x+cent_x
               x2 = (blc_x-mapsize/2+plotsize-1.5)*pix_x+cent_x
               y1 = (blc_y-mapsize/2-1.5)*pix_y+cent_y
               y2 = (blc_y-mapsize/2+plotsize-1.5)*pix_y+cent_y
               call pgenv(x1,x2,y1,y2,just,0)

! Zoom out
            else if ((ch=='X').or.(ch=='o').or.(ch=='O')) then
               blc_x = blc_x-plotsize/2
               blc_y = blc_y-plotsize/2
               BLCX = blc_x
               BLCY = blc_y
               plotsize = plotsize*2
               x1 = (blc_x-mapsize/2-1.5)*pix_x+cent_x
               x2 = (blc_x-mapsize/2+plotsize-1.5)*pix_x+cent_x
               y1 = (blc_y-mapsize/2-1.5)*pix_y+cent_y
               y2 = (blc_y-mapsize/2+plotsize-1.5)*pix_y+cent_y
               call pgenv(x1,x2,y1,y2,just,0)

! Draw circle
            elseif (ch=='b') then
               circle=.true.

! Get value
            else if ((ch=='V').or.(ch=='v')) then
               write(*,*)
               write(*,*) 'Pixel ',BLCX+int(x),' ',BLCY+int(y)
               write(*,*) 'Value at cursor is ', &
                           map(BLCX+int(x),BLCY+int(y))
               write(*,*)
               
! Help 
            else if ((ch=='h').or.(ch=='H')) then
               write(*,*) 'Help'
               write(*,*)
               write(*,*) 'r        reset plot'
               write(*,*) 'q        quit'
               write(*,*) 'g        grey scale'
               write(*,*) 'c        contours'
               write(*,*) 'v        get value at cursor'
               write(*,*) 'p        write postscript file'
               write(*,*) '<SPACE>  centre plot'
               write(*,*) 'Mouse 1  change transfer function'
               write(*,*) 'Mouse 2  zoom in'
               write(*,*) 'Mouse 3  zoom out'

! Select contours
            else if ((ch=='c').or.(ch=='C')) then
               call io_getc('Grey scale as well as contours (y/n):','n',&
                    & char1,status)
               if ((char1=='y').or.(char1=='Y')) then
                  switch = 'bth'
               else
                  switch = 'cnt'
               end if

! Choose which type of contours and allocate contour array
               get_contours : do
                  if (allocated(cont)) then
                     deallocate(cont)
                  end if
                  call io_getc('Standard, linear, exponential, or defined &
                       &contours (l,e,d)','s',char1,status) 

! Exponential contours
                  if ((char1=='e').or.(char1=='E')) then
                     if ((minlev <= 0.0).or.((maxlev/minlev) > 33000.0)) then
                        minlev = maxlev/33000.0
                     end if
                     ncont = int(2*log(maxlev/minlev)/log(2.0))
                     allocate (cont(ncont))
                     cont(1) = minlev
                     do i = 2, ncont
                        cont(i) = cont(i-1)*1.4142135
                     end do

! User defined contours
                  else if ((char1=='d').or.(char1=='D')) then
                     call io_geti('Number of contours','10',ncont,status)
                     allocate(cont(ncont)) 
                    do i= 1,ncont
                       call io_getr('Contour','1.0',cont(i),status)
                     end do

! Linear contours
                  else if ((char1=='l').or.(char1=='L')) then
                     call io_getr('Start level','1.0',startlev,status)
                     call io_getr('Increment level','1.0',inclev,status)
                     call io_geti('Number of contours','10',ncont,status)
                     allocate (cont(ncont))
                     cont(1) = startlev
                     do i = 2,ncont
                        cont(i) = cont(i-1)+inclev
                        if (cont(i)==0.0) then
                           cont(i) = cont(i)+inclev
                        end if
                     end do

! Standard contours
                  else
                     call io_geti('Number of contours','10',ncont,status)
                     allocate (cont(ncont))
                     cont(1) = minlev
                     inclev = (maxlev-minlev)/(ncont+1)
                     do i = 2,ncont
                        cont(i) = cont(i-1)+inclev
                        if (cont(i)==0.0) then
                           cont(i) = cont(i)+inclev
                        end if
                     end do
                  end if
                  exit get_contours
               end do get_contours
! Set max and min levels
            else if (ch=='s') then
               call io_getr('Max level:','0.3',maxlev,status)
               call io_getr('Min level:','0.0',minlev,status)

! Choose Grey scale plot
            else if ((ch=='g').or.(ch=='G')) then
               call pgpage
               switch = 'gre'

! Write postsrcipt file
            else if ((ch=='p').or.(ch=='P')) then
               prompt = 'Filename for postscript plot'
               which_dir = postscript_dir
               rootname = 'plot'
               extn = '.ps'
               call get_write_filename(ps_file)
               
               l1 = chr_lenb(ps_file)
               ps_file = ps_file(1:l1)//'/CPS'
               istat = pgopen(ps_file)
               if (istat.le.0) then
                  write(*,*) 'Error opening postscript file'
                  write(*,*) 'iostat:',istat
                  exit plot_loop
               end if
               call pgslct(istat)
               call pgvsize(2., 6.5, 1., 5.5)
               x1 = (blc_x-mapsize/2-1.5)*pix_x+cent_x
               x2 = (blc_x-mapsize/2+plotsize-1.5)*pix_x+cent_x
               y1 = (blc_y-mapsize/2-1.5)*pix_y+cent_y
               y2 = (blc_y-mapsize/2+plotsize-1.5)*pix_y+cent_y
               call pgenv(x1,x2,y1,y2,just,0)
               call pgask(.false.)
               call plot_map(map, mapsize, switch, blc_x, blc_y, plotsize, &
                    ncont,cont,maxlev,minlev,pix_x,pix_y,cent_x,cent_y)
               if (circle) then
                  call pgsfs(2)
                  call pgsci(0)
                  call pgcirc(0.,0.,600.)
                  call pgcirc(0.,0.,1200.)
               end if
               call pgclos
               call pgslct(pl_dev1)

            end if

! Replot map
            call plot_map(map, mapsize, switch, blc_x, blc_y, plotsize, &
                 ncont,cont,maxlev,minlev,pix_x,pix_y,cent_x,cent_y)
            if (circle) then
               call pgsfs(2)
               call pgsci(1)
               call pgcirc(0.,0.,600.)
               call pgcirc(0.,0.,1200.)
            end if

         end do plot_loop

! End of interactive adjustments
      end if

! Finish off
      call pgask(.true.)

      BLCX = blc_x
      BLCY = blc_y
      PLTSZ = plotsize

      if (allocated(cont)) then
         deallocate(cont)
      end if

   end subroutine display_map

! ************************************************************************

! Do the 2-d plot
   subroutine plot_map(map, mapsize, switch, blc_x, blc_y, plotsize, &
        ncont,cont,maxlev,minlev,pix_x,pix_y,cent_x,cent_y) 

      use kind_def

      implicit none

      real(kind=dp), dimension(mapsize,mapsize), intent(in) :: map
      real, dimension(ncont), intent(in) :: cont
      integer, intent(in) :: mapsize, plotsize, ncont, blc_x, blc_y
      character*3, intent(in) :: switch
      integer :: i, j, j1, i1
      real :: pix_x, pix_y,cent_x,cent_y
      real, intent(in) :: maxlev,minlev
      real, dimension(plotsize,plotsize) :: array
      logical, dimension(plotsize,plotsize) :: exist_array
      real, dimension(6) :: tr

      tr(1) = (blc_x-mapsize/2-2)*pix_x+cent_x
      tr(4) = (blc_y-mapsize/2-2)*pix_y+cent_y

      tr(2) = pix_x
      tr(3) = 0.
      tr(6) = pix_y
      tr(5) = 0.

! Find where array to be plotted overlaps map
      do i = 1, plotsize
         i1 = blc_y+i-1
         do j = 1,plotsize
            j1 = blc_x+j-1
            exist_array(j,i) = ((i1 <= mapsize).and.(i1 >= 1).and. &
                 (j1 <= mapsize).and.(j1 >= 1))
         end do
      end do

! Read in array to be plotted
      where (exist_array)
         array = map(blc_x:blc_x+plotsize-1,blc_y:blc_y+plotsize-1)
         elsewhere
         array = minlev
      end where

! Redraw from scratch if using contours
      if ((switch=='cnt').or.(switch=='bth')) then
         call pgpage
         call pgbox('BCNITS', 0., 0., 'BCNITS', 0., 0.)
      end if

! Draw grey scale
      if ((switch=='gre').or.(switch=='bth')) then
         if (verbose) then
            write(*,*) 'max = ', maxlev,' min = ',minlev
         end if
         call pggray(array,plotsize,plotsize,1,plotsize,1,plotsize,&
              maxlev,minlev,tr)
      end if

! Draw contours
      if ((switch=='cnt').or.(switch=='bth').or.(switch=='sup')) then
         write(*,*) 'contour levels '
         do i = 1,ncont
            write(*,*) cont(i)
         end do
         call pgcont(array,plotsize,plotsize,1,plotsize,1,plotsize,cont,&
              & ncont,tr)
      end if

   end subroutine plot_map

! ************************************************************************

   subroutine get_map_rectangle(map_rectangle,mapsize,pix_xd,pix_yd,&
                                cent_xd,cent_yd)

! read box co-ordinates from the pgplot device

      implicit none

      real, dimension(4), intent(out) :: map_rectangle
      real, dimension(4) :: orig_rect
      real(kind=dp), intent(in) :: pix_xd, pix_yd, cent_xd, cent_yd
      real :: temp_pt
      real :: pix_x, pix_y, cent_x, cent_y
      integer :: mapsize
      character*1 :: cursor

      pix_x = pix_xd
      pix_y = pix_yd
      cent_x = cent_xd
      cent_y = cent_yd

      cursor = ' '
      do while (cursor.ne.'A')
         call pgband(0,1,0,0,map_rectangle(1),map_rectangle(2),cursor)
         call pgband(2,0,map_rectangle(1),map_rectangle(2),&
              map_rectangle(3),map_rectangle(4),cursor)
      end do

      orig_rect = map_rectangle

!      write(*,*) 'x1:', orig_rect(1)
!      write(*,*) 'y1:', orig_rect(2)
!      write(*,*) 'x2:', orig_rect(3)
!      write(*,*) 'y2:', orig_rect(4)

!      write(*,*) 'cent_x:', cent_x
!      write(*,*) 'pix_x:',pix_x
!      write(*,*) 'mapsize:', mapsize

! change from world coordinates back to array positions
      pix_x = cellsize
      pix_y = cellsize
      map_rectangle(1) = (map_rectangle(1)-cent_x)/pix_x+mapsize/2
      map_rectangle(2) = (map_rectangle(2)-cent_y)/pix_y+mapsize/2
      map_rectangle(3) = (map_rectangle(3)-cent_x)/pix_x+mapsize/2
      map_rectangle(4) = (map_rectangle(4)-cent_y)/pix_y+mapsize/2

! if necessary rearrange map_rectangle so that it is blc_x, blc_y, trc_x, trc_y
      if (map_rectangle(1)>map_rectangle(3)) then
         temp_pt = map_rectangle(1)
         map_rectangle(1) = map_rectangle(3)
         map_rectangle(3) = temp_pt
      end if
      if (map_rectangle(2)>map_rectangle(4)) then
         temp_pt = map_rectangle(2)
         map_rectangle(2) = map_rectangle(4)
         map_rectangle(4) = temp_pt
      end if

!      write(*,*) 'blc = ',map_rectangle(1)+BLCX,map_rectangle(2)+BLCY
!      write(*,*) 'trc = ',map_rectangle(3)+BLCX,map_rectangle(4)+BLCY
      write(*,*) 'blc = ',map_rectangle(1),map_rectangle(2)
      write(*,*) 'trc = ',map_rectangle(3),map_rectangle(4)

      call pgsfs(2)
      call pgrect(orig_rect(1),orig_rect(3),&
           & orig_rect(2),orig_rect(4))
      call pgsfs(1)

!      map_rectangle(1) = map_rectangle(1)+BLCX
!      map_rectangle(2) = map_rectangle(2)+BLCY
!      map_rectangle(3) = map_rectangle(3)+BLCX
!      map_rectangle(4) = map_rectangle(4)+BLCY

      
   end subroutine get_map_rectangle

! ************************************************************************

! AMS 03/03/2004
! define a circular area on the map rather than a rectangle

   subroutine get_map_circle(circle_radius,circle_cent,mapsize,&
        pix_xd,pix_yd,cent_xd,cent_yd)
    
      implicit none

      real,dimension(4) :: map_line
      real, dimension(2) :: orig_cent
      real :: orig_radius
      real(kind=dp), intent(in):: pix_xd, pix_yd,cent_xd,cent_yd
      real :: temp_pt
      real :: pix_x,pix_y,cent_x,cent_y
      real, dimension(2), intent(out) :: circle_cent
      real, intent(out) :: circle_radius
      integer :: mapsize
      character*1 :: cursor

      pix_x = pix_xd
      pix_y = pix_yd
      cent_x = cent_xd
      cent_y = cent_yd

      cursor = '  '

      write(*,*) 'Draw diameter of circle - left to right'

      do while (cursor .ne. 'A')
         call pgband(0,1,0,0,map_line(1), map_line(2), cursor)
         call pgband(1,0,map_line(1), map_line(2), map_line(3),&
                     map_line(4), cursor)
      end do

! define circle in world co-ordinates

      circle_cent(1) = (map_line(3) - map_line(1))/2 + map_line(1)
      circle_cent(2) = (map_line(4) - map_line(2))/2 + map_line(2)
      circle_radius = ((map_line(3)-circle_cent(1))**2&
           + (map_line(4)-circle_cent(2))**2)**0.5

! debug

!      write(*,*) 'map line 3:', map_line(3), 'mapline 1', map_line(1)
!      write(*,*) 'x centre:', circle_cent(1)
!      write(*,*) 'y centre:', circle_cent(2)
      write(*,*) 'radius:', circle_radius

      orig_cent = circle_cent
      orig_radius = circle_radius

      write(*,*) 'central x co-ordinate:', orig_cent(1)
      write(*,*) 'central y co-ordinate', orig_cent(2)
      
      write(*,*) 'radius:', orig_radius


! change from world co-ords to array positions

      pix_x = cellsize
      pix_y = cellsize
      map_line(1) = (map_line(1)-cent_x)/pix_x+mapsize/2
      map_line(2) = (map_line(2)-cent_y)/pix_y+mapsize/2
      map_line(3) = (map_line(3)-cent_x)/pix_x+mapsize/2
      map_line(4) = (map_line(4)-cent_y)/pix_y+mapsize/2

! define circle in array

      circle_cent(1) = (map_line(3) - map_line(1))/2 + map_line(1)
      circle_cent(2) = (map_line(4) - map_line(2))/2 + map_line(2)
      circle_radius = ((map_line(3)-circle_cent(1))**2&
           + (map_line(4) - circle_cent(2))**2)**0.5

! debug

!      write(*,*) 'circle centre x:', circle_cent(1)
!      write(*,*) 'circle centre y:', circle_cent(2)
!      write(*,*) 'circle radius:', circle_radius

! draw circle

      call pgsfs(2)
      call pgcirc(orig_cent(1),orig_cent(2),orig_radius)
      call pgsfs(1)

   end subroutine get_map_circle


!**************************************************************************

! Display equally spaced array of points on pgplot device, determining best
! scales
   subroutine display_profile_lin(profile,profsize,x_step)

      use kind_def

      implicit none
      integer, intent(in) :: profsize
      real(kind=dp), dimension(profsize), intent(in) :: profile
      real(kind=dp), intent(in) :: x_step
      real, dimension(profsize) :: y_points, x_points
      real :: min_y, max_y, min_x, max_x
      integer :: i, status, pgopen
      character*1 :: ch, y, xlbl, ylbl, toplbl

      external pgopen

! Prepare x and y arrays
      x_points = (/(i,i=1,profsize)/)
      x_points = (x_points-1)*x_step
      y_points = profile

! Find max and mins
      min_y = minval(y_points)
      max_y = maxval(y_points)   
      min_x = minval(x_points)
      max_x = maxval(x_points)

! Adjust y scale if necessary
      if ((min_y>0.0).and.(max_y>2*min_y)) then 
         min_y = 0.0
      end if

      if ((max_y<0.0).and.(min_y<2*max_y)) then 
         max_y = 0.0
      end if

      if (max_y==min_y) then
         if (max_y==0.0) then
            max_y = 0.05
            min_y = -0.05
         else if (max_y>0.0) then
            max_y = 1.05*max_y
            min_y = 0.95*min_y
         else
            max_y = 0.95*max_y
            min_y = 1.05*min_y
         end if
      end if

      if (ch=='y') then
! Plot profile
         call pgvsize(1.5,8.5,1.5,6.0)
         call pgwindow(min_x,max_x,min_y,max_y)
         call pgbox('BCLNITS', 0., 0, 'BCNITS', 0., 0)
         call pglab('r(arcsec)','Jy/bm','Profile')
         call pgline(profsize,x_points,y_points)
         
      else
         
         call pgvsize(1.5,8.5,1.5,6.0)
         call pgwindow(min_x,max_x,min_y,max_y)
         call pgbox('BCNITS', 0., 0, 'BCNITS', 0., 0)
         call pgline(profsize,x_points,y_points)
         
! Plot 0 line if necessary
         if (min_y.lt.0.0) then
            call pgmove(min_x,0.0)
            call pgdraw(max_x,0.0)
         end if

      end if

   end subroutine display_profile_lin

! ***************************************************************************

! Display profile of 2 arrays of points on pgplot device, determining best 
! scale
   subroutine display_profile_2d(x_profile,y_profile,profsize,adj)

      use kind_def

      implicit none
      integer, intent(in) :: profsize
      real(kind=dp), dimension(profsize), intent(in) :: x_profile, y_profile
      logical, intent(in) :: adj
      real, dimension(profsize) :: y_points, x_points
      real :: min_y, max_y, min_x, max_x
      integer :: i

! Prepare x and y arrays
      x_points = x_profile
      y_points = y_profile

! Find max and mins
      min_y = minval(y_points)
      max_y = maxval(y_points)   
      min_x = minval(x_points)
      max_x = maxval(x_points)

! Adjust scales?
      if (adj) then 

! Adjust y scale if necessary
         if ((min_y>0.0).and.(max_y>2*min_y)) then 
            min_y = 0.0
         end if

         if ((max_y<0.0).and.(min_y<2*max_y)) then 
            max_y = 0.0
         end if

         if (max_y==min_y) then
            if (max_y==0.0) then
               max_y = 0.05
               min_y = -0.05
            else if (max_y>0.0) then
               max_y = 1.05*max_y
               min_y = 0.95*min_y
            else
               max_y = 0.95*max_y
               min_y = 1.05*min_y
            end if
         end if

! Adjust x scale if necessary
         if ((min_x>0.0).and.(max_x>2*min_x)) then 
            min_x = 0.0
         end if

         if ((max_x<0.0).and.(min_x<2*max_x)) then 
            max_x = 0.0
         end if

         if (max_x==min_x) then
            if (max_x==0.0) then
               max_x = 0.05
               min_x = -0.05
            else if (max_x>0.0) then
               max_x = 1.05*max_x
               min_x = 0.95*min_x
            else
               max_x = 0.95*max_x
               min_x = 1.05*min_x
            end if
         end if
      end if

! Plot profile
      call pgvsize(1.5,8.5,1.5,6.0)
      call pgwindow(min_x,max_x,min_y,max_y)
      call pgbox('BCNITS', 0., 0, 'BCNITS', 0., 0)
      call pgline(profsize,x_points,y_points)

! Plot 0 line if necessary
      if (min_y.lt.0.0) then
         call pgmove(min_x,0.0)
         call pgdraw(max_x,0.0)
      end if

   end subroutine display_profile_2d

! **************************************************************************

! Display 1-d plot of points on pgplot device
   subroutine display_def_profile_lin(profile,profsize,x_step,&
        & x1,x2,y1,y2)

      use kind_def

      implicit none
      integer, intent(in) :: profsize
      real(kind=dp), dimension(profsize), intent(in) :: profile
      real(kind=dp), intent(in) :: x_step
      real, dimension(profsize) :: y_points, x_points
      real(kind=dp), intent(in) :: x1,x2,y1,y2
      real :: min_y, max_y, min_x, max_x
      integer :: i

! Prepare x and y arrays
      x_points = (/(i,i=1,profsize)/)
      x_points = (x_points-1)*x_step
      y_points = profile

! Prepare max and mins
      min_y = y1
      max_y = y2   
      min_x = x1
      max_x = x2

! Plot profile
      call pgvsize(1.5,8.5,1.5,6.0)
      call pgwindow(min_x,max_x,min_y,max_y)
      call pgbox('BCNITS', 0., 0, 'BCNITS', 0., 0)
      call pgline(profsize,x_points,y_points)

! Plot 0 line if necessary
      if (min_y.lt.0.0) then
         call pgmove(min_x,0.0)
         call pgdraw(max_x,0.0)
      end if

   end subroutine display_def_profile_lin

! ***************************************************************************

! Display 2-d plot of points on pgplot device
   subroutine display_def_profile_2d(x_profile,y_profile,profsize,x1,x2,y1,y2)

      use kind_def

      implicit none
      integer, intent(in) :: profsize
      real(kind=dp), dimension(profsize), intent(in) :: x_profile, y_profile
      real(kind=dp), intent(in) :: x1,x2,y1,y2
      real, dimension(profsize) :: y_points, x_points
      real :: min_y, max_y, min_x, max_x
      integer :: i

! Prepare x and y arrays
      x_points = x_profile
      y_points = y_profile

! Prepare max and mins
      min_y = y1
      max_y = y2   
      min_x = x1
      max_x = x2

! Plot profile
      call pgvsize(1.5,8.5,1.5,6.0)
      call pgwindow(min_x,max_x,min_y,max_y)
      call pgbox('BCNITS', 0., 0, 'BCNITS', 0., 0)
      call pgline(profsize,x_points,y_points)

! Plot 0 line if necessary
      if (min_y.lt.0.0) then
         call pgmove(min_x,0.0)
         call pgdraw(max_x,0.0)
      end if

   end subroutine display_def_profile_2d

! **************************************************************************

   subroutine display_points(x_points,y_points_1,y_points_2,y_rms,use_point,&
        & do_err,n_points,n_use,col1,col2)

      use kind_def

      implicit none
      integer, intent(in) :: n_points, n_use
      logical, intent(in) :: do_err
      real(kind=dp), dimension(n_points), intent(in) :: x_points, &
           & y_points_1, y_points_2, y_rms
      logical, dimension(n_points), intent(in) :: use_point
      integer, intent(in) :: col1, col2
      real, dimension(n_use) :: x, y1, y2, err
      real :: max_x, min_x, max_y1, min_y1, max_y
      real :: min_y, max_y2, min_y2, max_err

! Write arrays inyo single precision ones for plotting
      x = pack(x_points,mask=use_point)
      y1 = pack(y_points_1,mask=use_point)
      y2 = pack(y_points_2,mask=use_point)
      err = pack(y_rms,mask=use_point)

! Determine what axes to use
!      max_x = maxval(x,mask=use_point)
!      min_x = minval(x,mask=use_point)
!      max_y1 = maxval(y1,mask=use_point)
!      min_y1 = minval(y1,mask=use_point)
!      max_y2 = maxval(y2,mask=use_point)
!      min_y2 = minval(y2,mask=use_point)
!      max_err = maxval(err,mask=use_point)
      max_x = maxval(x)
      min_x = minval(x)
      max_y1 = maxval(y1)
      min_y1 = minval(y1)
      if (do_err) then
         max_y2 = maxval(y2+err)
         min_y2 = minval(y2-err)
      else
         max_y2 = maxval(y2)
         min_y2 = minval(y2)
      end if

      if (min_y1.lt.min_y2) then
         min_y = min_y1
      else
         min_y = min_y2
      end if
      if (max_y1.gt.max_y2) then
         max_y = max_y1
      else
         max_y = max_y2
      end if

! Have graphs going to zero if all points are positive
      if (min_x.gt.0) min_x = 0.0
      if (min_y.gt.0) min_y = 0.0

! Have graphs going to zero if all points are negative
      if (max_x.lt.0) max_x = 0.0
      if (max_y.lt.0) max_y = 0.0

! Leave a little room to right of points
      max_x = max_x*1.05

      call pgask(.true.)
      call pgpage
      call pgvsize(1.5,8.5,1.5,6.0)
      call pgsci(plot_col)
      call pgwindow(min_x,max_x,min_y,max_y)
      call pgbox('BCNITS', 0., 0, 'BCNITS', 0., 0)

      call pgsci(col2)
      call pgpt(n_use,x,y2,-4)
      if (do_err) then
         call pgerry(n_use,x,y2+err,y2-err,0)
      end if
      call pgsci(col1)
      call pgpt(n_use,x,y1,9)
      call pgline(n_use,x,y1)
      call pgsci(plot_col)

   end subroutine display_points

! ***************************************************************************

   subroutine plot_uv(max_basel,n_samp,ha_start,ha_inc,dec,write_vis)

      implicit none

      integer :: n_samp, i, j, k, l, m, status
      real :: max_basel, x, y, z, minpos, maxpos, diam
      real(kind=dp), dimension(n_antennas,n_antennas,n_samp) :: u_arr,v_arr
      real(kind=dp), dimension(n_antennas,n_antennas,n_samp) :: w_arr
      integer, dimension(n_antennas,n_antennas,n_samp) :: col 
      real(kind=dp) :: ha_start, ha_inc
      real(kind=dp), dimension(n_antennas,n_antennas) :: u,v,w
      real    :: dum_re, dum_im, noise
      character(len=fname_len) :: filename
      integer, dimension(n_antennas,n_antennas) :: shadow
      integer :: n_shadow, pgopen, iunit
      real(kind=dp) :: ha, dec
      logical :: write_vis, ex
      logical :: io_yesno
      external pgopen, io_yesno

      maxpos = max_basel*obsfreq/const_c
      minpos = -max_basel*obsfreq/const_c

      dum_re = 0.0
      dum_im = 0.0
      noise = 0.01

      col = 1

      call pgsci(plot_col)
      call pgslw(line_width)
      call pgsch(char_hgt)
      call pgvsize(2., 6.5, 1., 5.5)
      call pgenv(minpos,maxpos,minpos,maxpos,1,1)
      call pgask(.false.)
      
! Calculate u, v and check for shadowing
      n_shadow = 0
      do k = 1, n_samp
         ha = (ha_start+(k-1)*ha_inc)
         call calc_uvw(ha,dec,u,v,w,shadow)
         u_arr(:,:,k) = u
         v_arr(:,:,k) = v
         w_arr(:,:,k) = w
         col(:,:,k) = shadow
         n_shadow = n_shadow+count(shadow(:,:).ne.1)
      end do

      write(*,*) n_shadow, ' out of ', n_samp*n_antennas*(n_antennas-1), &
                 ' visibilities to be flagged due to shadowing'

! Plot circle indictating geometrical shadowing
      diam = dish_diameter*obsfreq/const_c
      call pgsfs(2)
      call pgcirc(0.0,0.0,diam)

! Open .vis file
      if (write_vis) then

! Get filename
         prompt = 'vis template filename'
         if (wis) then
            extn = '.wis'
         else
            extn = '.vis'
         end if
         which_dir = profile_data_dir
         rootname = 'default'
         call get_write_filename(filename)

         call io_nxtlun(iunit,status)
         open (iunit,file=filename,form='formatted')
      end if

! Plot uv points
      l=0
      do i = 1, n_antennas
         do j = 1, n_antennas
            if (i.ne.j) then
               do k = 1, n_samp
                  do m = 1, nchan
                     l = l+1
                     x = u_arr(i,j,k)*nu(m)/obsfreq
                     y = v_arr(i,j,k)*nu(m)/obsfreq
                     z = w_arr(i,j,k)*nu(m)/obsfreq
                     if (write_vis) then
                        if (wis) then
                           write(iunit,*) l,x,y,z,dum_re,dum_im,noise
                        else
                           write(iunit,*) l,x,y,dum_re,dum_im,noise
                        end if
                     end if
                     call pgsci(col(i,j,k))
                     call pgpt1(x,y,1)
                  end do
               end do
            end if
         end do
      end do

! Close the file
      if (write_vis) then
         close(iunit)
      end if

      call pgsci(plot_col)

   end subroutine plot_uv

! ***************************************************************************

   subroutine superpose_cntr(map,mapsize,pix_xd,pix_yd,&
                             cent_xd,cent_yd)

      real(kind=dp), dimension(mapsize,mapsize), intent(in) :: map
      real, allocatable, dimension(:) :: cont
      integer, intent(in) :: mapsize
      integer :: blc_x, blc_y, plotsize, ncont, status, i
      real(kind=dp), intent(in) :: pix_xd, pix_yd, cent_xd, cent_yd
      real :: pix_x, pix_y, cent_x, cent_y
      real :: maxlev, minlev, centrelev, startlev, inclev
      real(kind=dp) :: max, min
      character*3 :: switch
      character*1 :: char1

! Set up as last call to plot_map
      plotsize = PLTSZ
      blc_x = BLCX
      blc_y = BLCY

      max = maxval(map)
      min = minval(map)
      maxlev = max
      minlev = min
      switch = 'sup'
      pix_x = pix_xd
      pix_y = pix_yd
      cent_x = cent_xd
      cent_y = cent_yd

! Choose which type of contours and allocate contour array
      get_contours : do
         if (allocated(cont)) then
            deallocate(cont)
         end if
         call io_getc('Standard, linear, exponential, or defined &
              &contours (l,e,d)','s',char1,status) 

! Exponential contours
         if ((char1=='e').or.(char1=='E')) then
            if ((minlev <= 0.0).or.((maxlev/minlev) > 33000.0)) then
               minlev = maxlev/33000.0
            end if
            ncont = int(2*log(maxlev/minlev)/log(2.0))
            allocate (cont(ncont))
            cont(1) = minlev
            do i = 2, ncont
               cont(i) = cont(i-1)*1.4142135
            end do

! User defined contours
         else if ((char1=='d').or.(char1=='D')) then
            call io_geti('Number of contours','10',ncont,status)
            allocate(cont(ncont)) 
            do i= 1,ncont
               call io_getr('Contour','1.0',cont(i),status)
            end do

! Linear contours
         else if ((char1=='l').or.(char1=='L')) then
            call io_getr('Start level','1.0',startlev,status)
            call io_getr('Increment level','1.0',inclev,status)
            call io_geti('Number of contours','10',ncont,status)
            allocate (cont(ncont))
            cont(1) = startlev
            do i = 2,ncont
               cont(i) = cont(i-1)+inclev
               if (cont(i)==0.0) then
                  cont(i) = cont(i)+inclev
               end if
            end do

! Standard contours
         else
            call io_geti('Number of contours','10',ncont,status)
            allocate (cont(ncont))
            cont(1) = minlev
            inclev = (maxlev-minlev)/(ncont+1)
            do i = 2,ncont
               cont(i) = cont(i-1)+inclev
               if (cont(i)==0.0) then
                  cont(i) = cont(i)+inclev
               end if
            end do
         end if
         exit get_contours
      end do get_contours

      call plot_map(map, mapsize, switch, blc_x, blc_y, plotsize, &
                 ncont,cont,maxlev,minlev,pix_x,pix_y,cent_x,cent_y)

      if (allocated(cont)) then
         deallocate(cont)
      end if

   end subroutine superpose_cntr

! ***************************************************************************

   subroutine plot_uvhist(n_samp,ha_start,ha_inc,dec,datmin,datmax,nbin)

      implicit none

      integer :: n_samp, nbin
      real(kind=dp) :: ha_start, ha_inc, dec
      real :: datmin, datmax
      real :: x1, x2, y1, y2
      integer :: i, j, k, l, n_hist, which_bin
      real(kind=dp), dimension(n_antennas,n_antennas,n_samp) :: u_arr,v_arr
      real(kind=dp), dimension(n_antennas,n_antennas,n_samp) :: w_arr
      real(kind=dp), dimension(n_antennas,n_antennas) :: u,v,w
      integer, dimension(n_antennas,n_antennas) :: shadow
      integer, dimension(n_antennas,n_antennas,n_samp) :: col 
      integer, dimension(nbin) :: hist_count
      real(kind=dp), dimension(nbin+1) :: x_bin
      real, dimension(n_antennas*(n_antennas-1)*n_samp/2) :: uv_dist
      real(kind=dp) :: ha

      x1 = 0.1
      x2 = 0.9
      y1 = 0.25
      y2 = 0.8

! Calculate u, v and check for shadowing
      do k = 1, n_samp
         ha = (ha_start+(k-1)*ha_inc)
         call calc_uvw(ha,dec,u,v,w,shadow)
         u_arr(:,:,k) = u
         v_arr(:,:,k) = v
         w_arr(:,:,k) = w
         col(:,:,k) = shadow
      end do

      n_hist = 0
      do i = 1, n_antennas-1
         do j = i+1, n_antennas
            do k = 1, n_samp

! Only include in histogram if unflagged 
               if (col(i,j,k).eq.1) then
                  n_hist = n_hist+1
                  uv_dist(n_hist) = sqrt(u_arr(i,j,k)**2+v_arr(i,j,k)**2)
               end if
            end do
         end do
      end do

! Write out numbers in bins in requested
      if (verbose) then
         hist_count = 0
         x_bin(1) = datmin
         do i = 2,nbin+1
            x_bin(i) = datmin+(i-1)*(datmax-datmin)/nbin
         end do
         n_hist = 0
         do i = 1, n_antennas-1
            do j = i+1, n_antennas
               do k = 1, n_samp

! Only include in histogram if unflagged 
                  if (col(i,j,k).eq.1) then
                     n_hist = n_hist+1
                     which_bin = (nbin*(uv_dist(n_hist)-datmin))/&
                                 (datmax-datmin)+1
                     if ((which_bin.gt.0).and.(which_bin.le.nbin)) then
                        hist_count(which_bin) = hist_count(which_bin)+1
                     end if
                  end if
               end do
            end do
         end do
         write(*,*) 'Histogram bins and count'
         do i = 1, nbin
            write(*,*) i,x_bin(i),x_bin(i+1),hist_count(i)
         end do
      end if

! Do plotting
      call pgvsize(2., 6.5, 1., 5.5)
      call pgsvp(x1,x2,y1,y2)
      call pghist(n_hist,uv_dist,datmin,datmax,nbin,0)
      call pglab('uv distance (lambda)','Number',&
                 'Histogram of radial uv-coverage')

   end subroutine plot_uvhist

! ***************************************************************************
   subroutine map_histogram(map,mapsize,datmin,datmax,nbin)

      implicit none

      integer, intent(in) :: mapsize, nbin
      real(kind=dp), dimension(mapsize,mapsize), intent(in) :: map
      real :: datmin, datmax
      real, dimension(mapsize*mapsize) :: flux_val
      integer :: i, j, count
      real :: x1, x2, y1, y2

      x1 = 0.1
      x2 = 0.9
      y1 = 0.25
      y2 = 0.8

      count = 1
      do i = 1, mapsize
         do j = 1, mapsize
            flux_val(count) = map(i,j)
            count = count+1
         end do
      end do
      count = mapsize*mapsize

! Do plotting
      call pgvsize(2., 6.5, 1., 5.5)
      call pgsvp(x1,x2,y1,y2)
      call pghist(count,flux_val,datmin,datmax,nbin,0)
      call pglab('Flux','Number','Histogram of fluxes')

   end subroutine map_histogram

! ***************************************************************************

end module kgraph





