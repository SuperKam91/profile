module Te_map_ams

contains


   subroutine make_Te_map(T)

      use kind_def
      use sz_globals
      use kgraph

      implicit none

      integer :: i,j,k,npoint,i1,j1
      real(dp), intent(out) :: T(maxsize,maxsize)
      real(dp) :: theta2,rmin2,rmax2
      real(dp), allocatable :: x1(:),x2(:),y(:)
      character*80 infile
      logical :: graph_adjust


      call io_getc('Input file:','A2142_temp.dat',infile,status)
      open(unit=11, file=infile)
      Call io_geti('How many temperature points?','8',npoint,status)
      allocate(x1(npoint))
      allocate(x2(npoint))
      allocate(y(npoint))

      do i = 1, npoint
         read(11,*) x1(i),x2(i),y(i)
      end do
      close(unit = 11)

      do i = 1, maxsize
         do j = 1, maxsize
            i1 = i-1
            j1 = j-1
            theta2 = float(i1**2+j1**2)
            do k = 1, npoint
               rmin2 = x1(k)**2
               rmax2 = x2(k)**2
               if ((theta2 .ge. rmin2) .and. (theta2 .lt. rmax2)) then
                  T(j,i) = y(k)
               else if (theta2 .ge. x2(npoint)**2) then
                  T(j,i) = y(npoint)
               end if
            end do
         end do
      end do

      graph_adj = .true.
      call display_map(T,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
      
      deallocate(x1)
      deallocate(x2)
      deallocate(y)

   end subroutine make_Te_map


end module Te_map_ams

