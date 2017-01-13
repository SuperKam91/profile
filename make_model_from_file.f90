! Routine which produces 2D n and T models from file. 
! Moved from make_skies 28/08/09 [KG] 
subroutine make_model_from_file(T,n)

   use kind_def
   use sz_globals
   use kgraph

   implicit none

   integer :: i,j, i1, j1, which_i
   real(kind=dp), dimension(maxsize,maxsize) :: n, T
   real(kind=dp) :: theta2, theta_cx2
   real(kind=dp) :: nt_x

   do i=1,maxsize
      do j=1,maxsize
         i1=i-1
         j1=j-1
         theta2=sqrt(float(i1**2+j1**2))
         nt_x=theta2/sqrt(theta_cx2)
         which_i = nint(nt_x/ntstep)
         if (which_i < num_nt) then
            T(j,i)=Tdat(which_i+1)
! Change 28/8/09 - remove multiplication by n0
!            n(j,i)=n0*ndat(which_i(nt_x,ntstep)+1)
            n(j,i)=ndat(which_i+1)
         else
            T(j,i) = 0
            n(j,i) = 0
         end if
      end do
   end do

   if (do_plot) then
      call display_map(T,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
   end if

   deallocate(Tdat)
   deallocate(ndat)
   from_file = .false.

end subroutine make_model_from_file
