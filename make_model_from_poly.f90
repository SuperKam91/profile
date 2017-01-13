! Routine which produces 2D T model with a polynomial temperature 
! distribution
! Moved from make_skies 28/08/09 [KG] 
subroutine make_model_from_poly(T)

   use kind_def
   use sz_globals

   implicit none

   integer :: i,j, i1, j1, which_i
   real(kind=dp), dimension(:), allocatable :: radius
   real(kind=dp), dimension(maxsize,maxsize) :: T
   real(kind=dp) :: top, bottom, factor, theta2, theta_cx2
   real(kind=dp) :: nt_x
   integer :: largenumber
 
   num_nt = 200
   ntstep = 0.05
   allocate(Tdat(num_nt))

   do i=1,num_nt
      Tdat(i) = (poly_a+poly_b*i+poly_c*i*i+poly_d*i*i*i)
      if (Tdat(i) < 0) then 
         Tdat(i) = 0.
      end if
   end do

! need to normalize this. 
   allocate(ndat(num_nt))
   allocate(radius(num_nt))
   do i=1,num_nt
      ndat(i) = (1+ (i*ntstep)**2)**(-beta*1.5)
      radius(i) = i*ntstep
   end do
! okay. That has produced a density profile of the required shape
   top = 0.
   bottom = 0.
   largenumber = 100
   do i=1, largenumber
      if (radius(i) < 3.) then 
         top = top + ndat(i)**2 * tdat(i)**1.3 * radius(i)**2 * (radius(i+1) - radius(i))
         bottom = bottom + ndat(i)**2 * tdat(i)**0.3 * radius(i)**2 * (radius(i+1) - radius(i))
      end if
   end do
   factor = top/bottom

   deallocate(radius)
   deallocate(ndat)

   do i=1,maxsize
      do j=1,maxsize
         i1=i-1
         j1=j-1
         theta2=sqrt(float(i1**2+j1**2))
         nt_x=theta2/sqrt(theta_cx2)
! nt_x is the dimensionless core radius, just as we've been using
! to do the fitting. Which is nice.
         which_i = nint(nt_x/ntstep)
         if (which_i< num_nt) then
            T(j,i)=Tdat(which_i+1)/factor
         else
            T(j,i) = 0
         end if
      end do
   end do
   deallocate(Tdat)

end subroutine make_model_from_poly
