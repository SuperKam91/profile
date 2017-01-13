! Downhill simplex method

module minimum

   use kind_def
   use structures
   use maths
   use stuff

   public :: amoeba
   public :: amotry

contains

!**************************************************************************

   subroutine amoeba(data_stuff,model_stuff,p,y,iter)

      integer, intent(out) :: iter
      type(data), intent(in) :: data_stuff
      type(model), intent(inout) :: model_stuff
      real(kind=dp), intent(inout), dimension(num_param+1,num_param) :: p
      real(kind=dp), intent(inout), dimension(num_param+1) :: y

      real(kind=dp), parameter :: ftol=1e-8
      integer, parameter :: nmax=20, itmax=5000
      integer :: i, ihi, ilo, inhi, j, m, n
      real(kind=dp) :: rtol, sum, swap, ysave, ytry
      real(kind=dp), dimension(nmax) :: psum
   
  
      iter=0
1     do n=1, num_param		!Enter here when starting or have just
         sum=0.0			!overall contracted.
         do m=1, num_param+1
            sum=sum+p(m,n)
         end do
         psum(n)=sum
      end do
2     ilo=1
      if(y(1).gt.y(2))then		!Determine which point is the highest
         ihi=1			!(worst), next-highest, and lowest
         inhi=2
      else
         ihi=2
         inhi=1
      endif

! DEBUG 
!        write (*,*) 'Got to point 2'

      do i=1,num_param+1
         if(y(i).le.y(ilo)) ilo=i
         if(y(i).gt.y(ihi)) then
            inhi=ihi
            ihi=i
         else if(y(i).gt.y(inhi))then
            if(i.ne.ihi) inhi=i
         endif
      end do
!        write(*,*) 'Got to point 3'
! Compute the fractional range from highest to lowest and return if 
! satisfactory
      rtol=2.0*abs(y(ihi)-y(ilo))/(abs(y(ihi))+abs(y(ilo)))

      if (loud_fitting) then
         write(*,*)'Iteration number ', iter
         write(*,*)'worst likelyhood = ',y(ihi),' best likelyhood=',y(ilo)
         write(*,*)'best parameter list:'
         write(*,*) p(ilo,:)
      end if

      if(rtol.lt.ftol)then		!If returning, put best point and 
         swap=y(1)		!value in slot 1.
         y(1)=y(ilo)
         y(ilo)=swap
         do n=1,num_param
            swap=p(1,n)
            p(1,n)=p(ilo,n)
            p(ilo,n)=swap
         end do
         return
      end if
      if(iter.ge.itmax) then 
!		pause'itmax exceeded in amoeba'
         return
      end if
      iter=iter+2
! Begin a new iteration. First extrapolate by a factor -1 through the face of 
! the simplex across from the high point, i.e., reflect the simplex from the 
! high point.
      ytry=amotry(data_stuff,model_stuff,p,y,psum,ihi,-1.0d0)
      if(ytry.le.y(ilo))then
         ytry=amotry(data_stuff,model_stuff,p,y,psum,ihi,2.0d0)
      else if(ytry.ge.y(inhi))then
         ysave=y(ihi)
         ytry=amotry(data_stuff,model_stuff,p,y,psum,ihi,0.5d0)
         if(ytry.ge.ysave)then
!DEBUG
!		write(*,*)'doing extra bit'
            do i=1,num_param+1
               if(i.ne.ilo)then
                  do j=1,num_param
                     psum(j)=0.5d0*(p(i,j)+p(ilo,j))
                     p(i,j)=psum(j)
                  end do
                  model_stuff%m_param%coeff=&
                       psum(1:num_param)
                  y(i)=-create_and_calculate&
                       (data_stuff,model_stuff)
               end if
            end do
            iter=iter+num_param
            goto 1
         end if
      else
         iter=iter-1
      end if
      goto 2

   end subroutine amoeba

!**************************************************************************

   function amotry(data_stuff,model_stuff,p,y,psum,ihi,fac)result(amotryres)

      type(data), intent(in) :: data_stuff
      type(model), intent(inout) :: model_stuff
      integer, intent(in) :: ihi
      real(kind=dp), intent(inout), dimension(num_param) :: psum
      real(kind=dp), intent(inout), dimension(num_param+1) :: y
      real(kind=dp), intent(inout), dimension(num_param+1, num_param) :: p
      real(kind=dp), intent(in) :: fac
      real(kind=dp) :: amotryres

      integer, parameter :: nmax=20
      integer :: j
      real(kind=dp) ::fac1,fac2,ytry
      real(kind=dp), dimension(nmax) :: ptry

      fac1=(1.0-fac)/num_param
      fac2=fac1-fac
      do j=1,num_param
         ptry(j)=psum(j)*fac1-p(ihi,j)*fac2
      end do
      model_stuff%m_param%coeff=ptry(1:num_param)
      ytry=-create_and_calculate(data_stuff,model_stuff)
      if(ytry.lt.y(ihi))then
         y(ihi)=ytry
         do j=1,num_param
            psum(j)=psum(j)-p(ihi,j)+ptry(j)
            p(ihi,j)=ptry(j)
         end do
      end if
      amotryres=ytry

   end function amotry

!**************************************************************************

end module minimum
