
subroutine two_d_tempmap(T)

   use sz_globals
   use kgraph

   implicit none
   
   integer :: i,j, x2,x3,y2,y3,irim,orim,inc
   real    :: T_dcont, minrad2, maxrad2, r2,r3,r1
   real    :: d1,d2,d3,theta_d
   real(dp), intent(out)    :: T(maxsize,maxsize)

   integer, PARAMETER :: nc=10
   REAL ::    c(nc),r(nc),g(nc),b(nc)
   REAL ::   min, max, immean, imsigma, trans(6)
   REAL ::    xmin,xmax,ymin,ymax
   
   real ::    array(maxsize,maxsize),prob(maxsize,maxsize),work(2*maxsize)
   complex :: carray(maxsize,maxsize),ckernel(maxsize,maxsize),cprob(maxsize,maxsize)
   integer :: icentre,jcentre,k,l,n,ninf,idum,job,iform,itype,inorm
   real    :: cell,dist,sigma,rdum,tiny,tmp,kernel,maxT_1,maxT_2,cft,cfr
   real    :: cf(maxsize,maxsize),dent,minT_1,vol

   tiny = 1e-32
   T = Te/keVtoK

! add a cooling flow:

   call io_getr('Extent of cooling flow:','30.',cfr,status)
   call io_getr('Minimum temperature of CF:','5.0',cft,status)

! make a Gaussian shaped CF:
   do i =1, maxsize
      do j = 1, maxsize
         dist = sqrt(real(i-maxsize/2)**2+real(j-maxsize/2)**2)*cellsize
         dent = exp(-(dist**2)/(2.*cfr*cfr))
         cf(i,j) = (T(i,j)-cft)*dent
         T(i,j) = T(i,j) - cf(i,j)
      end do
   end do

! add an arc shaped discontinuity:
   call io_getr('Temperature of discontinuity:','14.0',T_dcont,status)
   call io_geti('Inner edge of discontinuity:','165',irim,status)
   call io_geti('Outer edge of discontinuity:','220',orim,status)
   
   y2 = maxsize/2 - (80/cellsize)
   y3 = maxsize/2 + (80/cellsize)
   x2 = maxsize/2
   x3 = maxsize/2

   r2 = 0.
   r3 = 0.

   inc = 160/cellsize
   r2 = real((irim + inc)**2)
   r3 = real((orim - inc)**2)

   do i = 1, maxsize
      do j =1, maxsize
         minrad2 = 0.
         maxrad2 = 0.
         minrad2 = real((i-x2-1)**2 + (j-y2-1)**2)
         minrad2 = minrad2*cellsize**2
         maxrad2 = real((i-x3-1)**2 + (j-y3-1)**2)
         maxrad2 = maxrad2*cellsize**2
         
         if ((minrad2 .ge. r2) .and. (maxrad2 .le. r3)) then
            T(i,j) = T_dcont
         end if
      end do
   end do

   trans(1)=0.
   trans(2)=1.
   trans(3)=0.
   trans(4)=0.
   trans(5)=0.
   trans(6)=1.
   xmin=1.
   xmax=real(maxsize)
   ymin=1.
   ymax=real(maxsize)
   max = maxval(T)
   min = 0.

!   call pgenv(xmin,xmax,ymin,ymax,1,-2)
!   call ctab (c,r,g,b,nc)
!   call pgctab(c,r,g,b,nc,1.,.5)
!   call pgimag(real(T),maxsize,maxsize,1,maxsize,1,maxsize,min,max,trans) 
!!   call pggray(T,maxsize,maxsize,1,maxsize,1,maxsize,min,max,trans) 
!   call pgbox('BCTSN',0.0,0,'BCTSN',0.0,0)
!   call pglabel('x (pixels)','y (pixels)',' ')
!   call pgwedg('RI',1.,4.,min,max,' ')

   maxT_1 = maxval(T)
   minT_1 = minval(T)
   sigma = 1.5
! gaussian for convolution:
   do i=1,maxsize
      do j=1,maxsize
         dist=sqrt(real(i-maxsize/2)**2+real(j-maxsize/2)**2)
         kernel=exp(-(dist**2)/(2*sigma*sigma))
         ckernel(i,j)=kernel
      end do
   end do
   vol = sum(ckernel)
   ckernel = ckernel/vol   
   min = minval(real(ckernel))
   max = maxval(real(ckernel))
!   call pgimag(real(ckernel),maxsize,maxsize,1,maxsize,1,maxsize,min,max,trans)
!   call pgwedg('RI',1.,4.,min,max,' ')

   job=1                      ! forward transform (+ve exponential)
   iform=0                    ! data are real
   call makefft(maxsize,maxsize,ckernel,job,iform,work)

      do i=1,maxsize
        do j=1,maxsize
          rdum=T(i,j)
          carray(i,j)=cmplx(rdum,0.)
        end do
      end do

      job=1                      ! forward transform (+ve exponential)
      iform=0                    ! data are real
      call makefft(maxsize,maxsize,carray,job,iform,work)

      do i=1,maxsize
        do j=1,maxsize
           cprob(i,j)=carray(i,j)*ckernel(i,j)
        end do
      end do

      job=-1                     ! backward transform (-ve exponential)
      iform=1                    ! data are complex
      call makefft(maxsize,maxsize,cprob,job,iform,work)

      do i=1,maxsize
        do j=1,maxsize
          T(i,j)=real(cprob(i,j))
        end do
      end do

      maxT_2 = maxval(T)

   trans(1)=0.
   trans(2)=1.
   trans(3)=0.
   trans(4)=0.
   trans(5)=0.
   trans(6)=1.
   xmin=1.
   xmax=real(maxsize)
   ymin=1.
   ymax=real(maxsize)
   max = maxval(T)
   min = minval(T)

!   call pgenv(xmin,xmax,ymin,ymax,1,-2)
!   call ctab (c,r,g,b,nc)
!   call pgctab(c,r,g,b,nc,1.,.5)
!   call pgimag(real(T),maxsize,maxsize,1,maxsize,1,maxsize,min,max,trans) 
!!   call pggray(T,maxsize,maxsize,1,maxsize,1,maxsize,min,max,trans) 
!   call pgbox('BCTSN',0.0,0,'BCTSN',0.0,0)
!   call pglabel('x (pixels)','y (pixels)',' ')
!   call pgwedg('RI',1.,4.,min,max,' ')
!!   call pgcirc(x2,y2,r2)
!!   call pgcirc(x3,y3,r3)

!   call pgend
   
   T = T*keVtoK

end subroutine two_d_tempmap

! -------------------------------------------------------------------

subroutine ctab(x,r,g,b,nc)
 
   integer, intent(in) :: nc
   real, intent(out)  :: x(nc),r(nc),g(nc),b(nc)
   
   x(1)=0.00
   r(1)=0.00
   g(1)=0.00
   b(1)=0.00
   x(2)=0.25
   r(2)=0.5
   g(2)=0.0
   b(2)=0.0
   x(3)=0.5
   r(3)=1.0
   g(3)=0.5
   b(3)=0.0
   x(4)=0.75
   r(4)=1.0
   g(4)=1.0
   b(4)=0.5
   x(5)=1.0
   r(5)=1.0
   g(5)=1.0
   b(5)=1.0
   
 
end subroutine ctab

! -------------------------------------------------------------------
! -----*-----------------------------------------------------------------

      subroutine makefft(nx,ny,map,job,iform,work)
      implicit none
      integer  nx,ny,job,iform
      real     work(2*nx)
      complex  map(nx,ny) 

      integer  ndim,nn(2)

! ... initialise variables
      ndim=2
      nn(1)=nx
      nn(2)=ny

! ... calculate FFT
      call cheqboard(map,nx,ny)
      call fourt(map,nn,ndim,job,iform,work) 
      call cheqboard(map,nx,ny)
      call normfft(map,nx,ny,job)

      return
   end subroutine makefft

! -----*-----------------------------------------------------------------

      subroutine cheqboard(array,nx,ny)
      implicit none
      integer    nx,ny
      complex    array(nx,ny)
      
      integer    i,j,isign

      do i=1,nx
        do j=1,ny
          isign=(-1)**(i+j)
          array(i,j)=array(i,j)*float(isign)
        end do
      end do      

      return
   end subroutine cheqboard

! -----*-----------------------------------------------------------------

      subroutine normfft(array,nx,ny,job)
      implicit none
      integer  nx,ny,job
      complex  array(nx,ny)

      integer  i,j
      real     pi,factor

      pi=3.1415926535

      if (job.eq.1) then
        factor=1.                     ! forward transform (+ve exponential)
      else if (job.eq.-1) then
        factor=(1./real(nx))**2       ! backward transform (-ve exponential)
      end if

      do i=1,nx
        do j=1,ny
          array(i,j)=array(i,j)*factor
        end do
      end do

      return
   end subroutine normfft

! -----*-----------------------------------------------------------------

!      real function kernel(d,sigma,itype,inorm)
!      implicit none
!      integer  itype,inorm
!      real d,sigma
!      real pi
!
!      pi=3.1415926535
!
!! ... top-hat      
!      if (itype.eq.1) then
!        kernel=0.
!        if (d.le.sigma) kernel=1.
!        if (inorm.eq.1) kernel=kernel/(pi*sigma*sigma)
!
!! ... Gaussian
!      else if (itype.eq.2) then
!        kernel=exp(-d*d/(2*sigma*sigma))
!        if (inorm.eq.1) kernel=kernel/(2*pi*sigma*sigma) 
!
!      end if
!
!      return
!   end function kernel

! -----*-----------------------------------------------------------------
