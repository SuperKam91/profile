module maths

   use kind_def

   implicit none

   interface do_fft
      module procedure do_2d_fft
   end interface

   interface j0
      module procedure calc_j0
   end interface

   interface interp
      module procedure interp_1d_sg
      module procedure interp_2d_sg
      module procedure interp_1d_db
      module procedure interp_2d_db
   end interface

   interface poidev
      module procedure calc_poidev
   end interface

   interface gasdev
      module procedure calc_gasdev
   end interface

   interface factorial
      module procedure factorial
   end interface

   interface log_factorial
      module procedure log_factorial
   end interface

   interface polint
      module procedure polint
   end interface

   interface polin2
      module procedure polin2
   end interface

   interface ran2
      module procedure ran2
   end interface

   interface locate
      module procedure locate_dp
   end interface

   public :: coth
   public :: gammln
   public :: best_fit_line
   public :: Bfnc
   public :: qtrap
   public :: trapzd
   public :: funlog10
   public :: hunt
   public :: rombint
   public :: spline

contains

!**************************************************************************
   function calc_j0(x) result(res)

! Calculates the j0 Bessel function

      implicit none

      real(kind=dp), intent(in) :: x
      real(kind=dp) :: res
      real(kind=dp), parameter :: pi = 3.14159265359
      real(kind=dp), parameter :: cp0= 1.0
      real(kind=dp), parameter :: cp2=-1.0*9.0          /( 2.0*(8.0**2))
      real(kind=dp), parameter :: cp4= 1.0*9.0*25.0*49.0/(24.0*(8.0**4))
      real(kind=dp), parameter :: cq1=-1.0              /(      8.0    )
      real(kind=dp), parameter :: cq3= 1.0*9.0*25.0     /( 6.0*(8.0**3))
      real(kind=dp) :: tm, vl, z, chi, p, q

      if (abs(x) .lt. 10.0) then
         res=1.0
         tm=1.0
         vl=0.0
         z =0.25*(x**2)
         loop1 : do
            if (abs(tm) .lt. 1.0e-06) then
               exit loop1
            end if
            vl=vl+1.0
            tm=-tm*(z/(vl**2))
            res=res+tm
         enddo loop1
      else
         z  =abs(x)
         chi=z-0.25*pi
         p  =cp0+(cp2/(z**2))+(cp4/(z**4))
         q  =    (cq1/ z    )+(cq3/(z**3))
         res =sqrt(2.0/(pi*z))*(p*cos(chi)-q*sin(chi))
      endif

   end function calc_j0

!*****************************************************************************

   subroutine do_2d_fft(re,im,nx,ny,dirn)
! Do a fft on a 2d array
!
!     Input:
!       re(nx,ny)  Real part of two dimensional array to be ffted
!       im(nx,ny)  Imaginary part of two dimensional array to be ffted 
!       nx         First dimension of array to be ffted
!       ny         Second dimension of array to be ffted
!       dirn       Do forward (dirn=1) or back (dirn=-1) transform
!
!     Output:
!       re(nx,ny)  ffted array(real)
!       im(nx,ny)  ffted array(imag)
      implicit none
      integer, intent(in) :: nx, ny
      integer, intent(in) :: dirn
      integer :: i, j, k, data_size 
      integer, dimension(2) :: nn  
      real(kind=dp), dimension(2*nx*ny) :: data1
      real(kind=dp), dimension(nx,ny), intent(inout) :: re, im

! rearrange quadrants (zapping or smileying)
      call xzapit(re,nx,ny)
      call xzapit(im,nx,ny)

      nn(1) = nx
      nn(2) = ny
      data_size = 2*product(nn)
      k = 1
      do j = 1,ny
         do i = 1,nx
            data1(k) = re(i,j)
            k = k + 1
            data1(k) = im(i,j)
            k = k + 1
         end do
      end do
      call fourn(data1,data_size,nn,2,dirn)

      k = 1
      do j = 1,ny
         do i = 1,nx
            re(i,j) = data1(k)
            k = k + 1
            im(i,j) = data1(k)
            k = k + 1
         end do
      end do

! rearrange quadrants (zapping or smileying)

      call xzapit(re,nx,ny)
      call xzapit(im,nx,ny)

   end subroutine do_2d_fft

!*********************************************************************

   subroutine xzapit(x,nx,ny)
!     does zapping (smilying) required before and after fft on 2-d map
!
!     Input:
!       x(nx,ny) Two dimensional array of data to be zapped
!       nx       First dimension of array to be zapped
!       ny       Second dimension of array to be zapped
!
!     Output:
!       x(nx,ny) Zapped array
!
!     Function:
!       Zaps the Fourier output file into a proper display

      implicit none
      integer, intent(in) :: nx, ny
      real(kind=dp), dimension(nx,ny), intent(inout) :: x
      real(kind=dp) :: dummy
      integer i,j,k,m

      do i = 1,nx
         do j = 1,ny/2
            k = (i+nx/2)
            k = k-(k/nx)*nx
            if (k == 0) then
               k = nx
            end if
            m = (j+ny/2)
            m = m-(m/ny)*ny
            if (m == 0) then
               m = ny
            end if
            dummy = x(i,j)
            x(i,j) = x(k,m)
            x(k,m) = dummy
         end do
      end do
   end subroutine xzapit

!*************************************************************************

! 2-d FFT from Numerical Recipies. Faster than NAG.
   subroutine fourn(data,data_size,nn,ndim,isign)

      implicit none

      real(kind=dp) :: wr,wi,wpr,wpi,wtemp,theta,tempr,tempi
      integer, intent(in) :: data_size, ndim, isign
      integer, dimension(ndim), intent(in) :: nn
      real(kind=dp), dimension(data_size), intent(inout) :: data
      integer :: ntot, idim, nprev, n, nrem, ip1, ip2, ip3, i2rev, i2
      integer :: i1, i3, i3rev, ibit, ifp1, ifp2, k1, k2

      ntot=1
      do idim=1,ndim
         ntot=ntot*nn(idim)
      end do
      nprev=1
      do idim=1,ndim
         n=nn(idim)
         nrem=ntot/(n*nprev)
         ip1=2*nprev
         ip2=ip1*n
         ip3=ip2*nrem
         i2rev=1
         do i2=1,ip2,ip1
            if(i2.lt.i2rev)then
               do i1=i2,i2+ip1-2,2
                  do i3=i1,ip3,ip2
                     i3rev=i2rev+i3-i2
                     tempr=data(i3)
                     tempi=data(i3+1)
                     data(i3)=data(i3rev)
                     data(i3+1)=data(i3rev+1)
                     data(i3rev)=tempr
                     data(i3rev+1)=tempi
                  end do
               end do
            endif
            ibit=ip2/2
            cycle_1: do
               if ((ibit.ge.ip1).and.(i2rev.gt.ibit)) then
                  i2rev=i2rev-ibit
                  ibit=ibit/2
               else
                  exit cycle_1
               end if
            end do cycle_1
            i2rev=i2rev+ibit
         end do
         ifp1=ip1
         cycle_2 : do
            if(ifp1.lt.ip2)then
               ifp2=2*ifp1
               theta=isign*6.28318530717959d0/(ifp2/ip1)
               wpr=-2.d0*dsin(0.5d0*theta)**2
               wpi=dsin(theta)
               wr=1.d0
               wi=0.d0
               do  i3=1,ifp1,ip1
                  do i1=i3,i3+ip1-2,2
                     do i2=i1,ip3,ifp2
                        k1=i2
                        k2=k1+ifp1
                        if((data(k2).ne.0).or.(data(k2+1).ne.0)) then
                           tempr=sngl(wr)*data(k2)-sngl(wi)*data(k2+1)
                           tempi=sngl(wr)*data(k2+1)+sngl(wi)*data(k2)
                        else
                           tempr = 0.0
                           tempi = 0.0
                        end if
                        data(k2)=data(k1)-tempr
                        data(k2+1)=data(k1+1)-tempi
                        data(k1)=data(k1)+tempr
                        data(k1+1)=data(k1+1)+tempi
                     end do
                  end do
                  wtemp=wr
                  wr=wr*wpr-wi*wpi+wr
                  wi=wi*wpr+wtemp*wpi+wi
               end do
               ifp1=ifp2
            else
               exit cycle_2
            endif
         end do cycle_2
         nprev=n*nprev
      end do

   end subroutine fourn

!***************************************************************************

   subroutine interp_1d_sg(table,size,index,value)

      implicit none
      integer, intent(in) :: size
      real, dimension(size), intent(in) :: table
      real, intent(in) :: index
      real, intent(out) :: value
      real :: temp1, temp2
      integer :: i

      i = int(index)
      if ((i.le.0).or.(i.ge.size)) then
         value = 0.0
      else
         temp1 = table(i)
         temp2 = table(i+1)
         value = temp1+((temp2-temp1)*(index-float(i)))
      endif

   end subroutine interp_1d_sg

!***************************************************************************

   subroutine interp_1d_db(table,size,index,value)

      implicit none
      integer, intent(in) :: size
      real(kind=dp), dimension(size), intent(in) :: table
      real(kind=dp), intent(in) :: index
      real(kind=dp), intent(out) :: value
      real(kind=dp) :: temp1, temp2
      integer :: i

      i = int(index)
      if ((i.le.0).or.(i.ge.size)) then
         value = 0.0
      else
         temp1 = table(i)
         temp2 = table(i+1)
         value = temp1+((temp2-temp1)*(index-float(i)))
      end if

   end subroutine interp_1d_db

!***************************************************************************

   subroutine interp_2d_sg(map,size_x,size_y,x,y,value)

      implicit none
      integer, intent(in) :: size_x, size_y
      real, dimension(size_x,size_y), intent(in) :: map
      real, intent(in) :: x, y
      real, intent(inout) :: value
      real :: t, u, t1, u1
      integer :: i, j

      i = int(x)
      j = int(y)
      t = x-i
      u = y-j
      t1 = 1-t
      u1 = 1-u
      if ((i.le.0).or.(i.ge.size_x).or.(j.le.0).or.(j.ge.size_y)) then
         value = 0.0
      else
         value = t1*u1*map(i,j)+t*u1*map(i+1,j)+t*u*map(i+1,j+1)+&
              & t1*u*map(i,j+1)
      end if

   end subroutine interp_2d_sg

!****************************************************************************

   subroutine interp_2d_db(map,size_x,size_y,x,y,value)

      implicit none
      integer, intent(in) :: size_x, size_y
      real(kind=dp), dimension(size_x,size_y), intent(in) :: map
      real(kind=dp), intent(in) :: x, y
      real(kind=dp), intent(inout) :: value
      real(kind=dp) :: t, u, t1, u1
      integer :: i, j


      i = int(x)
      j = int(y)
      t = x-i
      u = y-j
      t1 = 1-t
      u1 = 1-u

      if ((i.le.0).or.(i.ge.size_x).or.(j.le.0).or.(j.ge.size_y)) then
         value = 0.0
      else
         value = t1*u1*map(i,j)+t*u1*map(i+1,j)+t*u*map(i+1,j+1)+t1*u*map(i,j+1)

      end if

   end subroutine interp_2d_db


! ---------------------------------------------------------------------------

!****************************************************************************

!  NR routine to give poisson number of mean xm
   function calc_poidev(xm,idum)

      use kind_def

      implicit none
      real(kind=dp) :: calc_poidev
      real(kind=dp) :: xm
      real(kind=dp) :: g, em, t, oldm, sq, alxm, y
      integer :: idum
      real(kind=dp), parameter :: pi=3.141592654
      data oldm /-1./
      if (xm.lt.12.)then
         if (xm.ne.oldm) then
            oldm=xm
            g=exp(-xm)
         endif
         em=-1
         t=1.
2        em=em+1.
         t=t*ran2(idum)
         if (t.gt.g) go to 2
      else
         if (xm.ne.oldm) then
            oldm=xm
            sq=sqrt(2.*xm)
            alxm=log(xm)
            g=xm*alxm-gammln(xm+1.)
         endif
1        y=tan(pi*ran2(idum))
         em=sq*y+xm
         if (em.lt.0.) go to 1
         em=int(em)
         t=0.9*(1.+y**2)*exp(em*alxm-gammln(em+1.)-g)
         if (ran2(idum).gt.t) go to 1
      endif
      calc_poidev=em
      return
   end function calc_poidev

! **************************************************************************
   function gammln(xx)

      real(kind=dp) ::  gammln,xx
      integer j
      real(kind=dp) :: ser,stp,tmp,x,y,cof(6)
      save cof,stp
      data cof,stp / 76.18009172947146d0,   -86.50532032941677d0,&
           &               24.01409824083091d0,   -1.231739572450155d0,&
           &               0.1208650973866179d-2, -0.5395239384953d-5,&
           &               2.5066282746310005d0 /

      x = xx
      y = x
      tmp = x + 5.5d0
      tmp = (x+0.5d0)*log(tmp)-tmp
      ser = 1.000000000190015d0
      do j=1,6
         y = y + 1.d0
         ser = ser + cof(j)/y
      end do
      gammln = tmp + log(stp*ser/x)

   end function gammln

! **************************************************************************

   function ran2 (idum)

      use kind_def

      implicit none

      integer       :: idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real(kind=dp) :: ran2,am,eps,rnmx
      parameter (im1=2147483563, im2=2147483399, am=1./im1, imm1=im1-1,&
           &           ia1=40014, ia2=40692, iq1=53668, iq2=52774, ir1=12211,&
           &           ir2=3791, ntab=32, ndiv=1+imm1/ntab, eps=1.2e-7,&
           &           rnmx=1.0-eps)
      integer idum2,j,k,iv(ntab),iy
      save    iv,iy,idum2
      data    idum2/123456789/,iv/ntab*0/,iy/0/

      if (idum.le.0) then
         idum = max(-idum,1)
         idum2 = idum
         do j=ntab+8,1,-1
            k = idum/iq1
            idum = ia1*(idum-k*iq1)-k*ir1
            if (idum.lt.0) idum = idum + im1
            if (j.le.ntab) iv(j) = idum
         end do
         iy = iv(1)
      endif
      k = idum/iq1
      idum = ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum = idum+im1
      k = idum2/iq2
      idum2 = ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2 = idum2 + im2
      j = 1 + iy/ndiv
      iy = iv(j) - idum2
      iv(j) = idum
      if (iy.lt.1) iy = iy + imm1
      ran2 = min(am*iy,rnmx) 
      return

   end function ran2

! **************************************************************************

   function calc_gasdev(idum)

      use kind_def

      implicit none
      real(kind=dp) calc_gasdev
      integer :: iset, idum
      real(kind=dp) :: fac, gset,rsq, v1, v2
      save iset, gset
      data iset/0/
      if (iset.eq.0) then
1        v1=2.*ran2(idum)-1.
         v2=2.*ran2(idum)-1.
         rsq=v1**2+v2**2
         if((rsq.ge.1.).or.(rsq==0.))go to 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         calc_gasdev=v2*fac
         iset=1
      else
         calc_gasdev=gset
         iset=0
      endif

   end function calc_gasdev

!***************************************************************************

   function factorial(n)result(facn)

      implicit none
      real(kind=dp), intent(in) :: n
      real(kind=dp) :: facn
      integer :: i, n_int

      facn=1.0
      n_int = n
      do i= 1,n_int
         facn=facn*i
      end do
   end function factorial



!*************************************************************************

   function log_factorial(n)result(logfacn)

      implicit none
      real(kind=dp), intent(in) :: n
      real(kind=dp) :: logfacn
      integer :: i, n_int

      if (n.le.1000.0) then

! Calculate exactly by summing the logs of the numbers 1 to n
         logfacn = 0.0
         n_int = n
         do i = 1,n_int
            logfacn = logfacn+log(1.0D0*i)
         end do
      else

! Use Stirling's approximation
         logfacn = n*log(n)-n
      end if

   end function log_factorial

! ************************************************************************

   subroutine polint(xa,ya,n,x,y,dy)
      integer n,NMAX
      real(kind=dp) dy,x,y,xa(n),ya(n)
      parameter (NMAX=25)
      integer i,m,ns
      real(kind=dp) den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.) then
             write(*,*) 'failure in polint'
             return
          end if
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return

   end subroutine polint

! ************************************************************************

   subroutine polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
   integer :: m,n
   integer :: NMAX,MMAX
   real(kind=dp) :: dy,x1,x2,y
   real(kind=dp), dimension(m) :: x1a
   real(kind=dp), dimension(n) :: x2a
   real(kind=dp), dimension(m,n) :: ya
   parameter (NMAX=25,MMAX=25)

! Uses polint

      integer j,k
      real(kind=dp) ymtmp(MMAX),yntmp(NMAX)
      do 12 j=1,m
        do 11 k=1,n
          yntmp(k)=ya(j,k)
11      continue
        call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
12    continue
      call polint(x1a,ymtmp,m,x1,y,dy)
      return

   end subroutine polin2

! ************************************************************************

   function coth(x)result(res)

      real(kind=dp), intent(in) :: x
      real(kind=dp) :: res

      res = cosh(x)/sinh(x)

   end function coth

! **************************************************************************

!Given an array xx and a given value x returns a value j such that x is 
!between xx(j) and xx(j+1). xx must be monotonic, either increasing or 
!decreasing. j = 0 or j = N is returned to indicate that x is out of range

! nr routine
function locate_dp(xx,x)

   use kind_def

   implicit none

   real(kind=dp), dimension(:), intent(in) :: xx
   real(kind=dp), intent(in) :: x
   integer :: locate_dp
   integer :: n, jl, jm, ju
   logical :: ascnd

   n = size(xx)
   ascnd = (xx(n) >= xx(1))
   jl = 0
   ju = n+1
   do 
      if (ju-jl <= 1) exit
      jm=(ju+jl)/2
      if (ascnd.eqv.(x>=xx(jm))) then
         jl=jm
      else
         ju = jm
      end if    
   end do
   if (x == xx(1)) then 
      locate_dp = 1
   else if (x == xx(n)) then
      locate_dp = n-1
   else
      locate_dp = jl
   end if

end function locate_dp

! **********************************************************************

! BETA FUNCTION
    function Bfnc(z,w)

       use kind_def

       real(dp) :: Bfnc, z,w

       Bfnc = exp(gammln(z)+gammln(w)-gammln(z+w))
       return

    end function Bfnc

! **********************************************************************

! Finds best fit of points to line y = mx+c
    subroutine best_fit_line(xpts,ypts,npts,fit_m,fit_c)

       use kind_def

       integer, intent(in) :: npts
       real(dp), dimension(npts), intent(in) :: xpts, ypts
       real(dp), intent(out) :: fit_m, fit_c
       real(dp), dimension(npts) :: xy, x2

       xy = xpts*ypts
       x2 = xpts*xpts
       fit_m = (npts*sum(xy)-sum(xpts)*sum(ypts))/(npts*sum(x2)-sum(xpts)**2)
       fit_c = (sum(ypts)-fit_m*sum(xpts))/npts

    end subroutine best_fit_line

! **************************************************************************

    subroutine qtrap(func,lim1,lim2,eps,S)

! Numerical recipes integration routine-the supplied real function
! func is integrated between limits A and B, and the result returned as
! the integral S. Progressively finer discretisations are used until 
! either
!
!     | S_{j}-S_{j-1} | 
!     ------------------- < eps
!          | S_{j-1} |
!
! or the range has to be split into more than 2^jmax sections. For
! smooth functions the latter shouldnt happen...
       
      implicit none
	
      real*8 S,lim1,lim2
      real*8 a1,b1,eps,olds,sign,t1,t2
      integer j,jmax
      parameter (jmax=25)
      
      INTERFACE
         function func(zz)
           real*8 zz,func
         end function func
      end INTERFACE

      a1 = lim1
      b1 = lim2
	
      ! First catch some stupidities:

      if( a1 == b1 ) then
         S = 0.d0
         goto 30
      elseif( a1 > b1 ) then        
         olds = b1
         b1 = a1
         a1 = olds
         sign = -1.d0
      else 
         sign = 1.d0
      endif
	 
      olds = -1.d30
      s = 0.d0
      
      do j = 1,jmax
         call trapzd(func,a1,b1,S,j)
         if ( j == 1 ) then
            t1 = func(a1)
            t2 = func(b1)
         endif
		 
         if( abs( s - olds ) < eps * abs(olds) .or. abs( s - olds ) == 0.d0 ) goto 20
         if( j == jmax ) goto 10
         
         olds = s
      enddo
      
10    write(*,*) 'QTRAP error: too many steps...'
      write(*,*) '   S = ',S
      write(*,*) '   oldS = ',oldS
      write(*,*) '   % difference = ',100.0*(S-oldS)/oldS
      write(*,*) '   limits were ',a1,b1
      write(*,*) '   function values at limits were ',t1,t2
      t1 = func(0.5*(a1+b1))
      write(*,*) '   function value at midpoint was ',t1
      goto 30
      
20    S = S * sign
30    return
    end subroutine qtrap

! **************************************************************************
      
    subroutine trapzd(func,a1,b1,S,n)
	
      implicit none
	
      real*8 a1,b1,S
      integer n
      integer it,j
      real*8 tnm,del,x,sum
      
      INTERFACE
         function func(zz)
           real*8 zz,func
         end function func
      end INTERFACE
      
      if (n.eq.1) then
         S=0.5*(b1-a1)*(func(a1)+func(b1))
         it=1
      else
         it=2**(n-2)
         tnm=it*1.0
         del=(b1-a1)/tnm
         x=a1+0.5*del
         sum=0.
         do 11 j=1,it
            sum=sum+func(x)
            x=x+del
11          continue
            s=0.5*(s+(b1-a1)*sum/tnm)
      endif
      
      return
    end subroutine trapzd
	
! **************************************************************************

    FUNCTION funlog10(x)

      IMPLICIT NONE

      REAL *8       ::x,funlog10

      IF(x.LT.0.0)THEN
         WRITE(*,*)'Error, tried to take log of: ',x
         STOP
      ELSEIF(x.LT.1.0d-45) THEN
         funlog10=-45.0
      ELSE
         funlog10=LOG10(x)
      ENDIF
      RETURN
    END FUNCTION funlog10

! **************************************************************************

    SUBROUTINE hunt(xx,n,x,jlo)

      implicit none

      INTEGER jlo,n
      REAL*8 x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd

      ascnd=xx(n).gt.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)return
      jm=(jhi+jlo)/2
      if(x.gt.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3

   end subroutine hunt

! **************************************************************************

   function rombint(f,a,b,tol)

!  Rombint returns the integral from a to b of using Romberg integration.
!  The method converges provided that f(x) is continuous in (a,b).
!  f must be double precision and must be declared external in the calling
!  routine.  tol indicates the desired relative accuracy in the integral.
!
        use kind_def

        real(dp) :: rombint
        real(dp) :: a, b, h, gmax, error, g0, g1, fourj, tol
        integer :: nint, i, jmax, k, j

        integer,parameter :: MAXITER=40,MAXJ=5

        real(dp), dimension(MAXJ+1) :: g
        real(dp), external :: f
!
        h=0.5d0*(b-a)
        gmax=h*(f(a)+f(b))
        g(1)=gmax
        nint=1
        error=1.0d20
        i=0
10        i=i+1 
          if (i.gt.MAXITER.or.(i.gt.5.and.abs(error).lt.tol)) go to 40
!  Calculate next trapezoidal rule approximation to integral.
          g0=0.0d0
            do 20 k=1,nint
            g0=g0+f(a+(k+k-1)*h)
20        continue
          g0=0.5d0*g(1)+h*g0
          h=0.5d0*h
          nint=nint+nint
          jmax=min(i,MAXJ)
          fourj=1.0d0
            do 30 j=1,jmax
!  Use Richardson extrapolation.
            fourj=4.0d0*fourj
            g1=g0+(g0-g(j))/(fourj-1.0d0)
            g(j)=g0
            g0=g1
30        continue         
          if (abs(g0).gt.tol) then
            error=1.0d0-gmax/g0
          else
            error=gmax
          end if
          gmax=g0
          g(jmax+1)=g0
        go to 10
40      rombint=g0
        if (i.gt.MAXITER.and.abs(error).gt.tol)&
         write(*,*) 'Rombint failed to converge; integral, error=',&
         rombint,error
        return

   end function rombint

! **************************************************************************

   subroutine spline(x,y,n,yp1,ypn,y2)

      implicit none

      INTEGER n,NMAX
      REAL*8 yp1,ypn,x(n),y(n),y2(n)
      PARAMETER (NMAX=10000)
      INTEGER i,k
      REAL*8 p,qn,sig,un,u(NMAX)
      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do 11 i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+&
      1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*&
      u(i-1))/p
11    continue
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do 12 k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue
      return

   end subroutine spline

! **************************************************************************

   FUNCTION Gammafun(x)
    
    IMPLICIT NONE
     REAL *8   Gammafun , x
     INTEGER    i,j,k,nn
     REAL *8  w, y , fun
     REAL*8   p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13
     PARAMETER(p0 = 0.999999999999999990d+00, p1 = -0.422784335098466784d+00)
     PARAMETER(p2 = -0.233093736421782878d+00,p3 = 0.191091101387638410d+00)
     PARAMETER(p4 = -0.024552490005641278d+00,p5 = -0.017645244547851414d+00)
     PARAMETER(p6 = 0.008023273027855346d+00 ,p7 = -0.000804329819255744d+00)
     PARAMETER(p8 = -0.000360837876648255d+00,p9 = 0.000145596568617526d+00)
     PARAMETER(p10 = -0.000017545539395205d+00,p11 = -0.000002591225267689d+00)
     PARAMETER(p12 = 0.000001337767384067d+00,p13 = -0.000000199542863674d+00)     
     
      nn = nint(x - 2)
      w = x - (nn + 2)
      y = ((((((((((((p13 * w + p12) * w + p11) * w + p10) * &
         w + p9) * w + p8) * w + p7) * w + p6) * w + p5) * &
         w + p4) * w + p3) * w + p2) * w + p1) * w + p0      
      if (nn .gt. 0) then
          w = x - 1
          do k = 2, nn
              w = w * (x - k)
          end do
      else
          w = 1
          do k = 0, -nn - 1
              y = y * (x + k)
          end do
      end if
!      write(*,*)w,y
        fun= w / y
	Gammafun =fun
      end function Gammafun     

! **************************************************************************

	!root funding using the Brent method
      !Numerical Recipes

      FUNCTION zbrent(funct,x1,x2,tol)
      INTEGER ITMAX
      REAL*8 zbrent,tol,x1,x2,func,EPS
      PARAMETER (ITMAX=100,EPS=3.d-8)
      INTEGER iter
      REAL*8 a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
      
      INTERFACE
            function funct(zz)
            	real*8 zz,funct
            end function funct
      end INTERFACE
      
      a=x1
      b=x2
      fa=funct(a)
      fb=funct(b)
      if((fa>0. .and. fb>0.) .or. (fa<0. .and. fb<0.))pause &
      'root must be bracketed for zbrent'
      c=b
      fc=fb
      do 11 iter=1,ITMAX
        if((fb>0. .and. fc>0.) .or. (fb<0. .and. fc<0.))then
          c=a
          fc=fa
          d=b-a
          e=d
        endif
        if(abs(fc)<abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
        endif
        tol1=2.*EPS*abs(b)+0.5*tol
        xm=.5*(c-b)
        if(abs(xm)<=tol1 .or. fb==0.)then
          zbrent=b
          return
        endif
        if(abs(e)>=tol1 .and. abs(fa)>abs(fb)) then
          s=fb/fa
          if(a==c) then
            p=2.*xm*s
            q=1.-s
          else
            q=fa/fc
            r=fb/fc
            p=s*(2.*xm*q*(q-r)-(b-a)*(r-1.))
            q=(q-1.)*(r-1.)*(s-1.)
          endif
          if(p>0.) q=-q
          p=abs(p)
          if(2.*p<min(3.*xm*q-abs(tol1*q),abs(e*q))) then
            e=d
            d=p/q
          else
            d=xm
            e=d
          endif
        else
          d=xm
          e=d
        endif
        a=b
        fa=fb
        if(abs(d)>tol1) then
          b=b+d
        else
          b=b+sign(tol1,xm)
        endif
        fb=funct(b)
11    continue
      pause 'zbrent exceeding maximum iterations'
      zbrent=b
      return
      END function zbrent
	
!=======================================================================

end module maths








