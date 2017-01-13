
module oneDFT

   public :: four1
   public :: realfft
   public :: twofft
   public :: convlv

   interface four1
      module procedure four1_real
      module procedure four1_comp
   end interface


contains
   
! ---------------------------------------------------------------

   subroutine four1_comp(data,nn,isign)
      

      integer, parameter     :: dp = kind(1.0D0)
      integer, intent(in)    :: isign, nn
      complex, intent(inout) :: data(nn)
      integer                :: i,istep,j,m,mmax,n
      real                   :: tempi, tempr
      real(dp)               :: theta,wi,wpi,wpr,wr,wtemp

!      write(*,*) 'got to four1'
      
      n = 2*nn
      j = 1
!      write(*,*) 'start of first loop'
      do i = 1,n,2
         if (j .gt. i) then
            tempr = data(j)
            tempi = date(j+1)
            data(j) = data(i)
            data(j+1) = data(i+1)
            data(i) = tempr
            data(i+1) = tempi
         end if
!         write(*,*) 'end of first if statement'
         m = nn
1        if ((m .ge. 2) .and. (j .gt. m)) then
            j = j-m
            m = m/2
            goto 1
         end if
!         write(*,*) 'end of second if statement',i
         j = j+m
      end do
!      write(*,*) 'end of first loop'
      mmax = 2
2     if (n .gt. mmax) then
         istep = 2*mmax
         theta = 6.28318530717959d0/(isign*mmax)
         wpr = -2.d0*sin(0.5d0*theta)**2
         wpi = sin(theta)
         wr = 1.d0
         wi = 0.d0
         do m = 1,mmax,2
            do i = m,n,istep
               j = i+mmax
               tempr = real(wr)*data(j)-real(wi)*data(j+1)
               tempi = real(wr)*data(j+1)+real(wi)*data(j)
               data(j) = data(i) - tempr
               data(j+1) = data(i+1) - tempi
               data(i) = data(i) + tempr
               data(i+1) = data(i+1) + tempi
            end do
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
         end do
         mmax = istep
         goto 2
      end if
      return
      
   end subroutine four1_comp

! ------------------------------------------------------------------

   subroutine four1_real(data,nn,isign)
      

      integer, parameter     :: dp = kind(1.0D0)
      integer, intent(in)    :: isign, nn
      real, intent(inout)    :: data(2*nn)
      integer                :: i,istep,j,m,mmax,n
      real                   :: tempi, tempr
      real(dp)               :: theta,wi,wpi,wpr,wr,wtemp

!      write(*,*) 'got to four1'
      
      n = 2*nn
      j = 1
!      write(*,*) 'start of first loop'
      do i = 1,n,2
         if (j .gt. i) then
            tempr = data(j)
            tempi = date(j+1)
            data(j) = data(i)
            data(j+1) = data(i+1)
            data(i) = tempr
            data(i+1) = tempi
         end if
!         write(*,*) 'end of first if statement'
         m = nn
1        if ((m .ge. 2) .and. (j .gt. m)) then
            j = j-m
            m = m/2
            goto 1
         end if
!         write(*,*) 'end of second if statement',i
         j = j+m
      end do
!      write(*,*) 'end of first loop'
      mmax = 2
2     if (n .gt. mmax) then
         istep = 2*mmax
         theta = 6.28318530717959d0/(isign*mmax)
         wpr = -2.d0*sin(0.5d0*theta)**2
         wpi = sin(theta)
         wr = 1.d0
         wi = 0.d0
         do m = 1,mmax,2
            do i = m,n,istep
               j = i+mmax
               tempr = real(wr)*data(j)-real(wi)*data(j+1)
               tempi = real(wr)*data(j+1)+real(wi)*data(j)
               data(j) = data(i) - tempr
               data(j+1) = data(i+1) - tempi
               data(i) = data(i) + tempr
               data(i+1) = data(i+1) + tempi
            end do
            wtemp = wr
            wr = wr*wpr - wi*wpi + wr
            wi = wi*wpr + wtemp*wpi + wi
         end do
         mmax = istep
         goto 2
      end if
      return
      
   end subroutine four1_real

! ------------------------------------------------------------------

   subroutine twofft(data1,data2,fft1,fft2,n)

      integer, intent(in)     :: n
      real, intent(in)        :: data1(n), data2(n)
      complex, intent(out)    :: fft1(n), fft2(n)
      integer                 :: j,n2
      real                    :: fft1_real(n)
      complex                 :: h1,h2,c1,c2

!      write(*,*) 'got to twofft'

      c1 = cmplx(0.5,0.0)
      c2 = cmplx(0.0,-0.5)
      do j = 1,n
         fft1(j) = cmplx(data1(j),data2(j))
      end do
!      fft1_real = real(fft1)
      call four1(fft1,n,1)
!      write(*,*) 'got back to twofft'
!      fft1 = cmplx(fft1_real)
!      write(*,*) 'debug'
      fft2(1) = cmplx(aimag(fft1(1)),0.0)
      fft1(1) = cmplx(real(fft1(1)),0.0)
      n2 = n+2
      do j = 2, n/2+1
         h1 = c1*(fft1(j)+conjg(fft1(n2-j)))
         h2 = c2*(fft1(j)-conjg(fft1(n2-j)))
         fft1(j) = h1
         fft1(n2-j) = conjg(h1)
         fft2(j) = h2
         fft2(n2-j) = conjg(h2)
      end do
      return
      
   end subroutine twofft

! -------------------------------------------------------------

   subroutine realfft(data,n,isign)

      use kind_def

      integer, intent(in)       :: isign,n
      complex, intent(inout)    :: data(n)
      integer                   :: i,i1,i2,i3,i4,n2p3
      real                      :: c1,c2,h1i,h2i,h2r,wis,wrs
      real(dp)                  :: theta,wi,wpi,wpr,wr,wtemp

!      write(*,*) 'got to realfft'

      theta = 3.141592653589793d0/dble(n/2)
      c1 = 0.5
      if (isign .eq. 1) then
         c2 =-0.5
         call four1(data,n/2,1)
      else
         c2 =0.5
         theta = -theta
      end if
      wpr = -2.0d0*sin(0.5d0*theta)**2
      wpi = sin(theta)
      wr = 1.0d0+wpr
      wi = wpi
      n2p3 = n+3
      do i = 2,n/4
         i1 = 2*i-1
         i2 = i1+1
         i3 = n2p3-i2
         i4 = i3+1
         wrs = real(wr)
         wis = real(wi)
         h1r = c1*(data(i1)+data(i3))
         h1i = c1*(data(i2)-data(i4))
         h2r = c2*(data(i2)+data(i4))
         h2i = c2*(data(i1)-data(i3))
         data(i1) = h1r+wrs*h2r-wis*h2i
         data(i2) = h1i+wrs*h2i+wis*h2r
         data(i3) = h1r-wrs*h2r+wis*h2i
         data(i4) = -h1i+wrs*h2i+wis*h2r
         wtemp=wr
         wr = wr*wpr-wi*wpi+wr
         wi = wi*wpr+wtemp*wpi+wi
      end do

      if (isign == 1) then
         h1r = data(1)
         data(1) = h1r+data(2)
         data(2) = h1r-data(2)
      else
         h1r = data(1)
         data(1) = c1*(h1r+data(2))
         data(2) = c1*(h1r-data(2))
         call four1(data,n/2,-1)
      end if
      return
   end subroutine realfft

! ---------------------------------------------------------------

   subroutine convlv(data,n,respns,m,isign,ans)

      integer, intent(in)        :: isign,m,n
      real, intent(in)           :: data(n)
      real, intent(inout)        :: respns(n)
      integer, parameter         :: nmax = 4096
      complex, intent(out)       :: ans(n)
      integer                    :: i,no2
      complex                    :: fft(nmax)
      real                       :: ans_real(n)

!      write(*,*) 'got to convlv'
      
      do i = 1,(m-1)/2
         respns(n+1-i) = respns(m+1-i)
      end do

      do i = (m+3)/2,n-(m-1)/2
         respns(i) = 0.0
      end do

      call twofft(data,respns,fft,ans,n)
      
      no2 = n/2
      do i = 1,no2+1
         if (isign .eq. 1) then
            ans(i) = fft(i)*ans(i)/no2
         else if (isign .eq. -1) then
            if (abs(ans(i)) .eq. 0.0) then
               write(*,*) 'deconvolving at response zero'
               return
            end if
            ans(i) = fft(i)/ans(i)/no2
         else
            write(*,*) 'no meaning for isign, that means isign is WRONG'
            return
         end if
      end do
      ans(1) = cmplx(real(ans(1)),real(ans(no2+1)))
      call realfft(ans,n,-1)

      return
   end subroutine convlv

! ---------------------------------------------------------------

















end module oneDFT


















