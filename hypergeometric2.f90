
module hypergeometric2
   
   public :: HYP2F1
   
contains
   
   subroutine HYP2F1(a, b, c, x, sum, n)

! Hypergeometric function 2F1 to fractional accuracy eps.
! Taken from Atlas for computing mathematical functions, William
! J. Thompson, Wiley (1997)
      
      implicit none
      
      integer, parameter :: dp = kind(1.0000D0)
      integer :: i, nn
      integer, intent(in) :: n
      real(dp) :: a,b,c,x,eps
      real(dp),intent(out) :: sum
      real(dp) :: ratio, aks, bks, kfact, cks, tk, xpow, eps10, term
      
      nn = n

200   sum = 1.0d0
      ratio = 10.0d0
      aks = a
      bks = b
      kfact = 1.0d0
      cks = c
      tk = 1.0d0
      xpow = 1.0d0
      eps10 = 0.1*eps
!      n = 1000000
      eps = 1.0e-6
      
!      write(*,*) 'from HYP2F1 n_it =',n
!      write(*,*) 'nn:', nn
      do i = 1, nn 
         if (ratio .gt. eps10) then            
            tk = tk*aks*bks/(kfact*cks)
            xpow = x*xpow
            term = tk*xpow
            sum = sum + term
            aks = aks + 1.0
            bks = bks + 1.0
            kfact = kfact + 1.0
            cks = cks + 1.0
            ratio = abs(term/sum)
         else
!            write(*,*) 'good accuracy'
            goto 400
         end if
      end do

!      write(*,*) 'bad accuracy'
!      write(*,*) 'hypgeo:',sum
      nn = 2*nn
      goto 200
!      sum = 1
400   continue
!write(*,*) 'done hypergeometric function'
!      write(*,*) 'series converged to sufficient accuracy'
      
   end subroutine HYP2F1

! **************************************************************
!
!function hypgeo(a,b,c,z)
!
!   use nrtype
!   use hypgeo_info
!   use nr, ONLY : bsstep,hypdrv,hypser,odeint
!   implicit none
!   
!   complex(spc), intent(in)   :: a,b,c,z
!   complex(spc)               :: hypgeo
!   real(sp), parameter        :: eps = 1.0e-6_sp
!   complex(spc), dimension(2) :: y
!   real(sp), dimension(4)     :: ry
!   
!   
!   if (real(z)**2 + aimag(z)**2 <= 0.25) then
!      call hypser(a,b,c,z,hypgeo,y(2))
!      return
!   else if (real(z) .lt. 0.0) then
!      hypgeo_z0 = cmplx(-0.5_sp,0.0_sp,kind=spc)
!   else if (real(z) .le. 1.0) then
!      hypgeo_z0 = cmplx(0.5_sp,0.0_sp,kind=spc)
!   else
!      hypgeo_z0 = cmplx(0.5_sp,sign(0.5_sp,aimag(z)),kind=spc)
!   end if
!
!   hypgeo_aa = a
!   hypgeo_bb = b
!   hypgeo_cc = c
!   hypgeo_dz =z-hypgeo_z0
!   call hypser(hypgeo_aa,hypgeo_bb,hypgeo_cc,hypgeo_z0,y(1),y(2))
!
!   ry(1:4:2) = real(y)
!   ry(2:4:2) = aimag(y)
!   call odeint(ry,0.0_sp,1.0_sp,eps,0.1_sp,0.0001_sp,hypdrv,bsstep)
!
!   y = cmplx(ry(1:4:2),ry(2:4:2),kind=spc)
!   hypgeo = y(1)
!
!end function hypgeo
!
!subroutine hypser(a,b,c,z,series,deriv)
!
!   use nrtype
!!   use nrutil, only : nrerror
!
!   implicit none
!
!   complex(spc), intent(in)       :: a,b,c,z
!   complex(spc), intent(out)      :: series,deriv
!   integer(I4B)                   :: n
!   integer(I4B), parameter        :: maxit = 1000
!   complex(spc)                   :: aa,bb,cc,fac,temp
!
!   deriv = cmplx(0.0_sp,0.0_sp,kind = spc)
!   fac = cmplx(1.0_sp,0.0_sp,kind = spc)
!   temp = fac
!   aa = a
!   bb = b
!   cc = c
!   do n = 1, maxit
!      fac = ((aa*bb)/cc)*fac
!      deriv = deriv+fac
!      fac = fac*z/n
!      series = temp + fac
!      if (series == temp) return
!      temp = series
!      aa = aa+1.0
!      bb = bb+1.0
!      cc = cc+1.0
!   end do
!
!   call nrerror('hypser: convergence failure')
!
!end subroutine hypser
!
!subroutine hypdrv(s,ry,rdyds)
!   
!   use nrtype
!   use hypgeo_info
!   
!   implicit none
!   
!   real(sp), intent(in)                :: s
!   real(sp), dimension(:), intent(in)  :: ry
!   real(sp), dimension(:), intent(out) :: rdyds
!   complex(spc), dimension(2)          :: y, dyds
!   complex(spc)                        :: z
!
!   y = cmplx(ry(1:4:2),ry(2:4:2),kind = spc)
!   z = hypgeo_z0 + s*hypgeo_dz
!   dyds(1) = y(2)*hypgeo_dz
!   dyds(2) = ((hypgeo_aa*hypgeo_bb)*y(1)-(hypgeo_cc-&
!        ((hypgeo_aa + hypgeo_bb) + 1.0_sp)*z)*y(2))*hypgeo_dz/(z*(1.0_sp - z))
!   rdyds(1:4:2) = real(dyds)
!   rdyds(2:4:2) = aimag(dyds)
!
!end subroutine hypdrv
!
!subroutine nrerror(string)
!
!   character(len=*), intent(in) :: string
!
!   write(*,*) 'nrerror:', string
!   stop 'program terminated by nrerror'
!
!end subroutine nrerror
!
! ---------------------------------------------------------------------

function gammln(xx)

          integer, parameter :: dp = kind(1.0D0)
          real(dp) :: gammln, xx
          integer :: j
          real(kind=dp) :: ser, stp, tmp, x,y,cof(6)
          save :: cof, stp
          data  cof, stp/76.18009172947146d0, -86.50532032941677d0,& 
               24.01409824083091d0, -1.231739572450155d0,&
               0.1208650973866179d-2, -0.5395239384953d-5,&
               2.5066282746310005d0/

          x=xx
          y=x
          tmp=x+5.5d0
          tmp=(x+0.5d0)*log(tmp)-tmp
          ser = 1.000000000190015d0
          do j = 1,6
             y = y+1.d0
             ser = ser +cof(j)/y
          end do
          gammln = tmp+log(stp*ser/x)
          return

       end function gammln

! ----------------------------------------------------------------------

!

































   
end module hypergeometric2



