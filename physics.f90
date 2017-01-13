module physics

   use kind_def
   use sz_globals
   use maths

   implicit none

   public :: planck_fn
   public :: rayleigh_jeans_fn
   public :: g_nu
   public :: g_nu_x2
   public :: g_nu_thermo
   public :: rel_corr_sz
   public :: dB_by_dT

contains

!***************************************************************************
   function planck_fn(Tmpr,freq)result(res)

! Calculates the conversion factor between temperature and surface brightness
! for the Planck Black Body spectrum
! B_nu = 2 h freq^3/c^2 (exp x -1) = I_0 x^3 / (exp x -1)

      implicit none

      real(kind=dp), intent(in) :: Tmpr, freq 
      real(kind=dp) :: x
      real(kind=dp) :: res

      x = const_h*freq/(k_b*Tmpr)
!      res = ((2.0*const_h*freq**3)/const_c2)/(exp(x)-1.0)
      res = I_0*x**3/(exp(x)-1.0)

   end function planck_fn

!***************************************************************************

   function rayleigh_jeans_fn(freq)result(res)

! Calculates the conversion factor between temperature and surface brightness
! for the Rayleigh-Jeans spectrum. Multiply this by Delta T to find the 
! intensity spectrum of a small temperature fluctuation.
! I = 2kT / lambda^2 
!   = I_0 x^2 (for T0)

      implicit none

      real(kind=dp), intent(in) :: freq 
      real(kind=dp) :: x
      real(kind=dp) :: res

!      x = const_h*freq/(k_b*Tmpr)
!      res = I_0*x**2

      res = 2.d0*k_b/(const_c/freq)**2

   end function rayleigh_jeans_fn

! **************************************************************************

   function g_nu(freq)result(res)

! Calculates the g_nu function used in SZ formula
! I_SZ = g_x * I_0 * y

      real(kind=dp), intent(in) :: freq 
      real(kind=dp) :: x
      real(kind=dp) :: res

      x = const_h*freq/(k_b*T0)
      res = ((x**4*exp(x)/((exp(x)-1)**2))*(x*coth(x/2.0)-4.0))

   end function g_nu

! **************************************************************************

   function g_nu_x2(freq)result(res)

! Calculates the g_nu/x^2 function used in determining the relationship between
! fractional BRIGHTNESS temperature decrement and the Comptonization y 
! parameter. This equals -2 in the R-J limit. 
! DeltaT_SZ_brightness = g_nu_x2 * y * T0
! (NB often DeltaT_SZ refers to the thermodynamic temperature decrement)
! DeltaT_SZ_thermo = (xcothx/2-4) * y * T0

      real(kind=dp), intent(in) :: freq 
      real(kind=dp) :: x
      real(kind=dp) :: res

      x = const_h*freq/(k_b*T0)
      res = ((x**4*exp(x)/((exp(x)-1)**2)))
      res = res*((x*coth(x/2.0)-4.0))/x**2

   end function g_nu_x2

! *************************************************************************

   function g_nu_thermo(freq)result(res)

! Calculates the g_nu/(dB/dT) function used in determining the relationship 
! between fractional THERMODYNAMIC temperature decrement and the
! Comptonization y parameter. This equals -2 in the R-J limit. 
! DeltaT_SZ_thermo = (xcothx/2-4) * y * T0

      real(kind=dp), intent(in) :: freq 
      real(kind=dp) :: x
      real(kind=dp) :: res

      x = const_h*freq/(k_b*T0)
      res = (x*coth(x/2.0)-4.0)

   end function g_nu_thermo

! *************************************************************************

   function rel_corr_sz(freq)result(res)

! Calculates the intensity of the S-Z effect making first and second order 
! corrections due to relativistic effect (taken from Anthony Challinor's code)
! Multiply by y to obtain answer in W m^(-2) Hz^(-1) sr^(-1)

      real(kind=dp), intent(in) :: freq 
      real(kind=dp) :: x, theta_e, zero_term, one_term, two_term, st, ct
      real(kind=dp) :: res

      if (freq.ne.0.0) then

         x = const_h*freq/(k_b*T0)
         theta_e = Te*k_b/(m_e*const_c2)
         ct=x*(exp(x)+1.d0)/(exp(x)-1.d0)
         st=2.d0*x/(exp(x/2.0d0)-exp(-1.0d0*x/2.0d0))

! The zero-order effect
         zero_term=ct-4.d0

! The first-order effect
         one_term=-10.d0+47.d0/2.d0*ct-42.d0/5.d0*ct**2+&
              7.d0/10.d0*ct**3+7.d0/5.d0*st**2*(ct-3.d0)

! The second-order effect
         two_term=-15.d0/2.d0+1023.d0/8.d0*ct-868.d0/5.d0*ct**2+&
              329.d0/5.d0*ct**3-44.d0/5.d0*ct**4+11.d0/30.d0*ct**5+&
              1.d0/30.0d0*st**2*(-2604.d0+3948.d0*ct-1452.d0*ct**2+&
              143.d0*ct**3)+1.d0/60.d0*st**4*(-528.d0+187.d0*ct)

         res = x**4*exp(x)/(exp(x)-1.d0)**2* &
              (2.0*(k_B*T0)**3/(const_h**2*const_c2))* &
              (zero_term+theta_e*one_term+theta_e**2*two_term)
      else
         res = 0.0
      end if

   end function rel_corr_sz

!***************************************************************************

   function dB_by_dT(Tmpr,freq)result(res)

! Calculates the differential of the Planck function with respect to 
! temperature. Multiply this by Delta T to find the intensity spectrum of
! a small temperature fluctuation.
! dB/dT = x exp x /(exp x - 1) B_nu (T) / T
!       = x^4 exp x /(exp x - 1)^2 I_0 / T0    for T=T0
      implicit none

      real(kind=dp), intent(in) :: Tmpr, freq 
      real(kind=dp) :: x
      real(kind=dp) :: res

      x = const_h*freq/(k_b*Tmpr)
      res = (x*exp(x)/(exp(x)-1.0))*planck_fn(Tmpr,freq)/Tmpr

   end function dB_by_dT

! *************************************************************************

end module physics
