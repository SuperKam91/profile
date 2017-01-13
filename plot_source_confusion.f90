! estimate the source confusion 'power spectrum' for ami
! Routine provided by Michael Brown

subroutine plot_source_confusion

   use kind_def
   use sz_globals
   use maths
   use kgraph

  implicit none
  integer :: i
  real(kind=dp), parameter :: eps = 1.d-3
  real(kind=dp) :: Slim, ell_step
  real(kind=dp) :: expterm, prefact, dBdT, lim1, lim2, Cl
  real(kind=dp) :: integrand, ell
  integer, parameter :: max_ell = 10000
  real(kind=dp), dimension(max_ell) :: p_confus
  external :: integrand
  
  Slim = 300.d-6
! dB/dT @ 15 GHz and 2.726 K
  expterm = exp(const_h*obsfreq/(k_b*T0))
  prefact = 2.d0*(const_h*obsfreq**2.d0)**2.d0/(const_c**2.d0*k_b*T0**2.d0)
  dBdT = prefact*expterm/(expterm-1.d0)**2.d0 !Joule s^-1 m^-2 sr^-1 Hz^-1 K^-1
  dBdT = 1.e26*dBdT  ! Jansky sr^-1 Hz^-1 K^-1 
  print*, 'dBdT: ', dBdT, '(Jy sr^-1 Hz^-1 K^-1)'

! do the integral
  lim1 = 0.d0
  lim2 = Slim
  print*, 'integration limits:', lim1, lim2, 'Jansky'
  call qtrap(integrand, lim1, lim2, eps, Cl)
  Cl = 1.d12 * Cl / dBdT / dBdT ! Jy^2 -> Delta_T^2 and then K^2 -> micro-K^2 

! print answer
  print*, 'Slim = ', int(1.e6 * Slim), 'micro-Jy' 
  print*, 'Power is', Cl, 'micro-K^2' 

  do i = 2, max_ell
     ell = real(i)
     p_confus(i) = ell*(ell+1.d0)*Cl/(2.d0*pi)
  end do
  
  ell_step = 1.d0
  call display_profile(p_confus,max_ell,ell_step)

  ! all done
end subroutine plot_source_confusion

!-----------------------------------------------------------------------------!

function integrand(S)

  implicit none
  real(8), parameter :: prefact = 376.d0, alpha = -1.80
  real(8) :: S, dNdS, integrand

  if (S .eq. 0.d0) then 
     integrand = 0.d0
  else
     dNdS = prefact * S**alpha
     integrand = S**2.d0 * dNdS 
  end if
     
end function integrand

