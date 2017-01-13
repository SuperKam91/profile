module cosmology

   use kind_def
   use sz_globals

   implicit none

   INTEGER, parameter :: nw=501
   DOUBLE PRECISION, dimension(nw) :: xa,za,z2a
 
   interface angdist
      module procedure angdist1
      module procedure angdist2
   end interface

   interface lumdist
      module procedure lumdist1
      module procedure lumdist2
   end interface

   public :: setup_da
   public :: one_over_h
   public :: angdist3
   public :: get_cosmology

contains

  subroutine angdist1

    use kind_def
    use sz_globals

    implicit none
    real(kind=dp) :: z1

    real(kind=dp) :: r,zstep,BigH,D_prop,GothicR,zcurrent,drdz
    integer i,noofpoints
    real(kind=dp), parameter :: limit = 1e-4
    
    z1 = 1.0+z

    if (abs(omegal)<limit) then 
       ! lambda is damn small, and so we'll treat it as zero.

!       write(*,*) 'No lambda, so doing it anayltically'

       q0 = omegam / 2.

       if (q0.ne.0.0) then

! D_theta in Mpc for general q_0. D_l = c/(H0 q^2) (qz+(q-1)((1+2qz)^0.5-1))
! D_theta = D_l/(1+z)^2. NB H0 is in km^-1 Mpc^-1 so need factor of 1000
          D_theta = (const_c/(H0*q0**2))*(q0*z+(q0-1.0)*&
               & (((1+2*q0*z)**0.5)-1.0))/(z1**2)
       else
          D_theta = (const_c/H0)*(z+(z**2)/2.0)/(z1**2)
       end if
       D_theta = D_theta/1000.0D0
    else

! lambda is large, so we're going to do this by integration...
! originally written by WFG. 

! the integration goes from z=0 to z=z

! Da = GothicR sin (r/GothicR)
! GothicR = c/(H0 *sqrt(omegalambda + omegamatter) -1 )
! r = int( dr/dz) wrt dz
! dr/dz = c/((1+z) H0 sqrt (omegam z +1 + omegal ((1/(1+z)^2)-1)))
! Dead straight forward. 

       noofpoints = 1000000*z
       zstep = (z-0.0)/(noofpoints*1.0)

       r = 0.0

       BigH = H0 * 1000.    ! BigH now in m/s/Mpc

       do i=0,(noofpoints-1)
          zcurrent = zstep * i
          z1 = zcurrent + 1.0
          drdz = const_c/(z1*BigH*sqrt(omegam*zcurrent+1+omegal*(z1**(-2.)-1)))
          if ((i==0).or.(i==noofpoints)) then 
             r = r + 0.5 * drdz
          else
             r = r + drdz
          end if
       end do
       r = r * zstep

! that's the integration done.

! Now to take account of the curvature of the universe.

       if ((omegal+omegam-1)>limit) then
          GothicR = const_c/(BigH*sqrt((omegam+omegal)-1))
          D_prop = GothicR * sin(r/GothicR)
       end if
       if ((omegal+omegam-1)<(-1.*limit)) then 
          GothicR = const_c/(BigH*sqrt(1-(omegam+omegal)))
          D_prop = GothicR * sinh(r/GothicR)
       end if
       if (abs(omegal+omegam-1)<limit) then
          ! this is just about flat....
          GothicR = 0
          D_prop = r    ! d'oh!
       end if

       D_theta = D_prop/z1

    end if

  end subroutine angdist1

! ***********************************************************************

  subroutine lumdist1

    use kind_def
    use sz_globals
    implicit none
    real(kind=dp) :: z1
    real(kind=dp), parameter :: limit = 1e-4

    ! For WFG's integration...
    real(kind=dp) :: r,zstep,BigH,D_prop,GothicR,zcurrent,drdz
    integer i,noofpoints

    z1 = 1.0+z

    if (abs(omegal)<limit) then 

! lambda is damn small, and so we'll treat it as zero.
       q0 = omegam / 2.

       if (q0.ne.0.0) then

! D_lum in Mpc for general q_0. D_l = c/(H0 q^2) (qz+(q-1)((1+2qz)^0.5-1))
! NB H0 is in km^-1 Mpc^-1 so need factor of 1000
          D_lum = (const_c/(H0*q0**2))*(q0*z+(q0-1.0)*(((1+2*q0*z)**0.5)-1.0))
       else
          D_lum = (const_c/H0)*(z+(z**2)/2.0)
       end if
       D_lum = D_lum/1000.0D0

    else
! lambda is large, so we're going to do this by integration...
! originally written by WFG. 


! the integration goes from z=0 to z=z

! Da = GothicR sin (r/GothicR)
! GothicR = c/(H0 *sqrt(omegalambda + omegamatter) -1 )
! r = int( dr/dz) wrt dz
! dr/dz = c/((1+z) H0 sqrt (omegam z +1 + omegal ((1/(1+z)^2)-1)))
! Dead straight forward. 

       noofpoints = 1000000*z
       zstep = (z-0.0)/(noofpoints*1.0)

       r = 0.0

       BigH = H0 * 1000.    ! BigH now in m/s/Mpc

       do i=0,(noofpoints-1)
          zcurrent = zstep * i
          z1 = zcurrent + 1.0
          drdz = const_c/(z1*BigH*sqrt(omegam*zcurrent+1+omegal*(z1**(-2.)-1)))
          if ((i==0).or.(i==noofpoints)) then 
             r = r + 0.5 * drdz
          else
             r = r + drdz
          end if
       end do
       r = r * zstep

! that's the integration done.

! Now to take account of the curvature of the universe.

       if ((omegal+omegam-1)>limit) then
          GothicR = const_c/(BigH*sqrt((omegam+omegal)-1))
          D_prop = GothicR * sin(r/GothicR)
       end if
       if ((omegal+omegam-1)<(-1.*limit)) then 
          GothicR = const_c/(BigH*sqrt(1-(omegam+omegal)))
          D_prop = GothicR * sinh(r/GothicR)
       end if
       if (abs(omegal+omegam-1)<limit) then
          ! this is just about flat....
          GothicR = 0
          D_prop = r    ! d'oh!
       end if

       D_lum = D_prop*z1

    end if
        
! (As D_lum is equal to D_theta with a factor of (1+z)^2 between them...


   end subroutine lumdist1

!***************************************************************************

   subroutine angdist2(z0,h,D_th)

      use kind_def
      use sz_globals

      implicit none
      real(kind=dp), intent(in) :: z0,h
      real(kind=dp), intent(out) :: D_th
      real(kind=dp) :: z1

      z1 = 1.0+z0

      if (q0.ne.0.0) then

! D_theta in Mpc for general q_0. D_l = c/(H0 q^2) (qz+(q-1)((1+2qz)^0.5-1))
! D_theta = D_l/(1+z)^2. NB H0 is in km^-1 Mpc^-1 so need factor of 1000
         D_th = (const_c/(h*q0**2))*(q0*z0+(q0-1.0)*&
              & (((1+2*q0*z0)**0.5)-1.0))/(z1**2)
      else
         D_th = (const_c/h)*(z0+(z0**2)/2.0)/(z1**2)
      end if
      D_th = D_th/1000.0D0

   end subroutine angdist2

! ***********************************************************************

   subroutine lumdist2(z0,h,D_l)

      use kind_def
      use sz_globals

      implicit none
      real(kind=dp), intent(in) :: z0, h
      real(kind=dp), intent(out) :: D_l
      real(kind=dp) :: z1

      z1 = 1.0+z0

      if (q0.ne.0.0) then

! D_lum in Mpc for general q_0. D_l = c/(H0 q^2) (qz+(q-1)((1+2qz)^0.5-1))
! NB H0 is in km^-1 Mpc^-1 so need factor of 1000
         D_l = (const_c/(h*q0**2))*(q0*z0+(q0-1.0)*(((1+2*q0*z0)**0.5)-1.0))
      else
         D_l = (const_c/h)*(z0+(z0**2)/2.0)
      end if
      D_l = D_l/1000.0D0

   end subroutine lumdist2

! **********************************************************************

  subroutine setup_da

! Routine from E.Komatsu to tabulate da(x), where x=ln(a)
! Proper angular diameter distance, returned in units of h^-1 Mpc.

    use kind_def
    use sz_globals
    use maths

    implicit none

    integer :: i
    real(dp) :: x,eta,tol=1d-7,ok0

    ok0=1d0-omegam-omegal
    do i=1,nw
       x = (i-nw)*0.01 ! x=ln(a)=-5.00,-4.99,...,0
       xa(i) = x
       eta=rombint(one_over_h,0d0,dexp(-x)-1d0,tol)
       if(abs(ok0)<1.d-5)then
          za(i)=2998d0*exp(x)*eta ! flat model
       else if(ok0>0)then
          za(i)=2998d0*exp(x)*sinh(sqrt(ok0)*eta)/sqrt(ok0)  ! open model
       else
          za(i)=2998d0*exp(x)*sin(sqrt(-ok0)*eta)/sqrt(-ok0) ! closed model
       endif
    enddo
    call spline(xa,za,nw,1.d30,1.d30,z2a) ! evaluate z2a(i)=d^2[da(i)]/dy^2
    return

  end subroutine setup_da

! **********************************************************************

  subroutine angdist3

! proper distance in units of h^-1 Mpc.
! Returns da(z) by interpolating the table generated by Setup_Da. 

     use kind_def
     use sz_globals
     use maths

     implicit none

     real(dp) :: a,b,hh,x
     integer :: jlo

     x=-dlog(1.+z) ! x=ln(a)=-ln(1+z)
     call hunt(xa,nw,x,jlo)
     hh=xa(jlo+1)-xa(jlo)
     a=(xa(jlo+1)-x)/hh
     b=(x-xa(jlo))/hh
     D_theta = (a*za(jlo)+b*za(jlo+1)+((a**3-a)*z2a(jlo)+&
               (b**3-b)*z2a(jlo+1))*(hh**2)/6.)*100.d0/H0
     return

  end subroutine angdist3

! **********************************************************************

  function one_over_h(z_dum)

! one_over_h(z) = H_0/H(z) (dimensionless)
! x=a
     use kind_def
     use sz_globals
     use maths

     implicit none
     real(dp) :: one_over_h
     real(dp), intent(in) :: z_dum
     real(dp) :: x,ok0,or0=0d0 ! radiation density is ignored!
     ok0=1d0-omegam-omegal
     x=1.d0/(1.d0+z_dum)
     one_over_h = 1d0/sqrt(omegam/x**3d0+ok0/x**2d0+or0/x**4d0+&
                           omegal/x**(3d0+3d0*w_de))

  end function one_over_h

! **************************************************************************

   subroutine get_cosmology

      use kind_def
      use sz_globals

      implicit none

      call io_getd('Hubble constant /km s^-1 Mpc^-1:','*',H0,status)
      call io_getd('Omega matter, :','*',omegam,status)
      call io_getd('Omega lambda, :','*',omegal,status) 
      call io_getd('Omega baryon, :','*',omegab,status) 

      call lumdist
      D_lum = D_lum*Mpc         ! D_lum in metres now....
      call angdist
      D_theta = D_theta*Mpc
      rhocrit0 = (3.0*(H0*1.D3/Mpc)**2)/(8.0*pi*G_c)
      write(*,*) 'Luminosity distance is : ',D_lum/Mpc,' Mpc'
      write(*,*) 'Angular distance is : ',D_theta/Mpc,' Mpc'

! Set up table with angular diameter distances
      call setup_da

   end subroutine get_cosmology

! **********************************************************************

end module cosmology
