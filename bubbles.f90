
subroutine bubbles

   use kind_def
   use sz_globals
   use maths
   use cosmology
   use physics

   implicit none

   integer  :: i,j,k,xb1,xb2,yb1,yb2
   real(dp)  :: rc(maxsize,maxsize),rb1(maxsize,maxsize)
   real(dp)  :: rb2(maxsize,maxsize)
   real(dp)  :: nb1,x,xx,x1,ab1,a1
   real(dp)  :: th,dl,k_sz

   real(dp)  :: B
   external  B

     
   call io_getd('Hubble constant /km s^-1 Mpc^-1:','*',H0,status)
   call io_getd('Omega matter, :','*',omegam,status)
   call io_getd('Omega lambda, :','*',omegal,status)
   call io_getd('Redshift:','*',z,status)
   
   call angdist

   n0 = 1d-2
   nb1 = 1d-3
   Te = 10.*keVtoK
   a1 = 180.
   dl = a1*sec2rad*D_theta*Mpc*1d4 ! integral in m
   th = n0*k_b*Te*sigma_T/(m_e*const_c**2)
   th = -g_nu_x2(nu(chan)*1d9)*th*dl*T0
   write(*,*) n0,k_b,Te,sigma_T,m_e,const_c,T0,-g_nu_x2(nu(chan)*1d9),dl
   k_sz = -g_nu_x2(nu(chan)*1d9)*T0*k_b*Te
   k_sz = k_sz*dl
   k_sz = k_sz*sigma_T/(m_e*const_c2)
   write(*,*) k_sz
   x = const_h*nu(chan)*1e9/(k_b*T0)
   x1 = 0.67 ! beta
   a1 = 180.
   ab1 = 25.
   xb1 = 25+maxsize/2+1
   yb1 = 25+maxsize/2+1
   xb2 = 0
   yb2 = 0
   

   I = 0.
   szsky = 0.
   do i = 1,maxsize
      do j = 1,maxsize
         rc(i,j) = sqrt(real((i-maxsize/2-1)**2+(j-maxsize/2-1)**2))*cellsize
         rb1(i,j) = sqrt(real((i-xb1)**2+(j-yb1)**2))*cellsize
         rb2(i,j) = sqrt(real((i-xb2)**2+(j-yb2)**2))*cellsize
         szsky(i,j,chan) = th*(1+(rc(i,j)**2/a1**2))**(-3*x1+0.5)*Bfnc(3.0d0*x1-0.5d0, 0.5d0)
         szsky(i,j,chan) = szsky(i,j,chan) &
              - sigma_T*nb1*dl*(x**3/(exp(x)-1))&
              *(1+(rb1(i,j)**2/ab1**2))**(-3*x1+0.5)*Bfnc(3.0d0*x1-0.5d0,0.5d0)
         
      end do
   end do

 

end subroutine bubbles

! **********************************************************************

