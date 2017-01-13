
subroutine plot_scaling

   use kind_def
   use sz_globals
   use kgraph
   use cosmology
   use map_handling

   implicit none

   real(kind=dp) :: z_l

   integer :: i,j, centre, ci, nloop, box500x, box500y
   character*1 :: chr1
   real(kind=dp) ::  z_s, z_ls, D_s, D_l, D_ls,crit,Hz
   real(kind=dp), dimension(4) :: mu
   real(kind=dp) k, delta, pix2metres, dsq
   real(kind=dp) :: av_mass, n_av, r_500
   real(kind=dp), dimension(2) :: alpha, alpha1, alpha2
   real(kind=dp) :: dt, M, radius, lim_rad,gauss
   real(kind=dp) :: im_x, im_y, tot_mass,ratio,omegaz,delz
   real(kind=dp), dimension(maxsize,maxsize) :: S,den
   real(kind=dp) :: r,M_200,r_200,rho_200,rho_200_dash,rho_500,rho_500_dash
   integer, parameter   :: n_it = 1000
   real(dp)             :: hypgeo,M_bar,n_200,n_500,M_500,A,al,T_new,h_70
   logical :: logg
   
   write(*,*) 'WElcome to PlOT-scALING.'

   call io_getc('Log plot?','y',chr1,status)
   if (chr1.eq.'y') logg = .true.

   if (logg) then
      call pgenv(log(0.01),log(2.),log(0.01),log(3.),0,1)
   else
      call pgenv(0.,2.0,0.,3.0,0,1)
   end if
   call pglab('redshift/z','R\d500\u [h\d70\u\u-1\d Mpc]','')

   write(*,*) 'COSMOLOGY:'
   call io_getd('Hubble constant:','70.',H0,status)
   call io_getd('Omega matter:','0.3',omegam,status)
100 call io_getc('(E)deS Universe or (S)pherical collapse?','E',chr1,status)

   al = 0.74
   A = 1.5
   WRITE(*,*) 'M-T SCALING RELATION:'
   call io_getd('A:','1.5',A,status)
   call io_getd('alpha:','0.74',al,status)      


   do i = 1,20
      z_l = 0.1+0.1*real(i-1)
      Hz = H0*(omegam*(1+z_l)**3 + 1. - omegam)**(0.5)
      omegaz = omegam*(1+z_l)**3*H0/Hz
      
      if (chr1.eq.'E') then 
         delz = 500.
      elseif (chr1.eq.'S') then
         delz = 500.*(1.+82.*(omegaz-1.)/(18*pi**2) - 39.*(omegaz-1.)**2/(18.*pi**2))
      else
         write(*,*) 'What kind of universe is that?'
         goto 100
      end if

! The critical surface density is c^2/(4 pi G D)
      Hz = Hz*(1.4e-18/70.) ! km s^-1 Mpc^-1 -> s^-1
      crit = 3*Hz**2/(8*pi*G_c)
            
      rho_200 = crit*200.
      rho_500 = crit*500.

! baryonic mass density - n_e*mu_e*m_p (Mason & Myers 2000)
      rho_200_dash = rho_200/(m_p*1.146*10**6)
      rho_500_dash = rho_500/(m_p*1.146*10**6)
      
      Hz = Hz/(1.4e-18/70.)
      h_70 = Hz/70. 

      open(74,file='scaling_10.dat')
! X-ray temperature from scaling relation (Ettori 2004):
      do j = 1,10
         T_new = real(j)
         M_500 = (T_new/(6.*(10**(-1.*al/A))))**A*(H0/Hz)
         M_500 = M_500*(h_70**(-1.))*m_0*10d0**14
         r_500 = (3.*M_500/(4*pi*crit*delz))**(1./3.)
         r_500 = r_500/(Mpc)

         if (j.eq.10) write(74,*) z_l,r_500
         gauss = 2.152*exp(-1.*z_l**0.971/(2.*0.805**2))
     
         call pgsci(j)
         if (logg) then
            call pgpt1(log(real(z_l)),log(real(r_500)),2)
         else
            call pgpt1(real(z_l),real(r_500),2)
!            call pgsci(2)
!            call pgpt1(real(z_l),real(gauss),2)
         end if
      end do
   end do
         
   close(unit=74)

end subroutine plot_scaling

