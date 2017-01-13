! Define cosmology
! 06/08/09 Removed from get_cluster_parameters [kg]
! Now moved to cosmology module

subroutine get_cosmology

   use kind_def
   use sz_globals
   use cosmology

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

end subroutine get_cosmology






