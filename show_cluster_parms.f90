subroutine show_cluster_parms

   use kind_def
   use sz_globals
   use define_telescope

   implicit none

   write(*,*) 'Current parameters'
   write(*,*) '------------------'
   write(*,*) 
   select case(model_type)
   case('d','D')
      write(*,*) 'Double beta model'
   case('b','B')
      write(*,*) 'Beta model - isothermal'
   case('p','P')
      write(*,*) 'Beta model - Polynomial T model'
   case('x','X')
      write(*,*) 'Beta model - X-ray profile T model'
   case('m','M') 
      write(*,*) 'Beta model - X-ray Map T model'
   case('f','F') 
      write(*,*) 'File with n and T profiles'
   case('h','H')
      write(*,*) 'Beta model - temperature (H)alo'
   case('g','G') 
      write(*,*) 'GNFW model'
   end select

   select case(model_type)
   case('g','G')
      write(*,*) 'GNFW alpha parameter:',a_GNFW
      write(*,*) 'GNFW beta parameter:',b_GNFW
      write(*,*) 'GNFW gamma parameter:',c_GNFW
      write(*,*) 'GNFW concentration parameter c_500:',c500_GNFW
      write(*,*) 'Y_tot (arcmin^2):',Ytot
      write(*,*) 'theta_s (arcmin):',theta_s
   case DEFAULT
      write(*,*) 'Central number density /cm^-3:     ', n0
      write(*,*) 'Beta parameter:                    ', beta
      write(*,*) 'Electron Temperature: (K):         ', Te
      if (ellipsoid) then
         write(*,*) 'Core Radius Axis 1:                ', theta_c(1)
         write(*,*) 'Core Radius Axis 2:                ', theta_c(2)
         write(*,*) 'Core Radius Axis 3:                ', theta_c(3)
!      write(*,*) 'Euler angle 1:                     ', ang(1)/deg2rad
!      write(*,*) 'Euler angle 2:                     ', ang(2)/deg2rad
!      write(*,*) 'Euler angle 3:                     ', ang(3)/deg2rad
         write(*,*) 'Angle of ellipsoid:                ', el_alpha/deg2rad
      else
         write(*,*) 'Core Radius:                       ', theta_cx
      end if
   end select
   write(*,*)
   write(*,*) 'Hubble constant /km s^-1 Mpc^-1:   ', H0
   write(*,*) 'Redshift:                          ', z
   write(*,*) 'Exposure Time (s):                 ', x_exp
   write(*,*) 'X-ray emission constant:           ', const
   write(*,*) 
   which_sz: select case (which_tel)
      case('c','C')
         write(*,*) 'Observing with the CAT telescope at ',obsfreq   
      case('v','V')
         write(*,*) 'Observing with the VLA telescope at ',obsfreq
      case('s','S')
         write(*,*) 'Observing with the VSA source subtractor at ',obsfreq
      case('m','M')
         write(*,*) 'Observing with the VSA main array at ',obsfreq
      case('e','E')
         write(*,*) 'Observing with the extended VSA main array at ',obsfreq 
      case('d','D') 
         call show_telescope
      case default
         write(*,*) 'Observing with the Ryle telescope at ',obsfreq,' Hz'
   end select which_sz
   write(*,*)
   which_xray: select case (x_tel)
   case('c','C')   
      write(*,*) 'Xray telescope = Chandra'
   case('p','P')
      write(*,*) 'Xray telescope = PSPC on ROSAT'
   case('h','H')
      write(*,*) 'Xray telescope = HRI on ROSAT'
   case('x1','X1')
      write(*,*) 'Xray telescope = mos1 on XMM'
   case('x2','X2')
      write(*,*) 'Xray telescope = mos2 on XMM'
   case('xp','XP')
      write(*,*) 'Xray telescope = pn on XMM'
   end select which_xray

   write(*,*) 'X-ray background level:            ', back 

end subroutine show_cluster_parms
