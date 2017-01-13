! This routine uses the SZ model to produce a 2-d grav. potential and hence 
! solve for the lensing.
! The relevent equations are 
!
! delsq psi(b) = 8 pi G S(b) 
!
! where psi is the potential, b is the distance from the lens centre in the 
! lens plane and S is the surface density of the lens, 
!
! alpha_ = grad psi / c^2
!
! where alpha_ is the deflection at the lens plane
!
! alpha(theta) = (D_ls/D_s) alpha_(D_l * theta)
!
! relates the deflection at the lens plane to the observed position alpha,
! and
!
! beta = theta - alpha(theta)
!
! relates the 'true' source position beta to the observed position theta
! (relative to the centre of the lens)
!
! The time delay is given by
!
!     t = (1+z_l)[ D alpha_^2 / 2c - psi(theta)/c^3]
!
! where D = D_l D_s / D_ls
!
! The magnification of the image is given by det(mu) where 
!
!     mu = d theta / d beta = inv[I - d alpha / d theta}
!
!        =  ( 1 - d ax/d tx       - d ax/d ty )^ -1
!           (   - d ay/d tx     1 - d ay/d ty )
!
subroutine lensing(z_l)

   use kind_def
   use sz_globals
   use kgraph
   use cosmology
   use map_handling
!   use hypergeometric2

   implicit none

   real(kind=dp) :: z_l

   integer :: i, j, centre, ci, nloop, box500x, box500y, n_lim_rad
   character*1 :: chr1
   real(kind=dp) ::  z_s, z_ls, D_s, D_l, D_ls,crit,Hz
   real(kind=dp), dimension(4) :: mu
   real(kind=dp) k, delta, pix2metres, dsq
   real(kind=dp) :: av_mass, n_av, r_500
   real(kind=dp), dimension(2) :: alpha, alpha1, alpha2
   real(kind=dp) :: dt, M, radius, lim_rad
   real(kind=dp) :: im_x, im_y, tot_mass,ratio,n_c
   real(kind=dp), dimension(maxsize,maxsize) :: S,den
   real(kind=dp) :: r,M_200,r_200,rho_200,rho_200_dash,rho_500,rho_500_dash
   integer, parameter   :: n_it = 1000
   real(dp)             :: hypgeo,M_bar,n_200,n_500,M_500,A,al,T_new,h_70
      
   centre = maxsize/2 + 1
   call angdist
   write(*,*) 'Angular diameter distance:', D_theta
      
   call io_getd('z of source?:','3.8',z_s,status)
   z_ls = ((1+z_s)/(1+z_l)) - 1
   write(*,*) 'z_ls = ', z_ls
   Hz = H0*(omegam*(1+z_l)**3 + 1. - omegam)**(0.5)
   call angdist(z_l, H0, D_l)
   call angdist(z_s, H0, D_s)
   call angdist(z_ls, Hz, D_ls)
   write(*,*) 'D_l, D_s, D_ls = ',D_l, D_s, D_ls,' Mpc'

   delta = cellsize*sec2rad*D_theta
   write(*,*) 'Pixelsize at lens = ',delta,' Mpc'
   pix2metres = delta*Mpc
   dsq = pix2metres**2
   delta = delta/cellsize
   write(*,*) ' = ', delta, ' Mpc per arcsec'

   call io_getd('Ratio of total mass to gas mass?:','10',ratio,status)

!     nsky is in cm^-3 x pixels: so x 10^6 to get m^-3, then x by pixel
!     size in m and the proton mass to S in get kg m^-2
   k = ratio*m_p*pix2metres*1.0D6

! nsky now defined differently in make_skies
!   S = nsky*k
!   den = nsky*k

   S = nsky*m_p*ratio/pix2metres**2
   den = nsky*m_p/pix2metres**2
   tot_mass = sum(S)
   write(*,*) 'Max and Min:',maxval(S),minval(S)

! Adjust total mass to take into account that lengths are in Mpc and then 
! divide by solar mass
   tot_mass = tot_mass*Mpc**2*(delta*cellsize)**2/m_0
   write(*,*) 'Total mass in lens = ',tot_mass, ' Mo'
   
! The critical surface density is c^2/(4 pi G D)
   crit = const_c2/(D_l*4.0D0*pi*G_c*Mpc)
   write(*,*) 'Critical surface density ',crit,' kg m^-2'
   
! r_200:
   Hz = Hz*(1.4e-18/70.) ! km s^-1 Mpc^-1 -> s^-1
   crit = 3*Hz**2/(8*pi*G_c)
   
   !crit = 3.45e-27
   write(*,*) 'Critical density at redshift', real(z),' = ',real(crit),' kg m^-3'
   
   write(*,*) '                                        = ',real(crit/(1.146*m_p*10**6)),' cm^-3'

   call io_getc('(c)alculate r200/r500 or from (f)ile?','f',chr1,status)
   if (chr1=='c') then
      rho_200 = crit*200.
      rho_500 = crit*500.
! baryonic mass density - n_e*mu_e*m_p (Mason & Myers 2000)
      rho_200_dash = rho_200/(m_p*1.146*10**6)
      rho_500_dash = rho_500/(m_p*1.146*10**6)
      do i = 1,maxsize
         do j = 1,maxsize
            if (nsky(i,j).ge.rho_200_dash) then
               r_200 = sqrt(real((i-(maxsize/2))**2+(j-maxsize/2)**2))
               exit
            end if
         end do
      end do
      write(*,*) 'R_200 = ',r_200*(delta*cellsize),'Mpc'
      do i = 1,maxsize
         do j = 1,maxsize
            if (nsky(i,j).ge.rho_500_dash) then
               r_500 = sqrt(real((i-(maxsize/2))**2+(j-maxsize/2)**2))
               exit
            end if
         end do
      end do
      write(*,*) 'R_500 = ',r_500*(delta*cellsize),'Mpc'
      
   elseif (chr1=='f') then

      call io_getd('R_200:','*',r_200,status)
      call io_getd('R_500:','*',r_500,status)
      
      write(*,*) 'R_200 = ',r_200,'Mpc =>',(r_200/D_theta)*(180./pi)*60.*60.,'arcsec.' 
      r_200 = (r_200/D_theta)*(180./pi)*60.*60. ! arcsec
      r_200 = r_200/cellsize ! pixels
      write(*,*) 'R_500 = ',r_500,'Mpc =>',(r_500/D_theta)*(180./pi)*60.*60.,'arcsec.' 
      r_500 = (r_500/D_theta)*(180./pi)*60.*60. ! arcsec
      r_500 = r_500/cellsize ! pixels
   end if

   
   M_200 = 0.
   n_200 = 0.
   M_500 = 0.
   n_500 = 0.
   do i = 1,maxsize
      do j = 1,maxsize
         r = sqrt(real((i-(maxsize/2))**2+(j-(maxsize/2))**2))
         if (r.lt.r_200) then
            n_200 = n_200+nsky(i,j)
         end if
         if (r.lt.r_500) then
            n_500 = n_500+nsky(i,j)
         end if
         if (r.lt.theta_c(1)) then
            n_c = n_c+nsky(i,j)
         end if
      end do
   end do
   
   !call HYP2F1(1.5d0,3.*beta/2.,2.5d0,-1.*(r_200/theta_c(1))**2,hypgeo, n_it) 
   !M_bar = (4./3.)*pi*n0*m_p*1.146*hypgeo
   !M_bar = M_bar*(r_200*Mpc*(delta*cellsize)*10**6)**3
   M_bar = n_200*1.146*m_p ! cm^-3 -> kg cm^-3
   M_bar = M_bar*dsq*1.0D6*pix2metres
!   M_bar = M_bar*pix2metres*1.0D6*Mpc**2*(delta*cellsize)**2
   M_500 = n_500*1.146*m_p*dsq*1.0D6*pix2metres
!   M_500 = n_500*1.146*m_p*pix2metres*1.0D6*Mpc**2*(delta*cellsize)**2
   
   write(*,*) 'r_200 Baryonic mass:',M_bar/m_0, ' Mo'
   write(*,*) '=> Total mass:',ratio*M_bar/m_0, ' Mo'
   write(*,*) 'r_500 Baryonic mass:',M_500/m_0, ' Mo'
   write(*,*) '=> Total mass:',ratio*M_500/m_0, ' Mo'
   write(*,*) 'r_c Baryonic mass:',n_c*1.146*m_p*dsq*1.0D6*pix2metres/m_0, ' Mo'
   write(*,*) '=> Total mass:',ratio*n_c*1.146*m_p*dsq*1.0D6*pix2metres/m_0, ' Mo'

   call io_getc('Temporary output?','y',chr1,status)
   if (chr1=='y') then
      open(13,file = 'cl_temp.dat')
      call io_getd('A:','1.5',A,status)
      call io_getd('alpha:','0.74',al,status)
      ! X-ray temperature from scaling relation (Ettori 2004):
      Hz = Hz/(1.4e-18/70.)
      h_70 = Hz/70.
      M_500 = M_500/(h_70)**(-1)
      write(*,*) 'M_500:',M_500,' 10^14 h_70^-1 Mo'

      T_new = 6.*(10**(-1.*al/A))*((Hz/H0)*M_500)**(1./A)
      write(13,*) T_new
   end if
   if (chr1=='n') then
      open(13,file='cl_mass.dat',POSITION='APPEND')
      write(13,*) M_bar/m_0, ratio*M_bar/m_0
   end if
   close(unit=13)
   return


! r_500 (Mohr 1999)
!   r_500 = 1.185*(100./H0)*(Te/(10.*keVtoK))**(0.5)
   
   write(*,*) 'Surface density in kg m^-2:'

   call display_map(S,maxsize,graph_adj,cellsize,cellsize,0D0,0D0)


!     now want delta in m/arcsec not Mpc/arcsec
   delta = delta*Mpc

!     now want (8 pi G M) / (2 pi c^2 ) = 4 G M/c^2
   k = 4.0D0*(delta*cellsize)**2*G_c/const_c2
   S = S*k

!     for a given position on the sky, find the true source position
   ci = 1
   source_loop: do
      call io_getd('Image position (RA, arcsec)','0',im_x,status)
      call io_getd('Image position (dec, arcsec)','0',im_y,status)

      ci = ci + 1
      call pgsci(ci)
      call pgpoint(1, real(centre+(im_x/cellsize)),&
           &        real(centre+(im_y/cellsize)), -3)

      call deflect(S,im_x,im_y,pix2metres,D_ls,D_s,alpha)

      write(*,*) 'Source position = ',im_x-alpha(1), im_y-alpha(2)
      call pgpoint(1, real(centre+((im_x-alpha(1))/cellsize)),&
           &        real(centre+((im_y-alpha(2))/cellsize)), -4)

!     now find the magnification
      dt = cellsize
      call deflect(S,im_x+dt,im_y,pix2metres,D_ls,D_s,alpha1)
      call deflect(S,im_x,im_y+dt,pix2metres,D_ls,D_s,alpha2)

      mu(1) = 1 - (alpha1(1) - alpha(1))/dt
      mu(2) = 0 - (alpha2(1) - alpha(1))/dt
      mu(3) = 0 - (alpha1(2) - alpha(2))/dt
      mu(4) = 1 - (alpha2(2) - alpha(2))/dt

      M = 1/(mu(1)*mu(4)-mu(3)*mu(2))
      if (M.lt.0) M = -M
      write(*,*) 'Magnification tensor = ( ',M*mu(4), -M*mu(3),' )'
      write(*,*) '                       ( ',-M*mu(2), M*mu(1),' )'
      write(*,*) 'Flux magnification is ', M

      call io_getc('Another image?:','y',chr1,status)
      if ((chr1=='n').or.(chr1=='N')) then
         exit source_loop
      end if
   end do source_loop
   
! reset pen colour
   call pgsci(plot_col)

   call io_getc('Determine mass in a region?:','y',chr1,status)
   if (chr1.eq.'y') then
      call io_getc('Determine mass within a radius or with a box:',&
           & 'b',chr1,status)
      if (chr1=='b') then
         write(*,*) 'Define region in which to find mass'

! Determine which area of the map to use to find mass
         call def_map_region

         tot_mass = sum(S,mask=ok_data)
         tot_mass = tot_mass*const_c2/(4.0D0*m_0*G_c)
         write(*,*) 'Total mass in region = ',tot_mass, ' Mo'
      else if (chr1=='r') then
         tot_mass = 0.0
         centre = maxsize/2+1
         call io_getd('Radius in which to integrate (arcsec)','300.0',&
              & lim_rad,status)
         lim_rad = lim_rad/cellsize
         n_lim_rad = lim_rad
         do nloop = 1, n_lim_rad
            tot_mass = 0.0
            av_mass = 0.0
            n_av = 0.
            do j = 1,maxsize
               do i = 1,maxsize
                  radius = sqrt(((j-centre)**2+(i-centre)**2)*1.0)
                  if (radius.le.nloop) then
                     tot_mass = tot_mass+S(i,j)
                     av_mass = av_mass+den(i,j)
                     n_av = n_av+1.
                  end if
               end do
            end do
            tot_mass = tot_mass*const_c2/(4.0D0*m_0*G_c)
            av_mass = av_mass/n_av
            !if (av_mass.le.200.*crit) then
            !   write(*,*) 'r_200 at',nloop*cellsize,' arcsec.'
            !end if
            write(*,*) 'Total mass in radius = ',nloop*cellsize,&
                        ' is ',tot_mass, ' Mo'
         end do
      end if

   end if

end subroutine lensing

!--------------------------------------------------------------------------
! calculate deflection due to mass distribution S, which is in (4G/c^2) kg
! im_x and im_y are the image position in arcsec
! pix2metres is the pixel size in metres
! cellsize is the pixel size in arcsec
! alpha is the image displacement in arcsec

subroutine deflect(S,im_x,im_y,pix2metres,D_ls,D_s,alpha)

   use kind_def
   use sz_globals

   implicit none

   real(kind=dp), dimension(maxsize,maxsize) :: S
   real(kind=dp) :: im_x, im_y, pix2metres
   real(kind=dp) :: D_ls, D_s
   real(kind=dp), dimension(2) :: alpha

   integer :: i, j, centre
   real(kind=dp) ::  x, y, k, m_x, m_y

   m_x = im_x/cellsize
   m_y = im_y/cellsize

   centre = maxsize/2+1
   alpha(1) = 0.0D0
   alpha(2) = 0.0D0

   do j = 1, maxsize
      do i = 1, maxsize
         x = dble(i-centre) 
         y = dble(j-centre)
         if ((x.ne.m_x).or.(y.ne.m_y)) then
            k = S(i,j)/((x-m_x)**2 + (y-m_y)**2)
         end if
         alpha(1) = alpha(1)+(k*(m_x-x))
         alpha(2) = alpha(2)+(k*(m_y-y))
      end do
   end do
   alpha(1) = alpha(1)/pix2metres
   alpha(2) = alpha(2)/pix2metres
   alpha(1) = alpha(1)*(D_ls/D_s)/sec2rad
   alpha(2) = alpha(2)*(D_ls/D_s)/sec2rad

end subroutine deflect

