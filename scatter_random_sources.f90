! Puts random sources onto S-Z sky
! 7/10/03 Source count normalisation changed by TLC to accommodate change 
! from Jy^-1 Sr^-1 as in papers to muJy^-1 arcmin^-2. Numbers entered in
! profile exactly as before, problem was in muJy conversion to Jy being 
! raised to (1-spectral index), so EMW and ACT counts didn't match before, 
! think its OK now. TLC has all this written down if anyone cares.
!
! 19/05/09 Following a great deal of confusion the code has been amended to
! deal with source counts in Jy^-1 Sr^-1 and removing all of Tom's hacks. [KG]
subroutine scatter_random_sources&
         (k_parm,gamma,s1,s2,alpha_mean,alpha_rms,to_file,src_file,iunit,idum)

   use kind_def
   use sz_globals
   use maths

   implicit none

   real(kind=dp),intent(in) :: k_parm, gamma, s1, s2, alpha_mean, alpha_rms
   integer,intent(in) :: iunit
   real(kind=dp) :: expected_n, field_size, dummy_dp, alpha, ra, dec
   real(kind=dp) :: s_density, ra2, dec2
   integer :: i, j, k, idum, actual_n
   integer :: tab_res, x1, src_x_pos, src_y_pos, n_levs
   integer :: write_time(3)
   character :: chra*16, chdec*16
   integer :: prra, prdec, lr, ld
   real(kind=dp), dimension(:),allocatable :: pd_s_p, pd_p_s
   real(kind=dp), dimension(:),allocatable :: s0, p0
   real(kind=dp) :: nsig, s_source, t_source, ran_index
   real(kind=dp), dimension(:),allocatable :: lim
   integer, dimension(:),allocatable :: n_lim
   logical :: report, to_file, src_file

   report = .false.
   prra = 1
   prdec = 0

! Check gamma choice won't cause a crash
   if (gamma.eq.1) then
      write(*,*) 'Inappropriate choice of gamma'
      return
   end if

! Allocate tables of probabilities
   tab_res = 10000
   allocate (pd_s_p(tab_res))
   allocate (pd_p_s(tab_res))
   allocate (s0(tab_res))
   allocate (p0(tab_res))

! Allocate arrays for binning
   n_levs = ((log(s2/s1))/log(2.0))+2
   allocate (lim(n_levs))
   allocate (n_lim(n_levs))

   do i = 1, n_levs
      n_lim(i) = 0
      lim(i) = s1*(2**(i-1))
   end do
   lim(n_levs) = s2

! Calculate field size in arcmin^2
   field_size = (maxsize*cellsize)**2

! Convert to steradians
   field_size = field_size*sec2rad**2

! Find expected number of sources in field and source density 
   expected_n = k_parm*field_size/(1-gamma)*(s2**(1-gamma)-s1**(1-gamma))
   s_density = expected_n/field_size*(deg2rad/60.0d0)**2

   if (verbose) write(*,*) 'Field size is ',field_size, ' sterad'
   if (verbose) write(*,*) 'Expected number of sources in field is ',expected_n
   if (verbose) write(*,*) 'Expected number of sources per arcmin^2',s_density

! Generate an actual number of sources in the field - assume that Poisson 
! distribution can be approximated by Gaussian
   nsig = gasdev(idum)
   actual_n = expected_n+nsig*sqrt(expected_n)

! Calculate probability distribution flux -> probability
   do i = 1, tab_res

! Following distribution of different flux levels gives approx uniform 
! spacing of probabilities for gamma = 2
      s0(i) = s1+(s2-s1)*((i-1.0)/(tab_res-1.0))**2

      pd_s_p(i) = (k_parm*field_size/(1-gamma)*&
                         (s0(i)**(1-gamma)-s1**(1-gamma)))/expected_n
   end do

! Make probability distribution probability -> flux by calling locate
   do i = 1, tab_res
     p0(i) = 1.0*(i-1)/(tab_res-1.0)
     x1 = locate(pd_s_p,p0(i))
     pd_p_s(i) = s0(x1)+&
                 ((p0(i)-pd_s_p(x1))/(pd_s_p(x1+1)-pd_s_p(x1)))*&
                 (s0(x1+1)-s0(x1))
   end do     

   ra2 = obs_raref*hr2rad
   call io_getra('RA','*',ra2,status)
   obs_raref = ra2/hr2rad
   dec2 = obs_decref*deg2rad
   call io_getdec('Declination','*',dec2,status)
   obs_decref = dec2/deg2rad

! Loop over all sources
   if (verbose) write(*,*) 'Placing ',actual_n, ' sources'
   do i = 1, actual_n

! Generate random flux from distribution
      ran_index = 1.0+ran2(idum)*(tab_res-1.0)
      call interp(pd_p_s,tab_res,ran_index,s_source)

! Bin fluxes
      do j = 1, n_levs-1
         if ((s_source.ge.lim(j)).and.(s_source.lt.lim(j+1))) then
            n_lim(j)= n_lim(j)+1
         end if
      end do

! Generate random position in field
      src_x_pos = 1+maxsize*ran2(idum)   
      src_y_pos = 1+maxsize*ran2(idum)

! Convert from flux density to brightness temperature
      t_source = s_source/((cellsize*sec2rad)**2*2.0*k_b*obsfreq**2&
               /const_c2*1.0D26)

! Generate spectral index
      alpha = alpha_mean+alpha_rms*gasdev(idum)

! Add source into sky
      do chan = 1, nchan
         szsky(src_x_pos,src_y_pos,chan) = &
                   szsky(src_x_pos,src_y_pos,chan)+&
                   t_source*(nu(chan)/obsfreq)**(-2.-alpha) 
      end do

      if (report) then
         write(*,*) 'Source of flux ',s_source,' temperature ',&
             t_source, ' at ', src_x_pos, src_y_pos,' spectrum ',alpha
      end if
      if (to_file) then
         if (src_file) then
! Calculate source ra and dec in degrees
            dec = obs_decref+(src_y_pos-maxsize/2)*cellsize/3600.0
            ra = obs_raref*15.-(src_x_pos-maxsize/2)*cellsize/&
                               (3600.0*cos(dec*deg2rad))
! Convert ra to hrs
            ra = ra/15.
            call chr_chdtos(ra,prra,chra,lr)
            call chr_chdtos(dec,prdec,chdec,ld)
            write(iunit,*) i, chra(1:lr),' ', chdec(1:ld),s_source,alpha
         else
            write(iunit,*) 'Source of flux ',s_source,' temperature ',&
                  t_source, ' at ', src_x_pos, src_y_pos,' spectrum ',alpha
         end if
      endif
   end do

   do i = 1, n_levs-1
      if (verbose) write(*,*) 'Number of sources in range ', &
                 lim(i),' - ', lim(i+1),' is ', n_lim(i)
   end do

! Deallocate probability tables
  deallocate (pd_s_p)
  deallocate (pd_p_s)
  deallocate (s0)
  deallocate (p0)

! Deallocate binning arrays
  deallocate (n_lim)
  deallocate (lim)

end subroutine scatter_random_sources



