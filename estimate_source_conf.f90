! Makes looped calls to scatter_random_sources and make_aperture to estimate
! source confusion from a given power law for source population
subroutine estimate_source_confusion(idum)

   use kind_def
   use sz_globals
   use kgraph

   implicit none

   real(kind=dp) :: k_parm,gamma,s1,s2
   real(kind=dp),dimension(:),allocatable :: m_r, m_i, m_r2, m_i2
   real(kind=dp),dimension(:),allocatable :: rad
   integer :: uv_resol, iunit, idum
   integer, parameter :: npoints = 50
   real (kind=dp), dimension(npoints) :: re_pt, im_pt
   integer :: n_it, i, j, k
   real(kind=dp) :: max_rad, rad_inc, ang_inc, psi, upt,vpt
   logical :: to_file

   k_parm = 44
   gamma = 2.06
   s1 = 1000
   s2 = 5000.0
   verbose = .false.
   n_it = 5
   uv_resol = 60
   max_rad = 440.0

   call io_getd('Value of K',&
        & '*',k_parm,status)
   call io_getd('Value of gamma',&
        & '*',gamma,status)
   call io_getd('Minimum flux, S1',&
        & '*',s1,status)
   call io_getd('Maximum flux, S2',&
        & '*',s2,status)

   call io_geti('Number of realisations','*',n_it,status)
   call io_geti('Number of uv bins','*',uv_resol,status)
   call io_getd('Maximum uv radius','*',max_rad,status)

   allocate (m_r(uv_resol))
   allocate (m_i(uv_resol))
   allocate (m_r2(uv_resol))
   allocate (m_i2(uv_resol))
   allocate (rad(uv_resol))

   rad_inc = max_rad/uv_resol
   do i = 1,uv_resol
      rad(i) = i*rad_inc
   end do

!ANNA HACK STARTS

m_r = 0
m_i = 0
m_r2 = 0
m_i2 = 0

!ANNA HACK ENDS

   do i = 1, n_it
      write(*,*) 'Doing realisation ',i
      szsky = 0.0
      iunit = 60
      to_file = .false.
      call scatter_random_sources(k_parm,gamma,s1,s2,to_file,iunit,idum)
      call make_aperture
      ang_inc = pi/npoints
      do j = 1, uv_resol
         psi = 0.0D0
         do k = 1, npoints
            upt = rad(j)*sin(psi)
            vpt = rad(j)*cos(psi)
!write(*,*) 'u and v', u,v 
!read(*,*)
           call extract_visibility&
                 & (re,im,maxsize,upt,vpt,cellsize,re_pt(k),im_pt(k))
!write(*,*) 'real', re_pt(k)
!write(*,*) 'imag', im_pt(k)
!read(*,*)
            psi = psi+ang_inc
         end do
!write(*,*) 'average of real', sum(re_pt)/npoints
!read(*,*)
  
         m_r(j) = m_r(j)+sum(re_pt)/npoints
         m_i(j) = m_i(j)+sum(im_pt)/npoints
         m_r2(j) = m_r2(j)+(sum(re_pt**2)/npoints)
         m_i2(j) = m_i2(j)+(sum(im_pt**2)/npoints)

!write(*,*) 'variance over npoints is', m_r2(j)
!read(*,*) 

      end do
!write(*,*) 'realisation', i
write(*,*) 'cumulative average over npoints', m_i(1)
read(*,*)
   end do

   do j = 1, uv_resol
      m_r(j) = m_r(j)/n_it
      m_i(j) = m_i(j)/n_it
      m_r2(j) = ((m_r2(j)/n_it)-m_r(j)**2)
!write(*,*) 'average is', m_r(j)
!write(*,*) 'variance is', m_r2(j)
      if(m_r2(j).gt.0.0) then
         m_r2(j) = sqrt(m_r2(j))
      else
         write(*,*) 'sqrt of neg'
      end if
      m_i2(j) = ((m_i2(j)/n_it)-m_i(j)**2)
      if(m_i2(j).gt.0.0) then
         m_i2(j) = sqrt(m_i2(j))
      else
         write(*,*) 'sqrt of neg'
      end if
   end do
   
   if (autoscale) then
      write(*,*) 'Displaying mean real part of confusing sources'
      call pgpage
      call display_profile(rad,m_r,uv_resol,.true.)
      call pglabel('Observing baseline ( \gl )','Predicted source contribution (&
           &Jy)','Source contribution (Real)')
      call pgpage
      call display_profile(rad,m_r2,uv_resol,.true.)
      call pglabel('Observing baseline ( \gl )','Predicted source confusion (&
           &Jy)','Source Confusion (Real)')
      call pgpage
      call display_profile(rad,m_i,uv_resol,.true.)
      call pglabel('Observing baseline ( \gl )','Predicted source contribution (&
           &Jy)','Source contribution (Imag)')
      call pgpage
      call display_profile(rad,m_i2,uv_resol,.true.)
      call pglabel('Observing baseline ( \gl )','Predicted source confusion (&
           &Jy)','Source Confusion (Imag)')
   else
      call pgpage
      call display_profile(rad,m_r,uv_resol,prof_x1,prof_x2,prof_y1,prof_y2)
      call pglabel('Observing baseline ( \gl )','Predicted flux (&
           &Jy)','Mean Real')
   end if


   deallocate (m_r)
   deallocate (m_i)
   deallocate (m_r2)
   deallocate (m_i2)
   deallocate (rad)


end subroutine estimate_source_confusion








