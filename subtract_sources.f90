

subroutine subtract_sources(bmaj,bmin)

   use kind_def
   use sz_globals
   use kgraph
   use mikes_fft
   use map_handling

   implicit none

   integer       :: i,j,k,l1
   real(dp)      :: temp,sum,norm,synth_beam,peak,norm2
   REAL          :: XFIT(1000),YFIT(1000),XX(1000),YY(1000)
   REAL          :: SIG(1000),AFIT(3),bmaj,bmin
   INTEGER       :: IAFIT(3),MA,NPC
   REAL          :: COVAR(3,3),CHISQ
   character*1   :: chr1
   character*80  :: plotfile, plotname

   external chr_lenb
   external pgopen
   external TWIST
   integer  chr_lenb
   integer  pgopen


   call def_map_region

   k = 0
   do i = 1,maxsize
      do j = 1,maxsize
         if (ok_data(i,j)) then
            if ((.not.ok_data(i+1,j)).or.(.not.ok_data(i,j+1)).or.&
                 (.not.ok_data(i-1,j)).or.(.not.ok_data(i,j-1))) then
               k = k+1
               xx(k) = i
               yy(k) = j
               yfit(k) = szsky(i,j,chan)
               xfit(k) = k
               sig(k) = 1.
            end if
         end if
      end do
   end do

! Fit to three parameters (AFIT(1)+AFIT(2)*XX+AFIT(3)*YY)

   MA=3
   NPC=3
   DO I=1,3
      AFIT(I)=0.0
      IAFIT(I)=1
   ENDDO

   CALL LFIT(XFIT,YFIT,SIG,k,AFIT,IAFIT,MA,COVAR,NPC,CHISQ,TWIST,XX,YY)
   
   WRITE(*,*)
   WRITE(*,*) 'Fitted twisted plane: A + B x + C y'
   WRITE(*,*) 'A:',AFIT(1)
   WRITE(*,*) 'B:',AFIT(2)
   WRITE(*,*) 'C:',AFIT(3)
   WRITE(*,*)

   write(*,*) 'Replacing sky with twisted plane.'

   norm = 0.
   norm = cellsize*sec2rad
   norm = norm*norm
   norm = norm*2.0*k_b*obsfreq**2/const_c2*1.0D26
   
   synth_beam = pi*bmaj*bmin*(sec2rad**2)/(4.*log(2.)) !!CHECK
   norm2 = 1.0D26*synth_beam*k_b*2.0*obsfreq/const_c2
   
   temp = 0.
   sum = 0.
   peak = 0.
   do i = 1,maxsize
      do j = 1,maxsize
         if (ok_data(i,j)) then
            temp = szsky(i,j,chan)
            szsky(i,j,chan) = AFIT(1) + AFIT(2)*real(i) + AFIT(3)*real(j)
            sum = sum + (temp - szsky(i,j,chan))
            if ((temp - szsky(i,j,chan))*norm2.gt.peak) then
               peak = (temp - szsky(i,j,chan))*norm2
            end if
         end if
      end do
   end do

   write(*,*) 'Displaying subtracted sky:'
   call display_map(szsky(:,:,chan),maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
   
   write(*,*) '--------------------------------'
   write(*,*) 'Peak flux:',peak,'Jy/beam'
   write(*,*) 'Integrated flux:',sum*(cellsize*sec2rad)**2*2.0*k_b*obsfreq**2/const_c2*1.0D26,'Jy.'
   write(*,*) '--------------------------------'
   
end subroutine subtract_sources


! ---------------------------------------------------------------------------------

!!uv style subtraction:


!
!   call make_aperture
!
!   call io_getc('Source list:','L860_src.dat',srclist,status)
!   call io_nxtlun(lunit,status)
!   open (lunit,file=srclist)
!   i = 1
!777 read(lunit,*,end=401) xpix(i),ypix(i),axx(i),axy(i),flux(i)
!   i = i+1
!   maxi = i-1
!   goto 777
!
!401 continue
!   
!   write(*,*) 'Read ',maxi,' sources from file.'
!! stupid scaling stuff from r-s-f and m-a:
!   norm = pi*bmaj*bmin*(sec2rad**2)/(4.*log(2.)) 
!   norm = 1.0D26*norm*k_b*2.0/obslambda**2
!   flux = flux*1.0/norm
!   norm = cellsize*sec2rad
!   norm = norm*norm
!   norm = norm*2.0*k_b/(obslambda**2)*1.0D26
!   flux = flux*norm
!
!   allocate(sourcesky(maxsize,maxsize))
!   allocate(beamsky(maxsize,maxsize))
!   allocate(work(2*maxsize*maxsize))
!   do k = 1,maxi
!      sourcesky = 0.
!      beamsky = 0.
!      sourcesky(xpix(k),ypix(k)) = cmplx(flux(k),0.)
!   
!      do i = 1,maxsize
!         do j = 1,maxsize
!            radx2 = cellsize**2*real(i-(maxsize/2+1))**2
!            rady2 = cellsize**2*real(j-(maxsize/2+1))**2
!            beamsky(i,j) = cmplx(exp(-radx2/(2.*axx(k)**2/(2.*log(2.))))&
!                 *exp(-rady2/(2.*axy(k)**2/(2.*log(2.)))),0.)
!         end do
!      end do
!      
!! ... fft of sky
!      job=1                     ! forward transform (-ve exponential)
!      iform=0                   ! data are real
!      call makefft(maxsize,maxsize,sourcesky,job,iform,work)   
!! ... fft of beam
!      job=1                     ! forward transform (-ve exponential)
!      iform=0                   ! data are real
!      call makefft(maxsize,maxsize,beamsky,job,iform,work)   
!      
!! to get convolution in image plane:
!      sourcesky = sourcesky*beamsky
!! subtract sources   
!      re(:,:,chan) = re(:,:,chan) - real(sourcesky)
!      im(:,:,chan) = im(:,:,chan) - aimag(sourcesky)
!
!   end do
!
!   deallocate(beamsky)
!
!   sourcesky = cmplx(re(:,:,chan),im(:,:,chan))
!
!! ... inverse fft of product
!   job=-1                     ! inverse transform (+ve exponential)
!   iform=1                    ! data are complex
!   call makefft(maxsize,maxsize,sourcesky,job,iform,work)   
! 
!   szsky(:,:,chan) = real(sourcesky)
!   re(:,:,chan) = 0.
!   im(:,:,chan) = 0.
!
!   deallocate(work)
!   deallocate(sourcesky)
!
!   write(*,*) 'Displaying subtracted sky:'
!   call display_map(szsky(:,:,chan),maxsize,graph_adj,cellsize,cellsize,0D0,0D0)
!
!
!end subroutine subtract_sources
