! Fits for n and H0 once a model aperture has been made.
subroutine fit_visibilities

   use kind_def
   use sz_globals
   use file_handling
   use kgraph

   implicit none

!   real(kind=dp), dimension(mvis) :: u, v, data_re, data_im, rms
   real(kind=dp), dimension(mvis) :: model_re, model_im, dummy1, dummy2, dummy3
   real(kind=dp), dimension(mvis) :: orig_model_re,orig_model_im,prev_model
   real(kind=dp), dimension(mvis) :: prev_model_re, prev_model_im
   real(kind=dp), dimension(:), allocatable :: bin_re, bin_im, bin_wt, bin_cent
   real(kind=dp), dimension(:), allocatable :: bin_mod_re, bin_mod_im
   real(kind=dp), dimension(:), allocatable :: bin_re2, bin_im2
   logical, dimension(:), allocatable :: use_bin
   real(kind=dp), dimension(mvis) :: n_array, theta, basel
   logical, dimension(mvis) :: use_vis, basel_plt
   real(kind=dp), dimension(:), allocatable :: chiarray, hub_c, likeli
   real(kind=dp) :: n, maxrad,minrad, binpos1, binpos2, binsize,fiddle,maxpost
   real(kind=dp) :: mean_re,mean_im,x1,x2,y1,y2,rescl,min_rms,lim,best_h0
   real(kind=dp) :: chisq,nstart,ninc,chisqmin,chisqmax,temp,sumweight,bestn
   integer :: i,dummy,iunit,nloop,j,k,nuse,nbin,nflag,plus,minus
   integer :: besti
   integer, dimension(1) :: posmin
   character*1 :: chr1, fit_imag

   dummy1 = 0.0
   dummy2 = 0.0
   dummy3 = 0.0
   u = 0.0
   v = 0.0
   rms = 0.0
   data_re = 0.0
   data_im = 0.0
   fit_imag = 'n'
   nflag = 0
   rescl = 1.0

   call read_vis_data
   write(*,*) nvis,' visibilities read in'
   min_rms = 10.0
   do i = 1,nvis
! Attempt to flag outlying points
!      if (abs(data_re(i)).gt.3.0*rms(i)) then
!         data_re(i) = 0.0
!         rms(i) = 0.0
!         nflag = nflag+1
!      end if
!      if (rms(i).gt.0.01) then
!         data_re(i) = 0.0
!         rms(i) = 0.0
!         nflag = nflag+1
!      end if 
!      if ((rms(i).ne.0.0).and.(rms(i).lt.min_rms)) min_rms=rms(i)
      if (rms(i).lt.min_rms) min_rms=rms(i)
   end do
   write(*,*) 'Number flagged is',nflag
   write(*,*) 'Minimum rms is',min_rms
   write(*,*) 'Maximum rms is',maxval(rms)
   write(*,*) 'Mean rms is',sum(rms)/(nvis-nflag)
   call io_getd('Factor for rescaling rms','*',rescl,status)
   rms = rms*rescl
!   rms = sum(rms)/nvis


! Need to flip over visibilities, reflecting in u axis of aperture plane
   u = -u
   v = v

! rms in visbin'ed files is on *amplitude*, so rms on real/imag is smaller
   rms = rms*0.7071
   where(rms.ne.0)
      dummy1 = data_re/rms**2
      dummy2 = data_im/rms**2
      dummy3 = 1./rms**2
   elsewhere 
      dummy1 = 0.0
      dummy2 = 0.0
      dummy3 = 0.0
   end where
   mean_re = sum(dummy1)
   mean_im = sum(dummy2)
   sumweight = sum(dummy3)

   mean_re = mean_re / sumweight
   mean_im = mean_im / sumweight
   sumweight = 1/sqrt(sumweight)

   write(*,*) 'mean real, imaginary',mean_re,mean_im,'+/-',sumweight

! Define which visibilities to use
   call io_getc('Use all visibilities (y/n):','y',chr1,status)
   if (chr1=='y') then
      use_vis(1:nvis) = .true.
      if (mvis.gt.nvis) then
         use_vis(nvis+1:mvis) = .false.
      end if
   else
      call io_getd('Minimum uv radius','*',minrad,status)
      call io_getd('Maximum uv radius','*',maxrad,status)
      use_vis = .false.
      where (((u**2+v**2).lt.maxrad**2).and.((u**2+v**2).gt.minrad**2))
         use_vis = .true.
      end where
   end if
   nuse = count(use_vis)

   call io_getc('Fit to imaginaries as well as reals (y/n):','n',fit_imag,&
                status)
 
! Get search parameters
   call io_geti('Number of searches over central density:','200',nloop,status)
   call io_getd('Start value of central density:','1.e-3',nstart,status)
   call io_getd('central density increment:','1.e-4',ninc,status)

! Allocate space to array for recording values of chi squared
   allocate(chiarray(nloop))
   allocate(hub_c(nloop))
   allocate(likeli(nloop))

! read in binned model visibilities
   write(*,*)'model zero-spacing flux = ',re(maxsize/2+1,maxsize/2+1,chan)

   do i = 1,nvis
      call extract_visibility(re,im,maxsize,u(i),v(i),cellsize,&
           & orig_model_re(i), orig_model_im(i))
   end do

! Changed 5/12/00 - SZ effect now decrement by default
! S-Z effect is a decrement in R-J region
   orig_model_re = -orig_model_re
   orig_model_im = -orig_model_im

! loop over various values of n for all other parameters remaining the same
   do j = 1,nloop
      n = nstart+(j-1)*ninc
      n_array(j) = n
      do k = 1,nvis
         model_re(k) = orig_model_re(k)*n/n0
         model_im(k) = orig_model_im(k)*n/n0
      end do

! calculated weighted chisquared
! For chisq need variance on *complex* number = 2*variance on sin or cos
      chisq = 0.
      do i = 1,nvis
         if (fit_imag=='y') then
            if ((rms(i).ne.0.0).and.(use_vis(i))) then
               temp = (((data_re(i)-model_re(i))**2)&
                     +(data_im(i)-model_im(i))**2)/(2*rms(i)**2)
               chisq = chisq + temp
            endif
         else
            if ((rms(i).ne.0.0).and.(use_vis(i))) then
               temp = (((data_re(i)-model_re(i))**2)/(rms(i)**2))
               chisq = chisq + temp
            endif
         end if
      end do
      chiarray(j) = chisq
   end do
   chisqmin = minval(chiarray)
   posmin = minloc(chiarray)
   bestn = nstart+(posmin(1)-1)*ninc

! Plot value of n against chi squared
   write(*,*) 'Assuming that H0 is ',H0 
   write(*,*) 'Plotting value of chi squared against n central'
   write(*,*) 
   write(*,*) 'Best value of n ',bestn
   if (plot_open) then
      call pgpage
!     call display_profile(chiarray,nloop,ninc)
!   There is a problem with display_profile (in kgraph). It doesn't allow for an inital point! D'oh!

      call wills_display_profile(chiarray,nstart,ninc,nloop)

      call pglabel("n central","Chi squared","")
   end if

! Is it a good fit? expect chisq/N = 1 +/- sqrt(2/N); the deviation from 1 is
! distributed as a gaussian, so use it to estimate w
   chisqmin = chisqmin/(2*nuse)

   write(*,*) 'Minimum (reduced) chisquared is ',chisqmin,&
        &    ' with ',2*nuse,' degrees of freedom'
   write(*,*) 'Expect chisquared = 1 +/-',sqrt(1/float(nuse))
   write(*,*) 'Chisq is ',abs(1-chisqmin)*sqrt(float(nuse)),&
        & ' sigma away from best possible fit'
   write(*,*) 
   write(*,*) 'When combined with x-ray data implies...'
   call io_getd('(with a fiddle factor of):','1.0', fiddle, status)
   best_h0 = H0*(n0/bestn)**2*fiddle
   write(*,*) 'Hubble constant = ',best_h0
   write(*,*) 'n0 = ',n0**2/bestn

!  Plot likelihoods
   if (plot_open) then

      do i = 1,nloop
         hub_c(i) = H0*(n0/(nstart+(i-1)*ninc))**2*fiddle
         likeli(i) = exp(-(chiarray(i)-minval(chiarray)))
      end do
      call pgpage
      x1 = 0.0
      x2 = 2*best_h0
      !   x2 = 150.0
      y1 = 0.0
      y2 = 1.0
      call display_profile(hub_c,likeli,nloop,x1,x2,y1,y2)
      call pglabel("Hubble Constant (km/s/Mpc)","Likelihood","")
      lim = 0.67
      call confidence(likeli, hub_c, nloop, plus, minus, lim)
      write(*,*) 'H0 error bars are +',hub_c(minus)-best_h0,hub_c(plus)-best_h0
      
      status = 0
      call io_getc('Multiply by Jefferies prior?:','y',chr1,status)
      if (chr1.eq.'y') then
         maxpost = 0.0
         do i = 1,nloop
            hub_c(i) = H0*(n0/(nstart+(i-1)*ninc))**2*fiddle
            likeli(i) = exp(-(chiarray(i)-minval(chiarray)))
            likeli(i) = likeli(i)/hub_c(i)
            if (likeli(i).gt.maxpost) then
               maxpost = likeli(i)
               besti = i
            end if
         end do
         do i = 1,nloop
            likeli(i) = likeli(i)/maxpost
         end do
         best_h0 = hub_c(besti)
         write(*,*) 'Best HO now = ',best_h0
         call pgpage
         x1 = 0.0
         x2 = 2*best_h0
         !   x2 = 150.0
         y1 = 0.0
         y2 = 1.0
         call display_profile(hub_c,likeli,nloop,x1,x2,y1,y2)
         call pglabel("Hubble Constant (km/s/Mpc)","Posterior Probability","")
         lim = 0.67
         call confidence(likeli, hub_c, nloop, plus, minus, lim)
         !      write(*,*) 'besti, plus, minus = ', besti, plus, minus
         write(*,*) 'H0 error bars are +',hub_c(minus)-best_h0,hub_c(plus)-best_h0
      end if


! Plot data and fit against either baseline or position angle 
      call io_getc('Plot against baseline or position angle or (b/a/r):','r',&
           chr1,status)
      if (chr1=='b') then
         where (use_vis)
            prev_model = orig_model_re*bestn/n0 
            basel = sqrt(u**2+v**2)
         elsewhere
            basel = 0.0
            data_re = 0.0
            data_im = 0.0
            prev_model = 0.0
            rms = 0.0
         end where
         call display_points(basel,prev_model,data_re,rms,use_vis,&
              & .true.,mvis,nuse,4,1)
          
         call pglabel("Baseline (lambda)","Re(data) Jy","")
                  
         where (use_vis)
            prev_model = orig_model_im*bestn/n0
         endwhere
         call display_points(basel,prev_model,data_im,rms,use_vis,&
              & .true.,mvis,nuse,4,1) 

           
         call pglabel("Baseline (lambda)","Im(data) Jy","")

      else if (chr1=='r') then

         nbin = 10
         binsize = 100.0
         nuse = 0
         where (use_vis)
            prev_model_re = orig_model_re*bestn/n0 
            prev_model_im = orig_model_im*bestn/n0 
            basel = sqrt(u**2+v**2)
         elsewhere
            basel = 0.0
            data_re = 0.0
            data_im = 0.0
            prev_model = 0.0
            rms = 0.0
         end where

         call io_geti('Number of bins?','*',nbin,status)
         call io_getd('Bin spacing?','*',binsize,status)
         allocate(bin_re(nbin))
         allocate(bin_im(nbin))
         allocate(bin_re2(nbin))
         allocate(bin_im2(nbin))
         allocate(bin_wt(nbin))
         allocate(bin_cent(nbin))
         allocate(use_bin(nbin))
         allocate(bin_mod_re(nbin))
         allocate(bin_mod_im(nbin))
         do i = 1,nbin
            bin_re(i) = 0.0
            bin_im(i) = 0.0
            bin_re2(i) = 0.0
            bin_im2(i) = 0.0
            bin_wt(i) = 0.0
            bin_cent(i) = 0.0
            use_bin(i) = .false.
            bin_mod_re(i) = 0.0
            bin_mod_im(i) = 0.0
         end do

         write(*,*) ' bin            Real                     Imaginary'
         do i = 1,nbin
            binpos1 = (i-1)*binsize
            binpos2 = i*binsize
            do j = 1, mvis
               if (use_vis(j)) then
                  if ((basel(j).lt.binpos2).and.(basel(j).gt.binpos1)) then
                     bin_re(i) = bin_re(i)+data_re(j)/rms(j)**2
                     bin_im(i) = bin_im(i)+data_im(j)/rms(j)**2
                     bin_re2(i) = bin_re2(i)+(data_re(j)**2)/rms(j)**2
                     bin_im2(i) = bin_im2(i)+(data_im(j)**2)/rms(j)**2
                     bin_wt(i) = bin_wt(i)+1.0/rms(j)**2
                     bin_cent(i) = bin_cent(i)+basel(j)/rms(j)**2
                     bin_mod_re(i) = bin_mod_re(i)+prev_model_re(j)/rms(j)**2
                     bin_mod_im(i) = bin_mod_im(i)+prev_model_im(j)/rms(j)**2
                  end if
               end if
            end do
            if (bin_wt(i).gt.0.0) then
               use_bin(i) = .true.
               bin_re(i) = bin_re(i)/bin_wt(i)
               bin_im(i) = bin_im(i)/bin_wt(i)
               bin_re2(i) = sqrt((bin_re2(i)/bin_wt(i)-bin_re(i)**2)/nvis)
               bin_im2(i) = sqrt((bin_im2(i)/bin_wt(i)-bin_im(i)**2)/nvis)
               bin_cent(i) = bin_cent(i)/bin_wt(i)
               bin_mod_re(i) = bin_mod_re(i)/bin_wt(i)
               bin_mod_im(i) = bin_mod_im(i)/bin_wt(i)
               bin_wt(i) = 1/sqrt(bin_wt(i))
               nuse = nuse+1
               write(*,100) &
               bin_cent(i),bin_re(i),'+/-',bin_wt(i),bin_im(i),'+/-',bin_wt(i)
100  format (ES10.3,4X,ES10.3,A3,ES10.3,4X,ES10.3,A3,ES10.3)
            else
               use_bin(i) = .false.
               bin_re(i) = 0.0
               bin_im(i) = 0.0
               bin_re2(i) = 0.0
               bin_im2(i) = 0.0
               bin_cent(i) = 0.0
               bin_mod_re(i) = 0.0
               bin_mod_im(i) = 0.0
               bin_wt(i) = 0.0
            end if
         end do
         call display_points(bin_cent,bin_mod_re,bin_re,bin_wt,use_bin,&
              & .true.,nbin,nuse,4,1)
          
         call pglabel("Baseline","Re(data) Jy","")

         call display_points(bin_cent,bin_mod_im,bin_im,bin_wt,use_bin,&
              & .true.,nbin,nuse,4,1)
          
         call pglabel("Baseline","Im(data) Jy","")


         deallocate(bin_re)
         deallocate(bin_im)
         deallocate(bin_re2)
         deallocate(bin_im2)
         deallocate(bin_wt)
         deallocate(use_bin)
         deallocate(bin_cent)
         deallocate(bin_mod_re)
         deallocate(bin_mod_im)

      else

         where (use_vis)
            theta = atan2(v,u)/deg2rad
            prev_model = orig_model_re*bestn/n0 
         elsewhere
            theta = 0.0
            data_re = 0.0
            data_im = 0.0
            prev_model = 0.0
            rms = 0.0
         end where


         call display_points(theta,prev_model,data_re,rms,use_vis,&
              & .true.,mvis,nuse,4,1)
          
         call pglabel("Baseline (lambda)","Re(data) Jy","")

         where (use_vis)
            prev_model = orig_model_im*bestn/n0
         endwhere

         call display_points(theta,prev_model,data_im,rms,use_vis,&
              & .true.,mvis,nuse,4,1) 

           
         call pglabel("pa /degrees","Im(data) Jy","")

      end if

   end if

   deallocate(chiarray)
   deallocate(hub_c)
   deallocate(likeli)

end subroutine fit_visibilities



subroutine wills_display_profile(profile,pstart,pinc,noofitems)

   use kind_def

   implicit none

   integer, intent(in) :: noofitems
   real(kind=dp), dimension(noofitems), intent(in) :: profile
   real(kind=dp), intent(in) :: pstart,pinc
   real, dimension(noofitems) :: x_points,y_points
   real :: minx,miny,maxx,maxy
   integer i,j,k     ! general counters

! Prepare x and y arrays
   do i=1,noofitems
      x_points(i) = pstart + (pinc * (i-1))
   end do
!      x_points = (/(i,i=1,noofitems)/) ! quick way to fill the array
!      x_points = (x_points-1)*x_step
   y_points = profile

! Find max and mins
   miny = minval(y_points)
   maxy = maxval(y_points)   
   minx = minval(x_points)
   maxx = maxval(x_points)

! Adjust y scale if necessary
   if ((miny>0.0).and.(maxy>2*miny)) then 
      miny = 0.0
   end if

   if ((maxy<0.0).and.(miny<2*maxy)) then 
      maxy = 0.0
   end if

   if (maxy==miny) then
      if (maxy==0.0) then
         maxy = 0.05
         miny = -0.05
      else if (maxy>0.0) then
         maxy = 1.05*maxy
         miny = 0.95*miny
      else
         maxy = 0.95*maxy
         miny = 1.05*miny
      end if
   end if

! Plot profile
   call pgvsize(1., 6.5, 1., 4.5)         ! in inches!
   call pgwindow(minx,maxx,miny,maxy)     ! defines a window in the viewport
   call pgbox('BCNITS', 0., 0, 'BCNITS', 0., 0)    ! draws a labelled frame
   call pgline(noofitems,x_points,y_points)           ! draws a line!!!!

end subroutine wills_display_profile

subroutine confidence(like, x, ndata, plus, minus, lim)

!  Subroutine to return % confidence intervals of a likelihood 
! array like(ndata) defined at points x(ndata) 
! lim is the confidence level required eg
! 0.67, indices of the intervals   returned in plus, minus

   use kind_def

   implicit none

   real (kind=dp) :: like(ndata), x(ndata), lim
   integer ndata, plus, minus
   
   real (kind=dp) :: conf, level, norm
   integer i
   
   character :: chr1
   integer :: status

   status = 0
   
   call io_getc('Attempt to find confidence limits?','n',chr1,status)
   if ((chr1=='y').or.(chr1=='Y')) then

      !     First find normalisation

      norm = 0
      do i = 1, ndata-1
         norm = norm + (((like(i) + like(i+1))/2) * (x(i+1)-x(i)))
      end do
      norm = 1 / norm
      
      level = 1
      conf = 0
      do while (conf.le.lim)
         level = level - 0.01
         conf = 0
         do i = 2, ndata-1
            if (like(i).gt.level) then 
               conf = conf + (like(i) * (x(i+1)-x(i-1))/2)
            end if
         end do
         conf = conf * norm
      end do
      write(*,*) 'level = ', level, 'conf = ',conf
      
      i = 0
      plus = 0.
      minus = 0.
      do while (minus.eq.0)
         i = i + 1
         if (like(i).gt.level) minus = i
      end do
      
      do while (plus.eq.0)
         i = i + 1
         if (like(i).lt.level) plus = i
      end do
   end if
   return
   
end subroutine confidence



