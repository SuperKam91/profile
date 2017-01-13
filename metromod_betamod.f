c-----*-----------------------------------------------------------------

      subroutine metro_betamod(mapin,okin)
      
      implicit none

      integer  m_ndim,m_nchains,m_neng,i,j,k
      parameter (m_ndim=5,m_nchains=10,m_neng=6)
      integer  m_iseed,m_nsteppl,m_verbose
      double precision m_lstart,m_lend,m_rate,m_evid
      double precision m_wmin,m_wmax,m_wfact
      double precision m_width(m_ndim),m_peng(m_neng)
      double precision m_xt(m_ndim,m_nchains)
      logical  m_varywidth
      
      integer  ic,id,istage,idum,dum1
      real     xmin,xrange,ymin,yrange,zmin,zrange,wmin,wrange
      real     rmin,rrange,xin,yin,zin,dum2,dum3
      
      integer pgopen, nchan
      external pgopen,poidev_ams

      double precision ran,poidev_ams
      double precision spec,gauss,ptsrc
      real xbr,ybr,zbr,wbr,rbr,xmax,ymax,zmax,wmax,rmax
      external ran,spec,gauss,ptsrc

      integer mapsize,status,maxsize,mchan,nc
      real cellsize,scale,incell,maxindex,minindex,maxsky,minsky
      parameter (mapsize = 512)
      parameter (cellsize = 15.)
      real index_map(mapsize,mapsize),inclev,c(20)
      real map_rectangle(4), tr(6), realdata(6)
      real counts(100),counts2(100),npix(100),r(100),r2,rtemp,ctsmax
      real cerr(100)
      double precision szsky(512,512),mapin(512,512),B
      integer nvis
      logical stat_fl, graph_adj, okpixel,ok_data(512,512),okin(512,512)

      character*80 fits_file,chr1,infile

      external B

      common /stage/ istage
      common /prior/ xmin,xrange,ymin,yrange,zmin,zrange,wmin,wrange,
     &     rmin,rrange
      common /outbar/ xbr,ybr,zbr,wbr,rbr,xmax,ymax,zmax,wmax,rmax
      common /data/ szsky,ok_data
      common /data2/ counts,cerr,r
      common /dimen/ mchan

      
      
c read in data:
      do i = 1,512
         do j = 1,512
            szsky(i,j)=mapin(i,j)
            ok_data(i,j) = okin(i,j)
         end do
      end do
c plot input data
      rtemp = 0.
      ctsmax = 0.
      do k = 1,100
         r(k) = real(k)*10.
         do i = 1,512
            do j = 1,512
               r2 = 16.*real(i-256)**2+real(j-257)**2
               if ((r2.le.r(k)**2).and.(r2.gt.rtemp**2)) then
                  counts(k) = counts(k)+real(szsky(i,j))
                  counts2(k) = counts2(k)+real(szsky(i,j))**2
                  npix(k) = npix(k)+1.
               end if
            end do
         end do
         rtemp = r(k)
         counts(k) = counts(k)/npix(k)
         counts2(k) = counts2(k)/npix(k)
         cerr(k) = (counts2(k) - counts(k)**2)**0.5
         if (counts(k).gt.ctsmax) ctsmax = counts(k)
         
      end do

      call pgenv(0.,1000.,0.,ctsmax,0,1)
      call pgpt(100,r,counts,2)
      
  
c     set likelihood
      
      idum = 9854
      
      
c set prior

      write(*,'(a,$)') ' Input [min, range] for beta: '
      read(*,*) xmin,xrange
      write(*,'(a,$)') ' Input [min, range] for theta1: '
      read(*,*) ymin,yrange
      write(*,'(a,$)') ' Input [min, range] for theta12: '
      read(*,*) zmin,zrange
      write(*,'(a,$)') ' Input [min, range] for PA: '
      read(*,*) rmin,rrange
      write(*,'(a,$)') ' Input [min, range] for norm: '
      read(*,*) wmin,wrange
      
      call pgenv(xmin,xmin+xrange,ymin,ymin+yrange,0,1)
    
      
c input general metro parameters

      write(*,'(a,$)') ' Input a random number seed integer: '
      read(*,*) m_iseed
      m_iseed=-iabs(m_iseed)
      
      write(*,*) 'You should go and make a cup of tea...'
      write(*,*) '               ...this may take some time.'

      do id=1,m_ndim
         m_width(id)=0.25d0
      end do
      m_wmin=0d0
      m_wmax=1d0
      m_wfact=1.414d0
      m_verbose=0
      m_peng(1)=1d0
      m_peng(2)=1d0
      m_peng(3)=0d0
      m_peng(4)=0d0
      m_peng(5)=0d0
      m_peng(6)=0d0
      
      
c     burn-in
      
      istage=1
      m_varywidth=.true.
      m_rate=0.1d0
      m_lstart=1d-6
      m_lend=1d0
      m_nsteppl=10
      do ic=1,m_nchains
         do id=1,m_ndim
            m_xt(id,ic)=ran(m_iseed)
         end do
      end do
      
      call metro(m_iseed,m_ndim,m_nchains,m_lstart,m_lend,m_rate,
     &     m_nsteppl,m_neng,m_peng,m_width,m_wmin,m_wmax,m_wfact,
     &     m_varywidth,m_evid,m_verbose,m_xt)
      
      write(*,*) 'log(evidence) = ',m_evid
      
c     sampling
      
      istage=2
      m_varywidth=.false.
      m_rate=0d0
      m_lstart=1d0
      m_lend=1d0
      m_nsteppl=10
      
      call metro(m_iseed,m_ndim,m_nchains,m_lstart,m_lend,m_rate,
     &     m_nsteppl,m_neng,m_peng,m_width,m_wmin,m_wmax,m_wfact,
     &     m_varywidth,m_evid,m_verbose,m_xt)
      
      call pgenv(0.,1000.,-0.,ctsmax,0,1)

      do i = 1,100
         npix(i) = wmax*(1+r(i)/ymax**2)**(-3.*xmax+0.5)*
     &        B(3.d0*xmax-0.5d0,0.5d0)
      end do
      call pgpt(100,r,counts,2)
      call pgsci(2)
      call pgline(100,r,npix)
      
      

      call pgend
      
      write(*,*) '-------------------------------------'
      write(*,*) '                DONE                   '
      write(*,*) '-------------------------------------'
      

      end

c-----*-----------------------------------------------------------------

      double precision function loglike(nd,cube)
      implicit none

      integer  nd
      double precision cube(nd)
      integer i,j,mchan,mapsize
      double precision x(nd)
      double precision y,pi,freq(100000,2)
      double precision cov(2,2),invcov(2,2),detcov,chisq
      real xmin,xrange,ymin,yrange,zmin,zrange,wmin,wrange
      real rmin,rrange,r2,xx,yy
      double precision szsky(512,512)
      real r(100),counts(100),cerr(100)
      integer nchan, nbin
      double precision l_min,l_max,spec,ptsrc,gauss,B
      logical ok_data(512,512)
      external spec,ptsrc,gauss,B
      
      common /prior/   xmin,xrange,ymin,yrange,zmin,zrange,wmin,wrange,
     &     rmin,rrange
      common /data/ szsky,ok_data
      common /data2/ counts,cerr,r
      common /dimen/ mchan

      pi = 3.14159d0

      x(1)=dble(xmin)+dble(xrange)*cube(1) !beta
      x(2)=dble(ymin)+dble(yrange)*cube(2) !theta1
      x(3)=dble(zmin)+dble(zrange)*cube(3) !theta2
      x(4)=dble(wmin)+dble(wrange)*cube(4) !norm
      x(5)=dble(rmin)+dble(rrange)*cube(5) !PA
      chisq = 0d0
      
c$$$      do i = 1,512
c$$$         do j = 1,512
c$$$            if(ok_data(i,j)) then
c$$$               xx = real(i-256)*4.
c$$$               yy = real(j-257)*4.
c$$$               r2 = (xx*cos(x(5)*pi/180.)+yy*sin(x(5)*pi/180.))**2
c$$$               y = (1+r2/x(2)**2)**(-3.*x(1)+0.5)*B(3.*x(1)-0.5,0.5)
c$$$               r2 = (yy*cos(x(5)*pi/180.)-xx*sin(x(5)*pi/180.))**2
c$$$               y = y*(1+r2/x(3)**2)**(-3.*x(1)+0.5)*B(3.*x(1)-0.5,0.5)
c$$$               y = y*x(4)
c$$$               chisq = chisq - szsky(i,j)*log(y) - y
c$$$            end if
c$$$         end do
c$$$      end do

      do i = 1,100
         y = x(4)*(1+r(i)/x(2)**2)**(-3.*x(1)+0.5)*
     :       B(3.d0*x(1)-0.5d0,0.5d0)
         chisq = chisq - (counts(i)-y)**2/cerr(i)**2
      end do

      loglike = chisq
      

      return
      end

c-----*-----------------------------------------------------------------

      subroutine metrobuild(cool,ndim,cube,lhood,lprior,retcode)
      implicit none

      integer ndim,retcode
      double precision cool,lhood,lprior,cube(*)

      double precision loglike
      external loglike

      lhood=loglike(ndim,cube)
      lprior=0d0
      retcode=0
      
      return
      end

c-----*-----------------------------------------------------------------

      subroutine metromonitor(iburn,cool,nchains,ic,nsteps,ndim,cube,
     &                        lhood,lprior,retcode)
      implicit none

      integer iburn,nchains,ic,nsteps,ndim,retcode,mchan
      double precision cool,lhood,lprior,cube(*),ltemp,llow

      real    xmin,xrange,ymin,yrange,zmin,zrange,wmin,wrange
      real    rmin,rrange
      real    xmax,ymax,zmax,wmax,rmax
      
      integer ns_afterburn
      save    ns_afterburn

      real    xbar,ybar,xsig,ysig,zbar,zsig,wbar,wsig,rbar,rsig

      
      double precision l_bin(100,100000),w_bin(100,100000)
      logical ok_bin(100,100000)
c      double precision freq(100),ell(10000),f1

      real    x,y,z,w,r
      integer istage

      common /stage/ istage
      common /prior/ xmin,xrange,ymin,yrange,zmin,zrange,wmin,wrange,
     &     rmin,rrange
      common /outbar/ xbar,ybar,zbar,wbar,rbar,xmax,ymax,zmax,wmax,rmax
      common /dimen/ mchan
      
      x=xmin+xrange*real(cube(1))
      y=ymin+yrange*real(cube(2))
      z=zmin+zrange*real(cube(3))
      w=wmin+wrange*real(cube(4))
      r=rmin+rrange*real(cube(5))
      
      if (istage.eq.1) then
        call pgsci(2)
        ns_afterburn=0
        ltemp = -10000000.
      else if (istage.eq.2) then
        call pgsci(1)
        ns_afterburn=ns_afterburn+1
        xbar=xbar+x
        xsig=xsig+x**2
        ybar=ybar+y
        ysig=ysig+y**2
        zbar=zbar+z
        zsig=zsig+z**2
        wbar=wbar+w
        wsig=wsig+w**2
        rbar=rbar+r
        rsig=rsig+r**2
        
      else if (istage.eq.3) then
        call pgsci(3)
      end if
      
      if (lhood.gt.ltemp) then
         ltemp = lhood
         xmax = x
         ymax = y
         wmax = w
      end if
      
      call pgpt1(x,y,1)

      retcode=0
      if ((istage.eq.1).and.(iburn.eq.0)) retcode=1
      if ((istage.eq.2).and.(ns_afterburn.gt.5000)) then
        retcode=1
c        write(*,*)
c        write(*,*) 'Total number of steps per chain = ',nsteps
        xbar=xbar/real(ns_afterburn)
        xsig=xsig/real(ns_afterburn)
        xsig=xsig-xbar**2
        xsig=sqrt(abs(xsig))
        ybar=ybar/real(ns_afterburn)
        ysig=ysig/real(ns_afterburn)
        ysig=ysig-ybar**2
        ysig=sqrt(abs(ysig))
        zbar=zbar/real(ns_afterburn)
        zsig=zsig/real(ns_afterburn)
        zsig=zsig-zbar**2
        zsig=sqrt(abs(zsig))
        wbar=wbar/real(ns_afterburn)
        wsig=wsig/real(ns_afterburn)
        wsig=wsig-wbar**2
        wsig=sqrt(abs(wsig))
        rbar=rbar/real(ns_afterburn)
        rsig=rsig/real(ns_afterburn)
        rsig=rsig-rbar**2
        rsig=sqrt(abs(rsig))
        


        write(*,*) 'xbar+/-xsig = ',xbar,'+/-',xsig
        write(*,*) 'ybar+/-ysig = ',ybar,'+/-',ysig
        write(*,*) 'zbar+/-zsig = ',zbar,'+/-',zsig
        write(*,*) 'wbar+/-wsig = ',wbar,'+/-',wsig
        write(*,*) 'rbar+/-rsig = ',rbar,'+/-',rsig
        write(*,*) 'minimum at:', xmax,ymax,wmax
        write(*,*) 'chisq:',ltemp
        write(*,*) 'reduced chisq:',ltemp/100.
      end if
      if ((istage.eq.3).and.(iburn.eq.0)) retcode=1
      

      return
      end

c-----*-----------------------------------------------------------------

! POISSON RANDOMNESS ADD-ER
      
      double precision function poidev_ams(xm, idum)
         
      integer idum
!     double precision poidev_ams
      double precision xm
      parameter ( pi = 3.141592 )
      real alxm, em, g, oldm, sq, t, y, ran1
      double precision gammln
      data oldm /-1./
         
      if (xm .lt. 12.) then
         if (xm .ne. oldm) then
            oldm = xm
            g = exp(-xm)
         end if
         
         em = -1
         t = 1.
 2       em = em + 1
         t = t*ran1(idum)
         
         if (t .gt. g) goto 2
      else
         if (xm .ne. oldm) then
            oldm = xm
            sq = sqrt(2.*xm)
            alxm = log(xm)
            g = xm*alxm-gammln(xm+1.d0)
         end if
         
 1       y = tan(pi*ran1(idum))
         em = sq*y + xm
         
         if (em .lt. 0.) goto 1
         
         em = int(em)
         t = 0.9*(1.+ y**2)*exp(em*alxm - gammln(em + 1.d0) - g)
         
         if (ran1(idum) .gt. t) goto 1
      end if
      
      poidev_ams = em
      
      return
      end
     
 ! RANDOM NUMBER GENERATOR
      
      real function ran1(idum)

      integer             idum
!     real    ran1
      parameter ( IA = 16807 )
      parameter ( IM = 2147483647)
      parameter    (  AM = 1./IM)
      parameter ( IQ = 127773)
      parameter ( IR = 2836)
      parameter (NTAB = 32)
      parameter  (NDIV = 1+(IM-1)/NTAB)
      parameter    ( EPS = 1.2e-7)
      parameter    ( RNMX = 1.-EPS)
      integer             j, k, iv(NTAB), iy
      
      save iv, iy
      data iv /NTAB*0/, iy /0/
      
      if (idum .le. 0 .or. iy .eq. 0) then
         idum = max(-idum, 1)
         do j = NTAB+8,1,-1
            k = idum/IQ
            idum = IA*(idum - k*IQ) - IR*k
            if (idum .lt. 0) idum = idum +IM
            if (j .le. NTAB) iv(j) = idum
         end do
         
         iy = iv(1)
         
      end if
      
      k = idum/IQ
      idum = IA*(idum-k*IQ) - IR*k
      
      if (idum .lt. 0) idum = idum+IM
      
      j = 1 + iy/ NDIV
      iy = iv(j)
      iv(j)  = idum
      ran1 = min(AM*iy, RNMX)
      
      return
      end
  
c-----*-----------------------------------------------------------
    
      double precision function gammln(xx)
         
!     double precision  gammln
      double precision  xx
      integer  j
      double precision ser, stp, tmp, x,y,cof(6)
      save  cof, stp
      data  cof, stp/76.18009172947146d0, -86.50532032941677d0, 
     &     24.01409824083091d0, -1.231739572450155d0,
     &     0.1208650973866179d-2, -0.5395239384953d-5,
     &     2.5066282746310005d0/
      
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser = 1.000000000190015d0
      do j = 1,6
         y = y+1.d0
         ser = ser +cof(j)/y
      end do
      gammln = tmp+log(stp*ser/x)
      
      return      
      end


c-----*----------------------------------------------------------------

c-----*----------------------------------------------------------------

      double precision function spec(freq,ay,bee,beta,alpha,cutoff)
      implicit none
      
      double precision freq,ay,bee,alpha,beta,cutoff
!     double precision spec

c      write(*,*) '***'
c      write(*,*) freq,ay,bee,beta,alpha,cutoff

      if (freq.lt.cutoff) then
         spec = ay*freq**(-1.*beta)
      else
         spec = bee*(freq - cutoff)**(-1.*alpha)
      end if

c      write(*,*) 'spec:',spec
c      read(*,*)
c      write(*,*) '***'
      

      return
      end

c----------------------------------------------------------------------

      double precision function gauss(freq,ell,ay,bee,sigma,alpha)
      implicit none

      double precision freq,ell,ay,bee,sigma,alpha,f1,norm
      double precision lambda,dee,pi
!     double precision gauss
      
      common /firstfreq/ f1


      pi = 3.14159
      lambda = 3.d8/(freq*1.d9)
      dee = 3.7/(2.*lambda)
c      dee = bee*dee

      norm = ay/(f1**(-1.*alpha))
      gauss = norm*exp(real(-1.*ell**2/(2.*(sigma**2))))
      gauss = (gauss)*freq**(-1.*alpha)
      gauss = gauss

      return
      end

c----------------------------------------------------------------------

      double precision function ptsrc(freq,ell,ay,bee)
      implicit none

      double precision freq,ell,ay,bee,f1
!     double precision ptsrc

      common /firstfreq/ f1

c      freq = freq*1d9
      ptsrc = ay*freq**(-1.*bee)
      ptsrc = ptsrc/(f1**(-1.*bee))

      return
      end

c----------------------------------------------------------------------

c BETA FUNCTION

      double precision function B(z,w)

       
!     double precision  B
      double precision  z,w
      double precision  gammln
      external gammln

      B = exp(gammln(z)+gammln(w)-gammln(z+w))
      return

      end 

c **********************************************************************


c XMM beam

c      subroutine xmm_beam


c      rc = 6.636
c      alpha = 1.525
c normalization of king profile:
c      sum = 0.
  
c      do i = 1,np
c         do j = 1,np
c            r = sqrt(real((np/2-(i-1))**2+(np/2-(j-1))**2))
c            r = r*cellsize      !asec
            
                                ! True beam:   
c            taper = ((1 + (r/rc)**2)**(alpha))**(-1)
                                ! Gaussian approx.:
c            taper = exp(-1.*(r**2)/(2*5.**2))
c            sum = sum + taper
c            beam(i,j)=taper
c         end do
c      end do
      
c      beam = beam/sum
      
c      do i=1,np
c         do j = 1,np
c            rdum=beam(i,j)
c            beam_map(i,j)=cmplx(rdum,0.)
c         end do
c      end do
      
c      job=1                     ! forward transform (+ve exponential)
c      iform=0                   ! data are real
c      call makefft(np,np,beam_map,job,iform,work)
      
      
      
      



