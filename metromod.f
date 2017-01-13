c-----*-----------------------------------------------------------------

      subroutine metro(m_iseed,m_ndim,m_nchains,m_lstart,m_lend,m_rate,
     &                 m_nsteppl,m_neng,m_peng,m_width,m_wmin,m_wmax,
     &                 m_wfact,m_varywidth,m_evid,m_verbose,m_xt)

      implicit none
      integer m_iseed,m_ndim,m_nchains,m_nsteppl,m_neng,m_verbose
      double precision m_lstart,m_lend,m_rate,m_evid
      double precision m_wmin,m_wmax,m_wfact
      double precision m_peng(m_neng),m_width(m_ndim)
      double precision m_xt(m_ndim,m_nchains)
      logical m_varywidth

      include 'metromod_params.inc'
      integer  netry(m_nengmax),neaccept(m_nengmax)
      integer  ntry(m_ndmax),naccept(m_ndmax)
      double precision xt(m_ndmax,m_ncmax),xtnew(m_ndmax,m_ncmax)
      double precision xmax(m_ndmax,m_ncmax),pmax(m_ncmax)
      double precision cube(m_ndmax)

      integer nsteps,nsamps,natoms,kpos,retcode
      integer ic,id,ie,nperl,iaccept,nlstep,iburn,ilambda
      double precision mu,mustart,mustep,lfactor,ratio
      double precision llhood,lprior,llmean,lpost,llmean0
      double precision lambda,lambdalo,lambdahi

      double precision ran,dmin
      external ran,dmin

      common /metro_state/ xt

c calculate derived control variables for geometric annealing schedule

c      if (m_rate.le.0d0) then   ! if m_rate =< 0 set lambda = 1
c        mustart=0d0
c        mustep=0d0
c      else
c        lfactor=1d0+m_rate
c        mustart=dlog(m_lstart)
c        mustep=dlog(lfactor)
c        mustart=-mustep*dint(-mustart/mustep) ! ensure lambdahi hits 1 exactly
c      end if
c      mu=mustart-mustep
      if (m_rate.eq.0d0) then
        mustart=0d0
        mustep=0d0
      else
        mustart=dlog(m_lstart)
        lfactor=1d0+m_rate
        nlstep=dint(dabs(dlog(m_lend)-dlog(m_lstart))/dlog(lfactor) 
     &              + 0.5d0)
        mustep=(dlog(m_lend)-dlog(m_lstart))/dble(nlstep)
      end if
      mu=mustart-mustep
      
c calculate likelihood and prior at input starting point for each chain
     
      do ic=1,m_nchains
        
        do id=1,m_ndim
          xt(id,ic)=m_xt(id,ic)
          cube(id)=xt(id,ic)
        end do
        
        call metrobuild(lambda,m_ndim,cube,llhood,lprior,retcode)
        pmax(ic)=llhood+lprior

      end do
     
c perform sampling

      iburn=1              ! start in burn-in state
      nsteps=0             ! number of steps taken by each chain
      nsamps=0             ! total number of samples
      natoms=1             ! metro does not support multiple atoms
      m_evid=0d0           ! evidence value accumulated
      ilambda=0            ! counter for lambda steps

c ... start loop over lambda 
 100  continue
      ilambda=ilambda+1

      call set_lambda(iburn,mu,mustep,lambda,lambdalo,lambdahi,m_lend,
     &                m_verbose)

      
      nperl=0              ! number of samples at present lambda
      llmean=0d0           ! mean loglikelihood at present lambda
      do id=1,m_ndim
        ntry(id)=0         ! engine 1 trials in id direction at lambda
        naccept(id)=0      ! engine 1 acceptances in id dirn at lambda
      end do
      do ie=1,m_neng
        netry(ie)=0        ! number of engine ie trails at lambda
        neaccept(ie)=0     ! number of engine ie acceptances at lambda
      end do
      
c ... start loop over steps
 200  continue
      nsteps=nsteps+1
 
      
c ... start loop over chains
      do ic=1,m_nchains
        nperl=nperl+1
        nsamps=nsamps+1
        call chain_step(ic,xt(1,ic),xtnew(1,ic),llhood,lprior,
     &                  lambda,m_ndim,m_nchains,m_iseed,ie,id,iaccept,
     &                  m_width,m_neng,m_peng)
        call metromonitor(iburn,lambda,m_nchains,ic,nsteps,m_ndim,
     &                    xtnew(1,ic),llhood,lprior,retcode)
        llmean=llmean+llhood
        lpost=llhood+lprior
  
        if (lpost.gt.pmax(ic)) then
          do id=1,m_ndmax
            xmax(id,ic)=xtnew(id,ic) ! record max point encountered 
          end do
          pmax(ic)=lpost
        end if

        if (ie.eq.1) then
          ntry(id)=ntry(id)+1
          naccept(id)=naccept(id)+iaccept
        end if
        netry(ie)=netry(ie)+1
        neaccept(ie)=neaccept(ie)+iaccept

        if (retcode.ne.0) goto 999 ! finished sampling
      end do  
c ... end of loop over chains 

      do ic=1,m_nchains
        do id=1,m_ndim
          xt(id,ic)=xtnew(id,ic)  ! update current position of chains
        end do
      end do
      if (nperl.lt.m_nsteppl*m_nchains) goto 200  ! check nperl value
c ... end of loop over steps

      llmean=llmean/dble(nperl)          

      if (iburn.eq.1) then
         m_evid=m_evid+llmean*dabs(lambdahi-lambdalo) ! update m_evid
         if ((ilambda.eq.1).and.(mustep.gt.0d0)) then
           m_evid=m_evid+llmean*lambdalo              ! correct end-point
         end if
         if ((ilambda.eq.nlstep).and.(mustep.lt.0d0)) then
           m_evid=m_evid+llmean*lambdahi              ! correct end-point
         end if
      end if

      if (m_varywidth) then ! vary proposal width
        do id=1,m_ndim
          if (ntry(id).gt.0) then
            ratio=dble(naccept(id))/dble(ntry(id))
            if (ratio.gt.0.3d0) m_width(id)=m_width(id)*m_wfact
            if (ratio.lt.0.3d0) m_width(id)=m_width(id)/m_wfact
            if (m_width(id).gt.m_wmax) m_width(id)=m_wmax
            if (m_width(id).lt.m_wmin) m_width(id)=m_wmin
          end if
        end do
      end if
 
      if (m_verbose.ge.1) then ! write out sampler diagnostics
        write(*,*) ' no. of sampler calls to each engine: ',
     &               (netry(ie),ie=1,m_neng)
        write(*,*) ' no. of acceptances  for each engine: ',
     &               (neaccept(ie),ie=1,m_neng)
        if (netry(1).gt.0) then
          do id=1,m_ndim
            write(*,*) ' engine 1:  calls, acceptances, width = ',
     &                   ntry(id),naccept(id),m_width(id)
          end do
        end if
      end if
      if (m_verbose.ge.2) then
        if (iburn.eq.1) then
          write(*,*) ' llmean, m_evid = ',llmean,m_evid
          write(*,*)
        end if
      end if

      goto 100
c ... end of loop over lambda

 999  continue

      if (iburn.eq.1) m_evid=m_evid


      do ic=1,m_nchains
        do id=1,m_ndim
          m_xt(id,ic)=xt(id,ic)  ! return final position of chains
        end do
      end do

      return
      end

c-----*-----------------------------------------------------------------

      subroutine chain_step(ic,cube,try,llhood,lprior,lambda,ndim,
     &                      nchains,idum,ie,id,iaccept,width,
     &                      ne,pe)
      implicit none

      include 'metromod_params.inc'
      double precision xt(m_ndmax,m_ncmax)

      integer  ic,ndim,nchains,idum,ie,id,iaccept,ne
      double precision lambda,llhood,lprior
      double precision cube(ndim),try(ndim),width(ndim),pe(ne)
      
      integer  n,ic2,ic3,retcode
      double precision u,alpha,llhoodtry,lpriortry,ptot

      integer irv
      external irv

      double precision proposal,proposal2,ran,dmin
      external proposal,proposal2,ran,dmin

      common /metro_state/ xt

c ... initialise iaccept
      iaccept=0

c ... scale to obtain probability for each engine
      ptot=0d0
      do ie=1,ne
        ptot=ptot+pe(ie)
      end do
      do ie=1,ne
        pe(ie)=pe(ie)/ptot
      end do

c ... choose MCMC engine at random
      ie=irv(idum,ne,pe)
      if ((nchains.lt.3).and.(ie.gt.2)) ie=1

c ... update single random parameter using (symmetric) proposal distribution
      if (ie.eq.1) then
        id=dint(ran(idum)*dble(ndim))+1  ! randomly pick parameter to update 
        do n=1,ndim
          try(n)=cube(n)
        end do
        try(id)=proposal(idum,id,ndim,cube,width)
        if (try(id).gt.1d0) try(id)=try(id)-dint(try(id))
        if (try(id).lt.0d0) try(id)=try(id)-dint(try(id)-1d0)
      end if

c ... or update all parameters using (symmetric) proposal distribution
      if (ie.eq.2) then
        do id=1,ndim
          try(id)=proposal2(idum,id,ndim,cube,width)
          if (try(id).gt.1d0) try(id)=try(id)-dint(try(id))
          if (try(id).lt.0d0) try(id)=try(id)-dint(try(id)-1d0)
        end do
      end if

c ... or update all parameters using the leapfrog method
      if (ie.eq.3) then
 100    ic2=dint(ran(idum)*dble(nchains))+1 ! randomly pick 'pivot' chain 
        if (ic2.eq.ic) goto 100
        do n=1,ndim
          try(n)=2d0*xt(n,ic2)-cube(n)
          if (try(n).gt.1d0) try(n)=try(n)-dint(try(n))
          if (try(n).lt.0d0) try(n)=try(n)-dint(try(n)-1d0)
        end do 
      end if

c ... or update all parameters using the cross walk method
      if (ie.eq.4) then
 200    ic2=dint(ran(idum)*dble(nchains))+1 ! randomly pick 'first' chain 
        if (ic2.eq.ic) goto 200
 300    ic3=dint(ran(idum)*dble(nchains))+1 ! randomly pick 'second' chain 
        if ((ic3.eq.ic2).or.(ic3.eq.ic)) goto 300
        do n=1,ndim
          try(n)=xt(n,ic2)+xt(n,ic3)-cube(n)
          if (try(n).gt.1d0) try(n)=try(n)-dint(try(n))
          if (try(n).lt.0d0) try(n)=try(n)-dint(try(n)-1d0)
        end do 
      end if

c ... or update all parameters using the guided walk method
      if (ie.eq.5) then
 400    ic2=dint(ran(idum)*dble(nchains))+1 ! randomly pick 'first' chain 
        if (ic2.eq.ic) goto 400
 500    ic3=dint(ran(idum)*dble(nchains))+1 ! randomly pick 'second' chain 
        if ((ic3.eq.ic2).or.(ic3.eq.ic)) goto 500
        do n=1,ndim
          try(n)=cube(n)+(xt(n,ic3)-xt(n,ic2))
          if (try(n).gt.1d0) try(n)=try(n)-dint(try(n))
          if (try(n).lt.0d0) try(n)=try(n)-dint(try(n)-1d0)
        end do 
      end if

c ... or update single random parameter by copy from a random chain (genetic)
      if (ie.eq.6) then
        ic2=dint(ran(idum)*dble(nchains))+1 ! randomly pick 'other' chain 
        id=dint(ran(idum)*dble(ndim))+1     ! randomly pick parameter to copy
        do n=1,ndim
          try(n)=cube(n)
        end do
        try(id)=xt(id,ic2)
      end if

c ... calculate ratio of posterior p(try)/p(cube) (Metropolis algorithm)
      call metrobuild(lambda,ndim,try,llhoodtry,lpriortry,retcode)
      call metrobuild(lambda,ndim,cube,llhood,lprior,retcode)
      alpha=lambda*(llhoodtry-llhood)+lpriortry-lprior

c ... decide whether to accept or reject candidate point
      u=dlog(ran(idum))
      if (u.le.alpha) then
        iaccept=1           ! if accepted, set iaccept=1   
        llhood=llhoodtry    ! if accepted, update llhood
        lprior=lpriortry    ! if accepted, update lprior
      else
        do n=1,ndim
          try(n)=cube(n)    ! if rejected, do not move chain
        end do
      end if

      return
      end

c-----*-----------------------------------------------------------------

      double precision function proposal(idum,n,ndim,cube,width)
      implicit none

      integer  idum,n,ndim,nc
      double precision cube(ndim),width(ndim)

      double precision ran,gdev
      external ran,gdev

c      proposal=ran(idum) ! uniform proposal distribution
      proposal=cube(n)+width(n)*gdev(idum) ! Gaussian prop. dist.

      return
      end

c-----*-----------------------------------------------------------------

      double precision function proposal2(idum,n,ndim,cube,width)
      implicit none

      integer  idum,n,ndim,nc
      double precision cube(ndim),width(ndim)

      double precision ran,gdev
      external ran,gdev

c      proposal=ran(idum) ! uniform proposal distribution
      proposal2=cube(n)+(width(n)/sqrt(real(ndim)))*gdev(idum) ! Gaussian 

      return
      end

c-----*-----------------------------------------------------------------

      subroutine set_lambda(iburn,mu,mustep,lambda,llo,lhi,lend,verbose)
      implicit none

      integer iburn,verbose
      double precision mu,mustep,lambda,llo,lhi,lend

      
      mu=mu+mustep
      llo=dexp(mu)
      lhi=dexp(mu+mustep)
      lambda=llo+0.5d0*(lhi-llo)
      if ((mustep.gt.0).and.(lambda.ge.lend)) then
        lambda=lend
        iburn=0
      end if
      if ((mustep.lt.0).and.(lambda.le.lend)) then
        lambda=lend
        iburn=0
      end if
c      if (verbose.ge.3) then
c         write(*,*) llo,lhi,lambda,iburn
c      end if

      return
      end

c-----*-----------------------------------------------------------------

      FUNCTION ran(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      DOUBLE PRECISION  ran,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1d0/IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-7,RNMX=1d0-EPS)
      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran=min(AM*iy,RNMX)
      return
      END

c-----*-----------------------------------------------------------------

      FUNCTION gdev(idum)
      INTEGER idum
      DOUBLE PRECISION gdev
CU    USES ran
      INTEGER iset
      DOUBLE PRECISION fac,gset,rsq,v1,v2,ran
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2d0*ran(idum)-1d0
        v2=2d0*ran(idum)-1d0
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=dsqrt(-2d0*dlog(rsq)/rsq)
        gset=v1*fac
        gdev=v2*fac
        iset=1
      else
        gdev=gset
        iset=0
      endif
      return
      END

c-----*-----------------------------------------------------------------

      double precision function dmin(x,y)
      implicit none
      double precision x,y

      if (x.le.y) then
        dmin=x
      else
        dmin=y
      end if

      return
      end

c-----*-----------------------------------------------------------------

c assigns a value between 1 and n to a discrete integer random variable 
c according to the probabilities p(n)

      integer function irv(idum,n,p)
      implicit none

      integer  idum,n
      double precision p(n)

      integer  i
      double precision dum,tmp,ran
      external ran

      dum=ran(idum)

      i=0
      tmp=0d0
 10   i=i+1
      tmp=tmp+p(i)
      if (dum.gt.tmp) goto 10
      irv=i

      return
      end

c-----*-----------------------------------------------------------------

