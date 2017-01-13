*+
      subroutine FITS_WRHDR (obs_name,obs_epoch,raref,decref,
     : rarot,decrot,ngroup, freq, nchan, status)
c
c
c  AMI Project - Reduce program, write FITS header.
c
c  Constructs entries in the FITS header and writes the complete header
c  to disc.
c
c  History:
c    12/07/99 - original VSA version, adapted from UV2FITS program [DJT].
c     1/11/00 - IF coordinate suppressed [DJT].
c    13/12/02 - use phase centre as reference if phase rotated [DJT].
c    24/04/03 - format for GCOUNT relaxed to allow more visibilities [DJT].
c    17/01/05 - AMI version has frequency channel data [djt].
c    20/01/06 - introduced argument nsrc to identify source [djt].
c    15/11/06 - uscale now set to 1.0 [djt].
c     1/06/07 - introduced argument noff to identify offset [djt].
c    22/06/07 - write raster pointing as RA,Dec coordinates [djt].
c     7/08/07 - cos(dec) term removed from hex raster offsets [djt].
c    10/08/07 - write observation details as comment [djt].
c    07/09/07 - TELESCOP and INSTRUME now distinguish between SA and LA [jz]
c    10/09/07 - Extracted RA/Dec calculation to enq_off_radec [jz]
c    13/11/09 - Frequencies now provided in Hz rather than GHz [KG]
c
*-
      real*4  freq(*),raref,decref,rarot,decrot
      integer nsrc, noff, ngroup, nchan,obs_epoch
      character*32  obs_name
      integer status

      include 'uvfits.inc'
      include 'reduce.inc'
      include 'history.inc'
c     include 'header.inc'
      include 'status.inc'
      include 'telescope.inc'
      include 'constants.inc'

      character chstr*80
      character scratch*32
      character source*16
      character units*8
      real*8    epoch
      real*8    fd
      integer   l,lc,ls
      integer   i,iy,id,im,j
c      integer   ia
      integer   iscr
      integer  idate(3)
      real*8   mjd1

      integer   nline
      parameter (nline=80)
      character line*(nline)
      real*4    rline(nline/4)
      equivalence (line,rline)

      real*8    ra, dec
      real*8    haoff, decoff, appra, appdec
      real*8    djutc, eq, dr, dd, rm, dm
      real*8    freq1, dfreq

      real*8    one, zero
      data      one, zero / 1.d0, 0.d0 /

      character stokes(7)*4
      data      stokes / 'BEAM', 'I', 'Q', 'U', 'V', 'I-Q', 'I+Q' /

      integer   chr_lenb
      external  chr_lenb

      if (status.ne.0) return

c  Get pointing centre, and phase centre, at standard epoch

      call util_enqdat(idate)
      call sla_cldj(idate(3),idate(2),idate(1),mjd1,j) 
      obs_mjd_start = mjd1
      source = obs_name
      if (obs_epoch.eq.1) epoch=1950.0d0
      if (obs_epoch.eq.2) epoch=2000.0d0

c RA and dec in degrees:
      raref = raref*15.  ! hours->degrees
      obs_rarot = rarot*15. ! hours->degrees
      obs_decrot = decrot

c  Entry, open scratch file
      call io_opescr(iscr,scratch,'write',0,status)

c  Initial items

      bitpix = 32
      naxis = 6
      naxisn(1) = 0
      naxisn(2) = 3
      naxisn(3) = 1
      naxisn(4) = nchan
      naxisn(5) = 1
      naxisn(6) = 1

      write(iscr,10) -bitpix,naxis,(naxisn(i),i=1,naxis)
 10   format(
     :  'SIMPLE  = ',19x,'T',1x,'/ Standard FITS format'/
     :  'BITPIX  = ',16x,i4,1x,'/ Bits per pixel'/
     :  'NAXIS   = ',17x,i3,1x,'/ Number of axes'/
     :  'NAXIS1  = ',16x,i4,1x,'/ No image: uv data'/
     :  'NAXIS2  = ',16x,i4,1x,'/ Complex: cos, sin, weight'/
     :  'NAXIS3  = ',16x,i4,1x,'/ Number of polarisations'/
     :  'NAXIS4  = ',16x,i4,1x,'/ Number of frequencies'/
     :  'NAXIS5  = ',16x,i4,1x,'/ RA'/
     :  'NAXIS6  = ',16x,i4,1x,'/ Dec'/
     :  'BLOCKED = ',19x,'T',1x,'/ Tape may be blocked')

c  Random groups format

      pcount = 6
      gcount = ngroup

      write(iscr,11) pcount,gcount
 11   format(
     :  'GROUPS  = ',19x,'T',1x,'/ Groups data structure'/
     :  'PCOUNT  = ',16x,i4,1x,'/ Parameters per group'/
     :  'GCOUNT  = ',12x,i8,1x,'/ Number of groups'/
     :  'EXTEND  = ',19x,'T',1x,'/ Extension is antenna table')

c  Object, telescope, date of observation

      call util_enqdat(idate)

      write(iscr,12) source(1:16),'Profile','Sim',idate(3),idate(2),
     :               idate(1)

 12   format(
     :  'OBJECT  = ',1h',a,1h',3x,'/ Source/field name'/
     :  'TELESCOP= ',1h',a,5x,1h',8x,'/ Profile'/
     :  'INSTRUME= ',1h',a,5x,1h',8x,'/ Sim'/
     :  'DATE-OBS= ',1h',i4,2('-',i2.2),1h',9x,
     :                             '/ Start date of observation')

c  Bunit, bzero, bscale, blank, epoch

      units = ' '
      bzero = zero
      bscale = one
      write(iscr,13) units,bzero,bscale
 13   format(
     :  'BUNIT   = ',1h',a,1h',11x,'/ Units of data'/
     :  'BZERO   = ',4x,1pe16.9,1x,'/ Data offset'/
     :  'BSCALE  = ',4x,1pe16.9,1x,'/ Data = FITS*BSCALE + BZERO')

      blank = zero
      if(blank.ne.0.d0) write(iscr,14) blank
 14   format(
     :  'BLANK   = ',4x,1pe16.9,1x,'/ Undefined values on tape')

      write(iscr,15) epoch
 15   format(
     :  'EPOCH   = ',10x,f10.4,1x,'/ Equinox of RA, Dec')

c  Pointing centre [degrees]
      write(iscr,16) raref, decref  
 16   format(
     :  'OBSRA   = ',4x,1pe16.9,1x,'/ Antenna pointing RA'/
     :  'OBSDEC  = ',4x,1pe16.9,1x,'/ Antenna pointing DEC')

c  Phase centre [degrees]
      raref = obs_rarot
      decref = obs_decrot
c      if (obs_rarot.ne.0.0 .or. obs_decrot.ne.0.0) then
c         if (obs_rarot.ne.obs_raobs .or. obs_decrot.ne.obs_decobs) then
c            eq = 2000.d0
c            ra = obs_rarot
c            dec = obs_decrot
c            djutc = obs_mjd_start
c            call sla_amp(ra,dec,djutc,eq,rm,dm)
c            raref = rm
c            decref =dm
c            if (obs_epoch.eq.1) then
c               eq = 1950.d0
c               call sla_fk54z(rm,dm,eq,raref,decref,dr,dd)
c            endif
c         endif
c      endif

c  Frequency channels

      freq1 = freq(1)
      dfreq = zero
      if (nchan.gt.1) dfreq = freq(2)-freq(1)

c  Details of group data

      write(iscr,17) one,one,one,zero
 17   format(
     :  'CTYPE2  = ',10h'COMPLEX ',11x,'/'/
     :  'CRVAL2  = ',4x,1pe16.9,1x,'/'/
     :  'CDELT2  = ',4x,1pe16.9,1x,'/'/
     :  'CRPIX2  = ',4x,1pe16.9,1x,'/'/
     :  'CROTA2  = ',4x,1pe16.9,1x,'/')

      write(iscr,18) stokes(1),one,one,one,zero
 18   format(
     :  'CTYPE3  = ',10h'STOKES  ',11x,'/ Stokes parameter : ',A/
     :  'CRVAL3  = ',4x,1pe16.9,1x,'/'/
     :  'CDELT3  = ',4x,1pe16.9,1x,'/'/
     :  'CRPIX3  = ',4x,1pe16.9,1x,'/'/
     :  'CROTA3  = ',4x,1pe16.9,1x,'/')

      write(iscr,19) freq1,dfreq,one,zero
 19   format(
     :  'CTYPE4  = ',10h'FREQ    ',11x,'/ Frequency, Hertz'/
     :  'CRVAL4  = ',4x,1pe16.9,1x,'/'/
     :  'CDELT4  = ',4x,1pe16.9,1x,'/'/
     :  'CRPIX4  = ',4x,1pe16.9,1x,'/'/
     :  'CROTA4  = ',4x,1pe16.9,1x,'/')

      write(iscr,20) raref,one,one,zero
 20   format(
     :  'CTYPE5  = ',10h'RA      ',11x,'/ Right Ascension, degrees'/
     :  'CRVAL5  = ',4x,1pe16.9,1x,'/'/
     :  'CDELT5  = ',4x,1pe16.9,1x,'/'/
     :  'CRPIX5  = ',4x,1pe16.9,1x,'/'/
     :  'CROTA5  = ',4x,1pe16.9,1x,'/')

      write(iscr,21) decref,one,one,zero
 21   format(
     :  'CTYPE6  = ',10h'DEC     ',11x,'/ Declination, degrees'/
     :  'CRVAL6  = ',4x,1pe16.9,1x,'/'/
     :  'CDELT6  = ',4x,1pe16.9,1x,'/'/
     :  'CRPIX6  = ',4x,1pe16.9,1x,'/'/
     :  'CROTA6  = ',4x,1pe16.9,1x,'/')

c  Details of group parameters

c     uscale = 1.0d-10
      uscale = one
      uzero = zero
      write(iscr,22) uscale,uzero
 22   format(
     :  'PTYPE1  = ',10h'UU      ',11x,'/ U coordinate, seconds'/
     :  'PSCAL1  = ',4x,1pe16.9,1x,'/ Real = FITS*PSCAL + PZERO'/
     :  'PZERO1  = ',4x,1pe16.9,1x,'/')

      vscale = uscale
      vzero = zero
      write(iscr,23) vscale,vzero
 23   format(
     :  'PTYPE2  = ',10h'VV      ',11x,'/ V coordinate, seconds'/
     :  'PSCAL2  = ',4x,1pe16.9,1x,'/ Real = FITS*PSCAL + PZERO'/
     :  'PZERO2  = ',4x,1pe16.9,1x,'/')

      wscale = uscale
      wzero = zero
      write(iscr,24) wscale,wzero
 24   format(
     :  'PTYPE3  = ',10h'WW      ',11x,'/ W coordinate, seconds'/
     :  'PSCAL3  = ',4x,1pe16.9,1x,'/ Real = FITS*PSCAL + PZERO'/
     :  'PZERO3  = ',4x,1pe16.9,1x,'/')

c  Date parameter, as two entries:
c   ptype1: pzero1 = intpt(jd at start-time of observations)
c           pscal1 = 1/128 days
c           value of type1 parameter = jd of sample to nearest
c              integral value of 1/128 days
c   ptype2: pzero2 = 0
c           pscal2 = 1/6e6; precision is about 0.01 second

      jdzero = int(obs_mjd_start+2400000.5d0)
      jdscal1 = 1.d0/128.d0
      jdscal2 = 1.d0/6.d6
      write(iscr,25) jdscal1,jdzero,jdscal2,zero
 25   format(
     :  'PTYPE4  = ',10h'DATE    ',11x,'/ Julian Date' /
     :  'PSCAL4  = ',4x,1pe16.9,1x,'/'/
     :  'PZERO4  = ',4x,1pe16.9,1x,'/'/
     :  'PTYPE5  = ',10h'DATE    ',11x,'/ Julian Date'/
     :  'PSCAL5  = ',4x,1pe16.9,1x,'/'/
     :  'PZERO5  = ',4x,1pe16.9,1x,'/')

      write(iscr,26) one,zero
 26   format(
     :  'PTYPE6  = ',10h'BASELINE',11x,'/ Ant1*256 + Ant2'/
     :  'PSCAL6  = ',4x,1pe16.9,1x,'/'/
     :  'PZERO6  = ',4x,1pe16.9,1x,'/')

c  Observation details

      write(iscr,27) 'Profile Simulation'

      chstr = 'simulated observation'
      ls = chr_lenb(chstr)
      l = chr_lenb(obs_name)
      write(iscr,27) obs_name(1:l),chstr(1:ls)
c     lc = chr_lenb(obs_comment)
c     if (lc.gt.0) write(iscr,27) obs_comment(1:lc)

      write(iscr,27) chstr(1:chr_lenb(chstr))
      
 27   format('COMMENT',2(X,A))

c  Data reduction history

      call set_history(status)
      do i = 1,lhist
         write(iscr,28) history(i)(1:chr_lenb(history(i)))
 28      format('HISTORY ',A)
      enddo

c  History records for AIPS sort-order etc

      write(iscr,29) one
 29   format(
     :  'HISTORY AIPS  SORT ORDER = ''TB'''/
     :  '                   / Where T is time (IAT), B is baseline num'/
     :  'HISTORY AIPS WTSCAL = ',4x,1pd16.9/
     :  '                   / Complex wts = WTSCAL*(FITS*BSCALE+BZERO)')

c  End of generating header to scratch file

      write(iscr,30)
 30   format('END')
      endfile(iscr)

c  Rewind scratch file, copy records to FITS output file

      rewind(iscr)
 31   read(iscr,'(A)',end=32) line
      call fits_write(rline,nline,status)
ccc      if (status .eq. -1) then
ccc         write(*,*) '1',line
ccc         write(*,*) '2',rline
ccc         write(*,*) '3',nline
ccc      endif
ccc      write(*,*) 'Copy into fits', status
      goto 31

c  End of header, fill current block with blanks

 32   line=' '
      call fits_fill(rline,status)
c      write(*,*) 'Fill with blanks', status

      close(iscr)

      end
