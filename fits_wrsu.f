*+
      subroutine FITS_WRSU (iscr, obs_name,obs_epoch,raref,decref,
     &     noff,ra_oset,dec_oset,sing_phase_cent,status)
c
c  AMI Project - Reduce program, write source (SU) table to FITS file.
c
c  Construct extension table containing details of offset sources
c  and appends to FITS file.
c
c  History:
c    07/09/07 - original version based on FITS_WRANT and FITS_WRHDR [jz]
c    17/09/07 - add offset index to source name [djt]
c    17/09/07 - change offset index to spiral pattern for hexagons
c                     and rings and standardise across obs [nhw]
c     9/05/08 - bug fix for nxm hex raster offset indexing [djt]
c
*-
      implicit none
      real*4  raref,decref
      integer obs_epoch
      character*32  obs_name
      integer noff
      real*4  ra_oset(noff),dec_oset(noff)
      logical sing_phase_cent
      integer status

      include 'uvfits.inc'
      include 'reduce.inc'
      include 'history.inc'
c     include 'header.inc'
      include 'status.inc'
      include 'telescope.inc'
      include 'constants.inc'
      include 'raster.inc'

      character chstr(noff)*32, chstr_temp*32
      character source*16
      character scratch*32

      integer   iscr, iunit

      integer   nline
      parameter (nline=240)
      integer   maxsize
      parameter (maxsize=256)
      character line*(nline)
      real*4    rline(nline/4)
      equivalence (line,rline)

      integer   ioff,l,offunit
c     integer   ioff,noff,l,offunit
      integer   off_list(MAX_NSTEP)
      integer   pt_list(MAX_NSTEP)
      real*8    off_ra(MAX_NSTEP),off_dec(MAX_NSTEP)
      real*8    off_app_ra(MAX_NSTEP),off_app_dec(MAX_NSTEP)
      real*8    epoch
      character cioff*1
      character offile*80

      integer  chr_lenb
      external chr_lenb

c Bin table vars
      character collabel(19)*10
      character colformat(19)*10
      character colunit(19)*10
      integer tfields
      character extname*10
      integer felem

c keyword varables
      character comment*72
      character value*16
      integer hdutype

      integer naxes(12)



c  Entry
      if(status.ne.0) return
 
      if (obs_epoch.eq.1) epoch = 1950.
      if (obs_epoch.eq.2) epoch = 2000.

c  Set up pointing name lists
      do ioff = 1, noff
         pt_list(ioff)=ioff
      enddo

      source = obs_name
      l = chr_lenb(source)
      do ioff = 1,noff 
         chstr_temp = source(1:l)
         write(chstr_temp(l+1:),'(I3.3)') pt_list(ioff)
         chstr(ioff) = chstr_temp(1:l+3)
      end do

c  Get pointing centre, and phase centre, at standard epoch
      raref = raref*15.
      do ioff = 1,noff
         off_dec(ioff) = decref+dec_oset(ioff)/3600.
         off_ra(ioff) = raref-ra_oset(ioff)/
     :                  (cos(off_dec(ioff)*const_d2r)*3600.)
         if (sing_phase_cent) then
            off_app_ra(ioff) = raref
            off_app_dec(ioff) = decref
         else
            off_app_ra(ioff) = off_ra(ioff)
            off_app_dec(ioff) = off_dec(ioff)
         end if
      end do

c Open fits
c      iunit = iscr
      call ftgiou(iunit, status)
      call FTREOPEN(iscr, iunit, status)
      if (status.ne.0) write(*,*) 'FTREOPEN in su failed', status

c Move to last existing HDU
      call FTMAHD(iunit,2, hdutype,status)
      if (status.ne.0) write(*,*) 'Moving to HDU', hdutype, status

c Setup Axis labels/properties for AN table

      collabel(1)  = 'ID. NO. '
      collabel(2)  = 'SOURCE  '
      collabel(3)  = 'QUAL    '
      collabel(4)  = 'CALCODE '
      collabel(5)  = 'IFLUX   '
      collabel(6)  = 'QFLUX   '
      collabel(7)  = 'UFLUX   '
      collabel(8)  = 'VFLUX   '
      collabel(9)  = 'FREQOFF '
      collabel(10) = 'BANDWIDTH'
      collabel(11) = 'RAEPO   '
      collabel(12) = 'DECEPO  '
      collabel(13) = 'EPOCH   '
      collabel(14) = 'RAAPP   '
      collabel(15) = 'DECAPP  '
      collabel(16) = 'LSRVEL  '
      collabel(17) = 'RESTFREQ'
      collabel(18) = 'PMRA    '
      collabel(19) = 'PMDEC   '

      colformat(1)  = 'J       '
      colformat(2)  = '16A     '
      colformat(3)  = 'J       '
      colformat(4)  = '4A      '
      colformat(5)  = 'D       '
      colformat(6)  = 'D       '
      colformat(7)  = 'D       '
      colformat(8)  = 'D       '
      colformat(9)  = 'D       '
      colformat(10) = 'D       '
      colformat(11) = 'D       '
      colformat(12) = 'D       '
      colformat(13) = 'D       '
      colformat(14) = 'D       '
      colformat(15) = 'D       '
      colformat(16) = 'D       '
      colformat(17) = 'D       '
      colformat(18) = 'D       '
      colformat(19) = 'D       '

      colunit(1)  = '        '
      colunit(2)  = '        '
      colunit(3)  = '        '
      colunit(4)  = '        '
      colunit(5)  = 'JY      '
      colunit(6)  = 'JY      '
      colunit(7)  = 'JY      '
      colunit(8)  = 'JY      '
      colunit(9)  = 'HZ      '
      colunit(10) = 'HZ      '
      colunit(11) = 'DEGREES '
      colunit(12) = 'DEGREES '
      colunit(13) = 'YEARS   '
      colunit(14) = 'DEGREES '
      colunit(15) = 'DEGREES '
      colunit(16) = 'M/SEC   '
      colunit(17) = 'HZ      '
      colunit(18) = 'DEG/DAY '
      colunit(19) = 'DEG/DAY '
      
      tfields = 19

c Create new BIN table for Antenna positions

      extname = 'AIPS SU '
      call ftibin(iunit, noff, tfields, collabel, 
     :     colformat,
     :     colunit, extname,0, status)
      
      felem = 1
      do ioff = 1,noff
         call ftpclj  (iunit,  1, ioff, felem, 1, ioff, status)
         call ftpcls  (iunit,  2, ioff, felem, 1, chstr(ioff), status)
         call ftpclj  (iunit,  3, ioff, felem, 1, 0, status)
         call ftpcls  (iunit,  4, ioff, felem, 1, '    ', status)
         call ftpcld  (iunit,  5, ioff, felem, 1, 0.0d0, status)
         call ftpcld  (iunit,  6, ioff, felem, 1, 0.0d0, status)
         call ftpcld  (iunit,  7, ioff, felem, 1, 0.0d0, status)
         call ftpcld  (iunit,  8, ioff, felem, 1, 0.0d0, status)
         call ftpcld  (iunit,  9, ioff, felem, 1, 0.0d0, status)
         call ftpcld  (iunit, 10, ioff, felem, 1, 0.0d0, status)
         call ftpcld  (iunit, 11, ioff, felem, 1, off_ra(ioff), status)
         call ftpcld  (iunit, 12, ioff, felem, 1, off_dec(ioff),status)
         call ftpcld  (iunit, 13, ioff, felem, 1, epoch, status)
         call ftpcld  (iunit, 14, ioff, felem, 1, off_app_ra(ioff)
     :        , status)
         call ftpcld  (iunit, 15, ioff, felem, 1, off_app_dec(ioff)
     :        , status)
         
         call ftpcld  (iunit, 16, ioff, felem, 1, 0.0d0, status)
         call ftpcld  (iunit, 17, ioff, felem, 1, 0.0d0, status)
         call ftpcld  (iunit, 18, ioff, felem, 1, 0.0d0, status)
         call ftpcld  (iunit, 19, ioff, felem, 1, 0.0d0, status)
      enddo


      call ftclos(iunit, status)
      call ftfiou(iunit, status)
      call FTGERR(status, comment)
      if (status.ne.0) write(*,*) 'BIN SU table', status, comment


      end
