!     *+
      subroutine FITS_WRANT 
     :     (iscr,n_antennas,x_pos,y_pos,z_pos,frq,status)
!c
!c  AMI Project - Reduce program, write antenna table to FITS file.
!c
!c  Construct extension table containing details of antenna geometry
!c  and appends to FITS file.
!c
!c  NB we employ the device of including two entries for each antenna,
!c  to provide AIPS with separate identifiers for '+' and '-' correlated
!c  spacings.  The '-' entries have station numbers offset by 10, i.e.
!c  ant 11->20 for the SA, ant 11->18 for the LA.
!c
!c  Coordinates are with respect to ground:
!c        x-axis  East -> West,
!c        y-axis  North -> South,
!c        z-axis  towards zenith.
!c
!c  History:
!c    17/01/05 - modified from VSA code [djt]
!c    18/07/07 - write out all antennas (not just working set) [djt]
!c    18/07/07 - add extra entries for 'minus' correlations [djt]
!c     7/08/07 - add header entries for geocentric array coordinates [djt]
!c    25/04/08 - add header entries for polarisation feed [djt/dag]
!c     6/06/08 - converted to binary table format [djt]
!c
!*-
      implicit none
      integer status, n_antennas
      real*4  x_pos(n_antennas),y_pos(n_antennas),z_pos(n_antennas)
      real*8  frq, freq
      integer   iscr, iunit

      include 'uvfits.inc'
      include 'reduce.inc'
      include 'history.inc'
      include 'header.inc'
      include 'status.inc'
      include 'telescope.inc'
      include 'constants.inc'

      character scratch*32
      character anname*8, station*8
      character poltya*4, poltyb*4
      real*8    stabxyz(3), gstia0, degpdy, fd, ut0
      real*4    polaa, polab, polcala(2), polcalb(2), staxof
      integer   nant, nosta, mant, mntsta
      integer   i, ii, id, im, iy, ls
      integer   idate(3)
      integer*8 ia, felem

      real*8    sla_gmst
      external  sla_gmst

      integer   nline
      parameter (nline=80)
      character line*(nline)
      real*4    rline(nline/4)
      equivalence (line,rline)

      character filename*120 

      integer   nbyte
      parameter (nbyte=76)
      real*4    record(nbyte/4)

      common /fits_ant/ anname,stabxyz,nosta,mntsta,staxof,
     :                  poltya,polaa,polcala,poltyb,polab,polcalb
      equivalence (record,anname)

c Bin table vars
      character collabel(12)*10
      character colformat(12)*10
      character colunit(12)*10
      integer tfields
      character extname*10
      integer firstrow, firstelem

c keyword varables
      character comment*72
      character value*16
      integer hdutype

      integer naxes(12)

c  Entry
      if(status.ne.0) return

      nant = n_antennas
      mant = n_antennas
      freq = frq

c Open fits

      call ftgiou(iunit, status)
      !if (status.ne.0) write(*,*) 'wrant,allocating new',iunit,status
      call FTREOPEN(iscr, iunit, status)
      !if (status.ne.0) write(*,*) 'wrant,reopen',iunit,iscr,status


c Setup Axis labels/properties for AN table

      collabel(1)  = 'ANNAME  '
      collabel(2)  = 'STABXYZ '
      collabel(3)  = 'ORBPARM '
      collabel(4)  = 'NOSTA   '
      collabel(5)  = 'MNTSTA  '
      collabel(6)  = 'STAXOF  '
      collabel(7)  = 'POLTYA  '
      collabel(8)  = 'POLAA   '
      collabel(9)  = 'POLCALA '
      collabel(10) = 'POLTYB  '
      collabel(11) = 'POLAB   '
      collabel(12) = 'POLCALB '


      colformat(1)  = '8A      '
      colformat(2)  = '3D      '
      colformat(3)  = '0D      '
      colformat(4)  = '1J      '
      colformat(5)  = '1J      '
      colformat(6)  = '1E      '
      colformat(7)  = '4A      '
      colformat(8)  = '1E      '
      colformat(9)  = '2E      '
      colformat(10) = '4A      '
      colformat(11) = '1E      '
      colformat(12) = '2E      '


      colunit(1)  = '        '
      colunit(2)  = 'METERS  '
      colunit(3)  = '        '
      colunit(4)  = '        '
      colunit(5)  = '        '
      colunit(6)  = 'METERS  '
      colunit(7)  = 'DEGREES '
      colunit(8)  = 'DEGREES '
      colunit(9)  = '        '
      colunit(10) = '        '
      colunit(11) = 'DEGREES '
      colunit(12) = '        '
      
      tfields = 12

c Create new BIN table for Antenna positions

      extname = 'AIPS AN '
      call FTGREC(iunit,1, filename,status)
      call FTMAHD(iunit,1, hdutype,status)
      if (status.ne.0) write(*,*) 'Moving', status
      call ftibin(iunit, n_antennas, tfields, collabel, 
     :     colformat,
     :     colunit, extname,0, status)
      if (status.ne.0) write(*,*) 'ftibin', status
      call FTGERR(status, filename)
      if (status.ne.0) write(*,*) 'FTGERR', status, filename


      ut0 = int(obs_mjd_start)
      gstia0 = sla_GMST(ut0)/const_d2r
      degpdy = 360.d0*const_sut2sst
      call util_enqdat(idate)
      call sla_CALDJ(idate(3), idate(2), idate(1),obs_mjd_start, i)
      call sla_djcl(obs_mjd_start,iy,im,id,fd,i)


      call ftpkyd(iunit,'ARRAYX',ARRAYX_SA,-10,
     :     'Array centre X coordinate',status)
      call ftpkyd(iunit,'ARRAYY',ARRAYY_SA,-10,
     :     'Array centre Y coordinate',status)
      call ftpkyd(iunit,'ARRAYZ',ARRAYZ_SA,-10,
     :     'Array centre Z coordinate',status)
      call ftpkyd(iunit,'GSTIA0',gstia0,-10,
     :     'GST at time=0 on ref date(deg)',status) 
      call ftpkyd(iunit,'DEGPDY',degpdy,-10,
     :     'Earth rotation rate (deg/day)',status) 
      call ftpkyd(iunit,'FREQ',freq,-10,
     :     'Obs. reference frequency (Hz)',status)
      write(value, '(I4 ,("-",I2.2),("-",I2.2))')  iy, im, id  
      call ftpkys
     :     (iunit,'RDATE', value,'Reference date',status)
      call ftpkyd (iunit,'POLARX',0.0d0,-10,
     :     'Polar position X (m) on ref date',status)
      call ftpkyd(iunit,'POLARY',0.0d0,-10,
     :     'Polar position Y (m) on ref date',status)
      !write(value,'(E16.9)'), 0.0 
      call ftpkyd(iunit,'UT1UTC',0.0d0,-10,'UT1-UTC (sec)',status)
      write(value,'(E16.9)'), 0.0 
      call ftpkyd(iunit,'DATUTC',0.0d0,-10,'Data time-UTC (sec)',status)
      call ftpkys
     :     (iunit,'TIMSYS','UTC','Time system',status)
      call ftpkys
     :     (iunit,'ARRNAM','Profile','Array name',status) 
      call ftpkyj
     :     (iunit,'NUMORB',0,'No of orbital parameters',status) 
      call ftpkyj
     :     (iunit,'NOPCAL',2,'No of polarisation constants',status)
      call ftpkys
     :     (iunit,'POLTYPE','X-Y LIN ','Feed polarisation',status)

c  Write comment

      write(comment, '(A)') 'Coordinates with respect to ground:'
      call ftpcom(iunit,comment, status)
      write(comment, '(A)') ' x-axis  East -> West'
      call ftpcom(iunit,comment, status)
      write(comment, '(A)') ' y-axis  North -> South'
      call ftpcom(iunit,comment, status)
      write(comment, '(A)') ' z-axis  towards zenith'
      call ftpcom(iunit,comment, status)


c  End of header records

      !call ftpkys(iunit,'END', '','',status)

c  Write antenna data records to FITS output file

c  Station identifier

      station = ' '
      station = 'Profile'

c  Polarisation parameters

      mntsta = 1
      staxof = 0.0
      poltya = 'X'
      poltyb = 'Y'
      polaa = 0.0
      polab = 90.0
      do i = 1,2
         polcala(i) = 0.0
         polcalb(i) = 0.0
      enddo


c  Convert to standard IEEE byte order

c$$$      call util_i4cvt(mntsta,1)
c$$$      call util_r4cvt(staxof,1)
c$$$      call util_r4cvt(polaa,1)
c$$$      call util_r4cvt(polab,1)
c$$$      call util_r4cvt(polcala,2)
c$$$      call util_r4cvt(polcalb,2)

c    Write antenna records 
      felem = 1
      station = 'SimTel'
      do ia = 1,n_antennas
         nosta = ia
         anname = ' '
         stabxyz(1) = x_pos(ia)/1.E3
         stabxyz(2) = y_pos(ia)/1.E3
         stabxyz(3) = z_pos(ia)/1.E3
         call util_i4cvt(nosta,1)
         call util_r8cvt(stabxyz,3)
         call chr_chitoc(ia,station(7:),ls)
         if (ia.le.nant) anname = station
         call ftpcls (iunit,  1, ia, felem, 1, anname, status)
         call ftpcld (iunit,  2, ia, felem, 3, stabxyz, status)
         call ftpclj (iunit,  4, ia, felem, 1, nosta, status)
         call ftpclj (iunit,  5, ia, felem, 1, mntsta, status)
         call ftpcle (iunit,  6, ia, felem, 1, staxof, status)
         call ftpcls (iunit,  7, ia, felem, 1, poltya, status)
         call ftpclj (iunit,  8, ia, felem, 1, polaa, status)
         call ftpcle (iunit,  9, ia, felem, 2, polcala, status)
         call ftpcls (iunit, 10, ia, felem, 1, poltyb, status)
         call ftpclj (iunit, 11, ia, felem, 1, polab, status)
         call ftpcle (iunit, 12, ia, felem, 2, polcalb, status)
         
      enddo

c  Fill last block with blanks

      call FTGERR(status, comment)
      !if (status.ne.0)write(*,*) 'wrant', status, comment
      call ftclos(iunit, status)
      call ftfiou(iunit, status)
      !if (status.ne.0)write(*,*) 'ftfiou', status


      end
