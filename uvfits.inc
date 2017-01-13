************************************************************************
*
*  AMI Project - Reduce program:  uvfits.inc
*
*  This file contains parameters associated with FITS output format.
*
*  History: 12/7/99  original version [DJT]
*
************************************************************************

c  FITS header items.

      integer bitpix
      integer naxis, naxisn(6)
      integer gcount, pcount

      real*8  blank
      real*8  bscale, bzero
      real*8  uscale, uzero
      real*8  vscale, vzero
      real*8  wscale, wzero

      real*8  jdscal1, jdscal2, jdzero

      common /fits_hdr/
     :      blank, bscale, bzero, uscale, uzero, vscale, vzero,
     :      wscale, wzero, jdscal1, jdscal2, jdzero,
     :      bitpix, naxis, naxisn, gcount, pcount

c   File control, output buffer for FITS records

      integer    ifits, blocks, bsize, mbuff
      parameter (blocks=1, bsize=2880)
      parameter (mbuff=blocks*bsize/4)
      real*4     buffer(mbuff)
      integer    ibuff, nblock

      common /fits_buff/ ifits, buffer, ibuff, nblock

