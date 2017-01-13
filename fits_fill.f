*+
      subroutine FITS_FILL (data, status)
c
c  AMI Project - Reduce program, fill FITS block
c
c  Pads the current FITS buffer using input data, and writes to file.
c
c  History:
c    12/07/99 - original VSA version [DJT].
c     7/05/06 - use io_fwrite for Linux compatibility [djt]
c
*-
       real*4   data
       integer  status

       include 'uvfits.inc'

       integer  i, block, bytes, ndata, offset

       if (status.ne.0) return

       if (ibuff.eq.0) return

c  Fill current buffer and write to file

       block = bsize/4
       ndata = ((ibuff-1)/block+1)*block

       if (ibuff.lt.ndata) then
          do i = ibuff+1, ndata
             buffer(i) = data
          enddo
       endif

       bytes = ndata*4
       offset = nblock*bsize
       call io_fwrite(ifits, offset, buffer, bytes, status)
       if (status.ne.0) call io_wrerr(status,' on write')
       nblock = nblock + ndata/block
       ibuff = 0

       end
