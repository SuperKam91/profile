*+
      subroutine FITS_WRITE (data, nbytes, status)
c
c  AMI Project - Reduce program, write data to FITS output file
c
c  Writes data to output buffer and thence to file when full.
c
c  History:
c    12/07/99 - original VSA version [DJT].
c     7/05/06 - use io_fwrite for Linux compatibility [djt]
c
*-
       real*4   data(*)
       integer  nbytes, status

       include 'uvfits.inc'

       integer  i, bytes, ndata, offset

       if (status.ne.0) return

       if (nbytes.le.0) return

c  Write data to output buffer, write to file when full

       ndata = (nbytes-1)/4+1

       do i = 1, ndata
          ibuff = ibuff+1
          buffer(ibuff) = data(i)
          if (ibuff.ge.mbuff) then
             bytes = mbuff*4
             offset = nblock*bsize
             call io_fwrite(ifits, offset, buffer, bytes, status)           
ccc             write(*,*) ifits, offset, bytes, status
             if (status.ne.0) call io_wrerr(status,' on write')
             nblock = nblock+blocks
             ibuff = 0
          endif
       enddo

       end
