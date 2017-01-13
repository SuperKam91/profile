*+
      subroutine FITS_OPEN (file, access, status)

c  AMI Project - Reduce program, open file for UVFITS output
c
c  Open output file, set file and output buffer byte pointers.
c
c  History:
c    12/07/99 - original VSA version [DJT].
c    25/04/01 - give warning if cnnot open file [KG].
c    23/08/06 - overwrite (delete) existing file [djt].
c
*-
       character file*(*), access*(*)
       integer   status

       include 'uvfits.inc'

       integer   lf
       logical   exists

       integer   chr_lenb
       external  chr_lenb

       if (status.ne.0) return

c  Check whether file already exists

       lf = chr_lenb(file)
       inquire ( file = file(1:lf), exist = exists, iostat = status )
       if (exists) then
          status = 0
          call io_delfil(file, 0, status)
       endif

c  Open output file

       call io_nxtlun(ifits, status)

       if (access.eq.'WRITE') then

          open(ifits, file=file, status='NEW', iostat=status)
          if (status.ne.0) then
             write(*,*) 'Cannot open FITS file'
             ifits = 0
          endif

       endif

       nblock = 0
       ibuff = 0

       end
