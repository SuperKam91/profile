*+
      subroutine RUN_CMD (filename, status)

c  AMI Project - Reduce program, run command file
c
c  History:
c    14/05/06 - original version [DJT].
c    23/11/06 - filename length increased to 64 chars [djt]
c    10/07/07 - filename length increased to 80 chars [djt/jz]
c    18/02/08 - filename length increased to 120 chars [djt/jz]
c
*-
       character filename*(*)
       integer   status

c      include  'reduce.inc'

       integer   icmd, iout

       if (status.ne.0) return

       call io_enqout(iout)

c  Open as command input

       call io_nxtlun(icmd,status)
       open (icmd, file=filename, status='OLD', iostat=status)
       if (status.ne.0) then
          write(iout,*) '*** problem opening ',filename
          call io_clrerr(status)
       else
	  write(iout,'(X,A)')filename
          call io_setin(icmd)
       endif

       end
