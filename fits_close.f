*+
      subroutine FITS_CLOSE (status)
c
c  AMI Project - Reduce program, close FITS output file
c
c  History:  12/7/99 - original VSA version [DJT].
c
*-
       integer  status

       include 'uvfits.inc'

       character chstr*8
       integer   iout, ls

       if (status.ne.0) return
      
       call io_enqout(iout)
       
c  Close output file

       call chr_chitoc(nblock,chstr,ls)
       write(iout,'(X,A)')
     :            chstr(1:ls)//' FITS blocks of 2880 bytes written'
       
       close(ifits)

       end
