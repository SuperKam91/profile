*+
       subroutine SET_HISTORY (status)
c
c  AMI Project - Reduce program, set processing history text
c
c  This routine sets up a text array with the current processing history
c  for display or output to reduced data file.  The history gives a text
c  record of the data reduction steps performed so far on the current
c  data file.
c
c  History:
c    20/02/02 - adapted from VSA reduce [DJT].
c    19/05/04 - added data reduction entries [DJT].
c    25/10/07 - added atmospheric absorption correction [DJT].
c
*-
       integer    status

       include 'header.inc'
       include 'version.inc'
       include 'reduce.inc'
       include 'history.inc'
       include 'status.inc'

       character chstr*16
       integer   l, ls, ns

       integer   chr_lenb
       external  chr_lenb

       if (status.ne.0) return

c  Report software version

       l = 1
       history(l) = 'AMI reduce, version '//versid

       if (histkey.eq.0) then
          l = l+1
          history(l) = 'Raw data, no reduction performed'
       endif

c  Data file

      l = l+1
      ls = chr_lenb(obs_file)
      history(l) = 'Sample file '//obs_file(1:ls)

c  Flagging operations

       if (and(flagkey,FLG_TIMES).gt.0) then
          l = l+1
          history(l) = 'Timestamp errors corrected'
       endif
       if (and(flagkey,FLG_CORR).gt.0) then
          l = l+1
          history(l) = 'Correlator errors flagged'
       endif
       if (and(flagkey,FLG_GAIN).gt.0) then
          l = l+1
          history(l) = 'Correlator gains flagged'
       endif
       if (and(flagkey,FLG_BASE).gt.0) then
          l = l+1
          history(l) = 'Faulty baselines flagged'
       endif
       if (and(flagkey,FLG_AGCS).gt.0) then
          l = l+1
          history(l) = 'AGC values flagged'
       endif
       if (and(flagkey,FLG_OFFSET).gt.0) then
          l = l+1
          history(l) = 'Offset changes flagged'
       endif
       if (and(flagkey,FLG_PCS).gt.0) then
          l = l+1
          history(l) = 'PC errors flagged'
       endif
       if (and(flagkey,FLG_POINT).gt.0) then
          l = l+1
          history(l) = 'Pointing errors flagged'
       endif
       if (and(flagkey,FLG_SHADOW).gt.0) then
          l = l+1
          history(l) = 'Shadowing effects flagged'
       endif
       if (and(flagkey,FLG_FRINGE).gt.0) then
          l = l+1
          history(l) = 'Slow fringes flagged'
       endif
       if (and(flagkey,FLG_RAIN).gt.0) then
          l = l+1
          history(l) = 'Rain gauge values flagged'
       endif
       if (and(flagkey,FLG_DATA).gt.0) then
          l = l+1
          history(l) = 'Data amplitudes flagged'
       endif

c  Data reduction operations

       if (and(histkey,RED_DEMOD).gt.0) then
          l = l+1
          history(l) = 'Rain gauge demodulation performed'
       endif
       if (and(histkey,RED_MEANS).gt.0) then
          l = l+1
          history(l) = 'Correlator lag means subtracted'
       endif
       if (and(histkey,RED_ZEROS).gt.0) then
          l = l+1
          history(l) = 'Correlator zero levels subtracted'
       endif
       if (and(histkey,RED_SMOOTH1).gt.0) then
          l = l+1
          ns = nsmooth1
          call chr_chitoc(ns,chstr,ls)
          history(l) = 'Integration x '//chstr(1:ls)//' performed'
       endif
       if (and(histkey,RED_LCAL).gt.0) then
          l = l+1
          history(l) = 'Correlator lag gains applied'
          ls = chr_lenb(lcal_file)
          if (ls.gt.0) then
             l = l+1
             history(l) = 'Calibration file '//lcal_file(1:ls)
          endif
       endif
       if (and(histkey,RED_FFT).gt.0) then
          l = l+1
          history(l) = 'Transform to frequency domain performed'
       endif
       if (and(histkey,RED_PCAL).gt.0) then
          l = l+1
          history(l) = 'Baseline gains applied'
          ls = chr_lenb(pcal_file)
          if (ls.gt.0) then
             l = l+1
             history(l) = 'Calibration file '//pcal_file(1:ls)
          endif
       endif
       if (and(histkey,RED_FILTER).gt.0) then
          l = l+1
          history(l) = 'Slow fringe rates filtered'
       endif
       if (and(histkey,RED_FRROT).eq.RED_FRROT) then
          l = l+1
          history(l) = 'Fringe rotation to phase centre performed'
       else
          if (and(histkey,RED_FROT_PATH).gt.0) then
             l = l+1
             history(l) = 'Fringe rotation for path performed'
          endif
          if (and(histkey,RED_FROT_PCS).gt.0) then
             l = l+1
             history(l) = 'Fringe rotation for PCs performed'
          endif
       endif
       if (and(histkey,RED_SMOOTH2).gt.0) then
          l = l+1
          ns = nsmooth2
          call chr_chitoc(ns,chstr,ls)
          history(l) =
     :        'Secondary integration x '//chstr(1:ls)//' performed'
       endif
       if (and(histkey,RED_RAIN).gt.0) then
          l = l+1
          history(l) = 'Rain gauge correction applied'
       endif
       if (and(histkey,RED_ABSORB).gt.0) then
          l = l+1
          history(l) = 'Atmospheric absorption correction applied'
       endif
       if (and(histkey,RED_SCAL).gt.0) then
          l = l+1
          history(l) = 'Secondary phase correction applied'
       endif
       if (and(histkey,RED_REWT).gt.0) then
          l = l+1
          history(l) = 'Data reweighted according to baseline rms'
       endif
       if (and(histkey,RED_COMB).gt.0) then
          l = l+1
          history(l) = 'Combination of +/- correlations performed'
       endif

       lhist = l

       end
