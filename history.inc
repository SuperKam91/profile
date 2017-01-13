************************************************************************
*
*  AMI Project - Reduce program:  history.inc
*
*  This file contains bit values for reduction history keyword.
*
*  History:
*    16/03/03 - adapted from VSA version [DJT]
*    19/05/04 - revised to include AMI reduction steps [DJT]
*    25/06/07 - flags for fringe rotation revised [djt]
*    31/03/08 - flagging operations moved to flagkey [djt]
*
************************************************************************

c  History text array.  This array is filled according to the values
c  in the history keyword, by the set_history routine.

      integer      MAX_HIST
      parameter   (MAX_HIST = 32)
      character*72 history(MAX_HIST)

c  Flagging keyword

      integer*4    flagkey

c  Data reduction history keyword

      integer*4    histkey

c  Length of current history

      integer      lhist

*  Bits within the flag keyword identify which flagging operations have
*  been performed.

      integer*4  FLG_TIMES, FLG_CORR, FLG_BASE, FLG_PCS, FLG_AGCS,
     :           FLG_OFFSET, FLG_POINT, FLG_RAIN, FLG_GAIN, FLG_DATA,
     :           FLG_WEATHER, FLG_SHADOW, FLG_FRINGE

      parameter (FLG_TIMES  =     x'01') !  flag timestamps
      parameter (FLG_CORR   =     x'02') !  flag correlator errors
      parameter (FLG_BASE   =     x'04') !  flag baselines
      parameter (FLG_PCS    =     x'08') !  flag PC spikes
      parameter (FLG_AGCS   =     x'10') !  flag AGC values
      parameter (FLG_POINT  =     x'20') !  flag pointing errors
      parameter (FLG_OFFSET =     x'40') !  flag offset changes
      parameter (FLG_RAIN   =     x'80') !  flag rain gauge values
      parameter (FLG_GAIN   =    x'100') !  flag correlator gains
      parameter (FLG_DATA   =    x'200') !  flag data values (interference)
      parameter (FLG_WEATHER=    x'400') !  flag data values (weather)
      parameter (FLG_SHADOW =    x'800') !  flag enclosure groundspill
      parameter (FLG_FRINGE =   x'1000') !  flag slow fringes

      integer*4  FLG_MASK                !  mask for flagging operations

      parameter (FLG_MASK = 
     :           FLG_TIMES+FLG_CORR+FLG_BASE+FLG_PCS+FLG_AGCS+FLG_POINT+
     :           FLG_OFFSET+FLG_RAIN+FLG_GAIN+FLG_DATA+FLG_WEATHER+
     :           FLG_SHADOW)

*  Bits within the history keyword identify which data reduction steps
*  have been performed on the raw data file.

      integer*4  RED_DEMOD, RED_MEANS, RED_ZEROS, RED_SMOOTH1,
     :           RED_LCAL,  RED_FFT, RED_PCAL, RED_FILTER,
     :           RED_FRROT, RED_FROT_PCS, RED_FROT_PATH,
     :           RED_SMOOTH2, RED_RAIN, RED_ABSORB, RED_SCAL,
     :           RED_REWT, RED_COMB

      parameter (RED_DEMOD  =     x'01') !  rain gauge demodulation
      parameter (RED_MEANS  =     x'02') !  subtract channel means
      parameter (RED_ZEROS  =     x'04') !  subtract zero levels
      parameter (RED_SMOOTH1=     x'08') !  initial smoothing (integration)
      parameter (RED_LCAL   =     x'10') !  apply correlator lag gains
      parameter (RED_FFT    =     x'20') !  perform FFT to frequency domain
      parameter (RED_PCAL   =     x'40') !  apply correlator baseline gains
      parameter (RED_FILTER =     x'80') !  filter out slow fringes
      parameter (RED_FROT_PATH=  x'100') !  fringe rotation for path
      parameter (RED_FROT_PCS=   x'200') !  fringe rotation for PCs
      parameter (RED_FRROT  =    x'300') !  fringe rotation (path+PCs)
      parameter (RED_SMOOTH2=    x'400') !  secondary smoothing (integration)
      parameter (RED_RAIN   =    x'800') !  atmospheric noise correction
      parameter (RED_ABSORB =   x'1000') !  atmospheric absorption correction
      parameter (RED_SCAL  =    x'2000') !  secondary phase calibration
      parameter (RED_REWT  =    x'4000') !  reweight data
      parameter (RED_COMB  =    x'8000') !  combine +/- correlations

      integer*4  RED_MASK            !  mask for data reduction operations
      integer*4  FROT_MASK           !  mask for fringe rotation operations

      parameter (RED_MASK = 
     :           RED_DEMOD+RED_MEANS+RED_ZEROS+RED_SMOOTH1+RED_LCAL+
     :           RED_FFT+RED_PCAL+RED_FILTER+RED_FRROT+RED_SMOOTH2+
     :           RED_RAIN+RED_ABSORB+RED_SCAL+RED_REWT+RED_COMB)

      parameter (FROT_MASK = RED_FROT_PATH+RED_FROT_PCS)

      common /history / history, flagkey, histkey, lhist

