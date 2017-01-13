*  AMI Project - telescope.inc
*
*  Include file containing fixed telescope parameters.
*
*  Last updated:
*    02/02/04 - adapted from VSA version [DJT]
*    17/11/04 - antenna and enclosure dimensions added [djt]
*    10/01/06 - site parameters included [djt]
*     2/07/06 - revised site parameters, antenna diameter for LA [djt]
*

* Site parameters (AMI-SA)  52 09 51.6 N  00 02 21.7 E

      real*8     SLATN_SA, SLONGE_SA, SHEIGHT_SA
      real*8     ARRAYX_SA, ARRAYY_SA, ARRAYZ_SA
      parameter (SLATN_SA   =  0.910439368775) ! geodetic latitude (radians)
      parameter (SLONGE_SA  =  0.000686980986) ! longitude East (radians)
      parameter (SHEIGHT_SA = 17.000000000000) ! height above sea level (m)
      parameter (ARRAYX_SA  =       3906330.0) ! array geocentric X (m)
      parameter (ARRAYY_SA  =          2950.0) ! array geocentric Y (m)
      parameter (ARRAYZ_SA  =       5030570.0) ! array geocentric Z (m)

* Site parameters (AMI-LA)  52 09 00.0 N  00 02 36.9 E

      real*8     SLATN_LA, SLONGE_LA, SHEIGHT_LA
      real*8     ARRAYX_LA, ARRAYY_LA, ARRAYZ_LA
      parameter (SLATN_LA   =  0.910189205000) ! geodetic latitude (radians)
      parameter (SLONGE_LA  =  0.000760818000) ! longitude East (radians)
      parameter (SHEIGHT_LA = 17.000000000000) ! height above sea level (m)
      parameter (ARRAYX_LA  =       3906330.0) ! array geocentric X (m)
      parameter (ARRAYY_LA  =          2950.0) ! array geocentric Y (m)
      parameter (ARRAYZ_LA  =       5030570.0) ! array geocentric Z (m)

* Antenna dimensions (AMI-SA)

      real       ANT_DIAM_SA, ANT_EXT_DIAM_SA
      real       PIVOT_HEIGHT_SA, ANT_SEC_DIAM_SA
      real       FOC_LEN_SA
      parameter (ANT_DIAM_SA     =  3.70)  ! Antenna diameter (metres)
      parameter (ANT_EXT_DIAM_SA =  3.96)  ! Extended diameter (metres)
      parameter (PIVOT_HEIGHT_SA =  2.122) ! Height of pivot from ground (m)
      parameter (ANT_SEC_DIAM_SA =  0.498) ! Secondary diameter (metres)
      parameter (FOC_LEN_SA      =  1.207) ! Focal length

* Antenna dimensions (AMI-LA)

      real       ANT_DIAM_LA
      real       FOC_LEN_LA
      parameter (ANT_DIAM_LA     = 12.80)  ! Antenna diameter (metres)
      parameter (FOC_LEN_LA      =  3.658) ! Focal length (metres)
       
* Enclosure dimensions (AMI-SA)

      real       ENC_EAST, ENC_WEST, ENC_NORTH, ENC_SOUTH, ENC_HEIGHT
      parameter (ENC_EAST   =   0.0)    ! Outer enclosure East limit (m)
      parameter (ENC_WEST   =  28.0)    ! Outer enclosure West limit (m)
      parameter (ENC_NORTH  =   0.0)    ! Outer enclosure North limit (m)
      parameter (ENC_SOUTH  =  40.0)    ! Outer enclosure South limit (m)
      parameter (ENC_HEIGHT =   4.5)    ! Enclosure height (m)

* Drive limits

      integer    HA_MIN, HA_MAX
      integer    DEC_MIN, DEC_MAX
      integer    EL_LIMIT
      parameter (HA_MIN = -89.00)       ! HA E limit (degrees)
      parameter (HA_MAX =  89.00)       ! HA W limit (degrees)
      parameter (DEC_MIN =  7.00)       ! Dec S limit (degrees)
      parameter (DEC_MAX = 90.00)       ! Dec N limit (degrees)
      parameter (EL_LIMIT= 20.00)       ! Elevation limit during track (degrees)

      integer    HA_MIN_ENC, HA_MAX_ENC
      integer    DEC_MIN_ENC, DEC_MAX_ENC
      parameter (HA_MIN_ENC =  32768)   ! HA E limit (encoder units)
      parameter (HA_MAX_ENC = 115000)   ! HA W limit (encoder units)
      parameter (DEC_MIN_ENC = 20000)   ! Dec minimum (encoder units)
      parameter (DEC_MAX_ENC =120000)   ! Dec maximum (encoder units)

* Encoder resolutions

      real*8     HA_BIT, DEC_BIT
      parameter (HA_BIT  = 9.887695)    ! HA encoder resolution (arcsec)
      parameter (DEC_BIT = 9.887695)    ! Dec encoder resolution (arcsec)

* Maximum permitted PC delay error

      real*4     MAX_PC_ERR
      parameter (MAX_PC_ERR = 40.0)     ! maximum permitted PC error (mm)

* Maximum permitted pointing error

      real*4     MAX_HA_ERR, MAX_DEC_ERR
      parameter (MAX_HA_ERR = 0.05)      ! maximum permitted HA error (degrees)
      parameter (MAX_DEC_ERR = 0.05)     ! maximum permitted Dec error (degrees)

* Threshholds for monitored voltages

      integer    VOLT_MIN, VOLT_MAX
      parameter (VOLT_MIN =  60)        ! minimum value
      parameter (VOLT_MAX = 150)        ! maximum value

* Threshholds for AGC values

      integer    AGC_MIN, AGC_MAX
      parameter (AGC_MIN = -500)        ! minimum value
      parameter (AGC_MAX =  500)        ! maximum value

* Threshhold for rain gauge

      real        RAIN_MIN
      parameter  (RAIN_MIN = 75.)       ! minimum value (% of nominal maximum)
      real        RAIN_MAX
      parameter  (RAIN_MAX = 102.)      ! maximum value (% of nominal maximum)

* Threshholds for cryostat temperature values (degrees K)

      integer    TCRYO_MIN, TCRYO_MAX
      parameter (TCRYO_MIN =  0)        ! minimum value
      parameter (TCRYO_MAX =  50)       ! maximum value

* Threshholds for vacuum levels (voltages)

      real*4     VAC_MIN, VAC_MAX
      parameter (VAC_MIN =  0.0)        ! minimum value
      parameter (VAC_MAX =  2.5)        ! maximum value

* Nominal values for ambient and CMB background temperatures (degrees K)

      real*4     T_ATMOS, T_CMB
      parameter (T_ATMOS = 270.0)
      parameter (T_CMB   = 2.728)

* The optical path length in 'perfect' atmospheric conditions.
      real*8     TAU_IDEAL
      parameter (TAU_IDEAL = 2.4d-2)

* Nominal value for rms noise on a 1 second sample from a correlator board
      real*4     RMS_CORR
      parameter (RMS_CORR = 2.5)  

* CHANGES 28/08/06 - ONLY LO FREQ AND CORR_STEP ARE INDEPENDENT PARAMETERS
* Frequency specifications (Hz)
* Nominal values for IF bandwidth, IF, LO and RF frequencies
      real*8     D_FREQ, IF_FREQ, LO_FREQ, RF_FREQ
c      parameter (D_FREQ = 6.0d9)
c      parameter (IF_FREQ = 9.0d9)
c      parameter (RF_FREQ = 15.0d9)
      parameter (LO_FREQ = 24.0D9)

* Nominal value for correlator lag spacing (m) [changed to 26mm 20/4/08]
      real*8     CORR_STEP
      parameter (CORR_STEP = 26.0d-3)
      parameter (D_FREQ = .2997925000000000D+09/(2*CORR_STEP))
      parameter (IF_FREQ = 1.5*D_FREQ)
      parameter (RF_FREQ = LO_FREQ-IF_FREQ)

* Threshold separation for reporting observation centre close to Sun, Moon
      real       MAX_SEPN
      parameter (MAX_SEPN = 50.0)

* Nominal rms of amplitude of visibility in Jy
      real       NOM_RMS_SA
      real       NOM_RMS_LA
      parameter (NOM_RMS_SA = 0.25)
      parameter (NOM_RMS_LA = 0.025)

