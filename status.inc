*  VSA Project - status.inc
*
*  Status values for telescope and micro communications.
*
*  Last updated: 29/01/02 this version for VSA tests [DJT]
 
* Telescope observing status values

      integer    SHUTDOWN, IDLE, STOWED, STARTING, OBSERVING, COMPLETE
      parameter (SHUTDOWN=0)            ! Control software shut down
      parameter (IDLE=1)                ! No observation in progress
      parameter (STOWED=2)              ! Telescope stowed
      parameter (STARTING=3)            ! Observation starting
      parameter (OBSERVING=4)           ! Observation in progress
      parameter (COMPLETE=5)            ! Observation complete

* Micro communication status values

      integer    DISCONNECT, NOCONNECT, CONNECTED
      parameter (DISCONNECT=-1)         ! Can't communicate with micro
      parameter (NOCONNECT=0)           ! Not connected with micro
      parameter (CONNECTED=1)           ! Connected with micro

* Correlator status values

      integer    RECORDING
      parameter (RECORDING=2)           ! Recording data to file

* Drive control values

      integer    DRIVE, TRACK
      parameter (DRIVE=2)               ! Antennas driving Az,El
      parameter (TRACK=3)               ! Antennas tracking RA,Dec

* Drive status values

      integer    DRIVING, ON_SOURCE, IN_LIMITS
      parameter (DRIVING=2)             ! Antennas driving
      parameter (ON_SOURCE=3)           ! Antennas on source
      parameter (IN_LIMITS=4)           ! Antenna drives in limits

* Offset regime values

      integer    NONE, UNSET, DEFOCUS
      integer    DRIFT, PC_DRIFT, PC_OFFST, PC_TEST, HA_OFFST, RA_OFFST
      integer    DEC_OFFST, FIVE_PT, SEVEN_PT, CAL_OFFST, DELAY_OFFST
      integer    RASTER_OFFST, HA_RASTER, RA_RASTER, HEX_RASTER
      integer    SRC_LIST
      parameter (UNSET=-1)              ! Observing mode unset
      parameter (NONE=0)                ! Default tracking mode
      parameter (DRIFT=1)               ! Drift observation mode
      parameter (PC_DRIFT=2)            ! PC drift observation
      parameter (PC_OFFST=3)            ! PC offset observation
      parameter (HA_OFFST=4)            ! Hour angle offset observation
      parameter (RA_OFFST=4)            ! Hour angle offset observation
      parameter (DEC_OFFST=5)           ! Declination offset observation
      parameter (FIVE_PT=6)             ! 5-point offset observation
      parameter (RASTER_OFFST=7)        ! Raster offset observation
      parameter (HA_RASTER=7)           ! Raster offset  observation
      parameter (RA_RASTER=7)           ! Raster offset observation
      parameter (PC_TEST=8)             ! PC bit sweep test
      parameter (CAL_OFFST=9)           ! Interleaved calibrator observation
      parameter (SRC_LIST=10)           ! List of observing centres
      parameter (DELAY_OFFST=11)        ! PC delay zero offset
      parameter (SEVEN_PT=12)           ! 7-point offset observation
      parameter (HEX_RASTER=13)         ! Hexagonal raster observation
      parameter (DEFOCUS=14)            ! Defocus mode (pointing offsets)

* Array identifiers

      integer    AMI_SA, AMI_LA
      parameter (AMI_SA=1)              ! AMI small array
      parameter (AMI_LA=2)              ! AMI large array

* Main array configuration

c     integer    COMPACT, EXTENDED
c     parameter (COMPACT=1)             ! Compact configuration
c     parameter (EXTENDED=2)            ! Extended configuration

* Observation type

      integer    POINT, FIELD
      parameter (POINT=1)               ! Point source observation
      parameter (FIELD=2)               ! Field observation

