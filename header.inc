*  AMI Project - header.inc
*
*  Data file header specification.
*
*  This Fortran include file must correspond exactly to the data structures
*  defined in the C header file /ami/include/data.h.
*
*  This format is common to both the Small and Large arrays.
*
*  History:
*     3/09/03 - for new data format [DJT]
*    21/07/05 - MAXBASE introduced [djt]
*     3/08/05 - MAXRACK increased to 16 [djt]
*     4/08/05 - format version 2, added src_list structure [djt]
*    15/01/07 - size of pc_delay table modified for LA [djt]
*
      integer MAXANT                    !  maximum no.of array antennas
      parameter (MAXANT=10)
      integer MAXRACK                   !  maximum no.of correlator racks
      parameter (MAXRACK=16)
      integer MAXBORD                   !  maximum no.of boards per rack
      parameter (MAXBORD=6)
      integer MAXCHAN                   !  maximum no.of data channels
      parameter (MAXCHAN=16)
      integer MAXBASE                   !  maximum no.of baselines
      parameter (MAXBASE=(MAXANT*(MAXANT-1))/2)
      integer MAXSP                     !  maximum no.of spacings
      parameter (MAXSP=2*MAXBASE)
      integer MAXSRC                    !  maximum size of source list
      parameter (MAXSRC=16)
      integer MAXPC                     !  maximum number of PC bits
      parameter (MAXPC=12)
      integer NTABPC                    !  size of PC delay lookup tables
      parameter (NTABPC=100)

      integer MAXANT_SA                 !  maximum no.of SA antennas
      parameter (MAXANT_SA=10)
      integer MAXANT_LA                 !  maximum no.of LA antennas
      parameter (MAXANT_LA=8)
      integer MAXBASE_SA                !  maximum no.of SA baselines
      parameter (MAXBASE_SA=(MAXANT_SA*(MAXANT_SA-1))/2)
      integer MAXBASE_LA                !  maximum no.of LA baselines
      parameter (MAXBASE_LA=(MAXANT_LA*(MAXANT_LA-1))/2)
      integer MAXRACK_SA                !  maximum no.of SA racks
      parameter (MAXRACK_SA=16)
      integer MAXRACK_LA                !  maximum no.of LA racks
      parameter (MAXRACK_LA=10)
      integer PCBITS                    !  no.of SA PC delay bits
      parameter (PCBITS=10)
      integer PCBITS_SA                 !  no.of SA PC delay bits
      parameter (PCBITS_SA=10)
      integer PCBITS_LA                 !  no.of LA PC delay bits
      parameter (PCBITS_LA=12)

*  Data file header

      integer hdr_len
      parameter (hdr_len=4608/4)
      integer*4 hdr(hdr_len)

      integer rt_ptr,rt_len
      integer obs_ptr,obs_len
      integer src_ptr,src_len
      integer tel_ptr,tel_len
      integer pc_ptr,pc_len
      parameter (rt_ptr=1, rt_len=24/4)
      parameter (obs_ptr=rt_ptr+rt_len, obs_len=568/4)
      parameter (src_ptr=obs_ptr+obs_len, src_len=1352/4)
      parameter (tel_ptr=src_ptr+src_len, tel_len=1976/4)
      parameter (pc_ptr=tel_ptr+tel_len, pc_len=400/4)

      integer*4 version                 !  format version
      integer*4 hdr_size                !  length of file header (bytes)
      integer*4 smp_size                !  length of sample data (bytes)
      integer*4 mon_size                !  length of monitor data (bytes)
      integer*4 mon_flag                !  key to monitored data
      integer*4 proc_flag               !  key to processing history
      integer*4 obs(obs_len)            !  description of observation
      integer*4 src(src_len)            !  description of source list
      integer*4 tel(tel_len)            !  description of telescope
      real*4 pc_delay(NTABPC)           !  pc delay table

      equivalence (hdr(1), version)
      equivalence (hdr(2), hdr_size)
      equivalence (hdr(3), smp_size)
      equivalence (hdr(4), mon_size)
      equivalence (hdr(5), mon_flag)
      equivalence (hdr(6), proc_flag)
      equivalence (hdr(obs_ptr), obs)
      equivalence (hdr(src_ptr), src)
      equivalence (hdr(tel_ptr), tel)
      equivalence (hdr(pc_ptr), pc_delay)

*  Description of observation

      integer*4 obs_version            !  format version
      integer*4 obs_status             !  observing status
      character*32 obs_name            !  observation title
      character*64 obs_file            !  data file name
      character*32 obs_observer        !  observer name
      character*256 obs_comment        !  text comment
      integer*4 obs_array              !  array flag (SA or LA)
      integer*4 obs_type               !  observation type (calibrator or field)
      integer*4 obs_planet             !  planet identifier
      integer*4 obs_drive              !  drive tracking flag
      integer*4 obs_gainsw             !  gain switch state (HI/LO)
      integer*4 obs_epoch              !  reference epoch (B1950,J2000)
      real*8 obs_raref, obs_decref     !  RA,Dec at reference epoch (radians)
      real*8 obs_raobs, obs_decobs     !  RA,Dec of observation date (radians)
      real*8 obs_rarot, obs_decrot     !  RA,Dec of phase centre (radians)
      real*8 obs_ha_start              !  HA at start of observation (deg)
      real*8 obs_dec_start             !  Dec at start of observation (deg)
      real*8 obs_mjd_start             !  mJd at start of observation
      real*8 obs_mjd_stop              !  mJd at start of observation
      real*8 obs_st_start              !  LAST of start of observation (rad)
      real*8 obs_st_stop               !  LAST of end of observation (rad)
      integer*4 obs_ut_start           !  UT of start of observation (sec)
      integer*4 obs_ut_stop            !  UT of end of observation (sec)
      integer*4 obs_off_mode           !  observing mode
      integer*4 obs_off_ant            !  key to offset antennas
      integer*4 obs_off_nstep          !  number of phases in offset cycle
      integer*4 obs_off_nsamp          !  number of samples at each phase
      integer*4 obs_off_step           !  offset step (arcsec or PC units)
      integer*4 obs_off_nstepx         !  number of steps in raster pattern
      integer*4 obs_off_nstepy         !
      integer*4 obs_off_stepx          !  offset step in raster (arcsec)
      integer*4 obs_off_stepy          !
      integer*4 obs_off_nsampc         !  number of samples at field centre
      integer*4 obs_intsam             !  number of samples per second
      integer*4 obs_nsamp              !  total number of samples

      equivalence (obs(1), obs_version)
      equivalence (obs(2), obs_status)
c      equivalence (obs(3), obs_name)
      equivalence (obs(11), obs_file)
      equivalence (obs(27), obs_observer)
      equivalence (obs(35), obs_comment)
      equivalence (obs(99), obs_array)
      equivalence (obs(100), obs_type)
      equivalence (obs(101), obs_planet)
      equivalence (obs(102), obs_drive)
      equivalence (obs(103), obs_gainsw)
c      equivalence (obs(104), obs_epoch)
      equivalence (obs(105), obs_raref)
      equivalence (obs(107), obs_decref)
      equivalence (obs(109), obs_raobs)
      equivalence (obs(111), obs_decobs)
      equivalence (obs(113), obs_rarot)
      equivalence (obs(115), obs_decrot)
      equivalence (obs(117), obs_ha_start)
      equivalence (obs(119), obs_dec_start)
      equivalence (obs(121), obs_mjd_start)
      equivalence (obs(123), obs_mjd_stop)
      equivalence (obs(125), obs_st_start)
      equivalence (obs(127), obs_st_stop)
      equivalence (obs(129), obs_ut_start)
      equivalence (obs(130), obs_ut_stop)
      equivalence (obs(131), obs_off_mode)
      equivalence (obs(132), obs_off_ant)
      equivalence (obs(133), obs_off_nstep)
      equivalence (obs(134), obs_off_nsamp)
      equivalence (obs(135), obs_off_step)
      equivalence (obs(136), obs_off_nstepx)
      equivalence (obs(137), obs_off_nstepy)
      equivalence (obs(138), obs_off_stepx)
      equivalence (obs(139), obs_off_stepy)
      equivalence (obs(140), obs_off_nsampc)
      equivalence (obs(141), obs_intsam)
      equivalence (obs(142), obs_nsamp)

*  Source list structure

      integer*4 src_version            !  format version
      integer*4 src_nsrc               !  number of sources
      character*64 src_file            !  source list file name
      character*16 src_name(MAXSRC)    !  source identifier
      integer*4 src_planet(MAXSRC)     !  planet identifier
      integer*4 src_nsamp(MAXSRC)      !  samples on this source
      integer*4 src_epoch(MAXSRC)      !  reference epoch (B1950/J2000)
      real*8 src_raref(MAXSRC)         !  RA at reference epoch
      real*8 src_decref(MAXSRC)        !  Dec at reference epoch
      real*8 src_raobs(MAXSRC)         !  RA at observation date
      real*8 src_decobs(MAXSRC)        !  Dec at observation date
      real*4 src_flux(MAXSRC)          !  nominal flux (Jy)
      character*16 src_flags(MAXSRC)   !  flags field

      equivalence (src(1), src_version)
      equivalence (src(2), src_nsrc)
      equivalence (src(3), src_file)
      equivalence (src(19), src_name)
      equivalence (src(19+4*MAXSRC), src_planet)
      equivalence (src(19+5*MAXSRC), src_nsamp)
      equivalence (src(19+6*MAXSRC), src_epoch)
      equivalence (src(19+7*MAXSRC), src_raref)
      equivalence (src(19+9*MAXSRC), src_decref)
      equivalence (src(19+11*MAXSRC), src_raobs)
      equivalence (src(19+13*MAXSRC), src_decobs)
      equivalence (src(19+15*MAXSRC), src_flux)
      equivalence (src(19+16*MAXSRC), src_flags)

*  Description of telescope

      integer*4 tel_version            !  format version
      integer*4 tel_array              !  array flag
      real*8 tel_along                 !  longitude east (radians)
      real*8 tel_glat                  !  geodetic latitude (radians)
      real*8 tel_height                !  elevation above sea level (metres)
      real*8 tel_ax(MAXANT)            !  antenna X table coordinate (mm)
      real*8 tel_ay(MAXANT)            !  antenna Y table coordinate (mm)
      real*8 tel_az(MAXANT)            !  antenna Z table coordinate (mm)
      real*8 tel_ha_axis(2,MAXANT)     !  HA axis alignment (radians)
      real*8 tel_dec_axis(2,MAXANT)    !  Dec axis alignment (radians)
      real*4 tel_freq                  !  nominal observing frequency (MHz)
      real*4 tel_freqIF(MAXANT)        !  nominal 1st IF frequency (MHz)
      real*4 tel_tsys(MAXANT)          !  nominal system temperatures (deg K)
      real*4 tel_rain(MAXANT)          !  nominal rain gauge maximum values
      real*4 tel_gain(MAXSP)           !  correlator amplitude corrections
      real*4 tel_phase(MAXSP)          !  correlator phase corrections (radians)
      integer*4 tel_ha_zero(MAXANT)    !  azimuth fixed offsets (encoder units)
      integer*4 tel_dec_zero(MAXANT)   !  elevation fixed offset (encoder units)
      integer*4 tel_pc_zero(MAXANT)    !  path compensation nominal delays (mm)
      integer*4 tel_ant(MAXANT)        !  antenna status (in use, etc)
      integer*4 tel_corr(2,MAXBASE)    !  correlator configuration map
      integer*4 tel_gainsw             !  gain switch state (HI/LO)
      integer*4 tel_nant               !  number of antennas in service
      integer*4 tel_nbase              !  number of baselines
      integer*4 tel_nspac              !  number of spacings 
      integer*4 tel_nchan              !  number of channels

      equivalence (tel(1), tel_version)
      equivalence (tel(2), tel_array)
      equivalence (tel(3), tel_along)
      equivalence (tel(5), tel_glat)
      equivalence (tel(7), tel_height)
      equivalence (tel(9), tel_ax)
      equivalence (tel(9+2*MAXANT), tel_ay)
      equivalence (tel(9+4*MAXANT), tel_az)
      equivalence (tel(9+6*MAXANT), tel_ha_axis)
      equivalence (tel(9+10*MAXANT), tel_dec_axis)
      equivalence (tel(9+14*MAXANT), tel_freq)
      equivalence (tel(10+14*MAXANT), tel_freqIF)
      equivalence (tel(10+15*MAXANT), tel_tsys)
      equivalence (tel(10+16*MAXANT), tel_rain)
      equivalence (tel(10+17*MAXANT), tel_gain)
      equivalence (tel(10+17*MAXANT+MAXSP), tel_phase)
      equivalence (tel(10+17*MAXANT+2*MAXSP), tel_ha_zero)
      equivalence (tel(10+18*MAXANT+2*MAXSP), tel_dec_zero)
      equivalence (tel(10+19*MAXANT+2*MAXSP), tel_pc_zero)
      equivalence (tel(10+20*MAXANT+2*MAXSP), tel_ant)
      equivalence (tel(10+21*MAXANT+2*MAXSP), tel_corr)
      equivalence (tel(10+21*MAXANT+3*MAXSP), tel_gainsw)
      equivalence (tel(11+21*MAXANT+3*MAXSP), tel_nant)
      equivalence (tel(12+21*MAXANT+3*MAXSP), tel_nbase)
      equivalence (tel(13+21*MAXANT+3*MAXSP), tel_nspac)
      equivalence (tel(14+21*MAXANT+3*MAXSP), tel_nchan)

      integer*4  file_lun              !  current open file logical unit

      common /hdr_struct/ hdr, file_lun

*  Index arrays for telescope configuration

      integer*2  ant_list(MAXANT)      !  list of active antennas
      integer*2  ant_index(MAXANT)     !  index of active antennas
      integer*2  tel_base(2,MAXBASE)   !  antenna pairs for each baseline
      integer    tel_nsp               !  number of spacings

      common /tel_config/ ant_list, ant_index, tel_base, tel_nsp

*  Path compensation delay lookup table

      real*4 delay_table(MAXPC,MAXANT) !  pc delay table (main array)

      common /pc_lookup/ delay_table

