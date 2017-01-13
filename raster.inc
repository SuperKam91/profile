*  AMI Project - raster.inc
*
*  This file contains lookup tables for centred hex-raster offsets.
*
*  History:
*     6/12/07 - original version [DJT/NHW]
*

*  Lookup table identifies offset index for each pointing.

c Default
      integer ptdf_list(100)
c Hexagons
      integer pt7_list(7)
      integer pt19_list(19)
c Rings
      integer pt18_list(18)
      integer pt24_list(24)
c Fat Rings - not worked this one out yet
c      integer 43pt_list43(43)

      data ptdf_list / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
     : 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     : 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
     : 31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
     : 41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
     : 51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
     : 61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
     : 71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
     : 81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
     : 91, 92, 93, 94, 95, 96, 97, 98, 99, 100 /
      data pt7_list / 6, 7, 5, 1, 2, 4, 3 /
      data pt19_list / 16, 17, 18, 15, 6, 7, 19, 14, 
     : 5, 1, 2, 8, 13, 4, 3, 9, 12, 11, 10 /
      data pt18_list / 32, 33, 34, 35, 36, 37, 
     : 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31 /
      data pt24_list / 54, 55, 56, 57, 58, 59, 60, 61, 
     : 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49,
     : 50, 51, 52, 53 /

