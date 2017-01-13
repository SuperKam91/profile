! This module contains structures and global variables for likelihood
! calculations

module structures

   use kind_def

! The following are global variables
   integer :: num_param, num_dat_pt

! The following are global structures
   type, public :: data
      real(kind=dp), dimension(:), pointer :: values
   end type data

   type, public :: param
      real(kind=dp), dimension(:), pointer :: coeff
   end type param

   type, public :: model
      type(param), pointer :: m_param
      type(data), pointer :: m_data
      real(kind=dp) :: likely
   end type model

end module structures
