module kind_def

  integer, parameter :: dp = kind(1.0D0)

  type array_pointer
     real (kind=dp), dimension(:,:), pointer :: pt
  end type array_pointer

end module kind_def
