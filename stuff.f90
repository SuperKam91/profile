! This module contains various things to start the poisson program going

module stuff

   use kind_def
   use structures
   use maths
   use sz_globals

   public :: create_param
   public :: create_data
   public :: set_model_param
   public :: create_model
   public :: calculate_likely
   public :: create_and_calculate
   public :: setup_simplex

contains

!**************************************************************************
! Sets up pointer to array of paramaters
   subroutine create_param(n,p,par)

      integer, intent(in) :: n
      real(kind=dp), intent(in), dimension(n), target :: p
      type(param) :: par

      par%coeff=>p
      num_param=n	

   end subroutine create_param

!**************************************************************************
! Sets up pointer to array of data points
   subroutine create_data(x,y,rect,dat)

      integer, intent(in) ::  x, y
      real(kind=dp), intent(in), dimension(x,y):: rect
      type(data), intent(out) :: dat
      integer :: allo_stat

      deallocate(dat%values,STAT=allo_stat)
      if (allo_stat.ne.0) then 
         write(*,*) 'Error deallocating in create data'
      else
         write(*,*) 'No problem deallocating in create data'
      end if
      allocate(dat%values(num_dat_pt),STAT=allo_stat)
      if (allo_stat.ne.0) then 
         write(*,*) 'Error allocating in create data'
      else
         write(*,*) 'No problem allocating in create data'
      end if
      dat%values=pack(rect,mask=ok_data) 

   end subroutine create_data

!**************************************************************************
! Sets the pointer to the parameters
   subroutine set_model_param(params,par_model)

      type(param), intent(in), target :: params
      type(model), intent(out) :: par_model

      par_model%m_param=>params

   end subroutine set_model_param

!**************************************************************************
! 
   subroutine create_model(par_model)

      type(model), intent(inout) :: par_model
      integer :: i,j
      real(kind=dp), dimension(maxsize,maxsize) :: model_values
      integer :: allo_stat

! Check that par_model does not already exist and deallocate if it does

      deallocate(par_model%m_data%values,stat=allo_stat)
      if (allo_stat.ne.0) then 
         write(*,*) 'Error deallocating in create_model'
      else
         write(*,*) 'No problem deallocating in create_model'
      end if

      allocate(par_model%m_data%values(num_dat_pt))
      if (allo_stat.ne.0) then 
         write(*,*) 'Error allocating in create_model'
      else
         write(*,*) 'No problem allocating in create_model'
      end if
! Calculate theoretical values according to the parameters and pack the
! observed ones into par_model
      call model_fn(par_model%m_param%coeff,num_param,model_values)
write(*,*) 'problem?'

      par_model%m_data%values = pack(model_values,mask=ok_data)

   end subroutine create_model

!**************************************************************************
! Calculate a log(probability) for each data point given a model assuming 
! Poisson statistics then sum these to get a likelihood

! *** HACKED 24/03/04 AMS ***

   subroutine calculate_likely(given,calculated)

      type(data), intent(in) :: given
!      real(kind=dp), intent(out) :: prob_sum
      type(model), intent(inout) :: calculated
      real(kind=dp), dimension(num_dat_pt) :: prob
      integer :: i
      real(kind=dp) :: dummy

! Poisson statistics: P(x=r) = exp(-lambda) lambda^r / r! 
! Find the log of the above
! Since given%values is actual data then it should have be an integer
      do i=1,num_dat_pt

! Don't pay too much attention to points where the model predicts virtually 
! no signal. NOT TOO SURE THAT THIS IS A VALID THING TO DO. 
         if(calculated%m_data%values(i).le.0.1)then
            prob(i)=0.0
         else
!        prob(i) = given%values(i)*log(calculated%m_data%values(i))&
!             & -1.0*calculated%m_data%values(i)-&
!             & log(1.0*factorial(given%values(i)))
            dummy=1.0*nint(given%values(i))
            prob(i) = given%values(i)*log(calculated%m_data%values(i))&
                 & -1.0*calculated%m_data%values(i)-&
                 & log_factorial(dummy)
         end if
      end do
      calculated%likely=sum(prob)
! added for new likelihood plot - AMS
!      prob_sum = sum(prob)

   end subroutine calculate_likely

!**************************************************************************
! Given some data and a model will calculate the likelihhod of the model given
! the data
   function create_and_calculate(data_stuff,model_stuff)result(likely)

      type(data), intent(in) :: data_stuff
      type(model), intent(inout) :: model_stuff
   
      real(kind=dp) :: likely

      call create_model(model_stuff)

      call calculate_likely(data_stuff,model_stuff)	
      likely=model_stuff%likely

   end function create_and_calculate

!**************************************************************************
! Sets up a simplex required for the amoeba minimiser
   subroutine setup_simplex(data_stuff, model_stuff,p,y)

      type(data), intent(in) :: data_stuff
      type(model), intent(inout) :: model_stuff
      real(kind=dp), dimension(num_param+1, num_param),&
           intent(out) :: p
      real(kind=dp), dimension(num_param+1), intent(out) :: y

      real(kind=dp), dimension(num_param) :: lambda
      
      integer :: i,j

      select case(num_param)

! spherical case with beta constant
      case(4)
         lambda(1)=1.0e-3
         lambda(2)=10.0d0
         lambda(3)=10.0d0
         lambda(4)=10.0d0 

! spherical case with beta variable
      case(5)
         lambda(1)=1.0e-3
         lambda(2)=0.1d0
         lambda(3)=10.0d0
         lambda(4)=10.0d0
         lambda(5)=10.0d0

! elliptical case with beta constant
      case(6)
         lambda(1)=1.0e-3
         lambda(2)=10.0d0
         lambda(3)=10.0d0
         lambda(4)=10.0d0
         lambda(5)=10.0d0
         lambda(6)=10.0d0

! Elliptical case with beta variable
      case(7)
         lambda(1)=1.0e-3
         lambda(2)=0.1d0
         lambda(3)=10.0d0
         lambda(4)=10.0d0
         lambda(5)=10.0d0
         lambda(6)=10.0d0
         lambda(7)=10.0d0
      end select

      do i=1,num_param
         p(1,i)=model_stuff%m_param%coeff(i)
      end do
      do i=2, num_param+1
         do j=1, num_param
            p(i,j)=p(1,j)
            if(j.eq.i-1)p(i,j)=p(i,j)+lambda(j)
         end do
      end do

! DEBUG
!     write(*,*)'the initialised simplex list is :'
!      write(*,*) p !  case(1)
      do i=1,num_param+1
         model_stuff%m_param%coeff=p(i,:)
         y(i)=-create_and_calculate(data_stuff,model_stuff)
      end do
!DEBUG
!      write(*,*)'the corresponding likelyhoods are :'
!      write(*,*) y
      model_stuff%m_param%coeff=p(1,:)

   end subroutine setup_simplex

!**************************************************************************

end module stuff




