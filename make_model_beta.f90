! Routine which produces a beta (or double beta) model n profile
! Moved from make_skies 28/08/09 [KG] 
subroutine make_model_beta(n)

   use kind_def
   use sz_globals
 
   implicit none

   integer :: i,j, i1, j1, which_i
   real(kind=dp), dimension(maxsize,maxsize) :: n
   real(kind=dp) :: co1, co2, pow, pow2, rad, ratio1, ratio2
   real(kind=dp) :: theta2, theta_cx2, theta_cx2b

   co2 = (cutoff2)**2
   co1 = (cutoff1)**2
   pow = -beta*1.5
   if (ellipsoid) then
      if (theta_c(1).gt.theta_c(2)) then 
         theta_cx = theta_c(1)
      else
         theta_cx = theta_c(2)
      end if
   end if

   theta_cx2 = theta_cx/cellsize
   theta_cx2 = theta_cx2**2

   theta_cx2b = theta_cx2a/cellsize
   theta_cx2b = theta_cx2b**2

  do i = 1, maxsize
     do j = 1, maxsize
        i1 = i-1
        j1 = j-1
        theta2 = float(i1**2+j1**2)

        select case(model_type)
! Double beta model
        case('d','D')
           ratio1 = theta2/theta_cx2
           ratio2 = theta2/theta_cx2b               
           pow2 = -beta2*1.5

           if ((ratio2 > co2)) then 
              n(j,i)=0.0
           else if (ratio2 < co1) then 
              n(j,i)=n0*((1+ratio1)**pow)+nx_2*((1+ratio2)**pow2)
           else 
              rad = sqrt(ratio2)
              if (hanning) then 
! Hanning Window
                 n(j,i) = (n0*((1.0+ratio2)**pow)+nx_2*((1+ratio2)**pow2)) * &
                      cos((rad-cutoff1)*pi/(2.d0*(cutoff2-cutoff1)))
              else
! Welch Window
                 n(j,i) = (n0*((1.0+ratio2)**pow)+nx_2*((1+ratio2)**pow2)) * &
                      (1-((rad-cutoff1)/(cutoff2-cutoff1)))**2
              end if
           end if

! Single beta model with cutoff
        case default 
           ratio1 = theta2/theta_cx2

           if (ratio1 > co2) then 
              n(j,i)=0.0
           else if (ratio1 < co1) then 
              n(j,i)=n0*((1+ratio1)**pow)
           else 
              rad = sqrt(ratio1)
              if (hanning) then
! Hanning Window
                 n(j,i) = n0*((1.0+ratio1)**pow) * &
                         cos((rad-cutoff1)*pi/(2.d0*(cutoff2-cutoff1)))
              else
! Welch Window
                 n(j,i) = n0*((1.0+ratio1)**pow) * &
                       (1-((rad-cutoff1)/(cutoff2-cutoff1)))**2     
              end if
           end if
        end select
     end do
  end do

end subroutine make_model_beta
