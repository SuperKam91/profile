
! Calculates the analytic contribution to a visibility from a point source with or without bandwidth smearing.


subroutine calc_vis(uu,vv,ww,dra,ddec,dec,bwith,ReIm,smear,vis)

   use kind_def
   use sz_globals

   implicit none


   real(dp)  :: uu,vv,ww,dra,ddec,dec,bwith,vis
   real(dp)  :: vis1,vis2,u_nu,v_nu,w_nu
   real(dp)  :: nu_upper,nu_lower
   logical   :: ReIm, smear


! no smearing:
! real part:   
   if ((.not.smear).or.((ddec==0.).and.(dra==0.))) then
      if (ReIm) then
         vis = cos(2.*pi*(uu*dra&
              +vv*ddec*dec &
              -ww*(dra**2+(ddec*dec)**2)/2.))
! imaginary part:
      else
         vis = sin(2.*pi*(uu*dra&
              +vv*ddec*dec &
              -ww*(dra**2+(ddec*dec)**2)/2.))
      end if

! with smearing:
! real part:
   else
      if (ReIm) then
         
         nu_upper = nu(chan)+bwith/2.
         nu_lower = nu(chan)-bwith/2.
         u_nu = (uu/nu(chan))*nu_upper
         v_nu = (vv/nu(chan))*nu_upper
         w_nu = (ww/nu(chan))*nu_upper
         
         vis1 = sin(2.*pi*(u_nu*dra+v_nu*ddec*dec&
              -w_nu*(dra**2+(ddec*dec)**2)/2.))
         vis1 = vis1/(2.*pi*((uu/nu(chan))*dra+(vv/nu(chan))*ddec*dec&
              -(ww/nu(chan))*((dra**2+(ddec*dec)**2)/2.)))
         vis1 = vis1/bwith ! normalization
         
         u_nu = (uu/nu(chan))*nu_lower
         v_nu = (vv/nu(chan))*nu_lower
         w_nu = (ww/nu(chan))*nu_lower
 
         vis2 = sin(2.*pi*(u_nu*dra+v_nu*ddec*dec&
              -w_nu*(dra**2+(ddec*dec)**2)/2.))
         vis2 = vis2/(2.*pi*((uu/nu(chan))*dra+(vv/nu(chan))*ddec*dec&
              -(ww/nu(chan))*((dra**2+(ddec*dec)**2)/2.)))
         vis2 = vis2/bwith

         vis = vis1 - vis2
  
! imaginary part:         
      else

         nu_upper = nu(chan)+bwith/2.
         nu_lower = nu(chan)-bwith/2.
         u_nu = (uu/nu(chan))*nu_upper
         v_nu = (vv/nu(chan))*nu_upper
         w_nu = (ww/nu(chan))*nu_upper
 
         vis1 = -1.*cos(2.*pi*(u_nu*dra+v_nu*ddec*dec&
              -w_nu*(dra**2+(ddec*dec)**2)/2.))
         vis1 = vis1/(2.*pi*((uu/nu(chan))*dra+(vv/nu(chan))*ddec*dec&
              -(ww/nu(chan))*((dra**2+(ddec*dec)**2)/2.)))
         vis1 = vis1/bwith ! normalization
         
         u_nu = (uu/nu(chan))*nu_lower
         v_nu = (vv/nu(chan))*nu_lower
         w_nu = (ww/nu(chan))*nu_lower

         vis2 = -1.*cos(2.*pi*(u_nu*dra+v_nu*ddec*dec&
              -w_nu*(dra**2+(ddec*dec)**2)/2.))
         vis2 = vis2/(2.*pi*((uu/nu(chan))*dra+(vv/nu(chan))*ddec*dec&
              -(ww/nu(chan))*((dra**2+(ddec*dec)**2)/2.)))
         vis2 = vis2/bwith 

         vis = vis1 - vis2

      end if
   end if

    


end subroutine calc_vis
