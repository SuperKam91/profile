! de-grids a visibility at specified u,v from an unconvolved gridded aperture
subroutine extract_visibility(re_map,im_map,mapsize,upt,vpt,cell,re_pt,im_pt)

   use kind_def
   use sz_globals
   use maths

   implicit none
   integer, intent(in) :: mapsize
   real(kind=dp), dimension(mapsize,mapsize), intent(in) :: re_map, im_map
   real(kind=dp), intent(in) :: upt,vpt,cell
   real(kind=dp), intent(out) :: re_pt,im_pt
   real(kind=dp) :: centre,x,y,uvcell

   uvcell = 1/(mapsize*cell*sec2rad)
   centre = (mapsize/2)+1
   x = (upt/uvcell)+centre
   y = (vpt/uvcell)+centre

   if ((x.lt.0.).or.(y.lt.0.).or.(x.gt.mapsize).or.(y.gt.mapsize)) then
      re_pt = 0.
      im_pt = 0.      
   else
      call interp(re_map,mapsize,mapsize,x,y,re_pt)
      call interp(im_map,mapsize,mapsize,x,y,im_pt)
   end if

end subroutine extract_visibility

