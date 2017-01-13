! defines graphics options
subroutine graphics_ops

   use kind_def
   use sz_globals

   implicit none

   logical :: io_yesno
   external io_yesno

   graph_adj = io_yesno &
        ('Allow graphics displays to be adjustable (y/n)','no',status)
   autoscale = io_yesno &
        ('Autoscale profiles instead of fixed axes (a/f)','yes',status)
   if (.not.autoscale) then
      call io_getd('X-axis start','*',prof_x1,status)
      call io_getd('X-axis finish','*',prof_x2,status)
      call io_getd('Y-axis start','*',prof_y1,status)
      call io_getd('Y-axis finish','*',prof_y2,status)
   end if

   if (io_yesno('Inverse video?','no',status)) then
      call pgscr(0,1.0,1.0,1.0)
      call pgscr(1,0.0,0.0,0.0)
   else
      call pgscr(0,0.0,0.0,0.0)
      call pgscr(1,1.0,1.0,1.0)
   end if
     
   write(*,*) 'Index   Colours'
   write(*,*)
   write(*,*) '  0      Black'
   write(*,*) '  1      White'
   write(*,*) '  2      Red'
   write(*,*) '  3      Green'
   write(*,*) '  4      Blue'
   write(*,*) '  5      Cyan'
   write(*,*) '  6      Magenta'
   write(*,*) '  7      Yellow'
   write(*,*) '  8      Orange'
   write(*,*) '  9      Light Green'
   write(*,*) '  10     Dark Green'
   write(*,*) '  11     Light Blue'
   write(*,*) '  12     Purple'
   write(*,*) '  13     Dark Red'
   write(*,*) '  14     Dark Grey'
   write(*,*) '  15     Light Grey'

   call io_geti('plot colour?','*',plot_col,status)

   call pgsci(plot_col)

   write(*,*)
   write(*,*) 'Index   Line Style'
   write(*,*)
   write(*,*) '  1     Full line'
   write(*,*) '  2     Dashed'
   write(*,*) '  3     Dot-Dash'
   write(*,*) '  4     Dotted'
   write(*,*) '  5     Dash-Dot-Dot-Dot'

   call io_geti('line style?','*',line_sty,status)

   call pgsls(line_sty)

   call io_geti('line width?','*',line_width,status)

   call pgslw(line_width)

   call io_getr('Character height?','*',char_hgt,status)

   call pgsch(char_hgt)

end subroutine graphics_ops



