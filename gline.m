function h = gline(fig)
% gline - Interactively draw milti-point line
%   gline(FIG) draws a line by left - clicking the mouse at the control
%   points of the spline in the figure FIG. Right clicking will end the
%   line.

%   H = gline(FIG) Returns the handle to the line.
%
%   gline with no input arguments draws in the current figure.


if nargin<1, 
  draw_fig = gcf;
  fig = draw_fig;
  s = 'start'; 
end

if isstr(fig), 
   s = fig;
   draw_fig = gcbf;
   ud = get(draw_fig,'UserData');
else
   s = 'start';
   draw_fig = fig; 
end

ax = get(draw_fig,'CurrentAxes');
if isempty(ax)
   ax = axes('Parent',draw_fig);
end

gcacolor = get(ax,'Color');

% Check to see if we're clicking to exit
mouse = get(draw_fig,'selectiontype');
% Normal:    Click left mouse button.
% Extend:    Shift - click left mouse button or click both left and right 
%            mouse buttons,  or click middle mouse button. 
% Alternate: Control - click left mouse button or click right mouse button.
% Open:      Double-click any mouse button.
if ~strcmp('start',s)
    switch lower(mouse(1:3))
        case 'nor'
        case 'ext'
            s = 'last';
        case 'alt'
            s = 'last';
        case 'ope'
            s = 'last';
    end
end

switch s
   case 'start'
       if ~isempty(get(gcf,'windowbuttonupfcn')); disp('gline:turn off zoom, then try again-- exiting'); return; end
   oldtag = get(draw_fig,'Tag');
   figure(draw_fig);
   if any(get(ax,'view')~=[0 90]), 
     set(draw_fig, 'Tag', oldtag);
     error('stats:gline:NotTwoDimensional','gline works only for 2-D plots.');
   end
   
   
   % Initialize line
   xlimits = get(ax,'Xlim');
   x = (xlimits + flipud(xlimits))./2;
   ylimits = get(ax,'Ylim');
   y = (ylimits + flipud(ylimits))./2;
   hline = line(x,y,'Parent',ax,'Visible','off','eraseMode','normal');

   % Save current window functions and data
   bdown = get(draw_fig,'WindowButtonDownFcn');
   bup = get(draw_fig,'WindowButtonUpFcn');
   bmotion = get(draw_fig,'WindowButtonMotionFcn');
   oldud = get(draw_fig,'UserData');
   
   % Create new window functions
   set(draw_fig,'WindowButtonDownFcn','gline(''first'')')
   set(draw_fig,'WindowButtonMotionFcn','gline(''motion'')')
   set(draw_fig,'WindowButtonupFcn','')
   
   set(draw_fig,'doublebuffer','on')
 
   % Save drawing data as in 'UserData'
   ud.hline = hline;
   ud.pts = [];
   ud.buttonfcn = {bdown; bup; bmotion};
   ud.oldud = oldud;
   ud.oldtag = oldtag;
   ud.npts = 1;
   ud.xlimits = xlimits;
   ud.ylimits = ylimits;
   set(draw_fig,'UserData',ud);
   
   if nargout == 1
      h = hline;
   end

case 'motion'
   set(draw_fig,'Pointer','crosshair');

   if isempty(ud.pts);
      return;
   end

   [xspline,yspline] = deal(ud.pts(:,1),ud.pts(:,2));
   
   set(ud.hline,'Xdata',xspline,'Ydata',yspline, ...
        'linestyle','-', 'marker','.','linewidth', 1.5, 'Color',1-gcacolor,'Visible','on');
   Pt2 = get(ax,'CurrentPoint'); 
   Pt2 = Pt2(1,1:2);    
   ud.pts(ud.npts+1,:) = Pt2;

   set(draw_fig,'UserData',ud);
   L=sum(abs(diff(ud.pts*[1; i])));
   if isappdata(ud.hline,'h_L_text'); 
       h_L_text=getappdata(ud.hline,'h_L_text'); 
       if ishandle(h_L_text);
           set(h_L_text,'string', ['L=', num2str(L)]);
       end
   end
   
case 'first'   
   Pt1 = get(ax,'CurrentPoint'); 
   ud.pts = [Pt1(1,1:2); Pt1(1,1:2)];
   win_pos=get(gcf,'position');
   h_win=msgbox('L=0','current line length');
   set(h_win,'position', [win_pos(1)+win_pos(3), win_pos(2)+.75*win_pos(4), .25*win_pos(3) .25*win_pos(4)]);
   h_L_text=findobj(h_win,'type','text');
   setappdata(ud.hline,'h_L_text',h_L_text);
   set(draw_fig,'WindowButtonDownFcn','gline(''down'')','UserData',ud)
     
case 'down'         
   Pt1 = get(ax,'CurrentPoint'); 
   ud.pts = [ud.pts; Pt1(1,1:2)]; 
   ud.npts = ud.npts + 1;   
   [xspline,yspline] = deal(ud.pts(:,1),ud.pts(:,2));
%    
%    set(ud.hline,'Xdata',xspline,'Ydata',yspline,'eraseMode','normal', ...
%         'linestyle','-', 'linewidth', 1.5, 'Color',1-gcacolor,'Visible','on');

   set(draw_fig,'WindowButtonDownFcn','gline(''down'')','UserData',ud)
      
case 'last'
   bfcns = ud.buttonfcn;
   set(draw_fig,'windowbuttondownfcn',bfcns{1},'windowbuttonupfcn',bfcns{2}, ...
         'windowbuttonmotionfcn',bfcns{3},'Pointer','arrow', ...
		 'Tag',ud.oldtag,'UserData',ud.oldud)
   set(ud.hline,'UserData',ud.pts)
   modify_path_points(ud.hline);
    otherwise
   error('stats:gline:BadCallBack','Invalid call-back.');
end

% Make sure the axis limits don't change
set(ax,'xlim',ud.xlimits,'ylim',ud.ylimits);


 

 