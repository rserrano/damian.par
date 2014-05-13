ALLPOINTS = [];
ALLMAGMAX = [];
ALLXMAX = [];
ALLYMAX = [];

surfs = strsplit(output, "\n");
mesh  = [];

first = 1;
Nrow = 0;
for i = 1:35
	str = sprintf('maxsim%d.m', i);
 	a = fopen(str, 'r');
	if ( a == -1 )
		continue
	else
		fclose(a);
	endif
	source(str)
	
	if ( first == 1 ) 
		first = 0;
	else
		cols = size(POINTS,2);
		for j = 1:cols-1
			Fn = j+Nrow;
			mesh = [mesh; Fn (Fn+1) (Fn+cols+1) (Fn+cols)];
		endfor
		Nrow = Nrow + cols;
	endif
	
	ALLPOINTS = [ALLPOINTS POINTS];
	ALLMAGMAX = [ALLMAGMAX MAGMAX];
	ALLXMAX   = [ALLXMAX XMAX];
	ALLYMAX   = [ALLYMAX YMAX];
endfor

savevtk('output.vtk', mesh, ALLPOINTS, [ALLMAGMAX; ALLXMAX; ALLYMAX]);



% trisurf(TRI(1:1,:),ALLPOINTS(1,1:1:end), ALLPOINTS(2,1:1:end), ALLPOINTS(3,1:1:end),ALLXMAX(1:1:end),'FaceColor','interp');
% view(2)
% print -dpng -color xmax.png
% print -dpdf -color xmax.pdf


% trisurf(TRI,ALLPOINTS(1,1:20:end), ALLPOINTS(2,1:20:end), ALLPOINTS(3,1:20:end),ALLYMAX(1:20:end),'EdgeAlpha',0,'FaceColor','interp');
% view(2)
% print -dpng -color ymax.png
% print -dpdf -color ymax.pdf


% trisurf(TRI,ALLPOINTS(1,1:20:end), ALLPOINTS(2,1:20:end), ALLPOINTS(3,1:20:end),ALLMAGMAX(1:20:end), 'EdgeAlpha',0,'FaceColor','interp');
% view(2)
% print -dpng -color magmax.png
% print -dpdf -color magmax.pdf

% trisurf(TRI,ALLPOINTS(1,1:20:end), ALLPOINTS(2,1:20:end), ALLPOINTS(3,1:20:end),ALLPOINTS(3,1:20:end), 'EdgeAlpha',0,'FaceColor','interp');
% view(2)
% print -dpng -color height.png
% print -dpdf -color height.pdf
% 
