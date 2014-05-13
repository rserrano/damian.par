InitY = 1169000;
InitX = 807200;
ChangeY = 1000;

[status, output] = system('ls -t zips/*.zip');
zips = strsplit(output, "\n");

for i = zips(1,1:size(zips,2)-1)
	
	name = strsplit(i{1},".");
	name = name(1);
	name = name{1};
	comm = sprintf('unzip -o %s', i{1});
	system(comm);
	source sheet.m
	UC = FIELD;
	XC = POINTS(1:2:end);
	XC = (XC-20000)/1000;
	YC = POINTS(2:2:end);
	YC = YC/1000;
	SC = SAMP;
	model = name(6:12);
	bemsh = strcat('bench/',model,'.m');
	source(bemsh);
	UO = FIELD;
	XO = POINTS(1:2:end);
	YO = POINTS(2:2:end);
	SO = SAMP;
	if ( model(7) == '1' )
		OC = 4;
	else
		OC = 3;
	end
	OO = 4;
	if ( model(1) == 'v' )
		OBS = [0.0 0.0 1.0 0.0 2.0 0.0];
	else
		OBS = [0.0 1.0 1.0 0.0 2.0 0.0];
	end

	name(1:4) = 'imgs';
	name(end) = 't';
	name = [name 'x'];
	createsintime(name,UC(:,1:2:end),UO(:,1:2:end),XC,YC,XO,YO,SC,SO,OC,OO,OBS(1:2:end),OBS(2:2:end));
	name(end) = 'y';
	createsintime(name,UC(:,2:2:end),UO(:,2:2:end),XC,YC,XO,YO,SC,SO,OC,OO,OBS(1:2:end),OBS(2:2:end));
endfor

