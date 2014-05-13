InitY = 1169000;
InitX = 807200;
ChangeY = 1000;
canyonsvx
canyonsvy
canyonpx
canyonpy
[status, output] = system('ls -t files/canyonsv1.zip');
zips = strsplit(output, "\n");
for i = zips(1,1:size(zips,2)-1)
	
	name = strsplit(i{1},".");
	name = name(1);
	name = name{1};
	comm = sprintf('unzip -o %s', i{1});
	system(comm);
	source sheet.m
	[a1,a2,a3] = transmex(1);
	name(end) = 's';
	name = [name 'x'];
	comparitrans(name, FIELD(:,1:2:end),SAMP,POINTS(1:2:end),a1,a2,a3,SVX);
	name(end) = 'y';
	comparitrans(name, FIELD(:,2:2:end),SAMP,POINTS(1:2:end),a1,a2,a3,SVY);
endfor
[status, output] = system('ls -t files/canyonp1.zip');
zips = strsplit(output, "\n");
for i = zips(1,1:size(zips,2)-1)
	
	name = strsplit(i{1},".");
	name = name(1);
	name = name{1};
	comm = sprintf('unzip -o %s', i{1});
	system(comm);
	source sheet.m
	[a1,a2,a3] = transmex(1);
	name(end) = 's';
	name = [name 'x'];
	comparitrans(name, FIELD(:,1:2:end),SAMP,POINTS(1:2:end),a1,a2,a3,PX);
	name(end) = 'y';
	comparitrans(name, FIELD(:,2:2:end),SAMP,POINTS(1:2:end),a1,a2,a3,PY);
endfor

