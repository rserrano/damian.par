InitY = 1169000;
InitX = 807200;
ChangeY = 1000;

set(gca, 'FontSize', 24);
[status, output] = system('ls zips/*.zip');
zips = strsplit(output, "\n");
for i = zips(1,1:size(zips,2)-1)
	
	name = strsplit(i{1},".");
	name = name(1);
	name = name{1};
	comm = sprintf('unzip -o %s', i{1});
	system(comm);
	source sheet.m
	name(1:4) = 'imgs';
	name(end) = 'f';
	name = [name 'x'];
	createtransfer(name, FIELD(:,1:2:end),SAMP,POINTS(1:2:end));
	name(end) = 'y';
	createtransfer(name, FIELD(:,2:2:end),SAMP,POINTS(1:2:end));
endfor

