InitY = 1169000;
InitX = 807200;
ChangeY = 1000;


set(gca, 'Position', get(gca, 'OuterPosition') - ...
get(gca, 'TightInset') * [-1 0 1 0; 0 -1 0 1; 0 0 1 0; 0 0 0 1]);
set(gca,'xtick',[],'ytick',[]);

[status, output] = system('ls -t zips/*.zip');
zips = strsplit(output, "\n");
for i = zips(1,1:size(zips,2)-1)
	
	name = strsplit(i{1},".");
	name = name(1);
	name = name{1};
	comm = sprintf('unzip -o %s', i{1});
	system(comm);
	source sheet.m
	name(1:4) = 'imgs';
	name(end) = 'x';
	createsheet(name, FIELD(:,1:2:end));
	name(end) = 'y';
	createsheet(name, FIELD(:,2:2:end));


endfor

