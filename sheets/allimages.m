InitY = 1169000;
InitX = 807200;
ChangeY = 1000;

[status, output] = system('ls zips/*.zip');
zips = strsplit(output, "\n");
for i = zips(1,1:size(zips,2)-1)
	
	name = strsplit(i{1},".");
	name = name(1);
	name = name{1};
	comm = sprintf('unzip -o %s', i{1});
	system(comm);
	source sheet.m
	n = min(mod(find(abs(FIELD)>0.01)-1,size(FIELD,1))+1);
	n = floor(n/2)-20
	name(1:4) = 'imgs';
	name(end) = 'x';
	createimage(name, FIELD(:,1:2:end),SAMP,n);
	name(end) = 'y';
	createimage(name, FIELD(:,2:2:end),SAMP,n);

endfor

