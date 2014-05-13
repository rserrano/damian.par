InitY = 1169000;
InitX = 807200;
ChangeY = 1000;

[status, output] = system('ls files/*.zip');
zips = strsplit(output, "\n");

for i = zips(1,1:size(zips,2)-1)
	
	name = strsplit(i{1},".");
	name = name(1);
	name = name{1};
	name = name(4:end-1);
	comm = sprintf('unzip -o %s', i{1});
	system(comm);
	source sheet.m

	magmax = [];
	xmax = [];
	ymax = [];
	for j = 1:2:size(FIELD,2)
		tmp = sqrt(FIELD(:,j).^2 + FIELD(:,j+1).^2);
		magmax = [magmax max(tmp)];
		xmax = [xmax max(abs(FIELD(:,j)))];
		ymax = [ymax max(abs(FIELD(:,j+1)))];
	endfor
	fname = sprintf('files/maxsim%s.m', name);
	file = fopen(fname, 'w');
	fprintf(file, 'NUM = %s;\n', name);
	Ycoord = InitY+ChangeY*str2num(name);
	PTS = [POINTS(1:2:end); POINTS(2:2:end)];
	fprintf(file, 'POINTS = [');
	fprintf(file, '%f ', PTS(1,1:30:end).+InitX);
	fprintf(file, '\n');
	Y = ones(1,size(PTS(:,1:30:end),2)).*Ycoord;
	fprintf(file, '%f ', Y);
	fprintf(file, '\n');
	fprintf(file, '%f ', 3100-PTS(2,1:30:end));
	fprintf(file, '];\n');
	fprintf(file, 'MAGMAX = [');
	fprintf(file, '%f ', magmax(1:30:end));
	fprintf(file, '];\n');
	fprintf(file, 'XMAX = [');
	fprintf(file, '%f ', xmax(1:30:end));
	fprintf(file, '];\n');
	fprintf(file, 'YMAX = [');
	fprintf(file, '%f ', ymax(1:30:end));
	fprintf(file, '];\n');
	fclose(file);
endfor

