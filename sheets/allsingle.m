InitY = 1169000;
InitX = 807200;
ChangeY = 1000;

set(gca, 'FontSize', 24);

[status, output] = system('ls -t zips/*.zip');
zips = strsplit(output, "\n");
for i = zips(1,1:size(zips,2)-1)
	
	name = strsplit(i{1},".");
	name = name(1);
	name = name{1};
	comm = sprintf('unzip -o %s', i{1});
	system(comm);

	model = name(6:12);
	
	if ( model(7) == '1' )
		OC = 4;
	else
		OC = 3;
	end
	OO = 4;

	[a] = transmex(1,1);
	
	source sheet.m
	disp('EXFEM')
	POINTS(1:2:end) = POINTS(1:2:end)-20000;
	POINTS = POINTS/1000;
	[ns, nf] = rangetime(OC,SAMP,size(FIELD,1),-2.5,10);
	ns
	nf

	FIELD=FIELD(ns:nf,:);
	
	UCX = FIELD(:,1:2:end);
	UCY = FIELD(:,2:2:end);
	
	disp(SAMP)
	
	[XC,FCX] = obtainsingle(UCX,SAMP,POINTS(1:2:end),a);
	[XC,FCY] = obtainsingle(UCY,SAMP,POINTS(1:2:end),a);

	bemsh = strcat('bench/',model,'.m');
	source(bemsh);
	if ( model(1) == 'v' )
		[POINTS, FIELD] = filterhappier(POINTS, FIELD);
	end
	[ns, nf] = rangetime(OO,SAMP,size(FIELD,1),-2.5,10);
	ns
	nf

	FIELD=FIELD(ns:nf,:);	
	UOX = FIELD(:,1:2:end);
	UOY = FIELD(:,2:2:end);
	disp('BEM')
	disp(SAMP)
	[XO,FOX] = obtainsingle(UOX,SAMP,POINTS(1:2:end),a);
	[XO,FOY] = obtainsingle(UOY,SAMP,POINTS(1:2:end),a);

	name(1:4) = 'imgs';
	name(end) = 's';
	name = [name 'x'];
	createsingle(name,XC,FCX,XO,FOX);
	name(end) = 'y';
	createsingle(name,XC,FCY,XO,FOY);
endfor

