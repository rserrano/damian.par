function [] = createtransfer(sheet,FIELD,SAMP,POINTS)
	Fs = 1/SAMP;
	L = size(FIELD,1);
	NFFT = 2^nextpow2(L);
	f = Fs/2*linspace(0,1,NFFT/2+1);

	% threat them like sheets.
	xax = POINTS(1:2:end);
	yax = f;
	Z = zeros(size(yax,2),size(xax,2));
	n = find(f>4,1,'first');
	Y = fft(FIELD,NFFT);
	Yabs = 2*abs(Y(1:NFFT/2+1,:));
	Yabs = Yabs(end:-1:1,:);
	yax  = yax(end:-1:1,:);
	imagesc((xax(1:end)-20000)/1000, yax(1:n), Yabs(1:n,1:end));
	colormap(flipud(gray));	
	ylabel('Frequency');
	xlabel('Surface Location');
	pdf = sprintf('%s.pdf',sheet);
	png = sprintf('%s.png',sheet);

	print(pdf, '-dpdfwrite','-tight');

	hold off
	
