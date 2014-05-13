function [Xloc,YF1] = obtainsingle(FIELD,SAMP,POINTS,a)
	Fs = 1/SAMP;
	L = size(FIELD,1);
	NFFT = 2^nextpow2(L);
	f = (Fs/2)*linspace(0,1,NFFT/2+1);
	
	
	% threat them like sheets.
	yax = f;
	xax = POINTS;
	norm = f;

	m = 0;
	
	
	n = find(norm>=1,1,'first');
	
	
	Y = fft(FIELD,NFFT);
	Yabs = (SAMP)*(2*abs(Y(1:NFFT/2+1,:)));
	Xloc = xax;
	YF1   = Yabs(n,:)/a;

