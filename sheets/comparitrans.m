function [] = comparitrans(sheet,FIELD,SAMP,POINTS,a1,a2,a3,KAWA)
	
	
	Fs = 1/SAMP;
	L = size(FIELD,1);
	NFFT = 2^nextpow2(L);
	f = Fs/2*linspace(0,1,NFFT/2+1);
	
	
	% threat them like sheets.
	yax = f;
	xax = (POINTS-6000)/1000;
	norm = f;

	n1 = find(norm>=1,1,'first');
	n2 = find(norm>=4,1,'first');
	n3 = find(norm>=16,1,'first');
	
	Y = fft(FIELD,NFFT);
	Yabs = 2*abs(Y(1:NFFT/2+1,:));
	plot(xax, Yabs(n1,:)/a1,'k',KAWA(:,1), KAWA(:,2), 'ko');


%	plot(xax, Yabs(n1,:),'r', xax, Yabs(n2,:), 'g' ,xax, Yabs(n3,:), 'b');
	pdf = sprintf('%s.pdf',sheet);
	print(pdf, '-dpdfwrite','-tight');
	hold off
	
