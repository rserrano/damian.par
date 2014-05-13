function [r] = transmex(fpeak,a)
	t = linspace(-5,5,1001);
	mexi = (1-2*(pi*fpeak*t).^2).*exp(-(pi*fpeak*t).^2);
	NFFT = 2^nextpow2(size(t,2));
	tran = fft(mexi,NFFT);
	tabs = (1/100)*(2*abs(tran(1:NFFT/2+1)));
	f = 100/2*linspace(0,1,NFFT/2+1);
	norm = f;
	n = find(norm>=a,1,'first');
	r = tabs(n);

