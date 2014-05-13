function [rinit, rfinish] = rangetime(start, samp, n, init, finish)
	TC        = (linspace(0,n-1,n)*samp)-start;
	rfinish   = find(TC>=finish,1,'first');
	rinit     = find(TC<=init,1,'last');
TC(rinit)
TC(rfinish)

