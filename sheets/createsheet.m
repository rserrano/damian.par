function [] = createsheet(sheet,FIELD)

	Sp = 5;
	Mx = size(FIELD,2);
	Sc = max(max(abs(FIELD)));
	
	Dist = 5*Sc/(size(FIELD,2)/Sp);
	Sep = 0;
	n = floor(size(FIELD,1)/4);
	for i = 1:Sp:Mx-1
		y = FIELD(1:n,i)'.+Sep;
		plot(1:n,y,'k');
		hold on
		Sep = Sep + Dist;
	endfor
	axis tight
	axis off	
	png = sprintf('%s.png',sheet);
	pdf = sprintf('%s.pdf',sheet);
	

	print(pdf, '-dpdfwrite','-tight');
	print(png, '-dpng');
	hold off	
	
