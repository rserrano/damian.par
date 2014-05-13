function [] = createimage(sheet,FIELD,SAMP,n)
	yax = linspace(-4500,4500,size(FIELD,2));
	xax = (0:(size(FIELD,1)-1)).*SAMP;
	imagesc(xax(n:end), yax, FIELD(n:end,:)');
	colormap(flipud(gray));	
	
	pdf = sprintf('%s.pdf',sheet);
	

	print(pdf, '-dpdfwrite','-tight');
	hold off	
	
