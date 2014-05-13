function [] = createsingle(sheet,XC,FC,XO,FO)
	clf;
	hold on;
	n1 = find(XO>=-2,1,'first');
	n2 = find(XO>=2,1,'first');
	set(gca, 'FontSize', 20);
	xlabel('X coordinate');
	ylabel('Displacements');
	plot(XC, FC,'k');
	plot(XO(n1:n2), FO(n1:n2),'ko');
	axis tight;
	pdf = sprintf('%s.pdf',sheet);
	print(pdf, '-dpdfwrite','-tight');
	hold off;
		
