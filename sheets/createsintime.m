function [] = createsintime(sheet,UC,UO,XC,YC,XO,YO,SC,SO,OC,OO,XOBS,YOBS)
	oobsInds = [];
	for i=1:size(XOBS,2)
		[S,I] = min((XO.-XOBS(i)).^2+(YO.-YOBS(i)).^2);
		oobsInds = [oobsInds, I];
	end
	
	obsInds = [];
	for i=oobsInds
		[S,I] = min((XC.-XO(i)).^2+(YC.-YO(i)).^2);
		obsInds = [obsInds, I];
	end

	% threat them like sheets.
	clf
	hold on;
	TC = linspace(0,size(UC,1)-1,size(UC,1)).*SC-OC;
	TO = linspace(0,size(UO,1)-1,size(UO,1)).*SO-OO;
	nc = find(TC>=20,1,'first');
	no = find(TO>=20,1,'first');
	
	set(gca, 'FontSize', 20);


	ylabel('Displacement');
	xlabel('Time');


	UO = -UO;
	for i = 1:size(oobsInds,2)
		UC(:,obsInds(i))  = UC(:,obsInds(i)).+XOBS(i);
		UO(:,oobsInds(i)) = UO(:,oobsInds(i)).+XOBS(i);
		plot(TC(1:nc),UC(1:nc,obsInds(i)),'k');
		plot(TO(1:2:no),UO(1:2:no,oobsInds(i)),'ko','MarkerSize',10);
	end
	pdf = sprintf('%s.pdf',sheet);
	print(pdf, '-dpdfwrite','-tight');
	hold off;
	
