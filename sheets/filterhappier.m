function [POINTS,FIELD] = filterhappier(POINTS,FIELD)

del = [];
for ind = 2:2:size(POINTS,2)
	if ( abs(POINTS(ind)) > 0.00001 ) 
		del = [del, ind-1];
		del = [del, ind];
	end
end

POINTS(del) = [];
FIELD(:,del) = [];
[B,I] = sort(POINTS(1:2:end));
IE = zeros(1,2*size(I,2));
IE(1:2:end) = 2*I-1;
IE(2:2:end) = 2*I;

POINTS(:) = POINTS(IE);
FIELD(:,:) = FIELD(:,IE);





