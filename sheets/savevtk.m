unction [] = savevtk(name, MESH, POINTS, FIELD)
	MESH = MESH.-1;
	file = fopen(name, 'w');


	fprintf(file, '# vtk DataFile Version 2.0\n');
	fprintf(file, 'example\n');
	fprintf(file, 'ASCII\n');
	fprintf(file, 'DATASET UNSTRUCTURED_GRID\n');
	fprintf(file, 'POINTS %d float\n', size(POINTS,2));
	for i = 1:size(POINTS,2)
		fprintf(file, '%f %f %f\n', POINTS(1,i), POINTS(2,i), 0.0);
	endfor
	fprintf(file, 'CELLS %d %d\n', size(MESH,1), size(MESH,1)*5);
	for i = 1:size(MESH,1)
		fprintf(file, '4 %d %d %d %d\n', MESH(i,1), MESH(i,2), MESH(i,3), MESH(i,4));
	endfor
	fprintf(file, 'CELL_TYPES %d\n', size(MESH,1));
	fprintf(file, '%d ', ones(1,size(MESH,1)).*7);
	fprintf(file, '\n');
	fprintf(file, 'POINT_DATA %d\n', size(POINTS,2));
	fprintf(file, 'VECTORS U float\n');
	for i = 1:size(POINTS,2)
		fprintf(file, '%f %f %f\n', FIELD(1,i), FIELD(2,i), FIELD(3,i));
	endfor
	fclose(file);
