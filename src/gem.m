function GeM(in_filename, out_filename, c = 0.85)
	fid = fopen(in_filename, 'r');
	[equipos, partidos] = fscanf(fid, '%u %u' ,"C");

	# Genero la matriz A
	A = zeros(equipos);
	while (partidos > 0)
		[f, e0, p0, e1, p1] = fscanf(fid, '%u %u %u %u %u', "C");

		if(p0 < p1)
			# si e0 perdió contra e1
			A(e0, e1) += p1 - p0;
		elseif(p0 > p1)
			# si e1 perdió contra e0
			A(e1, e0) += p0 - p1;
		endif

		partidos--;
	endwhile
	fclose(fid);

	# Genero la matriz H
	H = zeros(equipos);
	a = zeros(equipos, 1);
	for i = 1:equipos
		suma = sum(A(i, :));
		if(suma > 0)
			H(i, :) = A(i, :)/sum(A(i, :));
		else
			# vector para equipos invictos
			a(i) = 1;
		endif
	endfor

	# Genero la matriz G
	e = ones(equipos, 1);

	# Genero vectores de personalización
	# (por defecto distribución uniforme)

	# equipos invictos
	u = ones(equipos, 1)/equipos;
	# teletransportación
	v = u;

	G = c*[H + a*u'] + (1 - c)*e*v';

	# Calculo autovalores y autovectores
	[V, l] = eig(G');

	# Busco el autovalor 1
	i = 1;
	while(abs(l(i, i) - 1) > 0.0001)
		i++;
	endwhile
	x = V(:, i);

	# Normalizo el vector solución
	x = abs(x)/sum(abs(x));

	# Ordeno las soluciones
	S = zeros(equipos, 2);
	for i = 1:equipos
		S(i, 1) = i;
		S(i, 2) = x(i);
	endfor

	S = sortrows(S, 2);

	# Escribo la solución
	fid = fopen(out_filename, 'w');

	for i = 0:equipos - 1
		fprintf(fid, '%u %f\n',S(equipos - i,1) ,S(equipos - i, 2));
	endfor

	fclose(fid);
endfunction
