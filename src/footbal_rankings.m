1;

function GeM(in_filename, out_filename, team_codes_filename = '', c = 0.85, date_limit = 0, pres = 0.0001)
	has_date_limit = date_limit != 0;

	fid = fopen(in_filename, 'r');
	[equipos, partidos] = fscanf(fid, '%u %u' ,"C");

	# Genero la matriz A
	A = zeros(equipos);
	while (partidos > 0)
		[f, e0, p0, e1, p1] = fscanf(fid, '%u %u %u %u %u', "C");

		if(has_date_limit && f > date_limit)
			partidos = 0;
		else
			if(p0 < p1)
				# si e0 perdió contra e1
				A(e0, e1) += p1 - p0;
			elseif(p0 > p1)
				# si e1 perdió contra e0
				A(e1, e0) += p0 - p1;
			endif

			partidos--;
		endif
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
	while(abs(l(i, i) - 1) > pres)
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

	save_solution(equipos, S, out_filename, team_codes_filename);
endfunction

function AFA(in_filename, out_filename, team_codes_filename = '', date_limit = 0)
	has_date_limit = date_limit != 0;

	fid = fopen(in_filename, 'r');
	[equipos, partidos] = fscanf(fid, '%u %u' ,"C");

	# Genero la matriz S
	S = zeros(equipos, 2);
	S(:, 1) = 1:equipos;
	while (partidos > 0)
		[f, e0, p0, e1, p1] = fscanf(fid, '%u %u %u %u %u', "C");

		if(has_date_limit && f > date_limit)
			partidos = 0;
		else
			if(p0 == p1)
				# si e0 empató contra e1
				S(e0, 2) += 1;
				S(e1, 2) += 1;
			elseif(p0 > p1)
				# si ganó e0
				S(e0, 2) += 3;
			else
				# si ganó e1
				S(e1, 2) += 3;
			endif

			partidos--;
		endif
	endwhile
	fclose(fid);

	save_solution(equipos, S, out_filename, team_codes_filename);
endfunction

function save_solution(equipos, S, out_filename, team_codes_filename)
	has_team_codes = !strcmp(team_codes_filename, '');

	S = sortrows(S, 2);

	# Si tengo el nombre de cada equipo, lo cargo
	if(has_team_codes)
		i = equipos;
		fid = fopen(team_codes_filename, 'r');
		teamcodes = cell(1, equipos);
		while (i > 0)
			team_code = fgetl(fid);
			codearray = strsplit(team_code, ',');
			teamcodes{str2num(codearray{1})} = codearray{2};
			i--;
		endwhile
		fclose(fid);
	endif

	# Escribo la solución
	fid = fopen(out_filename, 'a');

	if(has_team_codes)
		for i = 0:equipos - 1
			fprintf(fid, '%u, %s %f\n', S(equipos - i,1), teamcodes{S(equipos - i,1)},S(equipos - i, 2));
		endfor
	else
		for i = 0:equipos - 1
			fprintf(fid, '%u %f\n', S(equipos - i,1), S(equipos - i, 2));
		endfor
	endif

	fclose(fid);
endfunction
