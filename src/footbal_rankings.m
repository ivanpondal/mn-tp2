1;

function res = GeM(in_filename, out_filename, team_codes_filename = '', c = 0.85, date_limit = 0, pres = 0.0001)
	has_date_limit = date_limit != 0;

	fid = fopen(in_filename, 'r');
	[teams, matches] = fscanf(fid, '%u %u' ,"C");

	# Genero la matriz A
	A = zeros(teams);
	while (matches > 0)
		[f, e0, p0, e1, p1] = fscanf(fid, '%u %u %u %u %u', "C");

		if(has_date_limit && f > date_limit)
			matches = 0;
		else
			if(p0 < p1)
				# si e0 perdió contra e1
				A(e0, e1) += p1 - p0;
			elseif(p0 > p1)
				# si e1 perdió contra e0
				A(e1, e0) += p0 - p1;
			endif

			matches--;
		endif
	endwhile
	fclose(fid);

	# Genero la matriz H
	H = zeros(teams);
	a = zeros(teams, 1);
	for i = 1:teams
		suma = sum(A(i, :));
		if(suma > 0)
			H(i, :) = A(i, :)/sum(A(i, :));
		else
			# vector para equipos invictos
			a(i) = 1;
		endif
	endfor

	# Genero la matriz G
	e = ones(teams, 1);

	# Genero vectores de personalización
	# (por defecto distribución uniforme)

	# Equipos invictos
	u = ones(teams, 1)/teams;
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
	S = zeros(teams, 2);
	for i = 1:teams
		S(i, 1) = i;
		S(i, 2) = x(i);
	endfor

	res = S;
	if(!strcmp(out_filename, ''))
		save_solution(teams, S, out_filename, team_codes_filename);
	endif
endfunction

function GeM_evolution(in_filename, out_filename, teams, number_of_dates)
	season_dates = cell(1, number_of_dates);
	for i = 1:number_of_dates
		season_dates{i} = GeM(in_filename, '', '', 0.85, i);
	endfor

	fid = fopen(out_filename, 'w');
	for i = 1:number_of_dates
		% Imprimo número de fecha
		fprintf(fid, '%u', i);
		S = season_dates{i};
		% Imprimo el ranking de cada equipo en esa fecha
		for j = 1:teams
			fprintf(fid, ' %f', S(j,2));
		endfor
		fprintf(fid, '\n');
	endfor
	fclose(fid);

endfunction

function AFA(in_filename, out_filename, team_codes_filename = '', date_limit = 0)
	has_date_limit = date_limit != 0;

	fid = fopen(in_filename, 'r');
	[teams, matches] = fscanf(fid, '%u %u' ,"C");

	# Genero la matriz S
	S = zeros(teams, 2);
	S(:, 1) = 1:teams;
	while (matches > 0)
		[f, e0, p0, e1, p1] = fscanf(fid, '%u %u %u %u %u', "C");

		if(has_date_limit && f > date_limit)
			matches = 0;
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

			matches--;
		endif
	endwhile
	fclose(fid);

	save_solution(teams, S, out_filename, team_codes_filename);
endfunction

function AFA_evolution(in_filename, out_filename, matches_per_date, normalize_score = 1)
	fid = fopen(in_filename, 'r');
	[teams, matches] = fscanf(fid, '%u %u' ,"C");

	# Genero la matriz S
	S = zeros(teams, 2);
	S(:, 1) = 1:teams;
	number_of_dates = matches/matches_per_date;
	season_dates = cell(1, number_of_dates);
	while (matches > 0)
		[f, e0, p0, e1, p1] = fscanf(fid, '%u %u %u %u %u', "C");

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

		matches--;
		if (mod(matches, matches_per_date) == 0)
			X = S;
			if (normalize_score == 1)
				X(:, 2) = X(:, 2)/sum(X(:, 2));
			endif
			season_dates{f} = X;
		endif
	endwhile
	fclose(fid);

	fid = fopen(out_filename, 'w');
	for i = 1:number_of_dates
		% Imprimo número de fecha
		fprintf(fid, '%u', i);
		S = season_dates{i};
		% Imprimo el ranking de cada equipo en esa fecha
		for j = 1:teams
			fprintf(fid, ' %f', S(j,2));
		endfor
		fprintf(fid, '\n');
	endfor
	fclose(fid);
endfunction

function save_solution(teams, S, out_filename, team_codes_filename)
	has_team_codes = !strcmp(team_codes_filename, '');

	S = sortrows(S, 2);

	# Si tengo el nombre de cada equipo, lo cargo
	if(has_team_codes)
		i = teams;
		fid = fopen(team_codes_filename, 'r');
		teamcodes = cell(1, teams);
		while (i > 0)
			team_code = fgetl(fid);
			codearray = strsplit(team_code, ',');
			teamcodes{str2num(codearray{1})} = codearray{2};
			i--;
		endwhile
		fclose(fid);
	endif

	# Escribo la solución
	fid = fopen(out_filename, 'w');

	if(has_team_codes)
		for i = 0:teams - 1
			fprintf(fid, '%u, %s %f\n', S(teams - i,1), teamcodes{S(teams - i,1)},S(teams - i, 2));
		endfor
	else
		for i = 0:teams - 1
			fprintf(fid, '%u %f\n', S(teams - i,1), S(teams - i, 2));
		endfor
	endif

	fclose(fid);
endfunction
