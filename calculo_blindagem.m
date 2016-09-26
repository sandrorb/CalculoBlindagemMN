# Autor: Sandro Roger Boschetti
#  Data: 25 de setembro de 2016

# Programa implementado para a realizacao de calculos de blindagem
# em medicina nuclear. O texto eh escrito sem acentos e cedilhas 
# para ficar compativel com alguns sistemas que nao utilizam a 
# encodificacao UTF-8, como o pacote listings do LaTeX utilizado
# para a listagem de codigos.

# Para executar esse codigo basta fazer um copiar e colar no site
# http://octave-online.net

# source('calculo_blindagem.m');

clc;

##########################################################################
function [x, y, z] = calculaEspessuras(mu, doseSemBlindagem, doseLimite)

	delta = 0.001;
	
	x = 0.0;
	doseComBlindagem = doseSemBlindagem .* exp (- mu(1,:) * x);
	doseInicial = sum(doseComBlindagem);
	
	while (sum(doseComBlindagem) > doseLimite) 
		x = x + delta;
		doseComBlindagem = doseSemBlindagem .* exp (- mu(1,:) * x);
	endwhile

	y = 0.0;
	doseComBlindagem = doseSemBlindagem .* exp (- mu(2,:) * y);
	doseInicial = sum(doseComBlindagem);

	while (sum(doseComBlindagem) > doseLimite) 
		y = y + delta;
		doseComBlindagem = doseSemBlindagem .* exp (- mu(2,:) * y);
	endwhile

	z = 0.0;
	doseComBlindagem = doseSemBlindagem .* exp (- mu(3,:) * z);
	doseInicial = sum(doseComBlindagem);

	while (sum(doseComBlindagem) > doseLimite) 
		z = z + delta;
		doseComBlindagem = doseSemBlindagem .* exp (- mu(3,:) * z);
	endwhile

endfunction
##########################################################################


################### Definicoes : Inicio ###################
sigla = cellstr(['Tc-99m'; 'I-131'; 'I-123'; 'Ga-67'; 'Tl-201'; 'Sm-153']);

# Meias-vidas fisicas
Tf = [6.02 192 13.2235 78.24 73.0104 46.284];

# Numero de pacientes por semana para cada radionuclideo
#N = [40 40 5 1 1 3];
N = [60 10 5 4 2 1];

# Atividade media administrada de cada radionuclideo em mCi
AmCi = [30 30 5 5 10 50];

# Atividade em MBq
A = AmCi .* 37;

# Gamao em (microSv m^2) / (Mbq h)
G = [0.0141 0.07647 0.07478 0.03004 0.02372 0.02440];

# Camadas semirredutoras em cm para o Pb
csrPb = [0.017 0.233 0.039 0.034 0.017 0.014];

# Camadas semirredutoras em cm para a barita
csrBarita = [0.272 1.776 0.580 0.508 0.149 0.056];

# Camadas semirredutoras em cm para a concreto
# CSR(concreto, Tc-99m) = 3.9 cm -> Facure 
#csrConcreto = [3.9 3.018 2.208 2.142 1.653 1.195];
csrConcreto = [1.911 3.018 2.208 2.142 1.653 1.195];


# Arranjo bidimensional dos dados das camadas semirredutoras. 
# Os ponto e virgulas separam as linhas. Espacos ou virgulas
# separam as colunas.
csr = [ csrPb; csrBarita; csrConcreto];

# Coeficiente de atenuacao linear em 1/cm
mu = log(2) ./ csr;


# Tempo em horas que cada paciente passa no recinto
t = 0.5;

# Tempo em horas de "uptake"
tu = 1.5;


############ Definicoes Relativas ao Ponto de Interesse ############
# Fator de ocupacao T
T = 1.0;

# Distancia em metros da fonte ao ponto de interesse a ser blindado.
d = 2.0;

# Dose Limite em microSv por semana
doseLimite = 20;
################### Definicoes : Fim ###################


####################### Calculo da Dose : Inicio ######################
Rt = Tf .* (1 - e.^( -(log(2) ./ Tf) * t)) / (log(2) * t);
Fu = e.^(- (log(2) ./ Tf) * tu);
doseSemBlindagem = (G .* A .* N .* Rt .* Fu .* T * t) / d^2;
#######################  Calculo da Dose : Fim   #######################



function doseSemBlindagem = calculaDose(G, A, N, t, tu, T, d, Tf)
	Rt = Tf .* (1 - e.^( -(log(2) ./ Tf) * t)) / (log(2) * t);
	Fu = e.^(- (log(2) ./ Tf) * tu);
	doseSemBlindagem = (G .* A .* N .* Rt .* Fu .* T * t) / d^2;
endfunction




##########################################################################
printf("Parede 1 - Fonte 1 - Ponto 1\n");
%             G       A    N    t   tU    T           d              Tf
%\calcDose{0.00705}{30*37}{60}{0.5}{1.5}{1/5}{5.9/2 + 0.15 + 0.30}{6.02}
G = 0.00705; T = 1/5; d = 2.75; doseLimite = 20;
doseSemBlindagem = calculaDose(G, A, N, t, tu, T, d, Tf);
[x, y, z] = calculaEspessuras(mu, doseSemBlindagem, doseLimite);
printf("A dose inicial de %.1f muSv pode ser blindada por:\n", sum(doseSemBlindagem));
printf("%6.3f cm de Pb\n", x);
printf("%6.3f cm de barita\n", y);
printf("%6.3f cm de concreto\n\n", z);
##########################################################################

##########################################################################
printf("Parede 2 - Fonte 1 - Ponto 2\n");
%             G       A    N    t   tU    T                   d              Tf
%\calcDose{0.00705}{30*37}{60}{0.5}{1.5}{1/5}{(3.55 - 0.15)/2 + 0.15 + 0.30}{6.02}
G = 0.00705; T = 1/20; d = 1.65; doseLimite = 20;
doseSemBlindagem = calculaDose(G, A, N, t, tu, T, d, Tf);
[x, y, z] = calculaEspessuras(mu, doseSemBlindagem, doseLimite);
printf("A dose inicial de %.1f muSv pode ser blindada por:\n", sum(doseSemBlindagem));
printf("%6.3f cm de Pb\n", x);
printf("%6.3f cm de barita\n", y);
printf("%6.3f cm de concreto\n\n", z);
##########################################################################

##########################################################################
printf("Parede 3 - Fonte 1 - Ponto 3\n");
%             G       A    N    t   tU    T                   d              Tf
%\calcDose{0.00705}{30*37}{60}{0.5}{1.5}{1/20}{2.98}{6.02}
G = 0.00705; T = 1/5; d = 2.36; doseLimite = 100;
doseSemBlindagem = calculaDose(G, A, N, t, tu, T, d, Tf);
[x, y, z] = calculaEspessuras(mu, doseSemBlindagem, doseLimite);
printf("A dose inicial de %.1f muSv pode ser blindada por:\n", sum(doseSemBlindagem));
printf("%6.3f cm de Pb\n", x);
printf("%6.3f cm de barita\n", y);
printf("%6.3f cm de concreto\n\n", z);
##########################################################################







##########################################################################
##########################################################################
##########################################################################
#                 O que vem abaixo sao meros testes                      #
##########################################################################
##########################################################################
##########################################################################


##########################################################################
printf("\n\nCalculo igual ao que esta em LaTeX somente para Tc-99m\n");
##########################################################################
printf("Parede 1 - Fonte 1 - Ponto 1\n");
A = [30 0 0 0 0 0] .* 37;
N = [60 0 0 0 0 0];
t = 0.5;
tu = 1.5;
T = 1/5;
d = 5.9/2 + 0.15 + 0.30;
Rt = Tf .* (1 - e.^( -(log(2) ./ Tf) * t)) / (log(2) * t);
Fu = e.^(- (log(2) ./ Tf) * tu);
doseSemBlindagem = (G .* A .* N .* Rt .* Fu .* T * t) / d^2;

# Valores sendo usados pelo calculo em LaTeX
mu(1,1) = log(2) / 0.03;
mu(2,1) = log(2) / 1.3863;
mu(3,1) = log(2) / 3.9;

[x,y,z] = calculaEspessuras(mu, doseSemBlindagem, doseLimite);

printf("   Dose sem blindagem: %7.3f muSv\n", doseSemBlindagem(1));
printf("      Espessura de Pb: %7.3f cm\n", x);
printf("  Espessura de Barita: %7.3f cm\n", y);
printf("Espessura de Concreto: %7.3f cm\n", z);
##########################################################################



################################################################################
printf("\n\nCalculo igual ao que esta em LaTeX somente para Tc-99m\n");
################################################################################
printf("Parede 8 - Fonte 2 - Ponto 8\n");
A = [1000 0 0 0 0 0] .* 37;
N = [1 0 0 0 0 0];
t = 24;
tu = 0.0;
T = 1/5;
d = 2.19/2 + 0.15 + 0.3;
Rt = Tf .* (1 - e.^( -(log(2) ./ Tf) * t)) / (log(2) * t);
Fu = e.^(- (log(2) ./ Tf) * tu);
doseSemBlindagem = (G .* A .* N .* Rt .* Fu .* T * t) / d^2;

# Valores sendo usados pelo calculo em LaTeX
mu(1,1) = log(2) / 0.03;
mu(2,1) = log(2) / 1.3863;
mu(3,1) = log(2) / 3.9;


[x,y,z] = calculaEspessuras(mu, doseSemBlindagem, doseLimite);

printf("   Dose sem blindagem: %7.3f muSv\n", doseSemBlindagem(1));
printf("      Espessura de Pb: %7.3f cm\n", x);
printf("  Espessura de Barita: %7.3f cm\n", y);
printf("Espessura de Concreto: %7.3f cm\n", z);
################################################################################

















