# Autor: Sandro Roger Boschetti
#  Data: 25 de setembro de 2016

# Programa implementado para a realizacao de calculos de blindagem
# em medicina nuclear. O texto eh escrito sem acentos e cedilhas 
# para ficar compativel com alguns sistemas que nao utilizam a 
# encodificacao UTF-8, como o pacote listings do LaTeX utilizado
# para a listagem de codigos.

# Para executar esse codigo basta fazer um copiar e colar no site
# http://octave-online.net

# Esse programa encontra-se por tempo indeterminado em:
# https://github.com/sandrorb/CalculoBlindagemMN
# https://sandrorb.github.io/CalculoBlindagemMN
# e pode ser retirado do ar a qualquer momento

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
# G(1) = 0.00705 para o Tc-99m quando a fonte eh o paciente e 0.0141 caso contrario
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
#Rt = Tf .* (1 - e.^( -(log(2) ./ Tf) * t)) / (log(2) * t);
#Fu = e.^(- (log(2) ./ Tf) * tu);
#doseSemBlindagem = (G .* A .* N .* Rt .* Fu .* T * t) / d^2;
#######################  Calculo da Dose : Fim   #######################



function doseSemBlindagem = calculaDose(G, A, N, t, tu, T, d, Tf)
	Rt = Tf .* (1 - e.^( -(log(2) ./ Tf) * t)) / (log(2) * t);
	Fu = e.^(- (log(2) ./ Tf) * tu);
	doseSemBlindagem = (G .* A .* N .* Rt .* Fu .* T * t) / d^2;
endfunction



function calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite)
	doseSemBlindagem = calculaDose(G, A, N, t, tu, T, d, Tf);
	[x, y, z] = calculaEspessuras(mu, doseSemBlindagem, doseLimite);
	printf("Dose: %.1f muSv (sem blindagem)\n", sum(doseSemBlindagem));
	printf("Limite: %.1f muSv\n", doseLimite);
	printf("Blindagem (cm de Pb, barita e concreto): %6.3f, %6.3f, %6.3f\n\n", x, y, z);
endfunction


printf("W significa Parede ou Porta. F, fonte e P, ponto de interesse.\n\n");


##########################################################################
printf("Calculos a partir da Sala de Exame:\n\n");

# parametros fixos para a Sala de Exames
G(1) = 0.00705;
t = 1.5;
tu = 1.5;

printf("W1, F1 e P1:\n");
T = 1/5; d = 1.82; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W2, F1 e P2: ");
T = 1/5; d = 2.89; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W3, F1 e P3: ");
T = 1/5; d = 2.51; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Porta da Sala de Exames
printf("W4, F1 e P4: ");
T = 1/5; d = 3.78; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W5, F1 e P5: ");
T = 1/5; d = 3.9; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W6, F1 e P6: ");
T = 1/20; d = 3.23; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################




##########################################################################
printf("Calculos a partir da Sala de Rejeitos (Atencao: Situacao Especial!!!):\n\n");

# parametros fixos para a Sala de Rejeitos
# apenas uma dose (N = 1)*5 de cada radionuclideo com atividade total da semana
# por 24 horas por dia decaindo fora da blindagem.
# Aproximacao: no caso dos radionuclideos exceto Tc-99m, ha um overlap de doses
# semanais que nao eh considerado.

G(1) = 0.0141;
t = 24.0;
tu = 0.0;
AmCi = [1000 50 10 10 10 50];
A = AmCi .* 37;
N = [1 1 1 1 1 1] * 5;

printf("W7, F2 e P7: ");
T = 1/5; d = 1.32; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
# doseSemBlindagem = calculaDose(G, A, N, t, tu, T, d, Tf)

printf("W8, F2 e P8: ");
T = 1/5; d = 1.47; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W9, F2 e P9: ");
T = 1/5; d = 1.44; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W10, F2 e P10: ");
T = 1/5; d = 1.62; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("Como tambem o Tc-99m e os outros radionuclideos estarao todos blindados,\n");
printf("baseado nos calculos acima, parece razoavel aplicacao de 2,5 cm de barita\n");
printf("em ambos os lados da parede\n\n");

printf("REVER TUDO ISSO POIS NAO ESTA BOM!!!!\n\n");
##########################################################################



##########################################################################
printf("Calculos a partir do Laboratorio:\n\n");

# suposicao que a maiorr dose usual de marcacao com Tc-99m (300 mCi) e as 
# doses usais de administracao aos pacientes dos outros radionuclideos 
# ficam expostas por um tempo aproximado de 2h por dia.

G(1) = 0.0141;
t = 2.0;
tu = 0.0;
AmCi = [300 30 5 5 10 50];
A = AmCi .* 37;
N = [1 1 1 1 1 1] * 5;

printf("W11, F3 e P11: ");
T = 1/5; d = 2.6; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W12, F3 e P12: ");
printf("ATENCAO: Na duvida sobre qual T usar\n");
T = 1; d = 2.09; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W13, F3 e P13: ");
T = 1/5; d = 1.89; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W14, F3 e P14: ");
T = 1/20; d = 2.13; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W15, F3 e P15: ");
T = 1/5; d = 1.89; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W16, F3 e P16: ");
T = 1/5; d = 1.73; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################



##########################################################################
printf("Calculos a partir da Sala de administracao de Radiofarmacos:\n\n");

# com relacao as consideracoes anteriores, aqui a dose de Tc-99m exposta
# e considerada de 30 mCi

printf("ATENCAO: Verificar se o tempo de exposicao (t)!!!\n");
G(1) = 0.0141;
t = 0.5;
tu = 0.0;
AmCi = [30 30 5 5 10 50];
A = AmCi .* 37;
N = [1 1 1 1 1 1] * 5;

printf("W17, F4 e P17: ");
printf("ATENCAO: Na duvida sobre qual T usar\n");
T = 1/20; d = 1.78; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W18, F4 e P18: ");
T = 1/20; d = 1.92; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W19, F4 e P19: ");
printf("ATENCAO: Na duvida sobre qual T usar\n");
T = 1; d = 1.77; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W20, F4 e P20: ");
T = 1/5; d = 1.91; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W21, F4 e P21: ");
T = 1/5; d = 1.43; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################


##########################################################################
printf("Calculos a partir da Sala de Espera de Pacientes Injetados:\n\n");

G(1) = 0.00705;
t = 0.5;
tu = 0.0;
AmCi = [30 30 5 5 10 50];
A = AmCi .* 37;
N = [60 10 5 4 2 1];

printf("W22, F5 e P22: ");
T = 1/20; d = 2.96; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W23, F5 e P23: ");
T = 1/20; d = 2.78; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W24, F5 e P24: ");
printf("ATENCAO: Na duvida sobre qual T usar\n");
printf("ATENCAO: PAREDE CRITICA!!!\n");
T = 1; d = 1.11; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W25, F5 e P25: ");
T = 1/20; d = 3.28; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W26, F5 e P26: ");
T = 1/5; d = 4.5; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W27, F5 e P27: ");
T = 1/5; d = 3.0; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W28, F5 e P28: ");
T = 1/5; d = 2.23; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################



##########################################################################
printf("Calculos a partir do Sanitario Exclusivo de Pacientes Injetados:\n\n");

# paciente fica em media 10 minutos no sanitario

G(1) = 0.00705;
t = 10/60;
tu = 0.0;
AmCi = [30 30 5 5 10 50];
A = AmCi .* 37;
N = [60 10 5 4 2 1];

printf("W29, F6 e P29: ");
T = 1/5; d = 1.16; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W30, F6 e P30: ");
T = 1; d = 1.82; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W31, F6 e P23: ");
T = 1/20; d = 1.85; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W32, F6 e P31: ");
T = 1/20; d = 1.26; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W33, F6 e P32: ");
T = 1/5; d = 1.55; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################


##########################################################################
printf("Calculos a partir da Ergometria:\n\n");

G(1) = 0.00705;
t = 0.5;
tu = 0.0;
AmCi = [30 30 5 5 10 50];
A = AmCi .* 37;
N = [60 10 5 4 2 1];

printf("W34, F7 e P33: ");
T = 1/20; d = 2.0; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W35, F7 e P34: ");
T = 1/20; d = 2.41; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W36, F7 e P35: ");
T = 1; d = 1.55; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W37, F7 e P30: ");
T = 1; d = 2.57; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W38, F7 e P36: ");
T = 1/20; d = 3.15; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W39, F7 e P37: ");
T = 1/5; d = 1.59; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("W40, F7 e P38: ");
T = 1/20; d = 2.12; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################


##########################################################################
# printf("Parede 1 - Fonte 1 - Ponto 1\n");
# %             G       A    N    t   tU    T           d              Tf
# %\calcDose{0.00705}{30*37}{60}{0.5}{1.5}{1/5}{5.9/2 + 0.15 + 0.30}{6.02}
# # isso aqui esta errado. G eh matriz e nao um escalar
# G = 0.00705; T = 1/5; d = 2.75; doseLimite = 20;
# doseSemBlindagem = calculaDose(G, A, N, t, tu, T, d, Tf);
# [x, y, z] = calculaEspessuras(mu, doseSemBlindagem, doseLimite);
# printf("A dose inicial de %.1f muSv pode ser blindada por:\n", sum(doseSemBlindagem));
# printf("%6.3f cm de Pb\n", x);
# printf("%6.3f cm de barita\n", y);
# printf("%6.3f cm de concreto\n\n", z);
##########################################################################

















