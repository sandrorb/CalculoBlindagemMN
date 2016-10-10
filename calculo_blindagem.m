# Autor: Sandro Roger Boschetti
#  Data: 10 de outubro de 2016 as 11h12min

# Programa implementado para a realizacao de calculos de blindagem
# em medicina nuclear. O texto eh escrito sem acentos e cedilhas 
# para ficar compativel com alguns sistemas que nao utilizam a 
# encodificacao UTF-8, como o pacote listings do LaTeX utilizado
# para a listagem de codigos.

# Para executar esse codigo basta fazer um copiar e colar do conteudo
# desse arquivo para o sitehttp://octave-online.net

# Esse programa encontra-se por tempo indeterminado em:
# https://github.com/sandrorb/CalculoBlindagemMN
# https://sandrorb.github.io/CalculoBlindagemMN
# e pode ser retirado do ar a qualquer momento

# source('calculo_blindagem.m');

# Variaveis globais usadas para acumular dados para impressao
# dentro de uma funcao. A variavel wfp significa Wall Fonte Ponto
global wfp;
global dadosParaImpressao; 

clc;

printf("Calculos realizados em 10 de outubro de 2016 as 11h12min\n\n");

########################### Definicoes : Inicio ###########################
sigla = cellstr(['Tc-99m'; 'I-131'; 'I-123'; 'Ga-67'; 'Tl-201'; 'Sm-153']);

# Meias-vidas fisicas
Tf = [6.02 192 13.2235 78.24 73.0104 46.284];

# Numero de pacientes por semana para cada radionuclideo
#N = [40 40 5 1 1 3];
N = [60 10 5 4 2 1];

# Atividade semanal total em Bq
AsemanalBq = [1500 100 25 20 20 50] .* 37;

# Atividade media administrada de cada radionuclideo em mCi
# Esse valores podem ser alterados em locais onde ha consideracoes especiais
# tais como Sala de Rejeitos, Laboratorio, Injecao e Ergometria
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
csrConcreto = [1.911 3.018 2.208 2.142 1.653 1.195];


# Arranjo bidimensional dos dados das camadas semirredutoras. 
# Os ponto-e-virgulas separam as linhas da matriz enquanto os espacos
# (ou virgulas) separam as colunas.
csr = [ csrPb; csrBarita; csrConcreto];

# Coeficiente de atenuacao linear em 1/cm
mu = log(2) ./ csr;
########################### Definicoes : Fim ###########################



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


##########################################################################
function doseSemBlindagem = calculaDose(G, A, N, t, tu, T, d, Tf)
	Rt = Tf .* (1 - e.^( -(log(2) ./ Tf) * t)) / (log(2) * t);
	Fu = e.^(- (log(2) ./ Tf) * tu);
	doseSemBlindagem = (G .* A .* N .* Rt .* Fu .* T * t) / d^2;
endfunction
##########################################################################


##########################################################################
function calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite)

	global wfp;
	global dadosParaImpressao;

	doseSemBlindagem = calculaDose(G, A, N, t, tu, T, d, Tf);
	[x, y, z] = calculaEspessuras(mu, doseSemBlindagem, doseLimite);
	
	printf("Parede: W%d, Fonte: F%d e Ponto: P%d: \n", wfp(1), wfp(2), wfp(3));
	printf("Dose: %.1f muSv | ", sum(doseSemBlindagem));
	printf("Limite: %.1f muSv | ", doseLimite);
	printf("Distancia: %.1f m\n", d);
	printf("Blindagem (cm de Pb, barita e concreto): %6.3f, %6.3f, %6.3f\n\n", x, y, z);
	
	aux = wfp;
	wfp = [wfp doseSemBlindagem doseLimite x y z];
 	dadosParaImpressao = vertcat(dadosParaImpressao, wfp);
 	wfp = aux;
	
endfunction
##########################################################################




printf("W significa Parede ou Porta. F, fonte e P, ponto de interesse.\n\n");


##########################################################################
printf("Sala de Exame:\n\n");

# parametros fixos para a Sala de Exames
G(1) = 0.00705;
t = 0.5;
tu = 1.5;

# A variavel wfp encerra os valores da parede (w), da fonte (f) e do ponto (p)
wfp = [1 1 1];
T = 1/5; d = 1.82; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2];
T = 1/5; d = 2.89; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3];
T = 1/5; d = 2.51; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Porta da Sala de Exames
wfp = [4 1 4];
T = 1/5; d = 3.78; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [5 1 5];
T = 1/5; d = 3.9; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [6 1 6];
T = 1/20; d = 3.23; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################




##########################################################################
printf("Sala de Rejeitos:\n\n");

# parametros fixos para a Sala de Rejeitos
# apenas uma dose (N = 1)*5 de cada radionuclideo com atividade total da semana
# por 24 horas por dia decaindo fora da blindagem.
# Aproximacao: no caso dos radionuclideos exceto Tc-99m, ha um overlap de doses
# semanais que nao eh considerado.

G(1) = 0.0141;
t = 24.0;
tu = 0.0;
# AmCi = [1000 50 10 10 10 50];
# Supondo que o I-131 seja muito bem blindado, entao admite-se a atividade
# dele como sendo 0 mCi
printf("Admitindo que o I-131 seja totalmente blindado (atividade 0)\n");
printf("Na verdade, todos os  radionuclideos sao blindados dentro da sala\n");
printf("Mesmo assim, admite-se, conservadoramente, que os demais radionuclideos\n");
printf("nao sao blindados e que toda atividade nao eh usada e guardada na sala.\n");
AmCi = [1000 0 10 10 10 50];
A = AmCi .* 37;
N = [1 1 1 1 1 1] * 5;
# Supondo que aproximadamente metade da atividade eh blindada 
# N = [1 1 1 1 1 1] * 2.5;

wfp = [7 2 7];
T = 1/5; d = 1.32; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [8 2 8];
T = 1/5; d = 1.47; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [9 2 9];
T = 1/5; d = 1.44; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [10 2 10];
T = 1/5; d = 1.62; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################



##########################################################################
printf("Laboratorio de Manipulacao e Armazenamento de Fontes em Uso:\n\n");

# Suposicao de que algum radionuclideo fica exposto sem blindagem por 2 h / dia
# Suposicao que cada um dos radionuclideos fica exposto um periodo igual nessas 2 h
# Suposicao de que toda atividade recebida fica exposta

G(1) = 0.0141;
t = 2.0  / numel(AmCi);
tu = 0.0;
#AmCi = [300 30 5 5 10 50];
AmCi = [1000 50 10 10 10 50];
A = AmCi .* 37;
N = [1 1 1 1 1 1] * 5;

wfp = [11 3 11]; T = 1/5; d = 2.6; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);


printf("ATENCAO: Na duvida sobre qual T usar\n");
wfp = [12 3 12]; T = 1; d = 2.09; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [13 3 13]; T = 1/5; d = 1.89; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [14 3 14]; T = 1/20; d = 2.13; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [15 3 15]; T = 1/5; d = 1.89; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [16 3 16]; T = 1/5; d = 1.73; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################



##########################################################################
printf("Sala de Administracao de Radiofarmacos:\n\n");

# Suposicao de que cada radionuclideo autorizado fica exposto por um
# periodo de 10 minutos em sua atividade tipica de "injecao".

G(1) = 0.0141;
t = 10 / 60;
tu = 0.0;
AmCi = [30 30 5 5 10 50];
A = AmCi .* 37;
N = [1 1 1 1 1 1] * 5;


printf("ATENCAO: Na duvida sobre qual T usar\n");
wfp = [17 4 17]; T = 1/20; d = 1.78; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [18 4 18]; T = 1/20; d = 1.92; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);


printf("ATENCAO: Na duvida sobre qual T usar\n");
wfp = [19 4 19]; T = 1; d = 1.77; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [20 4 20]; T = 1/5; d = 1.91; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [21 4 21]; T = 1/5; d = 1.43; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################


##########################################################################
printf("Sala de Espera de Pacientes Injetados:\n\n");

printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
printf("Rever o tempo de espera de cada paciente. Acho que nao eh 1/2 h e sim 1,5h.\n");
printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");
printf("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n");


G(1) = 0.00705;
t = 1.5;
tu = 0.0;
AmCi = [30 30 5 5 10 50];
A = AmCi .* 37;
N = [60 10 5 4 2 1];

wfp = [22 5 22]; T = 1/20; d = 2.96; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [23 5 23]; T = 1/20; d = 2.78; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("ATENCAO: Na duvida sobre qual T usar\n");
printf("ATENCAO: PAREDE CRITICA!!!\n");
wfp = [24 5 24]; T = 1; d = 1.11; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [25 5 25]; T = 1/20; d = 3.28; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [26 5 26]; T = 1/5; d = 4.5; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [27 5 27]; T = 1/5; d = 3.0; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [28 5 28]; T = 1/5; d = 2.23; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################



##########################################################################
printf("Sanitario Exclusivo de Pacientes Injetados:\n\n");

# paciente fica em media 10 minutos no sanitario

G(1) = 0.00705;
t = 10/60;
tu = 0.0;
AmCi = [30 30 5 5 10 50];
A = AmCi .* 37;
N = [60 10 5 4 2 1];

wfp = [29 6 29]; T = 1/5; d = 1.16; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [30 6 30]; T = 1; d = 1.82; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [31 6 23]; T = 1/20; d = 1.85; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [32 6 31]; T = 1/20; d = 1.26; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [33 6 32]; T = 1/5; d = 1.55; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################


##########################################################################
printf("Ergometria:\n\n");

# Consideracao de que apenas fontes de Tc-99m com dose tipica de 30 mCi
# ficam expostas por um periodo aproximado do procedimento de 30 minutos

G(1) = 0.00705;
t = 30 / 60;
tu = 0.0;
#AmCi = [30 30 5 5 10 50];
AmCi = [30 0 0 0 0 0];
A = AmCi .* 37;
N = [60 10 5 4 2 1];


wfp = [34 7 33]; T = 1/20; d = 2.0; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [35 7 34]; T = 1/20; d = 2.41; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [36 7 35]; T = 1; d = 1.55; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [37 7 30]; T = 1; d = 2.57; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [38 7 36]; T = 1/20; d = 3.15; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [39 7 37]; T = 1/5; d = 1.59; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [40 7 38]; T = 1/20; d = 2.12; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################



##########################################################################
printf("Sala Acima da Sala de Exames: \n\n");
# parametros fixos para a Sala de Exames
G(1) = 0.00705;
t = 0.5;
tu = 1.5;

wfp = [41 1 39]; 
peDireito = 3.0;
T = 1; d = peDireito - 1.5 + 1.5; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################



















