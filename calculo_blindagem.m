# Autor: Sandro Roger Boschetti
#  Data: 10 de outubro de 2016 as 17h13min

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

printf("Calculos realizados em 10 de outubro de 2016 as 17h13min\n\n");

########################### Definicoes : Inicio ###########################
sigla = cellstr(['Tc-99m'; 'I-131'; 'I-123'; 'Ga-67'; 'Tl-201'; 'Sm-153']);

# Meias-vidas fisicas
Tf = [6.02 192 13.2235 78.24 73.0104 46.284];

# Numero de pacientes por semana para cada radionuclideo
# Para cada area considerada esses valores podem mudar
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
	printf("Dose: %.2f muSv | ", sum(doseSemBlindagem));
	printf("Limite: %.2f muSv | ", doseLimite);
	printf("Distancia: %.2f m\n", d);
	printf("Blindagem (cm de Pb, barita e concreto): %6.3f, %6.3f, %6.3f\n\n", x, y, z);
	
	aux = wfp;
	wfp = [wfp doseLimite sum(doseSemBlindagem) x y z];
 	dadosParaImpressao = vertcat(dadosParaImpressao, wfp);
 	wfp = aux;
	
endfunction
##########################################################################



printf("W significa Parede ou Porta. F, fonte e P, ponto de interesse.\n\n\n");

##########################################################################
printf("Sala de Exame:\n\n");

# parametros fixos para a Sala de Exames

# o primeiro elemento da matriz dos gamas (Tc-99m) eh mudado aqui pois
# essa area consta de fontes que sao pacientes cujos fotons sao blindados
# em cerca de 50
G(1) = 0.00705;
AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [60 10 5 4 2 1];
t = 0.5; tu = 1.5;

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

printf("\n");
##########################################################################




##########################################################################
printf("Sala de Rejeitos:\n\n");

# parametros fixos para a Sala de Rejeitos
# apenas uma dose (N = 1)*5 de cada radionuclideo com atividade total da semana
# por 24 horas por dia decaindo fora da blindagem.
# Aproximacao: no caso dos radionuclideos exceto Tc-99m, ha um overlap de doses
# semanais que nao eh considerado.

# o Gamao do Tc-99m volta ao normal ja que nao ha auto blindagem
G(1) = 0.0141;

printf("Suposicao que somente 10%% no material adquiro semanalmente vai para o rejeito\n");

# atividade no rejeito eh 10%
A = AsemanalBq .* 0.1;

# printf("Admitindo que o I-131 seja totalmente blindado (atividade 0)\n");
# printf("Na verdade, todos os  radionuclideos sao blindados dentro da sala\n");
# printf("Mesmo assim, admite-se, conservadoramente, que os demais radionuclideos\n");
# printf("nao sao blindados e que toda atividade nao eh usada e guardada na sala.\n\n");
# atividade de iodo no rejeito eh zero
#A(2) = 0;

N = [1 1 1 1 1 1] * 5;
t = 24.0;
tu = 0.0;

wfp = [7 2 7];
T = 1/5; d = 1.32; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [8 2 8];
T = 1/5; d = 1.47; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# foi posto aqui a dose limite de publico no intuito de se proteger o
# detector da radiacao da sala de rejeitos
printf("Para esse conjunto de fonte/parede/ponto abaixo, foi estabelecido o limite\n");
printf("de publico para se ter uma blindagem boa para o detector e minimizar\n");
printf("a possibilidade de interferencia nos exames devido a radiacao da Sala\n");
printf("de Rejeitos\n\n");
wfp = [9 2 9];
T = 1/5; d = 1.44; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [10 2 10];
T = 1/5; d = 1.62; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");
##########################################################################



##########################################################################
printf("Laboratorio de Manipulacao e Armazenamento de Fontes em Uso:\n\n");

printf("Suposicao de que algum radionuclideo fica exposto sem blindagem por 2 h / dia\n");
printf("Suposicao que cada um dos radionuclideos fica exposto um periodo igual nessas 2 h\n");
printf("Suposicao de que toda atividade recebida fica exposta\n\n");

G(1) = 0.0141;
A = AsemanalBq;
N = [1 1 1 1 1 1] * 5;
t = 2.0  / numel(AmCi);
tu = 0.0;

wfp = [11 3 11]; T = 1/5; d = 2.6; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

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

printf("\n");
##########################################################################



##########################################################################
printf("Sala de Administracao de Radiofarmacos:\n\n");

# Suposicao de que cada radionuclideo autorizado fica exposto por um
# periodo de 10 minutos em sua atividade tipica de "injecao".

printf("Suposicao de que cada radionuclideo autorizado fica exposto por um\n");
printf("periodo de 10 minutos em sua atividade tipica de <<injecao>>\n");

G(1) = 0.0141;
A = AsemanalBq;
N = [1 1 1 1 1 1] * 5;
t = 10 / 60;
tu = 0.0;

wfp = [17 4 17]; T = 1/5; d = 1.78; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [18 4 18]; T = 1/20; d = 1.92; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# o T = 1 e o limite = 20 devem ser iguais aos do ponto 12
wfp = [19 4 19]; T = 1; d = 1.77; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [20 4 20]; T = 1/5; d = 1.91; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [21 4 21]; T = 1/5; d = 1.43; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");
##########################################################################


##########################################################################
printf("Sala de Espera de Pacientes Injetados:\n\n");

# o Gamao volta ao valor com a autoblindagem
G(1) = 0.00705;
AmCi = [30 30 5 5 10 50];
A = AmCi .* 37;
N = [60 10 5 4 2 1];
t = 1.5;
tu = 0.0;

wfp = [22 5 22]; T = 1/20; d = 2.96; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("Como eh ponto cruzado, o limite eh dividido por 2\n");
wfp = [23 5 23]; T = 1/20; d = 2.78; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("O T = 1/5 justifica-se pois o funcionario nao fica a parede o tempo todo\n");
wfp = [24 5 24]; T = 1/5; d = 1.11; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [25 5 25]; T = 1/20; d = 3.28; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [26 5 26]; T = 1/5; d = 4.5; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [27 5 27]; T = 1/5; d = 3.0; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [28 5 28]; T = 1/5; d = 2.23; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");
##########################################################################



##########################################################################
printf("Sanitario Exclusivo de Pacientes Injetados:\n\n");

printf("Suposicao de que cada paciente fica em media 10 minutos no sanitario\n\n");

G(1) = 0.00705;
AmCi = [30 30 5 5 10 50];
A = AmCi .* 37;
N = [60 10 5 4 2 1];
t = 10/60;
tu = 0.0;

wfp = [29 6 29]; T = 1/5; d = 1.16; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("Como eh ponto cruzado, o limite eh dividido por 2\n");
wfp = [30 6 30]; T = 1; d = 1.82; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("Como eh ponto cruzado, o limite eh dividido por 2\n");
wfp = [31 6 23]; T = 1/20; d = 1.85; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [32 6 31]; T = 1/20; d = 1.26; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [33 6 32]; T = 1/5; d = 1.55; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");
##########################################################################


##########################################################################
printf("Ergometria:\n\n");

# Consideracao de que apenas fontes de Tc-99m com dose tipica de 30 mCi
# ficam expostas por um periodo aproximado do procedimento de 30 minutos

printf("Suposicao de que essa sala eh usada apenas para pacientes com Tc-99m\n\n");

G(1) = 0.00705;
AmCi = [30 0 0 0 0 0];
A = AmCi .* 37;
N = [60 0 0 0 0 0];
t = 30 / 60;
tu = 0.0;

wfp = [34 7 33]; T = 1/20; d = 2.0; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [35 7 34]; T = 1/20; d = 2.41; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [36 7 35]; T = 1; d = 1.55; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("Como eh ponto cruzado, o limite eh dividido por 2\n");
wfp = [37 7 30]; T = 1; d = 2.57; doseLimite = 20/2;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [38 7 36]; T = 1/20; d = 3.15; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [39 7 37]; T = 1/5; d = 1.59; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [40 7 38]; T = 1/20; d = 2.12; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");
##########################################################################



##########################################################################
printf("Sala Acima da Sala de Exames: \n\n");

printf("Suposicao de que a dose parte de uma altura de 1,5 m. O ponto de\n");
printf("calculo eh 1,5 m acima do piso do andar de cima.\n");
printf("O pe-direito eh de 3,0 m.");

G(1) = 0.00705;
AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [60 10 5 4 2 1];
t = 0.5; tu = 1.5;

wfp = [41 1 39]; 
peDireito = 3.0;
T = 1; d = peDireito - 1.5 + 1.5; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");
##########################################################################

##########################################################################
printf("Muro/Parede ao longo do corredor: \n\n");

# Para calcular o tempo que um paciente gasta em um trecho de 1 m
# caminhando 4 km/h = 4000 m / h

printf("O calculo eh feito para cada metro de parede. Supoe-se que o paciente\n");
printf("anda a uma velocidade de 4 km/h (4000 m/h). O tempo que cada paciente\n");
printf("fica em cada metro de parede eh de (1 m) / (4000 m/h)\n\n");

wfp = [42 8 40]; 

G(1) = 0.00705;
AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [60 10 5 4 2 1];

t = 1 / 4000;
tu = 0.0;

T = 1; d = 1.2/2 + 0.15 + 0.3; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
##########################################################################



# Saida de dados para tabela LaTeX
# Estrutura do array dadosParaImpressao: [W F P limite dose Pb Barita Concreto]

fid = fopen("tabela_dados.tex", "w");
fprintf(fid, "\\textbf{W} & \\textbf{F} &  \\textbf{P} & ");
fprintf(fid, "\\textbf{Limite} & \\textbf{Dose} & \\textbf{Pb} & ");
fprintf(fid, "\\textbf{Barita} & \\textbf{Concreto}  \\\\ \\hline \n");
for i = 1:rows(dadosParaImpressao)
	fprintf(fid, "%3d & %3d & ", dadosParaImpressao(i,1), dadosParaImpressao(i,2));
	fprintf(fid, "%3d & %10.2f & ", dadosParaImpressao(i,3), dadosParaImpressao(i,4));
	fprintf(fid, "%10.2f & %10.2f & ", dadosParaImpressao(i,5), dadosParaImpressao(i,6));
	fprintf(fid, "%10.2f & ", dadosParaImpressao(i,7));
	fprintf(fid, "%10.2f \\\\ \n", dadosParaImpressao(i,8));
endfor
fclose(fid);

printf("\n");






