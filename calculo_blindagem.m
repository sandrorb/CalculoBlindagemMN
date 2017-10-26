#       Autor: Sandro Roger Boschetti
#        Data: 22 de novembro de 2016 às 11h09min
# Atualizacao: 26 de outubro de 2017 às 10h19min

# Programa implementado para a realização de cálculos de blindagem
# em medicina nuclear.

# Para executar esse código basta fazer um copiar e colar do conteúdo
# desse arquivo para o site http://octave-online.net

# Esse programa encontra-se por tempo indeterminado em:
# https://github.com/sandrorb/CalculoBlindagemMN
# https://sandrorb.github.io/CalculoBlindagemMN
# e pode ser retirado do ar a qualquer momento

# source('calculo_blindagem.m');

# Variáveis globais usadas para acumular dados para impressão
# dentro de uma função. A variável wfp significa Wall Fonte Ponto
global wfp;
global dadosParaImpressao;

clc;

printf("Calculos realizados em 21 de outubro de 2017 as 10h19min\n\n");

########################### Definicoes : Inicio ###########################
sigla = cellstr(['Tc-99m'; 'I-131'; 'I-123'; 'Ga-67'; 'Tl-201'; 'Sm-153']);

# Meias-vidas fisicas
Tf = [6.02 192 13.2235 78.24 73.0104 46.284];

# Numero de pacientes por semana para cada radionuclideo
# Para cada area considerada esses valores podem mudar
NumeroPacientesTc99m = 120;
N = [NumeroPacientesTc99m 10 5 4 2 1];

# Atividade media administrada de cada radionuclideo em mCi
# Esse valores podem ser alterados em locais onde ha consideracoes especiais
# tais como Sala de Rejeitos, Laboratorio, Injecao e Ergometria
AmCi = [30 30 5 5 10 50];

# Atividade adquirida semanalmente a ser usada na sala de rejeitos
# AmCi = [1500 100 25 20 20 50];

# Atividade em MBq
A = AmCi .* 37;

# Gamao em (microSv m^2) / (Mbq h)
# G(1) = 0.00705 para o Tc-99m quando a fonte eh o paciente e 0.0141 caso
# contrario. Mas a CNEN nao aceitou considerar a atenuacao.
G = [0.0141 0.07647 0.07478 0.03004 0.02372 0.02440];

# Camadas semirredutoras em cm para o Pb
#csrPb = [0.017 0.233 0.039 0.034 0.017 0.014];
csrPb = [0.025 0.233 0.039 0.034 0.017 0.014];

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
	
  # mu(1,:) Pb 
	while (sum(doseComBlindagem) > doseLimite) 
		x = x + delta;
		doseComBlindagem = doseSemBlindagem .* exp (- mu(1,:) * x);
	endwhile

	y = 0.0;
	doseComBlindagem = doseSemBlindagem .* exp (- mu(2,:) * y);
	doseInicial = sum(doseComBlindagem);

  # mu(2,:) Barita
	while (sum(doseComBlindagem) > doseLimite) 
		y = y + delta;
		doseComBlindagem = doseSemBlindagem .* exp (- mu(2,:) * y);
	endwhile

	z = 0.0;
	doseComBlindagem = doseSemBlindagem .* exp (- mu(3,:) * z);
	doseInicial = sum(doseComBlindagem);

  # mu(3,:) Concreto
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
  global contador;

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


##########################################################################
# Saida de dados para tabela LaTeX
# Estrutura do array dadosParaImpressao: [W F P limite dose Pb Barita Concreto]
function printLatex(fn)
  global dadosParaImpressao;
  fid = fopen(fn, "w");
  fprintf(fid, "\\textbf{W} & \\textbf{F} &  \\textbf{P} & ");
  fprintf(fid, "\\textbf{Limite} & \\textbf{Dose} & \\textbf{Pb} & ");
  fprintf(fid, "\\textbf{Barita} & \\textbf{Concreto}  \\\\ \\hline \n");
  for i = 1:rows(dadosParaImpressao)
	  fprintf(fid, "%3d & %3d & ", dadosParaImpressao(i,1), dadosParaImpressao(i,2));
	  fprintf(fid, "%3d & %10.2f & ", dadosParaImpressao(i,3), dadosParaImpressao(i,4));
	  fprintf(fid, "%10.2f & %10.3f & ", dadosParaImpressao(i,5), dadosParaImpressao(i,6));
	  if (dadosParaImpressao(i,7) > 2.5)	
      	fprintf(fid, "\\red{%10.3f} & ", dadosParaImpressao(i,7));
      else
    	  fprintf(fid, "%10.3f & ", dadosParaImpressao(i,7));
      endif
	  fprintf(fid, "%10.3f \\\\ \n", dadosParaImpressao(i,8));
  endfor
  fclose(fid);
endfunction
##########################################################################


printf("W significa Parede ou Porta. F, fonte e P, ponto de interesse.\n\n\n");

##########################################################################
printf("Ergometria:\n\n");

# Consideracao de que apenas fontes de Tc-99m com dose tipica de 30 mCi
# ficam expostas por um periodo aproximado do procedimento de 30 minutos

printf("Suposição de que essa sala é usada apenas para pacientes com Tc-99m\n\n");

AmCi = [30 0 0 0 0 0]; A = AmCi .* 37;
NumeroPacientesTc99m = 120;
N = [NumeroPacientesTc99m 0 0 0 0 0];
t = 30 / 60;
tu = 0.0;

wfp = [1 1 1]; T = 1/5; d = 2.0; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1; d = 1.8; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/40; d = 2.21; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/5; d = 2.2; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [5 1 5]; T = 1/5; d = 2.24; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [6 1 5]; T = 1/5; d = (1.5 + 0.1 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [7 1 5]; T = 1; d = (3.4 - 1.5 + 0.1 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_ergometria.tex");
##########################################################################



wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Laboratório:\n\n");

printf("Admite-se que todas as fontes estão blindadas a não ser durante\n");
printf("o preparo de cada dose. Conservadoramente o tempo de preparo é 40h\n");
printf("semanais dividido pelo número de pacientes atendidos por semana\n\n");

AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = (8 * 5) / sum(N);
tu = 0.0;

wfp = [1 1 1]; T = 1/5; d = 1.98; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/5; d = 2.12; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1; d = 1.64; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/40; d = 2.31; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [5 1 5]; T = 1/5; d = 1.93; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [6 1 5]; T = 1/5; d = (1.5 + 0.1 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [7 1 5]; T = 1; d = (3.4 - 1.5 + 0.1 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_laboratorio.tex");
##########################################################################


wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Administração de Radiofármacos:\n\n");

# Consideracao de que apenas fontes de Tc-99m com dose tipica de 30 mCi
# ficam expostas por um periodo aproximado do procedimento de 30 minutos

printf("Suposicao de que cada administração tem duração média de 10 minutos\n\n");

AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = 10 / 60;
tu = 0.0;

wfp = [1 1 1]; T = 1/5; d = 1.96; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/5; d = 1.41; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/40; d = 2.41; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/5; d = 1.44; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [5 1 5]; T = 1/5; d = 2.06; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [6 1 5]; T = 1/5; d = (1.5 + 0.1 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [7 1 5]; T = 1; d = (3.4 - 1.5 + 0.1 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_administracao.tex");
##########################################################################



#clear wfp dadosParaImpressao;
clear -all;


