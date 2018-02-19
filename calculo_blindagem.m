#       Autor: Sandro Roger Boschetti
#        Data: 22 de novembro de 2016 às 11h09min
# Atualizacao: 19 de fevereiro de 2018 às 11h48min

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

printf("Cálculos realizados em 19 de fevereiro de 2018 às 11h48min\n\n");

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
# AmCi = [1500 560 25 20 20 50];

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
csr = [csrPb; csrBarita; csrConcreto];

# Coeficiente de atenuacao linear em 1/cm
mu = log(2) ./ csr;
########################### Definicoes : Fim ###########################



##########################################################################
# x eh a espessura de Pb, y eh a de Barita e z a de Concreto
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




##########################################################################
# Algumas distâncias existentes na planta
peDireitoSMN = 3.4;
peDireitoAndarSuperior = 3.4;
peDireitoAndarInferior = 2.9
espessuraLaje = 0.09;

# distância do paciente (1.5 m do chão) a 30 cm abaixo da laje inferior
dPacAlvoAndarInferior = 1.5 + espessuraLaje + 0.3; # exceção sala de rejeitos
dPacAlvoAndarSuperior = peDireitoSMN - 1.5 + espessuraLaje + 0.3;

# distância conservadora de 30 cm cada lado mais 10 cm de espessura = 0.7 m
##########################################################################


printf("W significa Parede ou Porta. F, fonte e P, ponto de interesse.\n\n\n");



##########################################################################
printf("Ergometria:\n\n");

# Consideração de que apenas fontes de Tc-99m com dose típica de 30 mCi
# ficam expostas por um período aproximado do procedimento de 30 minutos.

printf("Suposição de que essa sala é usada apenas para pacientes com Tc-99m\n\n");

AmCi = [30 0 0 0 0 0]; A = AmCi .* 37;
NumeroPacientesTc99m = 120; # numero super estimado
N = [NumeroPacientesTc99m 0 0 0 0 0];
t = 30 / 60;
tu = 0.0;

wfp = [1 1 1]; T = 1/5; d = 0.7; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1; d = 0.7; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/40; d = 0.9; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/5; d = 0.7; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [5 1 5]; T = 1/5; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [6 1 6]; T = 1; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_ergometria.tex");
##########################################################################


wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Administração de Radiofármacos:\n\n");
printf("Suposicao de que cada administração tem duração média de 10 minutos\n\n");

AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = 10 / 60;
tu = 0.0;

wfp = [1 1 1]; T = 1/5; d = 0.7; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/5; d = 0.7; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/40; d = 0.9; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/5; d = 0.7; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [5 1 5]; T = 1/5; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [6 1 6]; T = 1; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_administracao.tex");
##########################################################################




wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Laboratório:\n\n");

#----------------------------------------------------------------
printf("Admite-se que todas as fontes estão blindadas a não ser durante\n");
printf("o preparo de cada dose (e eluições e marcações) e que cada dose\n");
printf("permanece típica permanece não blindada por todo o tempo de atividade.\n");
printf("do serviço de medicina nuclear, ou seja, 40h/sem.\n");

#AmCi = [1500 0 25 20 20 50]; # sem iodo:
#AmCi = [30 0 5 5 10 50];     # sem iodo
AmCi = [30 30 5 5 10 50];     # com iodo: solicitação do SPR e RT
A = AmCi .* 37;
N = [1 1 1 1 1 1]; # uma dose típica de cada radionuclídeo
t = 8*5;           # o tempo todo durante 40h por semana
tu = 0.0;
#----------------------------------------------------------------

wfp = [1 1 1]; T = 1/5; d = 0.70; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/5; d = 0.70; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/40; d = 0.90; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Devido ao posicionamento do rejeito, consideraremos 0.30 cm
# Os rejeitos estão bem blindados e sobra apenas o que se
# considerou anteriormente, ou seja, doses típicas sem blindagens.
wfp = [4 1 4]; T = 1/5; d = 0.70; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

#--------------------------------------------------
#AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
#N = [1 1 1 1 1 1]; #[NumeroPacientesTc99m 10 5 4 2 1];
#t = 8 * 5;
#tu = 0.0;
printf("Ponto especial com considerações especial para proteção da câmara\n");
wfp = [5 1 5]; T = 1; d = 1.90; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
#--------------------------------------------------

# Piso
wfp = [6 1 6]; T = 1/5; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [7 1 7]; T = 1; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_laboratorio.tex");
##########################################################################




wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Sala de Exame:\n\n");

AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = 30 / 60;
tu = 90 / 60;

wfp = [1 1 1]; T = 1/5; d = 1.65; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/5; d = 3.10; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/5; d = 2.35; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/5; d = 2.30; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [5 1 5]; T = 1/40; d = 2.11; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [6 1 6]; T = 1/40; d = 1.70; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [7 1 7]; T = 1/5; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [8 1 8]; T = 1; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_exame.tex");
##########################################################################


wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Sanitário Exclusivo de Pacientes Injetados (ao lado de I. S. Func.):\n\n");

AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = 5 / 60;
tu = 0;

wfp = [1 1 1]; T = 1/20; d = 0.90; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/5; d = 1.41; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/5; d = 1.40; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1; d = 0.85; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [5 1 5]; T = 1/5; d = (1.5 + 0.09 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [6 1 5]; T = 1; d = (3.4 - 1.5 + 0.09 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_sanitario_esq.tex");
##########################################################################



wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Sanitário Exclusivo de Pacientes Injetados (ao lado do corredor):\n\n");

AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = 5 / 60;
tu = 0;

wfp = [1 1 1]; T = 1/8; d = 1.48; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/5; d = 1.65; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/20; d = 0.90; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1; d = 0.85; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [5 1 5]; T = 1/5; d = (1.5 + 0.09 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [6 1 5]; T = 1; d = (3.4 - 1.5 + 0.09 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_sanitario_dir.tex");
##########################################################################


wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Sala Exclusiva de Pacientes Injetados:\n\n");

AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = 90 / 60;
tu = 0;

wfp = [1 1 1]; T = 1/8; d = 1.80; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/8; d = 1.85; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/8; d = 3.80; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/20; d = 3.10; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [5 1 5]; T = 1/20; d = 2.5; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [6 1 6]; T = 1/20; d = 2.5; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [7 1 7]; T = 1/5; d = (1.5 + 0.09 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [8 1 8]; T = 1; d = (3.4 - 1.5 + 0.09 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_espera_injetados.tex");
##########################################################################


wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Corredor de Saída para Área Interna do Hospital:\n\n");

AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = 2 / 60;
tu = 0;

wfp = [1 1 1]; T = 1/8; d = 3.10; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1; d = 1.30; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/5; d = 1.07; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [4 1 7]; T = 1/5; d = (1.5 + 0.09 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [5 1 8]; T = 1; d = (3.4 - 1.5 + 0.09 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_corredor.tex");
##########################################################################

wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Vestiário Ergometria:\n\n");

AmCi = [30 0 0 0 0 0]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 0 0 0 0 0];
t = 10 / 60;
tu = 0;

wfp = [1 1 1]; T = 1; d = 0.80; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1; d = 0.90; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/5; d = 1.20; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/5; d = 0.85; doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [5 1 7]; T = 1/5; d = (1.5 + 0.09 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [6 1 8]; T = 1; d = (3.4 - 1.5 + 0.09 + 0.3); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

printLatex("tabela_dados_vestiario_ergo.tex");
##########################################################################

#clear wfp dadosParaImpressao;
clear -all;


