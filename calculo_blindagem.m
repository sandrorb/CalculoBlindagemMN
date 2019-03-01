#       Autor: Sandro Roger Boschetti
#     Contato: linkedin.com/in/sandroboschetti
#        Data: 22 de novembro de 2016 às 11h09min
# Atualizacao: 27 de fevereiro de 2019 às 16h23min

# Programa implementado para a realização de cálculos de blindagem
# em medicina nuclear.

# Programas em Octave podem ser executados online em http://octave-online.net
# Para executar este programa sugere-se utilizar os programas gratuitos que
# podem ser encontrados em https://www.gnu.org/software/octave

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

printf("Cálculos realizados em 27 de fevereiro de 2019 às 16h23min\n\n");

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



########################################################################
# As variáveis x, y e z são as espessuras de Pb, Barita e
# Concreto respectivamente.
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
	#wfp = [wfp doseLimite sum(doseSemBlindagem) x y z];
  wfp = [wfp doseLimite sum(doseSemBlindagem) x y z t tu T d];
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
# Saida de dados para tabela LaTeX
# Estrutura do array dadosParaImpressao: [W F P limite dose Pb Barita Concreto]
function printLatex2(fn)
  global dadosParaImpressao;
  fid = fopen(fn, "w");

# Preparação da primeira linha do cobeçalho
  fprintf(fid, "\\\\ \\cline{8-12}\n");
  fprintf(fid, "\\multicolumn{7}{c|}{} & \n");
  fprintf(fid, "\\multicolumn{2}{c|}{\\textbf{Doses em $\\mu$Sv}}  & \n"); 
  fprintf(fid, "\\multicolumn{3}{c|}{\\textbf{Espessuras em cm}} \\\\ \\hline \n");

# Preparação da segunda linha do cabeçalho
  fprintf(fid, "\\textbf{W} & \\textbf{F} &  \\textbf{P} & ");
  fprintf(fid, "\\textbf{t} & \\textbf{t$_u$} & \\textbf{T} & \\textbf{d (m)} & ");
  fprintf(fid, "\\textbf{Limite} & \\textbf{Dose} & \\textbf{Pb} & ");
  fprintf(fid, "\\textbf{Barita} & \\textbf{Concreto} \\\\ \\hline \n");
  
  for i = 1:rows(dadosParaImpressao)
	  fprintf(fid,    "%3d & ", dadosParaImpressao(i,1)); # w (parede)
    fprintf(fid,    "%3d & ", dadosParaImpressao(i,2)); # f (fonte)
	  fprintf(fid,    "%3d & ", dadosParaImpressao(i,3)); # p (ponto de interesse)
    
    fprintf(fid, "%10.1f & ", dadosParaImpressao(i,9));  # t (tempo de permanência
    fprintf(fid, "%10.1f & ", dadosParaImpressao(i,10)); # tu (tempo de uptake
    fprintf(fid, "%10.2f & ", dadosParaImpressao(i,11)); # T (fator de ocupação)
    fprintf(fid, "%10.3f & ", dadosParaImpressao(i,12)); # d (distância e m)
    
    fprintf(fid, "%10.2f & ", dadosParaImpressao(i,4)); # dose limite
	  fprintf(fid, "%10.2f & ", dadosParaImpressao(i,5)); # dose sem blindagem
    
    fprintf(fid, "%10.3f & ", dadosParaImpressao(i,6)); # x espessura de Pb

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
# Os dados abaixo são valores padrão mas que podem ser modificados em
# outras partes do programa para atender aos requisitos de cada área do SMN.

# Algumas distâncias existentes na planta em metros
peDireitoSMN = 2.95;
peDireitoAndarSuperior = 2.95;
peDireitoAndarInferior = 2.95;
espessuraLaje = 0.09;

# distância do paciente (1.5 m do chão) a 30 cm abaixo da laje inferior
dPacAlvoAndarInferior = 1.5 + espessuraLaje + 0.3;

# distância do paciente a 30 cm da superfície da laje superior
dPacAlvoAndarSuperior = peDireitoSMN - 1.5 + espessuraLaje + 0.3;

fatorOcupAndarInf = 1;
fatorOcupAndarSup = 1;
##########################################################################


printf("W significa Parede ou Porta. F, fonte e P, ponto de interesse.\n\n\n");


# Esta auxiliar para calular a real distância em metros a partir do
# valor obtido na planta em PDF em milímetro reduzido pela escala 1/50.
function d = mm2m(mm)
  escala = 1/50;
  d = (mm / 1000) / escala;
endfunction




##########################################################################
##########################################################################
#                  AQUI INICIAM-SE DE FATO OS CÁLCULOS
##########################################################################
##########################################################################


##########################################################################
printf("Sala de Exame:\n\n");

AmCi = [30 5 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = 30 / 60;
tu = 90 / 60;

wfp = [1 1 1]; T = 1/5; d = mm2m(60.75); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/5; d = mm2m(55.25); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1; d = mm2m(52.33); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/5; d = mm2m(56.35); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [5 1 5]; T = 1; d = mm2m(54.73); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [6 1 6]; T = 1/5; d = mm2m(66.19); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [7 1 7]; T = fatorOcupAndarInf; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [8 1 8]; T = fatorOcupAndarSup; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

#printLatex("tabela_dados_exame.tex");
printLatex2("tabela_dados_exame.tex");
##########################################################################




wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Sanitário Exclusivo de Pacientes Injetados:\n\n");

AmCi = [30 5 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = 5 / 60;
tu = 0;

wfp = [1 1 1]; T = 1/40; d = mm2m(33.71); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/5; d = mm2m(24.71); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1; d = mm2m(36.95); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1; d = mm2m(25.17); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [5 1 5]; T = fatorOcupAndarInf; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [6 1 5]; T = fatorOcupAndarSup; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

#printLatex("tabela_dados_sanitario.tex");
printLatex2("tabela_dados_sanitario.tex");
##########################################################################




wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Sala Exclusiva de Pacientes Injetados:\n\n");

AmCi = [30 5 5 5 10 50]; A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1];
t = 90 / 60;
tu = 0;

wfp = [1 1 1]; T = 1/40; d = mm2m(36.25); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1; d = mm2m(62.50); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/5; d = mm2m(56.22); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1; d = mm2m(39.86); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [5 1 5]; T = 1; d = mm2m(45.98); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [6 1 6]; T = 1/20; d = mm2m(64.20); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [7 1 7]; T = fatorOcupAndarInf; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [8 1 8]; T = fatorOcupAndarSup; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

#printLatex("tabela_dados_espera_injetados.tex");
printLatex2("tabela_dados_espera_injetados.tex");
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

wfp = [1 1 1]; T = 1/40; d = mm2m(37.32); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1; d = mm2m(38.79); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1/5; d = mm2m(44.80); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/5; d = mm2m(34.44); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [5 1 5]; T = fatorOcupAndarInf; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [6 1 6]; T = fatorOcupAndarSup; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

#printLatex("tabela_dados_administracao.tex");
printLatex2("tabela_dados_administracao.tex");
##########################################################################



wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Laboratório:\n\n");

#----------------------------------------------------------------
printf("Admite-se que todas as fontes estão blindadas, a não ser durante\n");
printf("o preparo de cada dose (eluições e marcações), e que cada dose\n");
printf("típica permanece não blindada por todo o tempo médio de 10 minutos.\n");
printf("São preparadas tantas doses quantos pacientes atendidos semanalmente.\n");

AmCi = [30 30 5 5 10 50];     # com iodo: solicitação do SPR e RT
A = AmCi .* 37;
N = [NumeroPacientesTc99m 10 5 4 2 1]; # uma dose típica de cada radionuclídeo
t = 10/60;
tu = 0.0;
#----------------------------------------------------------------

wfp = [1 1 1]; T = 1/40; d = mm2m(35.68); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/20; d = mm2m(34.06); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1; d = mm2m(47.20); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1; d = mm2m(55.44); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [5 1 5]; T = 1; d = mm2m(68.35); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [6 1 6]; T = 1/5; d = mm2m(65.97); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [7 1 7]; T = 1; d = mm2m(29.53); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Devido ao posicionamento do rejeito, consideraremos 0.30 cm
# Os rejeitos estão bem blindados e sobra apenas o que se
# considerou anteriormente, ou seja, doses típicas sem blindagens.
#wfp = [4 1 4]; T = 1/5; d = 0.70; doseLimite = 100;
#calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

#--------------------------------------------------
#AmCi = [30 30 5 5 10 50]; A = AmCi .* 37;
#N = [1 1 1 1 1 1]; #[NumeroPacientesTc99m 10 5 4 2 1];
#t = 8 * 5;
#tu = 0.0;
#printf("Ponto especial com considerações especial para proteção da câmara\n");
#wfp = [5 1 5]; T = 1; d = 1.90; doseLimite = 20;
#calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);
#--------------------------------------------------

# Piso
wfp = [8 1 8]; T = fatorOcupAndarInf; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [8 1 9]; T = fatorOcupAndarSup; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

#tabela_dados_laboratorio.tex");
printLatex2("tabela_dados_laboratorio.tex");
##########################################################################



wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Sala de Rejeitos:\n\n");

printf("Suposicao de que 10%% de cada dose semanal é descartada e de que são\n");
printf("armazenadas em blindagens de 5 mm de espessura de chumbo.\n\n");


AmCi = 0.1 .* [1000 50 15 20 10 100] .* exp(-log(2) * 0.5 ./ csrPb);
A = AmCi .* 37;
#A = A .* exp(-log(2) * 0.5 ./ csrPb);

N = [1 1 1 1 1 1];
t = 40;
tu = 0.0;

wfp = [1 1 1]; T = 1/40; d = mm2m(34.64); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/20; d = mm2m(25.64); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1; d = mm2m(39.14); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/8; d = mm2m(29.77); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [5 1 5]; T = fatorOcupAndarInf; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [6 1 6]; T = fatorOcupAndarSup; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

#printLatex("tabela_dados_rejeitos.tex");
printLatex2("tabela_dados_rejeitos.tex");
##########################################################################



wfp = [];
dadosParaImpressao = [];

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

wfp = [1 1 1]; T = 1; d = mm2m(39.72); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1; d = mm2m(32.67); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1; d = mm2m(29.08); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1/5; d = mm2m(43.88); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [5 1 5]; T = 1/5; d = mm2m(56.86); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [6 1 6]; T = 1; d = mm2m(51.63); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [7 1 7]; T = 1/5; d = mm2m(50.94); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [8 1 8]; T = fatorOcupAndarInf; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [9 1 9]; T = fatorOcupAndarSup; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

#printLatex("tabela_dados_ergometria.tex");
printLatex2("tabela_dados_ergometria.tex");
##########################################################################


wfp = [];
dadosParaImpressao = [];

##########################################################################
printf("Inalação:\n\n");

printf("Suposição de que essa sala é usada apenas para pacientes com Tc-99m\n\n");

AmCi = [30 0 0 0 0 0]; A = AmCi .* 37;
NumeroPacientesTc99m = 20; # numero super estimado
N = [NumeroPacientesTc99m 0 0 0 0 0];
t = 30 / 60;
tu = 0.0;

wfp = [1 1 1]; T = 1/5; d = mm2m(34.94); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [2 1 2]; T = 1/5; d = mm2m(31.03); doseLimite = 100;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [3 1 3]; T = 1; d = mm2m(29.40); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

wfp = [4 1 4]; T = 1; d = mm2m(17.80); doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Piso
wfp = [5 1 5]; T = fatorOcupAndarInf; d = dPacAlvoAndarInferior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

# Teto
wfp = [6 1 6]; T = fatorOcupAndarSup; d = dPacAlvoAndarSuperior; doseLimite = 20;
calculoParede(G, A, N, t, tu, T, d, Tf, mu, doseLimite);

printf("\n");

#printLatex("tabela_dados_inalacao.tex");
printLatex2("tabela_dados_inalacao.tex");
##########################################################################



#clear wfp dadosParaImpressao;
clear -all;


