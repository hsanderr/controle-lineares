%% Definicao da matriz de incertezas (delta A)

clear; % limpa variaveis do espaco de trabalho
options=timeoptions;
options.Grid='on'; % usa grade por default nos graficos

% matriz de estado
A = [0          376.9911    0           0        ;
     -0.15685    0           -0.0787    4        ;
     -0.16725    0           -0.46296   0.166667 ;
     1572.825    0           -5416.98   -100    ];

% matriz delta A
DA = A;
DA(1, 2) = 0;
DA(4, 4) = 0;
DA = DA * 0.05 .* (rand(size(DA, 1), size(DA, 2)) * 2 - 1);

save("deltaA.mat", "A", "DA");

%% Calculo do ganho otimo K

load("deltaA.mat", "A", "DA")

% Q = eye(4) .* (rand(4, 4) * 1e-3); % matriz Q do funcional LQR

Q = [8.78e-04   0           0           0;
     0          7.85e-04    0           0;
     0          0           3.34e-04    0;
     0          0           0           0.94e-04]; % matriz Q do funcional LQR

R = [1]; % matriz R do funcional LQR

B = [0 0 0 10000]'; % matriz de entrada

% determinacao do ganho otimo e solucao da equacao algebrica Riccati
% K: ganho otimo
% P: solucao da equacao algebrica Riccati associada
% Po: polos do sistema em malha fechada
[K, P, Po] = lqr(A + DA, B, Q, R);

%% Calculo do funcional J

% definicao das matrizes C e D para criar um espaco de estados com ss
C = eye(4);
D = zeros(4, 1);

sys = ss((A + DA - B * K), B, C, D); % criacao do espaco de estado
t = 0:1e-3:0.5; % vetor tempo (de 0 a 0,5s com passos de 1ms)
x0 = [0 0 0 0.1]'; % condicoes iniciais
[y, t, x] = initial(sys, x0, t); % resposta do sistema para condicoes iniciais
x = x'; % initial retorna o vetor de estado transposto

u = -K * x; % lei de controle por realimentacao de estados

% calculo de J pelo metodo 1
JD = x' * Q * x + u' * R * u; % integrando
J = zeros(1, length(JD)); % inicializacao de J
for i = 1:(length(JD))
    J(1, i) = JD(i, i);
end
J1 = trapz(t, J); % integracao trapezoidal

% calculo de J pelo metodo 2
J2 = x0' * P * x0;
   
%% Gerando o grafico da resposta para as condicoes iniciais

figure(1); 

subplot(4,1,1); plotx1=plot(t,x(1,:)); title('x_1');
ylabel('x'); xlabel('tempo (s)'); set(plotx1,'Linewidth',1.5); grid on;

subplot(4,1,2); plotx2=plot(t,x(2,:)); title('x_2');
ylabel('x'); xlabel('tempo (s)'); set(plotx2,'Linewidth',1.5); grid on;

subplot(4,1,3); plotx3=plot(t,x(3,:)); title('x_3');
ylabel('x'); xlabel('tempo (s)'); set(plotx3,'Linewidth',1.5); grid on;

subplot(4,1,4); plotx4=plot(t,x(4,:)); title('x_4');
ylabel('x'); xlabel('tempo (s)'); set(plotx4,'Linewidth',1.5); grid on;

%% Calculando os 100 funcionais

% Inicializacao dos vetores
J1n = zeros(1, 100);
J2n = zeros(1, 100);
n = 0:99;

% inicializacao dos erros quadraticos medios
mse1 = 0;
mse2 = 0;

% laco para criacao dos 100 modelos nao nominais e calculo dos funcionais
for i = 1:100

    % matriz delta A (incertezas)
    DAn = A;
    DAn(1, 2) = 0;
    DAn(4, 4) = 0;
    DAn = DAn * 0.05 .* (rand(size(DAn, 1), size(DAn, 2)) * 2 - 1);
    
    % calculo do ganho otimo e solucao da equacao algebrica Riccati
    [Kn, Pn, En] = lqr(A + DAn, B, Q, R);
    
    sysn = ss((A + DAn - B * K), B, C, D); % criaco do espaco de estados
    [yn, t, xn] = initial(sysn, x0, t); % resposta do sistema para condicoes iniciais
    xn = xn'; % initial retorna o vetor de estado transposto
    
    un= - K * xn; % lei de controle por realimentacao de estados
    
    % calculo do funcional pelo metodo 1
    JDn = xn' * Q * xn + un' * R * un;
    Jn = zeros(1,length(JDn));
    for k = 1:(length(JDn))
        Jn(1, k) = JDn(k, k);
    end
    J1n(1, i) = trapz(t, Jn);
    mse1 = mse1 + (J1n(1, i) - J1) ^ 2;

    % calculo do funcional pelo metodo 2
    J2n(1, i) = x0' * Pn * x0;
    mse2 = mse2 + J2n(1, i) - J2;
end

mse1 = mse1 / 100; % erro quadratico medio pelo metodo 1
mse2 = mse2 / 100; % erro quadratico medio pelo metodo 2

%% Plotando os resultados

figure(2);

subplot(2, 1, 1); hold on; grid on;
scatter(n, J1n);
plotJ1n = plot(n, ones(size(n)) * J1);
title('Funcional J por integração - comparação com o inicial');
ylabel('Funcional'); xlabel('Iteração'); 
set(plotJ1n, 'Linewidth', 1.5); 

subplot(2, 1, 2); hold on; grid on;
scatter(n, J2n);
plotJ2n = plot(n, ones(size(n)) * J2);
title('Funcional J pela equação - comparação com o inicial');
ylabel('Funcional'); xlabel('Iteração');
set(plotJ2n, 'Linewidth', 1.5);