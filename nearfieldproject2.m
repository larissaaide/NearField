clear all 
close all 
clc 

%% DEFINIÇÃO DE PARÂMETROS E ENTRADAS DO NEAR FIEL - ONDA ESFÉRICA 

c=3e8; %velocidade da onda eletromagnética 
fo = 60e9;
lambda = c/fo; 
deltaf = 120e3;
T = 1024; %número de amostras
f = fo + linspace(-(T-1)/2,(T-1)/2, T)*deltaf;
N = 65; %número de sensores (ímpares)
dx = lambda/2; %distância entre os sensores. 
D = dx*N;
rmax = 0.62*sqrt(D^3/lambda);


%% CÁLCULO DO RUÍDO
snr_db = 30; %snr em dB
SNR = 10^(snr_db/10); %dB para linear
sigma = 1/sqrt(SNR); 
ruido = sigma*(randn(N,T) + 1j*randn(N,T))/sqrt(2); % geração do ruído

rkvec = 0.5:0.1:8;
thetakvec = 0:0.1:pi;

% GERANDO O CANAL EM PATH LOS
%rk = rand()*(7.5) + 0.5; % distância da fonte ao sensor central (qualificado para o campo próximo)
%thetak = rand()*pi; % ângulo de chegada da fonte no sensor central
rk = 1;
thetak = 0.5;
posUser = [cos(pi/6),sin(pi/6)]*rk;
tau_los = delay_los(rk, c); % atraso refente a linha de visada 
path_loss = path_los(rk,lambda); % caminho referente a linha de visada
strVec_sphh = strVec_sph(lambda,thetak,rk, dx, N); %steering vector onda esférica
chanell_los = chanellos(strVec_sphh,fo,f,rk,c); %canal com atraso referente a linha de visada

s = sign(randi([0,1],1,T) - 0.5); % gerando o sinal da fonte

S = diag(s);

s_los = chanell_los*S; 
s_los = s_los/(path_loss) + ruido; %snr = 1/N0B

% s_los = s.*chanell_los; 
% s_los = s_los/norm(s_los); %snr = 1/N0B
%y = strVec_sphh*s_los + ruido; % sinal recebido

r_seq = (s_los/S); % sequencia piloto para cálculo do OMP

%r_seq = (y*s')/norm(s); % sequencia piloto para cálculo do OMP

idx = zeros(length(rkvec)*length(thetakvec),2);
A = zeros(N,length(rkvec)*length(thetakvec));
l = 1;


for i = 1:length(rkvec)
    for j = 1:length(thetakvec)
       idx(l,:) = [rkvec(i),thetakvec(j)]; 
       A(:,l)   = strVec_sph(lambda,thetakvec(j),rkvec(i),dx,N);
       l = l+1;
    end
end


[coeff,dictatom,atomidx,errnorm] = ompdecomp(r_seq,A,'MaxSparsity',1);
val_estimado = idx(atomidx,:);
angulo_estimado = val_estimado(1,2);
rk_estimado = val_estimado(1,1);
erro = norm(rk_estimado-rk);
posUser_est = [cos(angulo_estimado),sin(angulo_estimado)]*rk_estimado;
error_dist = norm(posUser - posUser_est);


%% Funções 
function A = strVec_sph(lambda,thetak,rk, dx, N)
    
   
    for n=1:N
        cent_dist = dx*(-(N-1)/2+(n-1));
        rmk(n)=sqrt(rk + cent_dist^2 - 2*cent_dist*rk*sin(thetak));
        phase(n) = 2*pi*dx/lambda*(rmk(n) - rk);
    end 
    
    A = exp(1j*phase).';
    
end

%Função de calcula o atraso do sinal Tau = d/c, em que d é a distancia do usuário ate a antena
function tau = delay_los(rk, c)
    tau = rk/c; 
end

%função que calcula a perda de caminho (redução na densidade de potência)
function pl = path_los(rk,lambda)
    pl = lambda/(4*pi*rk);
end

function chanel_los = chanellos(strVec_sph,fo,f,rk,c)
    
    delaylos = delay_los(rk, c);
    h_f = path_los(fo,rk);
    phase = 2*pi*f*delaylos + rand(1)*2*pi; 
    chanel_los = strVec_sph*h_f*exp(-1j*phase); 
    
end

function chanel_2raios = chanell2raios(strVec_sph,strVec_sph_nlos,ratio_db,fo,f,rk,c, delta_tau)
    
    delaylos = delay_los(rk, c);
    delaynlos = delaylos + delta_tau;
    ratio = 10^(ratio_db/10);
    
    h_f = path_los(fo,rk);
    h_fnlos = sqrt(ratio)*h_f;
    phase = 2*pi*f*delaylos + rand(1)*2*pi; 
    chanel_los = strVec_sph*h_f*exp(-1j*phase); 
    phase_nlos = 2*pi*f*delaynlos + rand(1)*2*pi; 
    chanel_nlos = strVec_sph_nlos*h_fnlos*exp(-1j*phase_nlos); 
    
    chanel_2raios = chanel_los + chanel_nlos;
    
    
end

