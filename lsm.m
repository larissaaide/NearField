clc
clear all

c=3e8; %velocidade da onda eletromagnética 
fo = 60e9;
lambda = c/fo; 
deltaf = 120e3;
T = 16; %número de amostras
f = fo + linspace(-(T-1)/2,(T-1)/2, T)*deltaf;
N = 65; %número de sensores (ímpares)
dx = lambda/2; %distância entre os sensores. 
D = dx*N;
rmax = 0.62*sqrt(D^3/lambda);
rkvec = 0.5:0.5:8;
thetakvec = 0:1:180;
rk = 1; %distância da fonte ao sensor central (qualificado para o campo próximo)
thetak = pi/6; %ângulo de chegada da fonte no sensor central

posuser = rk*[cos(pi/6),sin(pi/6)];

snr_db = 20;%variação da SNR
snr_vetor= length(snr_db); 


SNR = 10^(snr_db/10); %dB para linear
sigma = 1/SNR; 
ruido = sigma*(randn(N,T) + 1j*randn(N,T))/sqrt(2); % geração do ruído

tau_los = delay_los(rk, c); %atraso refente a linha de visada 
path_loss = path_los(rk,lambda); %caminho referente a linha de visada
strVec_sphh = strVec_sph(lambda,thetak,rk, dx, N); %steering vector onda esférica
chanell_los = chanellos(N,strVec_sphh,fo,f,rk,c); %canal com atraso referente a linha de visada

s = sign(randi([0,1],1,T) - 0.5); %gerando o sinal da fonte
S = diag(s);
s_los = chanell_los*S; 
s_los_matriz = s_los/(path_loss) + ruido; %snr = 1/N0B
r_seq_matriz = (s_los_matriz/S); %sequencia piloto para cálculo do OMP 
 

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
  
dt = ones(T);
X = r_seq_matriz; 

mu = 1e-3;

nIterMax = length(dt);
w = complex(zeros(N,1),0);
e = complex(zeros(nIterMax,1),0);

    for ii = 1:nIterMax
        y = X(:,ii);
        e(ii) = dt(ii) - conj(w)'*y;
        w = w + mu*conj(e(ii))*y*1/norm(y);  %NLMS - normalized LMS
    end

w_lms = conj(w); 
y_lms = w_lms.*A; 



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

function chanel_los = chanellos(N,strVec_sph,fo,f,rk,c)
    
    delaylos = delay_los(rk, c);
    h_f = path_los(fo,rk);
    phase = 2*pi*f*delaylos + rand(1)*2*pi; 
    chanel_los = h_f*(exp(-1j*phase).*strVec_sph); 
end

% function lms = f_lms(X,dt,N,mu)
% 
%     nIterMax = length(dt);
%     w = zeros(N);
%     e = zeros(nIterMax);
% 
%     for ii = 1:nIterMax
%         y = X(:,ii);
%         e(ii) = dt(ii) - conj(w)'*y;
%         w = w + mu*conj(e(ii))*y*1/norm(y);  %NLMS - normalized LMS
% 
%     end
% end
