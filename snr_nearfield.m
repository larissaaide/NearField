clear all 
close all 
clc 

%% DEFINIÇÃO DE PARÂMETROS E ENTRADAS DO NEAR FIEL - ONDA ESFÉRICA 

c=3e8; % velocidade da onda eletromagnética 
fo = 60e9;
lambda = c/fo; 
deltaf = 120e3;
T = 1024; %número de amostras
f = fo + linspace(-(T-1)/2,(T-1)/2, T)*deltaf;
N = 65; % número de sensores (ímpares)
dx = lambda/2; % distância entre os sensores. 
D = dx*N;
rmax = 0.62*sqrt(D^3/lambda);
rkvec = 0.5:0.5:8;
thetakvec = 0:1:180;
rk = 1; % distância da fonte ao sensor central (qualificado para o campo próximo)
thetak = pi/6; % ângulo de chegada da fonte no sensor central


%% CÁLCULO DO RUÍDO

snr_db = 0:1:20;%variação da SNR
snr_vetor= length(snr_db); 
MC = 1:1:1; %número de simulações monte carlo
MC_vetor= length(MC); 

s_los_matriz = zeros(N,T,snr_vetor,MC_vetor); 
r_seq_matriz = zeros(N,T,snr_vetor,MC_vetor);
val_estimado_matriz = zeros(snr_vetor,2,MC_vetor);
atomidx_matriz = zeros(1,snr_vetor,MC_vetor,MC_vetor);
angulo_estimado_matriz = zeros(snr_vetor,MC_vetor);
rk_estimado_matriz = zeros(snr_vetor,MC_vetor);
erro_matriz = zeros(snr_vetor,MC_vetor);
%erro_matriz = norm(rk_estimado-rk);

var_mc = 1;

for mc = 1:MC_vetor

    var_snr = 1;

    for snr = snr_db

        SNR = 10^(snr/10); %dB para linear
        sigma = 1/SNR; 
        ruido = sigma*(randn(N,T) + 1j*randn(N,T))/sqrt(2); % geração do ruído

        tau_los = delay_los(rk, c); % atraso refente a linha de visada 
        path_loss = path_los(rk,lambda); % caminho referente a linha de visada
        strVec_sphh = strVec_sph(lambda,thetak,rk, dx, N); %steering vector onda esférica
        chanell_los = chanellos(N,strVec_sphh,fo,f,rk,c); %canal com atraso referente a linha de visada

        s = sign(randi([0,1],1,T) - 0.5); % gerando o sinal da fonte
    %    s_los = s.*chanell_los; 
    %    s_los = s_los/norm(s_los); %snr = 1/N0B
    %    y = strVec_sphh.*s_los + ruido; % sinal recebido
    %    r_seq = (y*s_los')/norm(s_los); % sequencia piloto para cálculo do OMP

        S = diag(s);
        s_los = chanell_los*S; 
        s_los_matriz(:,:,var_snr,var_mc) = s_los/(path_loss) + ruido; %snr = 1/N0B
        r_seq_matriz(:,:,var_snr,var_mc) = (s_los_matriz(:,:,var_snr,var_mc)/S); % sequencia piloto para cálculo do OMP 

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
        
  %Esparsidade através do OMP      
        [coeff,dictatom,atomidx_matriz(:,var_snr,var_mc),errnorm] = ompdecomp(r_seq_matriz(:,:,var_snr,var_mc),A,'MaxSparsity',1);
        val_estimado_matriz(var_snr,:,var_mc) = idx(atomidx_matriz(:,var_snr,var_mc),:);
        angulo_estimado_matriz(var_snr,var_mc) = val_estimado_matriz(var_snr,2,var_mc);
        rk_estimado_matriz(var_snr,var_mc) = val_estimado_matriz(var_snr,1,var_mc);
        erro_matriz(var_snr,var_mc) = norm(rk_estimado_matriz(var_snr,var_mc)-rk);


        var_snr=var_snr+1;
        
%Esparsidade através do CVX        
    
    
%     cvx_begin
%         variable x(n)
%         minimize(norm( A*b - r_seq_matriz(:,:,var_snr,var_mc), 2))
%         subject to
%             norm(b,1 <=) * x == d
%             norm( x, Inf ) <= e
%     cvx_end
%     
    end

    var_mc=var_mc+1;
end

erro_estimadorr = zeros(snr_vetor,1);
  
    for coluna_estimadorr=1:snr_vetor   
        for linhar=1:MC_vetor 
           erro_estimadorr(coluna_estimadorr) = erro_estimadorr(coluna_estimadorr) + erro_matriz(coluna_estimadorr,linhar);     
        end
    end
  
 media_erro_estimado = zeros(snr_vetor,1);  
 media_erro_estimado = erro_estimadorr./MC_vetor;
 
figure(5)
semilogy(snr_db, media_erro_estimado, 'r','linewidth',1);
hold on;
title('MSE vs Eb/N0');
ylabel('ERRO MÉDIO QUADRÁTICO NORMALIZADO - MSE (dB)');
xlabel('SNR(dB)');
legend ('Erro estimado do canal resultante', ...
'location', 'southwest');
grid on 


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


%função do canal NLOS

