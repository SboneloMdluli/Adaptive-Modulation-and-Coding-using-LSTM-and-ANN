%% AMC System simulation
% Date: 21 August 2020
% Autho: Sbonelo Mdluli
% Student Number: 1101772
%%
clear
clc
warning('off')

n = 15;
k = 5;
Tx = 2;              %  number of Tx antennas

fig = figure;
grid on;
ax = fig.CurrentAxes;
hold(ax,'on');


gf_square = 3 ;%0:1:10;
a = 0.96; %0:0.1:1;
frame_size = 2000;
no_frames = 1000;
SNR = 0:2:40;
figure

enc = comm.BCHEncoder(n,k);

dec = comm.BCHDecoder(n,k);


errorRate = comm.ErrorRate;
errorCalc2 = comm.ErrorRate;
errorCalc3 = comm.ErrorRate;

qam16_error = zeros(3,21);
qam8_error = zeros(3,21);
qam4_error = zeros(3,21);

prefectPred_error = zeros(3,21);


zeta16qam = zeros(1,no_frames);
states = zeros(1,no_frames);


% keep track of how many times a system is in use
sy1num =0;
sy2num =0;
sy3num =0;

for n= 1:length(a)
    
    for m = 1:length(gf_square)
    
    for pos = 1:length(SNR)
        reset(errorRate);
        reset(errorCalc2);
        reset(errorCalc3);
        H = complex(eye(Tx,Tx))/sqrt(2); % normalised
        for bits = 1:no_frames
            H_herrm = ctranspose(H); % H^H matrix
            
            data = randi([0 1], frame_size, 1); % generate random bits
            encodedData = enc(data);
            
            QAM16 = qammod(encodedData,16,'InputType','bit');
            QAM8 = qammod(encodedData,8,'InputType','bit');
            QAM_4 = qammod(data,4,'InputType','bit');
            
            Es16 = mean(abs(QAM16).^2); % calculate energy
            Es8 = mean(abs(QAM8).^2);
            EsQAM4 = mean(abs(QAM_4).^2);
            
            X = (reshape( QAM16, Tx, [])); % convert to serial parrallel stream, 16QAM
            xq = (reshape( QAM8, Tx, [])); % 8QAM
            xpsk = (reshape( QAM_4, Tx, [])); % 4QAM
            
            
            sigSquared = (Es16*n*10^(-SNR(pos)/10))/(k*log2(16));
            sigSquared8qam = (Es8*n*10^(-SNR(pos)/10))/(k*log2(8));
            sigSquared4QAM = (EsQAM4*n*10^(-SNR(pos)/10))/(k*log2(4));
            
            sig = sqrt(sigSquared);
            sig8qam = sqrt(sigSquared8qam);
            sig4QAM = sqrt(sigSquared4QAM);
            
            N16qam = (sig/sqrt(2))*(normrnd(0,1,[Tx,1]) + 1i*normrnd(0,1,[Tx,1])); % noise matrix
            N8qam = (sig8qam/sqrt(2))*(normrnd(0,1,[Tx,1]) + 1i*normrnd(0,1,[Tx,1])); % noise matrix
            N4qam = (sig4QAM/sqrt(2))*(normrnd(0,1,[Tx,1]) + 1i*normrnd(0,1,[Tx,1])); % noise matrix
            
            W16QAM = ((H_herrm*H+sigSquared*eye(Tx,Tx))^-1)*(H_herrm);
            W8qam = ((H_herrm*H+sigSquared8qam*eye(Tx,Tx))^-1)*(H_herrm);
            W4QAM = ((H_herrm*H+sigSquared4QAM*eye(Tx,Tx))^-1)*(H_herrm);
            
            Y_hat = (X +N16qam);%16QAM
            Yhatqam = (xq +N8qam);%8QAM
            Y_hatPsk = (xpsk +N4qam); % 8PSK
            
            outRx16 = Y_hat(:); % 16QAM
            outRx8 = Yhatqam(:);% 8QAM
            outRx4 = Y_hatPsk(:);% 8QAM
            
            
            demodSignal = qamdemod(outRx16,16,'OutputType','bit'); % 16QAM
            qamdemodout = qamdemod(outRx8,8,'OutputType','bit'); % 8QAM
            pskmodout = qamdemod(outRx4,4,'OutputType','bit');
            
            received16QAM = dec(demodSignal);
            received8QAM = dec(qamdemodout);
            
            
            qam4_error(:,pos) = errorCalc3(data, pskmodout); % 4QAM
            qam16_error(:,pos) = errorRate(data, received16QAM); % 16QAM
            qam8_error(:,pos) = errorCalc2(data, received8QAM); % 8QAM
            
            GF = normrnd(0,gf_square(m),Tx,Tx) + 1i*normrnd(0,gf_square(m),Tx,Tx);
            H = a(n)*H+(1-a(n))*GF; %;
            
%             
%             if SNR(pos) == 30 %assuming perfect channel prediction
%                 zeta16qam(1,bits) = 10*log(abs(sum((H*W16QAM).^2,'all')/sum(W16QAM.^2,'all')));
%                 
%                 if zeta16qam(1,bits) < -15
%                     states(:,bits) = 1;
%                     sy1num = sy1num +1;
%                     prefectPred_error (:,pos) = errorCalc3(data, pskmodout);
%                 elseif zeta16qam(1,bits) > -5
%                     states(:,bits) = 3;
%                     sy3num = sy3num +1;
%                     prefectPred_error (:,pos) = errorCalc3(data, pskmodout);
%                 else
%                     states(:,bits) = 2;
%                     sy2num = sy2num +1;
%                     prefectPred_error (:,pos) = errorRate(data, received16QAM);
%                 end
%             end
%             
%             
%             zeta16qam(1,bits) = 10*log(abs(sum((H*W16QAM).^2,'all')/sum(W16QAM.^2,'all')));
%             
%             if zeta16qam(1,bits) < -15
%                 
%                 prefectPred_error (:,pos) = errorCalc3(data, pskmodout);
%             elseif zeta16qam(1,bits) > 0
%                 
%                 prefectPred_error (:,pos) = errorCalc3(data, pskmodout);
%             else
%                 
%                 prefectPred_error (:,pos) = errorRate(data, received16QAM);
%             end
%             
%             
            
        end
        
    end
    
    %BER curve fitting
   % QAM16 = berfit(SNR, qam16_error(1,:));
  %  QAM8 = berfit(SNR, qam8_error(1,:));
    %QAM4 = berfit(SNR, qam4_error(1,:));
   % perfect = berfit(SNR, prefectPred_error(1,:));
    
    tiledlayout(4,1)
    
    %SNR1 = SNR(:,1:length(qam16_error));
    %SNR2 = SNR(:,1:length(qam16_error));
    %SNR3 = SNR(:,1:length(qam16_error));
 %   SNR4 = SNR(:,1:length(perfect));
    
    
    semilogy(SNR,  qam4_error(1,:), 'r', ...
        SNR, qam8_error(1,:), 'g',...
        SNR, qam16_error(1,:), 'b')
    
    legend('4QAM',...
        '8QAM', ...
        '16QAM');
    xlabel('SNR (dB)');
    ylabel('BER');
    
    title('System performance over AWGN with no fading');
    xlim([SNR(1), SNR(end)]);
    
    nexttile
    
    x = 0:1:no_frames-1;
    plot(x, zeta16qam);
    title('\zeta VS frames');
    ylabel('10log(\zeta)') ;
    xlabel('frame number') ;
    
    avg = mean(zeta16qam);
    
    
    dataRate = (10^-3)*(sy1num*log2(4) + sy2num*log2(8)*(1/3) + sy3num*log2(16)*(1/3))*(2);
    % 2 antenna
    % 1/3 encoder
    % threshold is 10^-3
     %a(n)
   % scatter3(a(n),gf_square(m),dataRate,'*')
   % view(-30,10)
   % xlabel('\alpha')
   % ylabel('\sigma')
   % zlabel('Data rate')
 

    % title('Data rate as a function of \alpha and \sigma')
    % hold on;
     
     %   drawnow;
        
 
    end

end

data = vertcat(zeta16qam,states);
writematrix(data)
writematrix(states)
