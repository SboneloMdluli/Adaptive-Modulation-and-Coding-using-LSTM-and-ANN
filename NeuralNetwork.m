%% NN AMC simulation
% Date: 21 August 2020
% Autho: Sbonelo Mdluli
% Student Number: 1101772
%%
data = readmatrix('data.txt');

frame = 0:1000;

y = num2cell(data(1:900));

ftdnn_net = timedelaynet([1:4],10);
ftdnn_net.trainParam.epochs = 1000;
ftdnn_net.divideFcn = '';

p = y(5:end);
t = y(5:end);
Pi=y(1:4);
ftdnn_net = train(ftdnn_net,p,t,Pi);

yp = ftdnn_net(p,Pi);
e = gsubtract(yp,t);
rmse = sqrt(mse(e))

% keep track of how many times a system is in use
sy1num =0;
sy2num =0;
sy3num =0;


pos =1;
frameSize = 1000;


yp = cell2mat(yp);
NNstates = zeros(1,length(yp));

for i = 1:length(yp)
            if yp(1,i) < -15
              
                NNstates (:,pos) = 1;
                pos = pos+1;
                sy1num = sy1num +1;
            elseif yp(1,i) > 0
                
                NNstates (:,pos) = 3;
                pos = pos+1;
                sy3num = sy3num +1;
            else
                
                NNstates (:,pos) = 2;
                pos = pos+1;
                sy2num = sy2num +1;
            end
end

dataRate1 = (10^-3)*(sy1num*log2(4) + sy2num*log2(8)*(1/3) + sy3num*log2(16)*(1/3))*(2);

writematrix(NNstates)


