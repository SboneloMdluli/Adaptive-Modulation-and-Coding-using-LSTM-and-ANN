
%% LSTM AMC simulation
% Date: 21 August 2020
% Autho: Sbonelo Mdluli
% Student Number: 1101772
%%
data = readmatrix('data.txt');

numTimeStepsTrain = floor(0.9*numel(data));

dataTrain = data(1:numTimeStepsTrain+1);
dataTest = data(numTimeStepsTrain+1:end);

no_frames = 1000;

LSTMstates = zeros(1,no_frames);

mu = mean(dataTrain);
sig = std(dataTrain);

dataTrainStandardized = (dataTrain - mu) / sig;

XTrain = dataTrainStandardized(1:end-1);
YTrain = dataTrainStandardized(2:end);

numFeatures = 1;
numResponses = 1;
numHiddenUnits = 250;

layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    fullyConnectedLayer(numResponses)
    fullyConnectedLayer(numResponses)
    regressionLayer];

options = trainingOptions('adam', ...
    'MaxEpochs',250, ...
    'GradientThreshold',1, ...
    'InitialLearnRate',0.005, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropPeriod',110, ...
    'LearnRateDropFactor',0.09, ...
    'Verbose',0, ...
    'Plots','training-progress');

net = trainNetwork(XTrain,YTrain,layers,options);

mu = mean(dataTrain);
sig = std(dataTrain);

dataTestStandardized = (dataTest - mu) / sig;
XTest = dataTestStandardized(1:end-1);
YTest = dataTest(2:end);
YPred = [];
numTimeStepsTest = numel(XTest);
for i = 1:numTimeStepsTest
    [net,YPred(:,i)] = predictAndUpdateState(net,XTest(:,i),'ExecutionEnvironment','cpu');
end

% keep track of how many times a system is in use
sy1num =0;
sy2num =0;
sy3num =0;

YPred = sig*YPred + mu;

rmse = sqrt(mean((YPred-YTest).^2))
pos =1;


for i = 1:numTimeStepsTest
            if YPred(1,i) < -15
              
                LSTMstates (:,pos) = 1;
                pos = pos+1;
                sy1num = sy1num +1;
            elseif YPred(1,i) > 0
                
                LSTMstates (:,pos) = 3;
                pos = pos+1;
                sy3num = sy3num +1;
            else
                
                LSTMstates (:,pos) = 2;
                pos = pos+1;
                sy2num = sy2num +1;
            end
end


dataRate2 = (10^-3)*(sy1num*log2(4) + sy2num*log2(8)*(1/3) + sy3num*log2(16)*(1/3))*(2);

figure
subplot(2,1,1)
plot(YTest)
hold on
plot(YPred,'.-')
hold off
legend(["Observed" "Predicted"])
ylabel("zeta")
xlabel("frames")
title("LSTM zeta predictor")

writematrix(LSTMstates)
