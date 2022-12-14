% Solve an Autoregression Problem with External Input with a NARX Neural Network
% Script generated by Neural Time Series app and then modified to construct
% a loop for training the best performed network for each compartment.
% Dependencies:
% Assessment.m      - to assess the simulation by RMSE and NSE
% Input.m,          - to load input data (tidal timeseries)
% Discharge.m       - to load discharge of 5 compartments 
%                     calculate and load cumulative volume
% Baichuan Yang, UCL

clc
clear
close all
%% load input and output data

[a,Y] = Discharge(1);     % output of training
[b,Y1] = Discharge(2);    % output of test
P = Input(1,size(Y,2));   % input of training
P1 = Input(2,size(Y1,2)); % input of test
Y = Y(5,:);               % choose compartment range from 1 to 5
Y1 = Y1(5,:);

%% loop for selecting best performed ANN
temp = -999;
for i = 1:50

disp([num2str(i),'/50']);
%%
X = tonndata(P,true,false);
T = tonndata(Y,true,false);

% Choose a Training Function
trainFcn = 'trainlm';  % Levenberg-Marquardt backpropagation

% Create a Nonlinear Autoregressive Network with External Input
inputDelays = 1:4;
feedbackDelays = 1:1;
hiddenLayerSize = 10;
net = narxnet(inputDelays,feedbackDelays,hiddenLayerSize,'open',trainFcn);

% Choose Input and Feedback Pre/Post-Processing Functions
% Settings for feedback input are automatically applied to feedback output
% For a list of all processing functions type: help nnprocess
% Customize input parameters at: net.inputs{i}.processParam
% Customize output parameters at: net.outputs{i}.processParam
net.inputs{1}.processFcns = {'mapminmax'};
net.inputs{2}.processFcns = {'mapminmax'};

% Prepare the Data for Training and Simulation
% The function PREPARETS prepares timeseries data for a particular network,
% shifting time by the minimum amount to fill input states and layer
% states. Using PREPARETS allows you to keep your original time series data
% unchanged, while easily customizing it for networks with differing
% numbers of delays, with open loop or closed loop feedback modes.
[x,xi,ai,t] = preparets(net,X,{},T);

% Setup Division of Data for Training, Validation, Testing
% For a list of all data division functions type: help nndivision

net.divideFcn = 'dividerand';  % Divide data randomly
net.divideMode = 'time';  % Divide up every sample
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;



% setting parameters
net.trainParam.epochs=1000;
ney.trainParam.goal=1e-3;
net.trainParam.showWindow = true;

% Choose a Performance Function
% For a list of all performance functions type: help nnperformance
net.performFcn = 'mse';  % Mean Squared Error

% Choose Plot Functions
% For a list of all plot functions type: help nnplot
net.plotFcns = {'plotperform','plottrainstate', 'ploterrhist', ...
    'plotregression', 'plotresponse', 'ploterrcorr', 'plotinerrcorr'};

% Train the Network

% Use parallel computing
[net,tr] = train(net,x,t,xi,ai,'useParallel','yes');

% Test the Network
y = net(x,xi,ai);
e = gsubtract(t,y);
performance = perform(net,t,y)

% Recalculate Training, Validation and Test Performance
trainTargets = gmultiply(t,tr.trainMask);
valTargets = gmultiply(t,tr.valMask);
testTargets = gmultiply(t,tr.testMask);
trainPerformance = perform(net,trainTargets,y);
valPerformance = perform(net,valTargets,y);
testPerformance = perform(net,testTargets,y);

% View the Network
view(net)

% Closed Loop Network
% Use this network to do multi-step prediction.
% The function CLOSELOOP replaces the feedback input with a direct
% connection from the output layer.
netc = closeloop(net);
netc.name = [net.name ' - Closed Loop'];
view(netc)
[xc,xic,aic,tc] = preparets(netc,X,{},T);
yc = netc(xc,xic,aic);
closedLoopPerformance = perform(net,tc,yc)


% Test
X1 = tonndata(P1,true,false);
T1 = tonndata(Y1,true,false);
[x1,xi1,ai1,t1] = preparets(netc,X1,{},T1);
y_simc = netc(x1,xi1,ai1);
[x1,xi1,ai1,t1] = preparets(net,X1,{},T1);
y_sim = net(x1,xi1,ai1);
error = gsubtract(t1,y_simc);

testperformanceC = perform(netc,t1,y_simc)
testperformance = perform(net,t1,y_sim)

%% Assessment by RMSE and NSE

[rmse,nse] = Assessment(Y1,y_simc); % assess the ANN

if max(nse) > 0.75
    ModelSavePath='.\02_SavedModel\';
    filename = [ModelSavePath,num2str(min(nse)),'model.mat'];
    model = struct('net',net,'closed_net',netc);
    save(filename,'model'); % save ANN model
    disp('Networks saved!')
end

if nse > temp
    temp = nse;
    % store the ANN model with biggest NSE
    model2 = struct('open_net',net,'closed_net',netc,'nse',nse); 
end

end





%% Visuallzation
%net = model.net;
%netc = model.closed_net;

% Training dataset
X = tonndata(P,true,false);
T = tonndata(Y,true,false);
[x,xi,ai,t] = preparets(net,X,{},T);
y = net(x,xi,ai);
[xc,xic,aic,tc] = preparets(netc,X,{},T);
yc = netc(xc,xic,aic);
closedLoopPerformance = perform(net,tc,yc)

% Test dataset
X1 = tonndata(P1,true,false);
T1 = tonndata(Y1,true,false);
[x1,xi1,ai1,t1] = preparets(netc,X1,{},T1);
y_simc = netc(x1,xi1,ai1);
[x1,xi1,ai1,t1] = preparets(net,X1,{},T1);
y_sim = net(x1,xi1,ai1);
testperformanceC = perform(netc,t1,y_simc)

yi = zeros(size(Y));    % add missing time steps
yi(:,(size(Y,2)-size(y,2))+1:end) = cell2mat(y);
yci = zeros(size(Y));   % add missing time steps
yci(:,(size(Y,2)-size(yc,2))+1:end) = cell2mat(yc);

figure()
time = 1:length(Y);
time = time*2*60/3600;   % convert to hour
plot(time,Y)
hold on
plot(time,yi)
plot(time,yci)
hold off
ylabel('Volume (m^3)');
xlabel('Time (hours)');
title('Training');
legend('boxoff')
legend('Telemac','open-ANN','closed-ANN','Fontsize',8);


yi1 = zeros(size(Y1));    % add missing time steps
yi1(:,(size(Y1,2)-size(y_sim,2))+1:end) = cell2mat(y_sim);
yci1 = zeros(size(Y1));   % add missing time steps
yci1(:,(size(Y1,2)-size(y_simc,2))+1:end) = cell2mat(y_simc);

figure()
time = 1:length(Y1);
time = time*2*60/3600;   % convert to hour
plot(time,Y1)
hold on
plot(time,yi1)
plot(time,yci1)
hold off
ylabel('Volume (m^3)');
xlabel('Time (hours)');
title('Test');
legend('boxoff')
legend('Telemac','open-ANN','closed-ANN','Fontsize',8);











