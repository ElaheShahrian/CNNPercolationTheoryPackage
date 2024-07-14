clear; close all; clc;
%% Loading Data
load CNN64.mat

CVModelPart = cvpartition(length(yt),'Holdout',0.1) ; % cross validation
indexTrain = training(CVModelPart) ;
indexTest = test(CVModelPart) ;

% Train
XTrain = xt(:,:,:,indexTrain) ; % 4D (64, 66, 1, 5040)
YTrain = yt(indexTrain,1) ; % 1D (5040, 1)
PTrain = pt(indexTrain,1) ; % 1D (5040, 1)

% Test
XTest = xt(:,:,:,indexTest) ; % 4D (64, 66, 1, 560)
YTest = yt(indexTest,1) ; % 1D (560, 1)
PTest = pt(indexTest,1) ; % 1D (560, 1)

% for i = 1:size(XT,4)
%     xt(:,:,1,i) = XT(:,:,1,i);
%     xt(:,:,2,i) = MyChannalGeneration(XT(:,:,1,i),4);
%     xt(:,:,3,i) = MyChannalGeneration(XT(:,:,1,i),8);
%     Conf(i,1) = mean(xt(:,:,:,i),"all");
% end
% XT = xt ;
% load NewYT.mat
%YT = NewYT ;
%NewYT = 26.03*YT - 3.082 ;

% plot(YT,NewYT,"*r")

% XT = cat(4,XT,XT) ;
% YT = [YT;YT];

% for i = 1:size(XT,4)
%     Poro(i) = nnz(XT(:,:,1,i))/1088 ;
% end
% YT = Poro' ;

% XTrain = XT(:,:,:,1:800) ;
% YTrain = YT(1:800);
% 
% 
% XTest = XT(:,:,:,50:100) ;
% YTest = YT(50:100);


%%
layers = [
    imageInputLayer([64 66 1])
    convolution2dLayer(3,16,'Padding','same')
    batchNormalizationLayer
    reluLayer
    averagePooling2dLayer(2,'Stride',2)
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    averagePooling2dLayer(2,'Stride',2)
    convolution2dLayer(3,64,'Padding','same')
    batchNormalizationLayer
    reluLayer
    convolution2dLayer(3,32,'Padding','same')
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];



% filterSize_init = 10 ;
% numberFilter_init = 32 ;
% 
% filterSize_2nd = 1 ;
% 
% W = ones([filterSize_init,filterSize_init,1,numberFilter_init]) ;
% 
% 
% layers = [
%     imageInputLayer([32 34 1])
%     convolution2dLayer(3,8,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
%     averagePooling2dLayer(2,'Stride',2)
%     convolution2dLayer(3,16,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
%     averagePooling2dLayer(2,'Stride',2)
%     convolution2dLayer(3,32,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
%     convolution2dLayer(3,32,'Padding','same')
%     batchNormalizationLayer
%     reluLayer
%     fullyConnectedLayer(1)
%     dropoutLayer(0.2)
%     regressionLayer];



% layers = [imageInputLayer([32 34 1])
% 		convolution2dLayer(10,10,'Padding',0,'Stride',1)
% 		reluLayer
% 		batchNormalizationLayer
% 		maxPooling2dLayer(2)
% 		
% 		convolution2dLayer(7,20,'Padding',0,'Stride',1)
% 		reluLayer
% 		batchNormalizationLayer
% 		maxPooling2dLayer(2)
% 		
% 		convolution2dLayer(5,40,'Padding',0,'Stride',1)
% 		reluLayer
% 		batchNormalizationLayer
% 		maxPooling2dLayer(2)
% 		
% 		convolution2dLayer(3,80,'Padding',0,'Stride',1)
% 		reluLayer
% 		batchNormalizationLayer
% 		maxPooling2dLayer(2)
% 		
% 		convolution2dLayer(2,160,'Padding',0,'Stride',1)
% 		reluLayer
% 		batchNormalizationLayer
% 		maxPooling2dLayer(2)
% 		
% 		convolution2dLayer(2,400,'Padding',0,'Stride',1)
% 		reluLayer
% 		batchNormalizationLayer
% 		maxPooling2dLayer(2)
% 		
% 		
% 		fullyConnectedLayer(10)
% 		tanhLayer
% 		
% 		fullyConnectedLayer(1)
% 		...transposedConv2dLayer ... or ???
% 		
% 		regressionLayer];

%% Iterations per epoch = Number of training samples ÷ MiniBatchSize
Iterations_Per_Epoch = 1 ;
miniBatchSize  = ceil((numel(YTrain)/Iterations_Per_Epoch)/5) ;
validationFrequency = 1;
options = trainingOptions('sgdm', ...
    'MiniBatchSize',miniBatchSize, ...
    'MaxEpochs',20, ...
    'InitialLearnRate',1e-3, ...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',15, ...
    ...'Shuffle','every-epoch', ...
    'ValidationData',{XTest,YTest}, ...
    'ValidationFrequency',validationFrequency, ...
    'Plots','training-progress', ...
    ...'Verbose',true,...
    'OutputNetwork','best-validation-loss');
%%
net = trainNetwork(XTrain,YTrain,layers,options);
%%
YPredicted = double(predict(net,xt));
Error = sqrt(mean((yt - YPredicted).^2)) ; % MAPE
R2 = 1 - (sum((yt-YPredicted).^2)/sum((yt-mean(yt)).^2)) ;
%
figure
scatter(yt,YPredicted,'*b')
hold on
plot([-1 10], [-1 10],'r--') ; title(['RMSE= ' num2str(Error) ' & R^{2}= ' num2str(R2)])
axis tight
xlabel("Real Target")
ylabel("Predict Target")
%%
