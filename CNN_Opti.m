function  [bestIdx,fileName] = CNN_Opti(AnimalID,trainData,trainLabel,testData,testLabel)
%% optimization
optimVars = [
    optimizableVariable('InitialLearnRate',[1e-2 0.5],'Transform','log')
    optimizableVariable('Momentum',[0.8 0.9])
    optimizableVariable('L2Regularization',[1e-10 1e-2],'Transform','log')];
ObjFcn = makeObjFcn(AnimalID,trainData,trainLabel,testData,testLabel);
BayesObject = bayesopt(ObjFcn,optimVars, ...
            'MaxObjectiveEvaluations',30,...
            'MaxTime',2*60*60, ...
            'IsObjectiveDeterministic',false, ...
            'UseParallel',false);

function ObjFcn = makeObjFcn(AnimalID,trainData,trainLabel,testData,testLabel)

ObjFcn = @valErrorFun; 
    function [valError,cons,fileName] = valErrorFun(optVars)
    
        layers = [
        sequenceInputLayer(size(trainData{1, 1},1),"MinLength",size(trainData{1, 1},2))
        convolution1dLayer(32,64,"Padding","same")
        maxPooling1dLayer(2)
        reluLayer
    
        convolution1dLayer(32,128,"Padding","same")
             maxPooling1dLayer(2)
    
        reluLayer
    
        convolution1dLayer(32,128,"Padding","same")
        maxPooling1dLayer(2)
    
        reluLayer
        
        lstmLayer(128)
        dropoutLayer(0.5)
        lstmLayer(128,"OutputMode","last")

        fullyConnectedLayer(100)
        reluLayer
        dropoutLayer(0.5)
        fullyConnectedLayer(3)
        softmaxLayer
        classificationLayer];
    
            miniBatchSize = 256;
            validationFrequency = floor(numel(trainLabel)/miniBatchSize);
            options = trainingOptions('sgdm', ...
                'ExecutionEnvironment','gpu', ...
                'InitialLearnRate',optVars.InitialLearnRate, ...
                'Momentum',optVars.Momentum, ...
                'MaxEpochs',60, ...
                'LearnRateSchedule','piecewise', ...
                'LearnRateDropPeriod',40, ...
                'LearnRateDropFactor',0.1, ...
                'MiniBatchSize',miniBatchSize, ...
                'L2Regularization',optVars.L2Regularization, ...
                'Shuffle','every-epoch', ...
                'Verbose',false, ...
                'Plots','training-progress', ...
                'ValidationData',{testData,testLabel}, ...
                'ValidationFrequency',validationFrequency);
    
            trainedNet = trainNetwork(trainData,trainLabel,layers,options);
            close(findall(groot,'Tag','NNET_CNN_TRAININGPLOT_UIFIGURE'))
    
            YPredicted = classify(trainedNet,testData);
            valError = 1 - mean(YPredicted == testLabel);
    
            Nparts = [char(AnimalID)  num2str(valError) + ".mat"];
            fileName = Nparts(1)+"_"+Nparts(2);
            save(fileName,'trainedNet','valError','options')
            cons = [];
       end
end
bestIdx = BayesObject.IndexOfMinimumTrace(end);
fileName = BayesObject.UserDataTrace{bestIdx};
savedStruct = load(fileName);
valError = savedStruct.valError
end
