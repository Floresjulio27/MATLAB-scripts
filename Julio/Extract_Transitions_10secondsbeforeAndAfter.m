clear all;

load AnalysisResults.mat;

%% TRANSITIONS EXTRACTION

%PFC
dataJF037_AwaketoNREM_PFC = AnalysisResults.JF037.Transitions.AWAKEtoNREM.GRAB_NE;
dataJF038_AwaketoNREM_PFC = AnalysisResults.JF038.Transitions.AWAKEtoNREM.GRAB_NE;
dataJF039_AwaketoNREM_PFC = AnalysisResults.JF039.Transitions.AWAKEtoNREM.GRAB_NE;
dataJF040_AwaketoNREM_PFC= AnalysisResults.JF040.Transitions.AWAKEtoNREM.GRAB_NE;
dataJF048_AwaketoNREM_PFC = AnalysisResults.JF048.Transitions.AWAKEtoNREM.GRAB_NE;
dataJF050_AwaketoNREM_PFC = AnalysisResults.JF050.Transitions.AWAKEtoNREM.GRAB_NE;

dataJF037_NREMtoAwake_PFC = AnalysisResults.JF037.Transitions.NREMtoAWAKE.GRAB_NE;
dataJF038_NREMtoAwake_PFC = AnalysisResults.JF038.Transitions.NREMtoAWAKE.GRAB_NE;
dataJF039_NREMtoAwake_PFC = AnalysisResults.JF039.Transitions.NREMtoAWAKE.GRAB_NE;
dataJF040_NREMtoAwake_PFC = AnalysisResults.JF040.Transitions.NREMtoAWAKE.GRAB_NE;
dataJF048_NREMtoAwake_PFC = AnalysisResults.JF048.Transitions.NREMtoAWAKE.GRAB_NE;
dataJF050_NREMtoAwake_PFC = AnalysisResults.JF050.Transitions.NREMtoAWAKE.GRAB_NE;

dataJF037_NREMtoREM_PFC = AnalysisResults.JF037.Transitions.NREMtoREM.GRAB_NE;
dataJF038_NREMtoREM_PFC = AnalysisResults.JF038.Transitions.NREMtoREM.GRAB_NE;
dataJF039_NREMtoREM_PFC = AnalysisResults.JF039.Transitions.NREMtoREM.GRAB_NE;
dataJF040_NREMtoREM_PFC = AnalysisResults.JF040.Transitions.NREMtoREM.GRAB_NE;
dataJF048_NREMtoREM_PFC = AnalysisResults.JF048.Transitions.NREMtoREM.GRAB_NE;
dataJF050_NREMtoREM_PFC = AnalysisResults.JF050.Transitions.NREMtoREM.GRAB_NE;

dataJF037_REMtoAwake_PFC = AnalysisResults.JF037.Transitions.REMtoAWAKE.GRAB_NE;
dataJF038_REMtoAwake_PFC = AnalysisResults.JF038.Transitions.REMtoAWAKE.GRAB_NE;
dataJF039_REMtoAwake_PFC = AnalysisResults.JF039.Transitions.REMtoAWAKE.GRAB_NE;
dataJF040_REMtoAwake_PFC = AnalysisResults.JF040.Transitions.REMtoAWAKE.GRAB_NE;
dataJF048_REMtoAwake_PFC = AnalysisResults.JF048.Transitions.REMtoAWAKE.GRAB_NE;
dataJF050_REMtoAwake_PFC = AnalysisResults.JF050.Transitions.REMtoAWAKE.GRAB_NE;

%S1BF

dataJF037_AwaketoNREM_S1BF = AnalysisResults.JF037.Transitions.AWAKEtoNREM.GRAB_ACh;
dataJF038_AwaketoNREM_S1BF = AnalysisResults.JF038.Transitions.AWAKEtoNREM.GRAB_ACh;
dataJF039_AwaketoNREM_S1BF = AnalysisResults.JF039.Transitions.AWAKEtoNREM.GRAB_ACh;
dataJF040_AwaketoNREM_S1BF = AnalysisResults.JF040.Transitions.AWAKEtoNREM.GRAB_ACh;
dataJF048_AwaketoNREM_S1BF = AnalysisResults.JF048.Transitions.AWAKEtoNREM.GRAB_ACh;
dataJF050_AwaketoNREM_S1BF = AnalysisResults.JF050.Transitions.AWAKEtoNREM.GRAB_ACh;

dataJF037_NREMtoAwake_S1BF = AnalysisResults.JF037.Transitions.NREMtoAWAKE.GRAB_ACh;
dataJF038_NREMtoAwake_S1BF = AnalysisResults.JF038.Transitions.NREMtoAWAKE.GRAB_ACh;
dataJF039_NREMtoAwake_S1BF = AnalysisResults.JF039.Transitions.NREMtoAWAKE.GRAB_ACh;
dataJF040_NREMtoAwake_S1BF = AnalysisResults.JF040.Transitions.NREMtoAWAKE.GRAB_ACh;
dataJF048_NREMtoAwake_S1BF = AnalysisResults.JF048.Transitions.NREMtoAWAKE.GRAB_ACh;
dataJF050_NREMtoAwake_S1BF = AnalysisResults.JF050.Transitions.NREMtoAWAKE.GRAB_ACh;

dataJF037_NREMtoREM_S1BF = AnalysisResults.JF037.Transitions.NREMtoREM.GRAB_ACh;
dataJF038_NREMtoREM_S1BF = AnalysisResults.JF038.Transitions.NREMtoREM.GRAB_ACh;
dataJF039_NREMtoREM_S1BF = AnalysisResults.JF039.Transitions.NREMtoREM.GRAB_ACh;
dataJF040_NREMtoREM_S1BF = AnalysisResults.JF040.Transitions.NREMtoREM.GRAB_ACh;
dataJF048_NREMtoREM_S1BF = AnalysisResults.JF048.Transitions.NREMtoREM.GRAB_ACh;
dataJF050_NREMtoREM_S1BF = AnalysisResults.JF050.Transitions.NREMtoREM.GRAB_ACh;

dataJF037_REMtoAwake_S1BF = AnalysisResults.JF037.Transitions.REMtoAWAKE.GRAB_ACh;
dataJF038_REMtoAwake_S1BF = AnalysisResults.JF038.Transitions.REMtoAWAKE.GRAB_ACh;
dataJF039_REMtoAwake_S1BF = AnalysisResults.JF039.Transitions.REMtoAWAKE.GRAB_ACh;
dataJF040_REMtoAwake_S1BF = AnalysisResults.JF040.Transitions.REMtoAWAKE.GRAB_ACh;
dataJF048_REMtoAwake_S1BF = AnalysisResults.JF048.Transitions.REMtoAWAKE.GRAB_ACh;
dataJF050_REMtoAwake_S1BF = AnalysisResults.JF050.Transitions.REMtoAWAKE.GRAB_ACh;

% Calculating the mean for columns 601 to 900. This corresponds to the 10s
% the transition 

%PFC

AwaketoNREM_PFC_mean_before_JF037 = mean(dataJF037_AwaketoNREM_PFC(1:900)); 
AwaketoNREM_PFC_mean_before_JF038 = mean(dataJF038_AwaketoNREM_PFC(1:900)); 
AwaketoNREM_PFC_mean_before_JF039 = mean(dataJF039_AwaketoNREM_PFC(1:900)); 
AwaketoNREM_PFC_mean_before_JF040 = mean(dataJF040_AwaketoNREM_PFC(1:900)); 
AwaketoNREM_PFC_mean_before_JF048 = mean(dataJF048_AwaketoNREM_PFC(1:900)); 
AwaketoNREM_PFC_mean_before_JF050 = mean(dataJF050_AwaketoNREM_PFC(1:900)); 

NREMtoAwake_PFC_mean_before_JF037 = mean(dataJF037_NREMtoAwake_PFC(1:900)); 
NREMtoAwake_PFC_mean_before_JF038 = mean(dataJF038_NREMtoAwake_PFC(1:900)); 
NREMtoAwake_PFC_mean_before_JF039 = mean(dataJF039_NREMtoAwake_PFC(1:900)); 
NREMtoAwake_PFC_mean_before_JF040 = mean(dataJF040_NREMtoAwake_PFC(1:900)); 
NREMtoAwake_PFC_mean_before_JF048 = mean(dataJF048_NREMtoAwake_PFC(1:900)); 
NREMtoAwake_PFC_mean_before_JF050 = mean(dataJF050_NREMtoAwake_PFC(1:900));

NREMtoREM_PFC_mean_before_JF037 = mean(dataJF037_NREMtoREM_PFC(1:900));
NREMtoREM_PFC_mean_before_JF038 = mean(dataJF038_NREMtoREM_PFC(1:900));
NREMtoREM_PFC_mean_before_JF039 = mean(dataJF039_NREMtoREM_PFC(1:900));
NREMtoREM_PFC_mean_before_JF040 = mean(dataJF040_NREMtoREM_PFC(1:900));
NREMtoREM_PFC_mean_before_JF048 = mean(dataJF048_NREMtoREM_PFC(1:900));
NREMtoREM_PFC_mean_before_JF050 = mean(dataJF050_NREMtoREM_PFC(1:900));

REMtoAwake_PFC_mean_before_JF037 = mean(dataJF037_REMtoAwake_PFC(1:900));
REMtoAwake_PFC_mean_before_JF038 = mean(dataJF038_REMtoAwake_PFC(1:900));
REMtoAwake_PFC_mean_before_JF039 = mean(dataJF039_REMtoAwake_PFC(1:900));
REMtoAwake_PFC_mean_before_JF040 = mean(dataJF040_REMtoAwake_PFC(1:900));
REMtoAwake_PFC_mean_before_JF048 = mean(dataJF048_REMtoAwake_PFC(1:900));
REMtoAwake_PFC_mean_before_JF050 = mean(dataJF050_REMtoAwake_PFC(1:900));

%S1BF

AwaketoNREM_S1BF_mean_before_JF037 = mean(dataJF037_AwaketoNREM_S1BF(1:900)); 
AwaketoNREM_S1BF_mean_before_JF038 = mean(dataJF038_AwaketoNREM_S1BF(1:900)); 
AwaketoNREM_S1BF_mean_before_JF039 = mean(dataJF039_AwaketoNREM_S1BF(1:900)); 
AwaketoNREM_S1BF_mean_before_JF040 = mean(dataJF040_AwaketoNREM_S1BF(1:900)); 
AwaketoNREM_S1BF_mean_before_JF048 = mean(dataJF048_AwaketoNREM_S1BF(1:900)); 
AwaketoNREM_S1BF_mean_before_JF050 = mean(dataJF050_AwaketoNREM_S1BF(1:900)); 

NREMtoAwake_S1BF_mean_before_JF037 = mean(dataJF037_NREMtoAwake_S1BF(1:900)); 
NREMtoAwake_S1BF_mean_before_JF038 = mean(dataJF038_NREMtoAwake_S1BF(1:900)); 
NREMtoAwake_S1BF_mean_before_JF039 = mean(dataJF039_NREMtoAwake_S1BF(1:900)); 
NREMtoAwake_S1BF_mean_before_JF040 = mean(dataJF040_NREMtoAwake_S1BF(1:900)); 
NREMtoAwake_S1BF_mean_before_JF048 = mean(dataJF048_NREMtoAwake_S1BF(1:900)); 
NREMtoAwake_S1BF_mean_before_JF050 = mean(dataJF050_NREMtoAwake_S1BF(1:900));

NREMtoREM_S1BF_mean_before_JF037 = mean(dataJF037_NREMtoREM_S1BF(1:900));
NREMtoREM_S1BF_mean_before_JF038 = mean(dataJF038_NREMtoREM_S1BF(1:900));
NREMtoREM_S1BF_mean_before_JF039 = mean(dataJF039_NREMtoREM_S1BF(1:900));
NREMtoREM_S1BF_mean_before_JF040 = mean(dataJF040_NREMtoREM_S1BF(1:900));
NREMtoREM_S1BF_mean_before_JF048 = mean(dataJF048_NREMtoREM_S1BF(1:900));
NREMtoREM_S1BF_mean_before_JF050 = mean(dataJF050_NREMtoREM_S1BF(1:900));

REMtoAwake_S1BF_mean_before_JF037 = mean(dataJF037_REMtoAwake_S1BF(1:900));
REMtoAwake_S1BF_mean_before_JF038 = mean(dataJF038_REMtoAwake_S1BF(1:900));
REMtoAwake_S1BF_mean_before_JF039 = mean(dataJF039_REMtoAwake_S1BF(1:900));
REMtoAwake_S1BF_mean_before_JF040 = mean(dataJF040_REMtoAwake_S1BF(1:900));
REMtoAwake_S1BF_mean_before_JF048 = mean(dataJF048_REMtoAwake_S1BF(1:900));
REMtoAwake_S1BF_mean_before_JF050 = mean(dataJF050_REMtoAwake_S1BF(1:900));

% Calculating the mean for columns 901 to 1200. This corresponds to the 10
% seconds after the transition
%PFC

dataJF037_AwaketoNREM_PFC = AnalysisResults.JF037.Transitions.AWAKEtoNREM.GRAB_NE;
dataJF038_AwaketoNREM_PFC = AnalysisResults.JF038.Transitions.AWAKEtoNREM.GRAB_NE;
dataJF039_AwaketoNREM_PFC = AnalysisResults.JF039.Transitions.AWAKEtoNREM.GRAB_NE;
dataJF040_AwaketoNREM_PFC= AnalysisResults.JF040.Transitions.AWAKEtoNREM.GRAB_NE;
dataJF048_AwaketoNREM_PFC = AnalysisResults.JF048.Transitions.AWAKEtoNREM.GRAB_NE;
dataJF050_AwaketoNREM_PFC = AnalysisResults.JF050.Transitions.AWAKEtoNREM.GRAB_NE;

dataJF037_NREMtoAwake_PFC = AnalysisResults.JF037.Transitions.NREMtoAWAKE.GRAB_NE;
dataJF038_NREMtoAwake_PFC = AnalysisResults.JF038.Transitions.NREMtoAWAKE.GRAB_NE;
dataJF039_NREMtoAwake_PFC = AnalysisResults.JF039.Transitions.NREMtoAWAKE.GRAB_NE;
dataJF040_NREMtoAwake_PFC = AnalysisResults.JF040.Transitions.NREMtoAWAKE.GRAB_NE;
dataJF048_NREMtoAwake_PFC = AnalysisResults.JF048.Transitions.NREMtoAWAKE.GRAB_NE;
dataJF050_NREMtoAwake_PFC = AnalysisResults.JF050.Transitions.NREMtoAWAKE.GRAB_NE;

dataJF037_NREMtoREM_PFC = AnalysisResults.JF037.Transitions.NREMtoREM.GRAB_NE;
dataJF038_NREMtoREM_PFC = AnalysisResults.JF038.Transitions.NREMtoREM.GRAB_NE;
dataJF039_NREMtoREM_PFC = AnalysisResults.JF039.Transitions.NREMtoREM.GRAB_NE;
dataJF040_NREMtoREM_PFC = AnalysisResults.JF040.Transitions.NREMtoREM.GRAB_NE;
dataJF048_NREMtoREM_PFC = AnalysisResults.JF048.Transitions.NREMtoREM.GRAB_NE;
dataJF050_NREMtoREM_PFC = AnalysisResults.JF050.Transitions.NREMtoREM.GRAB_NE;

dataJF037_REMtoAwake_PFC = AnalysisResults.JF037.Transitions.REMtoAWAKE.GRAB_NE;
dataJF038_REMtoAwake_PFC = AnalysisResults.JF038.Transitions.REMtoAWAKE.GRAB_NE;
dataJF039_REMtoAwake_PFC = AnalysisResults.JF039.Transitions.REMtoAWAKE.GRAB_NE;
dataJF040_REMtoAwake_PFC = AnalysisResults.JF040.Transitions.REMtoAWAKE.GRAB_NE;
dataJF048_REMtoAwake_PFC = AnalysisResults.JF048.Transitions.REMtoAWAKE.GRAB_NE;
dataJF050_REMtoAwake_PFC = AnalysisResults.JF050.Transitions.REMtoAWAKE.GRAB_NE;

%S1BF

dataJF037_AwaketoNREM_S1BF = AnalysisResults.JF037.Transitions.AWAKEtoNREM.GRAB_ACh;
dataJF038_AwaketoNREM_S1BF = AnalysisResults.JF038.Transitions.AWAKEtoNREM.GRAB_ACh;
dataJF039_AwaketoNREM_S1BF = AnalysisResults.JF039.Transitions.AWAKEtoNREM.GRAB_ACh;
dataJF040_AwaketoNREM_S1BF = AnalysisResults.JF040.Transitions.AWAKEtoNREM.GRAB_ACh;
dataJF048_AwaketoNREM_S1BF = AnalysisResults.JF048.Transitions.AWAKEtoNREM.GRAB_ACh;
dataJF050_AwaketoNREM_S1BF = AnalysisResults.JF050.Transitions.AWAKEtoNREM.GRAB_ACh;

dataJF037_NREMtoAwake_S1BF = AnalysisResults.JF037.Transitions.NREMtoAWAKE.GRAB_ACh;
dataJF038_NREMtoAwake_S1BF = AnalysisResults.JF038.Transitions.NREMtoAWAKE.GRAB_ACh;
dataJF039_NREMtoAwake_S1BF = AnalysisResults.JF039.Transitions.NREMtoAWAKE.GRAB_ACh;
dataJF040_NREMtoAwake_S1BF = AnalysisResults.JF040.Transitions.NREMtoAWAKE.GRAB_ACh;
dataJF048_NREMtoAwake_S1BF = AnalysisResults.JF048.Transitions.NREMtoAWAKE.GRAB_ACh;
dataJF050_NREMtoAwake_S1BF = AnalysisResults.JF050.Transitions.NREMtoAWAKE.GRAB_ACh;

dataJF037_NREMtoREM_S1BF = AnalysisResults.JF037.Transitions.NREMtoREM.GRAB_ACh;
dataJF038_NREMtoREM_S1BF = AnalysisResults.JF038.Transitions.NREMtoREM.GRAB_ACh;
dataJF039_NREMtoREM_S1BF = AnalysisResults.JF039.Transitions.NREMtoREM.GRAB_ACh;
dataJF040_NREMtoREM_S1BF = AnalysisResults.JF040.Transitions.NREMtoREM.GRAB_ACh;
dataJF048_NREMtoREM_S1BF = AnalysisResults.JF048.Transitions.NREMtoREM.GRAB_ACh;
dataJF050_NREMtoREM_S1BF = AnalysisResults.JF050.Transitions.NREMtoREM.GRAB_ACh;

dataJF037_REMtoAwake_S1BF = AnalysisResults.JF037.Transitions.REMtoAWAKE.GRAB_ACh;
dataJF038_REMtoAwake_S1BF = AnalysisResults.JF038.Transitions.REMtoAWAKE.GRAB_ACh;
dataJF039_REMtoAwake_S1BF = AnalysisResults.JF039.Transitions.REMtoAWAKE.GRAB_ACh;
dataJF040_REMtoAwake_S1BF = AnalysisResults.JF040.Transitions.REMtoAWAKE.GRAB_ACh;
dataJF048_REMtoAwake_S1BF = AnalysisResults.JF048.Transitions.REMtoAWAKE.GRAB_ACh;
dataJF050_REMtoAwake_S1BF = AnalysisResults.JF050.Transitions.REMtoAWAKE.GRAB_ACh;

% Calculating the mean for columns 601 to 900. This corresponds to the 10s
% the transition 

% Calculating the mean for columns 901 to 1800. This corresponds to the 10s
% after the transition 

%PFC

AwaketoNREM_PFC_mean_after_JF037 = mean(dataJF037_AwaketoNREM_PFC(901:1800)); 
AwaketoNREM_PFC_mean_after_JF038 = mean(dataJF038_AwaketoNREM_PFC(901:1800)); 
AwaketoNREM_PFC_mean_after_JF039 = mean(dataJF039_AwaketoNREM_PFC(901:1800)); 
AwaketoNREM_PFC_mean_after_JF040 = mean(dataJF040_AwaketoNREM_PFC(901:1800)); 
AwaketoNREM_PFC_mean_after_JF048 = mean(dataJF048_AwaketoNREM_PFC(901:1800)); 
AwaketoNREM_PFC_mean_after_JF050 = mean(dataJF050_AwaketoNREM_PFC(901:1800)); 

NREMtoAwake_PFC_mean_after_JF037 = mean(dataJF037_NREMtoAwake_PFC(901:1800)); 
NREMtoAwake_PFC_mean_after_JF038 = mean(dataJF038_NREMtoAwake_PFC(901:1800)); 
NREMtoAwake_PFC_mean_after_JF039 = mean(dataJF039_NREMtoAwake_PFC(901:1800)); 
NREMtoAwake_PFC_mean_after_JF040 = mean(dataJF040_NREMtoAwake_PFC(901:1800)); 
NREMtoAwake_PFC_mean_after_JF048 = mean(dataJF048_NREMtoAwake_PFC(901:1800)); 
NREMtoAwake_PFC_mean_after_JF050 = mean(dataJF050_NREMtoAwake_PFC(901:1800));

NREMtoREM_PFC_mean_after_JF037 = mean(dataJF037_NREMtoREM_PFC(901:1800));
NREMtoREM_PFC_mean_after_JF038 = mean(dataJF038_NREMtoREM_PFC(901:1800));
NREMtoREM_PFC_mean_after_JF039 = mean(dataJF039_NREMtoREM_PFC(901:1800));
NREMtoREM_PFC_mean_after_JF040 = mean(dataJF040_NREMtoREM_PFC(901:1800));
NREMtoREM_PFC_mean_after_JF048 = mean(dataJF048_NREMtoREM_PFC(901:1800));
NREMtoREM_PFC_mean_after_JF050 = mean(dataJF050_NREMtoREM_PFC(901:1800));

REMtoAwake_PFC_mean_after_JF037 = mean(dataJF037_REMtoAwake_PFC(901:1800));
REMtoAwake_PFC_mean_after_JF038 = mean(dataJF038_REMtoAwake_PFC(901:1800));
REMtoAwake_PFC_mean_after_JF039 = mean(dataJF039_REMtoAwake_PFC(901:1800));
REMtoAwake_PFC_mean_after_JF040 = mean(dataJF040_REMtoAwake_PFC(901:1800));
REMtoAwake_PFC_mean_after_JF048 = mean(dataJF048_REMtoAwake_PFC(901:1800));
REMtoAwake_PFC_mean_after_JF050 = mean(dataJF050_REMtoAwake_PFC(901:1800));

%S1BF

AwaketoNREM_S1BF_mean_after_JF037 = mean(dataJF037_AwaketoNREM_S1BF(901:1800)); 
AwaketoNREM_S1BF_mean_after_JF038 = mean(dataJF038_AwaketoNREM_S1BF(901:1800)); 
AwaketoNREM_S1BF_mean_after_JF039 = mean(dataJF039_AwaketoNREM_S1BF(901:1800)); 
AwaketoNREM_S1BF_mean_after_JF040 = mean(dataJF040_AwaketoNREM_S1BF(901:1800)); 
AwaketoNREM_S1BF_mean_after_JF048 = mean(dataJF048_AwaketoNREM_S1BF(901:1800)); 
AwaketoNREM_S1BF_mean_after_JF050 = mean(dataJF050_AwaketoNREM_S1BF(901:1800)); 

NREMtoAwake_S1BF_mean_after_JF037 = mean(dataJF037_NREMtoAwake_S1BF(901:1800)); 
NREMtoAwake_S1BF_mean_after_JF038 = mean(dataJF038_NREMtoAwake_S1BF(901:1800)); 
NREMtoAwake_S1BF_mean_after_JF039 = mean(dataJF039_NREMtoAwake_S1BF(901:1800)); 
NREMtoAwake_S1BF_mean_after_JF040 = mean(dataJF040_NREMtoAwake_S1BF(901:1800)); 
NREMtoAwake_S1BF_mean_after_JF048 = mean(dataJF048_NREMtoAwake_S1BF(901:1800)); 
NREMtoAwake_S1BF_mean_after_JF050 = mean(dataJF050_NREMtoAwake_S1BF(901:1800));

NREMtoREM_S1BF_mean_after_JF037 = mean(dataJF037_NREMtoREM_S1BF(901:1800));
NREMtoREM_S1BF_mean_after_JF038 = mean(dataJF038_NREMtoREM_S1BF(901:1800));
NREMtoREM_S1BF_mean_after_JF039 = mean(dataJF039_NREMtoREM_S1BF(901:1800));
NREMtoREM_S1BF_mean_after_JF040 = mean(dataJF040_NREMtoREM_S1BF(901:1800));
NREMtoREM_S1BF_mean_after_JF048 = mean(dataJF048_NREMtoREM_S1BF(901:1800));
NREMtoREM_S1BF_mean_after_JF050 = mean(dataJF050_NREMtoREM_S1BF(901:1800));

REMtoAwake_S1BF_mean_after_JF037 = mean(dataJF037_REMtoAwake_S1BF(901:1800));
REMtoAwake_S1BF_mean_after_JF038 = mean(dataJF038_REMtoAwake_S1BF(901:1800));
REMtoAwake_S1BF_mean_after_JF039 = mean(dataJF039_REMtoAwake_S1BF(901:1800));
REMtoAwake_S1BF_mean_after_JF040 = mean(dataJF040_REMtoAwake_S1BF(901:1800));
REMtoAwake_S1BF_mean_after_JF048 = mean(dataJF048_REMtoAwake_S1BF(901:1800));
REMtoAwake_S1BF_mean_after_JF050 = mean(dataJF050_REMtoAwake_S1BF(901:1800));


%% Creating a matrix with 'before' and 'after' means

%PFC

combined_means = [
AwaketoNREM_PFC_mean_before_JF037, AwaketoNREM_PFC_mean_after_JF037;
AwaketoNREM_PFC_mean_before_JF038, AwaketoNREM_PFC_mean_after_JF038;
AwaketoNREM_PFC_mean_before_JF039, AwaketoNREM_PFC_mean_after_JF039;
AwaketoNREM_PFC_mean_before_JF040, AwaketoNREM_PFC_mean_after_JF040;
AwaketoNREM_PFC_mean_before_JF048, AwaketoNREM_PFC_mean_after_JF048;
AwaketoNREM_PFC_mean_before_JF050, AwaketoNREM_PFC_mean_after_JF050;

NREMtoAwake_PFC_mean_before_JF037, NREMtoAwake_PFC_mean_after_JF037;
NREMtoAwake_PFC_mean_before_JF038, NREMtoAwake_PFC_mean_after_JF038;
NREMtoAwake_PFC_mean_before_JF039, NREMtoAwake_PFC_mean_after_JF039;
NREMtoAwake_PFC_mean_before_JF040, NREMtoAwake_PFC_mean_after_JF040;
NREMtoAwake_PFC_mean_before_JF048, NREMtoAwake_PFC_mean_after_JF048;
NREMtoAwake_PFC_mean_before_JF050, NREMtoAwake_PFC_mean_after_JF050;

NREMtoREM_PFC_mean_before_JF037, NREMtoREM_PFC_mean_after_JF037;
NREMtoREM_PFC_mean_before_JF038, NREMtoREM_PFC_mean_after_JF038;
NREMtoREM_PFC_mean_before_JF039, NREMtoREM_PFC_mean_after_JF039;
NREMtoREM_PFC_mean_before_JF040, NREMtoREM_PFC_mean_after_JF040;
NREMtoREM_PFC_mean_before_JF048, NREMtoREM_PFC_mean_after_JF048;
NREMtoREM_PFC_mean_before_JF050, NREMtoREM_PFC_mean_after_JF050;

REMtoAwake_PFC_mean_before_JF037, REMtoAwake_PFC_mean_after_JF037;
REMtoAwake_PFC_mean_before_JF038, REMtoAwake_PFC_mean_after_JF038;
REMtoAwake_PFC_mean_before_JF039, REMtoAwake_PFC_mean_after_JF039;
REMtoAwake_PFC_mean_before_JF040, REMtoAwake_PFC_mean_after_JF040;
REMtoAwake_PFC_mean_before_JF048, REMtoAwake_PFC_mean_after_JF048;
REMtoAwake_PFC_mean_before_JF050, REMtoAwake_PFC_mean_after_JF050;

AwaketoNREM_S1BF_mean_before_JF037, AwaketoNREM_S1BF_mean_after_JF037;
AwaketoNREM_S1BF_mean_before_JF038, AwaketoNREM_S1BF_mean_after_JF038;
AwaketoNREM_S1BF_mean_before_JF039, AwaketoNREM_S1BF_mean_after_JF039;
AwaketoNREM_S1BF_mean_before_JF040, AwaketoNREM_S1BF_mean_after_JF040;
AwaketoNREM_S1BF_mean_before_JF048, AwaketoNREM_S1BF_mean_after_JF048;
AwaketoNREM_S1BF_mean_before_JF050, AwaketoNREM_S1BF_mean_after_JF050;

NREMtoAwake_S1BF_mean_before_JF037, NREMtoAwake_S1BF_mean_after_JF037;
NREMtoAwake_S1BF_mean_before_JF038, NREMtoAwake_S1BF_mean_after_JF038;
NREMtoAwake_S1BF_mean_before_JF039, NREMtoAwake_S1BF_mean_after_JF039;
NREMtoAwake_S1BF_mean_before_JF040, NREMtoAwake_S1BF_mean_after_JF040;
NREMtoAwake_S1BF_mean_before_JF048, NREMtoAwake_S1BF_mean_after_JF048;
NREMtoAwake_S1BF_mean_before_JF050, NREMtoAwake_S1BF_mean_after_JF050;

NREMtoREM_S1BF_mean_before_JF037, NREMtoREM_S1BF_mean_after_JF037;
NREMtoREM_S1BF_mean_before_JF038, NREMtoREM_S1BF_mean_after_JF038;
NREMtoREM_S1BF_mean_before_JF039, NREMtoREM_S1BF_mean_after_JF039;
NREMtoREM_S1BF_mean_before_JF040, NREMtoREM_S1BF_mean_after_JF040;
NREMtoREM_S1BF_mean_before_JF048, NREMtoREM_S1BF_mean_after_JF048;
NREMtoREM_S1BF_mean_before_JF050, NREMtoREM_S1BF_mean_after_JF050;

REMtoAwake_S1BF_mean_before_JF037, REMtoAwake_S1BF_mean_after_JF037;
REMtoAwake_S1BF_mean_before_JF038, REMtoAwake_S1BF_mean_after_JF038;
REMtoAwake_S1BF_mean_before_JF039, REMtoAwake_S1BF_mean_after_JF039;
REMtoAwake_S1BF_mean_before_JF040, REMtoAwake_S1BF_mean_after_JF040;
REMtoAwake_S1BF_mean_before_JF048, REMtoAwake_S1BF_mean_after_JF048;
REMtoAwake_S1BF_mean_before_JF050, REMtoAwake_S1BF_mean_after_JF050;
];

% combined_means = [
%     mean_before_JF037, mean_after_JF037;
%     mean_before_JF038, mean_after_JF038;
%     mean_before_JF039, mean_after_JF039;
%     mean_before_JF040, mean_after_JF040;
%     mean_before_JF048, mean_after_JF048;
%     mean_before_JF050, mean_after_JF050;
% ];

combined_substracted = [
AwaketoNREM_PFC_mean_after_JF037 - AwaketoNREM_PFC_mean_before_JF037;
AwaketoNREM_PFC_mean_after_JF038 - AwaketoNREM_PFC_mean_before_JF038;
AwaketoNREM_PFC_mean_after_JF039 - AwaketoNREM_PFC_mean_before_JF039;
AwaketoNREM_PFC_mean_after_JF040 - AwaketoNREM_PFC_mean_before_JF040;
AwaketoNREM_PFC_mean_after_JF048 - AwaketoNREM_PFC_mean_before_JF048;
AwaketoNREM_PFC_mean_after_JF050 - AwaketoNREM_PFC_mean_before_JF050;

NREMtoAwake_PFC_mean_after_JF037 - NREMtoAwake_PFC_mean_before_JF037;
NREMtoAwake_PFC_mean_after_JF038 - NREMtoAwake_PFC_mean_before_JF038;
NREMtoAwake_PFC_mean_after_JF039 - NREMtoAwake_PFC_mean_before_JF039;
NREMtoAwake_PFC_mean_after_JF040 - NREMtoAwake_PFC_mean_before_JF040;
NREMtoAwake_PFC_mean_after_JF048 - NREMtoAwake_PFC_mean_before_JF048;
NREMtoAwake_PFC_mean_after_JF050 - NREMtoAwake_PFC_mean_before_JF050;

NREMtoREM_PFC_mean_after_JF037 - NREMtoREM_PFC_mean_before_JF037;
NREMtoREM_PFC_mean_after_JF038 - NREMtoREM_PFC_mean_before_JF038;
NREMtoREM_PFC_mean_after_JF039 - NREMtoREM_PFC_mean_before_JF039;
NREMtoREM_PFC_mean_after_JF040 - NREMtoREM_PFC_mean_before_JF040;
NREMtoREM_PFC_mean_after_JF048 - NREMtoREM_PFC_mean_before_JF048;
NREMtoREM_PFC_mean_after_JF050 - NREMtoREM_PFC_mean_before_JF050;

REMtoAwake_PFC_mean_after_JF037 - REMtoAwake_PFC_mean_before_JF037;
REMtoAwake_PFC_mean_after_JF038 - REMtoAwake_PFC_mean_before_JF038;
REMtoAwake_PFC_mean_after_JF039 - REMtoAwake_PFC_mean_before_JF039;
REMtoAwake_PFC_mean_after_JF040 - REMtoAwake_PFC_mean_before_JF040;
REMtoAwake_PFC_mean_after_JF048 - REMtoAwake_PFC_mean_before_JF048;
REMtoAwake_PFC_mean_after_JF050 - REMtoAwake_PFC_mean_before_JF050;

AwaketoNREM_S1BF_mean_after_JF037 - AwaketoNREM_S1BF_mean_before_JF037;
AwaketoNREM_S1BF_mean_after_JF038 - AwaketoNREM_S1BF_mean_before_JF038;
AwaketoNREM_S1BF_mean_after_JF039 - AwaketoNREM_S1BF_mean_before_JF039;
AwaketoNREM_S1BF_mean_after_JF040 - AwaketoNREM_S1BF_mean_before_JF040;
AwaketoNREM_S1BF_mean_after_JF048 - AwaketoNREM_S1BF_mean_before_JF048;
AwaketoNREM_S1BF_mean_after_JF050 - AwaketoNREM_S1BF_mean_before_JF050;

NREMtoAwake_S1BF_mean_after_JF037 - NREMtoAwake_S1BF_mean_before_JF037;
NREMtoAwake_S1BF_mean_after_JF038 - NREMtoAwake_S1BF_mean_before_JF038;
NREMtoAwake_S1BF_mean_after_JF039 - NREMtoAwake_S1BF_mean_before_JF039;
NREMtoAwake_S1BF_mean_after_JF040 - NREMtoAwake_S1BF_mean_before_JF040;
NREMtoAwake_S1BF_mean_after_JF048 - NREMtoAwake_S1BF_mean_before_JF048;
NREMtoAwake_S1BF_mean_after_JF050 - NREMtoAwake_S1BF_mean_before_JF050;

NREMtoREM_S1BF_mean_after_JF037 - NREMtoREM_S1BF_mean_before_JF037;
NREMtoREM_S1BF_mean_after_JF038 -NREMtoREM_S1BF_mean_before_JF038;
NREMtoREM_S1BF_mean_after_JF039 - NREMtoREM_S1BF_mean_before_JF039;
NREMtoREM_S1BF_mean_after_JF040 - NREMtoREM_S1BF_mean_before_JF040;
NREMtoREM_S1BF_mean_after_JF048 - NREMtoREM_S1BF_mean_before_JF048;
NREMtoREM_S1BF_mean_after_JF050 - NREMtoREM_S1BF_mean_before_JF050;

REMtoAwake_S1BF_mean_after_JF037 - REMtoAwake_S1BF_mean_before_JF037;
REMtoAwake_S1BF_mean_after_JF038 - REMtoAwake_S1BF_mean_before_JF038;
REMtoAwake_S1BF_mean_after_JF039 - REMtoAwake_S1BF_mean_before_JF039;
REMtoAwake_S1BF_mean_after_JF040 - REMtoAwake_S1BF_mean_before_JF040;
REMtoAwake_S1BF_mean_after_JF048 - REMtoAwake_S1BF_mean_before_JF048;
REMtoAwake_S1BF_mean_after_JF050 - REMtoAwake_S1BF_mean_before_JF050;
];

% combined_substracted = [
%     mean_after_JF037 - mean_before_JF037;
%     mean_after_JF038 - mean_before_JF038;
%     mean_after_JF039 - mean_before_JF039;
%     mean_after_JF040 - mean_before_JF040;
%     mean_after_JF048 - mean_before_JF048;
%     mean_after_JF050 - mean_before_JF050;
% ];

% % Converting the matrix into a table for better readability
% combined_table = array2table(combined_means, ...
%     'VariableNames', {'Mean_Before', 'Mean_After'}, ...
%     'RowNames', {'JF037', 'JF038', 'JF039', 'JF040', 'JF048', 'JF050'});

% % Display the table
% disp(combined_table);
