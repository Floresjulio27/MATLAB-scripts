function [] = Fig_Stim_SST_ResponseIndividual_Alcohol_Stats(rootFolder,saveState,delim)
%________________________________________________________________________________________________________________________
% Written by Kevin L. Turner
% The Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner
%
% This code was originally written by Dr. Kevin L. Turner, modified by MD
% Shakhawat Hossain and adapted by Julio Flores-Cuadra
%________________________________________________________________________________________________________________________
path = [rootFolder delim 'Results_SST_project'];
cd(path)
resultsStruct = 'Results_Evoked';
load(resultsStruct)

% Set up group, day, and solenoid info
groups = {'Water', 'Alcohol'};
days = {'Day1', 'Day2', 'Day3', 'Day4', 'Day5', 'Day6'};
solenoidNames = {'LPadSol','RPadSol','AudSol'};
compDataTypes = {'Ipsi','Contra','Auditory'};

% Loop through groups
for ee = 1:length(groups)
    currentGroup = groups{1,ee};

    % Loop through animals
    for  xx = 1:length(days)
        currentDay = days{1,xx};
        % Skip if this day doesn't exist
        if ~isfield(Results_Evoked.(currentGroup), currentDay)
            continue;
        end
        FPanimalIDs = fieldnames(Results_Evoked.(currentGroup).(currentDay));
      
         if isempty(FPanimalIDs)
            continue
        end

        % Loop through days
        for aa = 1:length(FPanimalIDs)
            animalID = FPanimalIDs{aa,1}; %it is organizing it as a column vector

            if ~isfield(Results_Evoked.(currentGroup).(currentDay).(animalID), 'Stim')
                continue
            end
           
            % Loop through solenoids
            for dd = 1:length(solenoidNames)
                solenoidName = solenoidNames{1,dd};

                data.(currentGroup).(currentDay).(animalID).P_NE.(solenoidName).count = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_NE.(solenoidName).count;
                data.(currentGroup).(currentDay).(animalID).P_NE.(solenoidName).CBV = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_NE.(solenoidName).CBV.CBV;
                data.(currentGroup).(currentDay).(animalID).P_ACh.(solenoidName).CBV = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_ACh.(solenoidName).CBV.CBV;
                data.(currentGroup).(currentDay).(animalID).P_NE.(solenoidName).GFP = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_NE.(solenoidName).GFP.GFP;
                data.(currentGroup).(currentDay).(animalID).P_ACh.(solenoidName).GFP = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_ACh.(solenoidName).GFP.GFP;

                % Standard deviations
                data.(currentGroup).(currentDay).(animalID).P_NE.(solenoidName).CBV_std = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_NE.(solenoidName).CBV.CBVStD;
                data.(currentGroup).(currentDay).(animalID).P_ACh.(solenoidName).CBV_std = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_ACh.(solenoidName).CBV.CBVStD;
                data.(currentGroup).(currentDay).(animalID).P_NE.(solenoidName).GFP_std = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_NE.(solenoidName).GFP.GFPStD;
                data.(currentGroup).(currentDay).(animalID).P_ACh.(solenoidName).GFP_std = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_ACh.(solenoidName).GFP.GFPStD;

                % Cortical activity
                data.(currentGroup).(currentDay).(animalID).cortical.(solenoidName).cortMUA = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.cortical.(solenoidName).MUA.corticalData;
                data.(currentGroup).(currentDay).(animalID).cortical.(solenoidName).cortGam = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.cortical.(solenoidName).Gam.corticalData;
                data.(currentGroup).(currentDay).(animalID).cortical.(solenoidName).timeVector = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.cortical.(solenoidName).timeVector;
                data.(currentGroup).(currentDay).(animalID).cortical.(solenoidName).cortS = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.cortical.(solenoidName).LFP.corticalS;
                data.(currentGroup).(currentDay).(animalID).cortical.(solenoidName).cortS_Gam = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.cortical.(solenoidName).LFP.corticalS(49:end,20:23);
                data.(currentGroup).(currentDay).(animalID).cortical.(solenoidName).T = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.cortical.(solenoidName).LFP.T;
                data.(currentGroup).(currentDay).(animalID).cortical.(solenoidName).F = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.cortical.(solenoidName).LFP.F;
            end 
               
                % concatenate the data from the contra and ipsi data
                % contra
                data.(currentGroup).(currentDay).(animalID).Contra.count = data.(currentGroup).(currentDay).(animalID).P_NE.LPadSol.count;
                data.(currentGroup).(currentDay).(animalID).Contra.P_AChCBV = data.(currentGroup).(currentDay).(animalID).P_ACh.RPadSol.CBV;
                data.(currentGroup).(currentDay).(animalID).Contra.P_AChGFP = data.(currentGroup).(currentDay).(animalID).P_ACh.RPadSol.GFP;
                data.(currentGroup).(currentDay).(animalID).Contra.P_NECBV = data.(currentGroup).(currentDay).(animalID).P_NE.LPadSol.CBV;
                data.(currentGroup).(currentDay).(animalID).Contra.P_NEGFP = data.(currentGroup).(currentDay).(animalID).P_NE.LPadSol.GFP;

                data.(currentGroup).(currentDay).(animalID).Contra.P_AChCBV_std = data.(currentGroup).(currentDay).(animalID).P_ACh.RPadSol.CBV_std;
                data.(currentGroup).(currentDay).(animalID).Contra.P_AChGFP_std = data.(currentGroup).(currentDay).(animalID).P_ACh.RPadSol.GFP_std;
                data.(currentGroup).(currentDay).(animalID).Contra.P_NECBV_std = data.(currentGroup).(currentDay).(animalID).P_NE.LPadSol.CBV_std;
                data.(currentGroup).(currentDay).(animalID).Contra.P_NEGFP_std = data.(currentGroup).(currentDay).(animalID).P_NE.LPadSol.GFP_std;

                data.(currentGroup).(currentDay).(animalID).Contra.cortMUA = data.(currentGroup).(currentDay).(animalID).cortical.LPadSol.cortMUA;
                data.(currentGroup).(currentDay).(animalID).Contra.cortGam = data.(currentGroup).(currentDay).(animalID).cortical.LPadSol.cortGam;
                data.(currentGroup).(currentDay).(animalID).Contra.timeVector = data.(currentGroup).(currentDay).(animalID).cortical.LPadSol.timeVector;
                data.(currentGroup).(currentDay).(animalID).Contra.cortS = data.(currentGroup).(currentDay).(animalID).cortical.LPadSol.cortS;
                data.(currentGroup).(currentDay).(animalID).Contra.cortS_Gam = data.(currentGroup).(currentDay).(animalID).cortical.LPadSol.cortS_Gam;
                data.(currentGroup).(currentDay).(animalID).Contra.T = data.(currentGroup).(currentDay).(animalID).cortical.LPadSol.T;
                data.(currentGroup).(currentDay).(animalID).Contra.F = data.(currentGroup).(currentDay).(animalID).cortical.LPadSol.F;

                %ipsi
                data.(currentGroup).(currentDay).(animalID).Ipsi.count = data.(currentGroup).(currentDay).(animalID).P_NE.RPadSol.count;
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_AChCBV_std = data.(currentGroup).(currentDay).(animalID).P_ACh.LPadSol.CBV_std;
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_AChGFP_std = data.(currentGroup).(currentDay).(animalID).P_ACh.LPadSol.GFP_std;
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_NECBV_std = data.(currentGroup).(currentDay).(animalID).P_NE.RPadSol.CBV_std;
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_NEGFP_std = data.(currentGroup).(currentDay).(animalID).P_NE.RPadSol.GFP_std;

                data.(currentGroup).(currentDay).(animalID).Ipsi.P_AChCBV = data.(currentGroup).(currentDay).(animalID).P_ACh.LPadSol.CBV;
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_AChGFP = data.(currentGroup).(currentDay).(animalID).P_ACh.LPadSol.GFP;
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_NECBV = data.(currentGroup).(currentDay).(animalID).P_NE.RPadSol.CBV;
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_NEGFP = data.(currentGroup).(currentDay).(animalID).P_NE.RPadSol.GFP;

                data.(currentGroup).(currentDay).(animalID).Ipsi.cortMUA = data.(currentGroup).(currentDay).(animalID).cortical.RPadSol.cortMUA;
                data.(currentGroup).(currentDay).(animalID).Ipsi.cortGam = data.(currentGroup).(currentDay).(animalID).cortical.RPadSol.cortGam;
                data.(currentGroup).(currentDay).(animalID).Ipsi.timeVector = data.(currentGroup).(currentDay).(animalID).cortical.RPadSol.timeVector;
                data.(currentGroup).(currentDay).(animalID).Ipsi.cortS = data.(currentGroup).(currentDay).(animalID).cortical.RPadSol.cortS;
                data.(currentGroup).(currentDay).(animalID).Ipsi.cortS_Gam = data.(currentGroup).(currentDay).(animalID).cortical.RPadSol.cortS_Gam;
                data.(currentGroup).(currentDay).(animalID).Ipsi.T = data.(currentGroup).(currentDay).(animalID).cortical.RPadSol.T;
                data.(currentGroup).(currentDay).(animalID).Ipsi.F = data.(currentGroup).(currentDay).(animalID).cortical.RPadSol.F;

                %auditory
                data.(currentGroup).(currentDay).(animalID).Auditory.count = data.(currentGroup).(currentDay).(animalID).P_NE.AudSol.count;

                data.(currentGroup).(currentDay).(animalID).Auditory.P_NECBV = data.(currentGroup).(currentDay).(animalID).P_NE.AudSol.CBV;
                data.(currentGroup).(currentDay).(animalID).Auditory.P_AChCBV = data.(currentGroup).(currentDay).(animalID).P_ACh.AudSol.CBV;
                data.(currentGroup).(currentDay).(animalID).Auditory.P_NEGFP = data.(currentGroup).(currentDay).(animalID).P_NE.AudSol.GFP;
                data.(currentGroup).(currentDay).(animalID).Auditory.P_AChGFP = data.(currentGroup).(currentDay).(animalID).P_ACh.AudSol.GFP;

                data.(currentGroup).(currentDay).(animalID).Auditory.P_NECBV_std = data.(currentGroup).(currentDay).(animalID).P_NE.AudSol.CBV_std;
                data.(currentGroup).(currentDay).(animalID).Auditory.P_AChCBV_std = data.(currentGroup).(currentDay).(animalID).P_ACh.AudSol.CBV_std;
                data.(currentGroup).(currentDay).(animalID).Auditory.P_NEGFP_std = data.(currentGroup).(currentDay).(animalID).P_NE.AudSol.GFP_std;
                data.(currentGroup).(currentDay).(animalID).Auditory.P_AChGFP_std = data.(currentGroup).(currentDay).(animalID).P_ACh.AudSol.GFP_std;

                data.(currentGroup).(currentDay).(animalID).Auditory.cortMUA = data.(currentGroup).(currentDay).(animalID).cortical.AudSol.cortMUA;
                data.(currentGroup).(currentDay).(animalID).Auditory.cortGam = data.(currentGroup).(currentDay).(animalID).cortical.AudSol.cortGam;
                data.(currentGroup).(currentDay).(animalID).Auditory.timeVector = data.(currentGroup).(currentDay).(animalID).cortical.AudSol.timeVector;
                data.(currentGroup).(currentDay).(animalID).Auditory.cortS = data.(currentGroup).(currentDay).(animalID).cortical.AudSol.cortS;
                data.(currentGroup).(currentDay).(animalID).Auditory.cortS_Gam = data.(currentGroup).(currentDay).(animalID).cortical.AudSol.cortS_Gam;
                data.(currentGroup).(currentDay).(animalID).Auditory.T = data.(currentGroup).(currentDay).(animalID).cortical.AudSol.T;
                data.(currentGroup).(currentDay).(animalID).Auditory.F = data.(currentGroup).(currentDay).(animalID).cortical.AudSol.F;
        end
    end
end

% Loop through groups
for ee = 1:length(groups)
    currentGroup = groups{1,ee};

    % Loop through animals
    for  xx = 1:length(days)
        currentDay = days{1,xx};
        % Skip if this day doesn't exist
        if ~isfield(Results_Evoked.(currentGroup), currentDay)
            continue;
        end
        FPanimalIDs = fieldnames(Results_Evoked.(currentGroup).(currentDay));
      
         if isempty(FPanimalIDs)
            continue
        end

        % Loop through days
        for aa = 1:length(FPanimalIDs)
            animalID = FPanimalIDs{aa,1}; %it is organizing it as a column vector

            if ~isfield(Results_Evoked.(currentGroup).(currentDay).(animalID), 'Stim')
                continue
            end

            for dd = 1:length(solenoidNames)
                solenoidName = solenoidNames{1,dd};
                
                % Determine matrix size
                MatLength = size(Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_NE.(solenoidName).CBV.CBVRaw, 1);
                
                if ~isfield(data.(currentGroup).(currentDay).(animalID).P_NE.(solenoidName), 'CBVRaw') || isempty(data.(currentGroup).(currentDay).(animalID).P_NE.(solenoidName).CBVRaw)
                    DataLength = 0;
                else
                    DataLength = size(data.(currentGroup).(currentDay).(animalID).P_NE.(solenoidName).CBVRaw, 1);
                end
                
                % Concatenate raw data
                data.(currentGroup).(currentDay).(animalID).P_NE.(solenoidName).CBVRaw(DataLength+1:DataLength+MatLength,:)  = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_NE.(solenoidName).CBV.CBVRaw;
                data.(currentGroup).(currentDay).(animalID).P_ACh.(solenoidName).CBVRaw(DataLength+1:DataLength+MatLength,:) = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_ACh.(solenoidName).CBV.CBVRaw;
                data.(currentGroup).(currentDay).(animalID).P_NE.(solenoidName).GFPRaw(DataLength+1:DataLength+MatLength,:)  = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_NE.(solenoidName).GFP.GFPRaw;
                data.(currentGroup).(currentDay).(animalID).P_ACh.(solenoidName).GFPRaw(DataLength+1:DataLength+MatLength,:) = Results_Evoked.(currentGroup).(currentDay).(animalID).Stim.P_ACh.(solenoidName).GFP.GFPRaw;
            end 

                % concatenate the data from the contra and ipsi data
                % contra
                data.(currentGroup).(currentDay).(animalID).Contra.P_AChCBVRaw = data.(currentGroup).(currentDay).(animalID).P_ACh.RPadSol.CBVRaw;
                data.(currentGroup).(currentDay).(animalID).Contra.P_AChGFPRaw = data.(currentGroup).(currentDay).(animalID).P_ACh.RPadSol.GFPRaw;
                data.(currentGroup).(currentDay).(animalID).Contra.P_NECBVRaw = data.(currentGroup).(currentDay).(animalID).P_NE.LPadSol.CBVRaw;
                data.(currentGroup).(currentDay).(animalID).Contra.P_NEGFPRaw = data.(currentGroup).(currentDay).(animalID).P_NE.LPadSol.GFPRaw;

                %ipsi
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_AChCBVRaw = data.(currentGroup).(currentDay).(animalID).P_ACh.LPadSol.CBVRaw;
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_AChGFPRaw = data.(currentGroup).(currentDay).(animalID).P_ACh.LPadSol.GFPRaw;
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_NECBVRaw = data.(currentGroup).(currentDay).(animalID).P_NE.RPadSol.CBVRaw;
                data.(currentGroup).(currentDay).(animalID).Ipsi.P_NEGFPRaw = data.(currentGroup).(currentDay).(animalID).P_NE.RPadSol.GFPRaw;

                %auditory
                data.(currentGroup).(currentDay).(animalID).Auditory.P_NECBVRaw = data.(currentGroup).(currentDay).(animalID).P_NE.AudSol.CBVRaw;
                data.(currentGroup).(currentDay).(animalID).Auditory.P_AChCBVRaw = data.(currentGroup).(currentDay).(animalID).P_ACh.AudSol.CBVRaw;
                data.(currentGroup).(currentDay).(animalID).Auditory.P_NEGFPRaw = data.(currentGroup).(currentDay).(animalID).P_NE.AudSol.GFPRaw;
                data.(currentGroup).(currentDay).(animalID).Auditory.P_AChGFPRaw = data.(currentGroup).(currentDay).(animalID).P_ACh.AudSol.GFPRaw;
        end
    end
end

%Calculate mean time vectors just for plotting 


for ee = 1:length(groups)
    currentGroup = groups{ee};
    for dd = 1:length(days)
        currentDay = days{dd};
        % skip missing days
        if ~isfield(data.(currentGroup), currentDay)
            continue;
        end
        animalIDs = fieldnames( data.(currentGroup).(currentDay) );
        if isempty(animalIDs)
            continue;
        end
        
        for ff = 1:length(compDataTypes)
            comp = compDataTypes{ff};
            
            % assume the field you want to average is .count,
            % and that each animal has a vector of same length N
            firstID = animalIDs{1};
            N = size( data.(currentGroup).(currentDay).(firstID).(comp).count, 1 );
            nAnimals = length(animalIDs);
            
            % build matrix [time × animals]
            allCounts = zeros(N, nAnimals);
            for ii = 1:nAnimals
                id = animalIDs{ii};
                allCounts(:,ii) = data.(currentGroup).(currentDay).(id).(comp).count;
            end
            
            % now average (and std) across the 2nd dim → one trace per currentDay
            data.(currentGroup).(currentDay).GetTimeVector.(comp).count = mean( allCounts, 2 );
            data.(currentGroup).(currentDay).GetTimeVector.(comp).count_std = std( allCounts, 0, 2 );
            
            % if you also want to propagate the timeVector:
            data.(currentGroup).(currentDay).GetTimeVector.(comp).mean_timeVector = data.(currentGroup).(currentDay).(firstID).(comp).timeVector;
        end
    end
end


% Compute Day-by-Day averages (and SEM) across animals

% ————————————————————————————————————————————————————————————————
% Compute Day-by-Day averages (and SEM) across animals. This is other way
% to organize it. 
% ————————————————————————————————————————————————————————————————
% Loop through groups  
for ee = 1:length(groups)
    currentGroup = groups{1,ee};

    % Loop through animals
    for  xx = 1:length(days)
        currentDay = days{1,xx};
        % Skip if this day doesn't exist
        if ~isfield(Results_Evoked.(currentGroup), currentDay)
            continue;
        end
        FPanimalIDs = fieldnames(Results_Evoked.(currentGroup).(currentDay));

        if isempty(FPanimalIDs)
            continue
        end
        
        for ff = 1:length(compDataTypes)
            comp = compDataTypes{ff};

            % gather all animals' raw traces for this day
            allAChCBV = [];
            allNECBV  = [];
            allAChGFP = [];
            allNEGFP  = [];
            for aa = 1:length(FPanimalIDs)
                animalID = FPanimalIDs{aa};
                if ~isfield(data.(currentGroup).(currentDay), animalID)
                    continue
                end
                S = data.(currentGroup).(currentDay).(animalID).(comp);
                allAChCBV = [allAChCBV; S.P_AChCBVRaw];   %#ok<AGROW>
                allNECBV  = [allNECBV;  S.P_NECBVRaw];    %#ok<AGROW>
                allAChGFP = [allAChGFP; S.P_AChGFPRaw];   %#ok<AGROW>
                allNEGFP  = [allNEGFP;  S.P_NEGFPRaw];    %#ok<AGROW>
            end

            % compute and store stats
            N = size(allAChCBV,1);
            data.(currentGroup).(currentDay).(comp).P_AChCBV_N_Exp  = N;
            data.(currentGroup).(currentDay).(comp).P_AChCBV_Mean   = mean(allAChCBV,   1);
            data.(currentGroup).(currentDay).(comp).P_AChCBV_SEM    = std(allAChCBV,0,1) / sqrt(N);

            N = size(allNECBV,1);
            data.(currentGroup).(currentDay).(comp).P_NECBV_N_Exp   = N;
            data.(currentGroup).(currentDay).(comp).P_NECBV_Mean    = mean(allNECBV,    1);
            data.(currentGroup).(currentDay).(comp).P_NECBV_SEM     = std(allNECBV,0,1)  / sqrt(N);

            N = size(allAChGFP,1);
            data.(currentGroup).(currentDay).(comp).P_AChGFP_N_Exp  = N;
            data.(currentGroup).(currentDay).(comp).P_AChGFP_Mean   = mean(allAChGFP,   1);
            data.(currentGroup).(currentDay).(comp).P_AChGFP_SEM    = std(allAChGFP,0,1) / sqrt(N);

            N = size(allNEGFP,1);
            data.(currentGroup).(currentDay).(comp).P_NEGFP_N_Exp   = N;
            data.(currentGroup).(currentDay).(comp).P_NEGFP_Mean    = mean(allNEGFP,    1);
            data.(currentGroup).(currentDay).(comp).P_NEGFP_SEM     = std(allNEGFP,0,1)  / sqrt(N);
        end
    end
end

%% Get the data snip 1 sec - 5sec after puff

for ee = 1:length(groups)
    currentGroup = groups{ee};

    for xx = 1:length(days)
        currentDay = days{xx};
        if ~isfield(Results_Evoked.(currentGroup), currentDay), continue; end
        FPanimalIDs = fieldnames(Results_Evoked.(currentGroup).(currentDay));
        if isempty(FPanimalIDs), continue; end

        for ff = 1:length(compDataTypes)
            comp = compDataTypes{ff};

            % ——— initialize your concat-arrays for this day/comp ———
            allSnip_CBV_S1BF = [];
            allSnip_GFP_S1BF = [];
            allSnip_CBV_PFC  = [];
            allSnip_GFP_PFC  = [];

            % get your day-level time vector & window just once
            Tday = data.(currentGroup).(currentDay).GetTimeVector.(comp).mean_timeVector;
            i0   = find(Tday >= 1, 1, 'first');
            i1   = find(Tday <= 5, 1, 'last');
            i2_5   = find(Tday >= 2.5, 1, 'first'); %for CBV S1BF
            i5   = find(Tday <= 5, 1, 'last'); %For CBV S1BF

            for aa = 1:length(FPanimalIDs)
                animalID = FPanimalIDs{aa};
                if ~isfield(data.(currentGroup).(currentDay), animalID)
                    continue;
                end

                % pull out each animal’s trial-by-trial raw AUCs
                CBV_S1BF = data.(currentGroup).(currentDay).(animalID).(comp).P_AChCBVRaw;
                GFP_S1BF = data.(currentGroup).(currentDay).(animalID).(comp).P_AChGFPRaw;
                CBV_PFC  = data.(currentGroup).(currentDay).(animalID).(comp).P_NECBVRaw;
                GFP_PFC  = data.(currentGroup).(currentDay).(animalID).(comp).P_NEGFPRaw;

                % compute the per-trial snip means for getting the
                % individual animal value 
                nTrials = size(CBV_S1BF,1);
                DataSnip_CBV_S1BF_mouse = nan(nTrials,1);
                DataSnip_GFP_S1BF_mouse = nan(nTrials,1);
                DataSnip_CBV_PFC_mouse  = nan(nTrials,1);
                DataSnip_GFP_PFC_mouse  = nan(nTrials,1);
                for tt = 1:nTrials
                    DataSnip_CBV_S1BF_mouse(tt) = mean( CBV_S1BF(tt, i2_5:i5) );
                    DataSnip_GFP_S1BF_mouse(tt) = mean( GFP_S1BF(tt, i0:i1) );
                    DataSnip_CBV_PFC_mouse(tt)  = mean( CBV_PFC(tt,  i0:i1) );
                    DataSnip_GFP_PFC_mouse(tt)  = mean( GFP_PFC(tt,  i0:i1) );
                end

                % save back if you still need per-animal stats
                data.(currentGroup).(currentDay).(animalID).(comp).DataSnip_CBV_S1BF_mouse = mean(DataSnip_CBV_S1BF_mouse);
                data.(currentGroup).(currentDay).(animalID).(comp).DataSnip_GFP_S1BF_mouse = mean(DataSnip_GFP_S1BF_mouse);
                data.(currentGroup).(currentDay).(animalID).(comp).DataSnip_CBV_PFC_mouse  = mean(DataSnip_CBV_PFC_mouse);
                data.(currentGroup).(currentDay).(animalID).(comp).DataSnip_GFP_PFC_mouse  = mean(DataSnip_GFP_PFC_mouse);

                %Just get the time window that you want
                DataSnip_CBV_S1BF = CBV_S1BF(:, i2_5:i5);
                DataSnip_GFP_S1BF = GFP_S1BF(:, i0:i1);
                DataSnip_CBV_PFC  = CBV_PFC(:,  i0:i1);
                DataSnip_GFP_PFC  = GFP_PFC(:,  i0:i1);

                % Concatonate animals according to the window of interest
                allSnip_CBV_S1BF = [allSnip_CBV_S1BF; DataSnip_CBV_S1BF];
                allSnip_GFP_S1BF = [allSnip_GFP_S1BF;  DataSnip_GFP_S1BF];
                allSnip_CBV_PFC  = [allSnip_CBV_PFC;  DataSnip_CBV_PFC];
                allSnip_GFP_PFC  = [allSnip_GFP_PFC; DataSnip_GFP_PFC];
            end

            % get the mean and SEM
            data.(currentGroup).(currentDay).(comp).allSnip_CBV_S1BF = allSnip_CBV_S1BF;
            data.(currentGroup).(currentDay).(comp).meanSnip_CBV_S1BF = mean(allSnip_CBV_S1BF);
            data.(currentGroup).(currentDay).(comp).semSnip_CBV_S1BF  = std(allSnip_CBV_S1BF)/sqrt(numel(allSnip_CBV_S1BF));
            data.(currentGroup).(currentDay).(comp).AUCSnip_CBV_S1BF = mean(data.(currentGroup).(currentDay).(comp).meanSnip_CBV_S1BF);

            data.(currentGroup).(currentDay).(comp).allSnip_GFP_S1BF = allSnip_GFP_S1BF;
            data.(currentGroup).(currentDay).(comp).meanSnip_GFP_S1BF = mean(allSnip_GFP_S1BF);
            data.(currentGroup).(currentDay).(comp).semSnip_GFP_S1BF  = std(allSnip_GFP_S1BF)/sqrt(numel(allSnip_GFP_S1BF));
             data.(currentGroup).(currentDay).(comp).AUCSnip_GFP_S1BF = mean(data.(currentGroup).(currentDay).(comp).meanSnip_GFP_S1BF);

            data.(currentGroup).(currentDay).(comp).allSnip_CBV_PFC  = allSnip_CBV_PFC;
            data.(currentGroup).(currentDay).(comp).meanSnip_CBV_PFC  = mean(allSnip_CBV_PFC);
            data.(currentGroup).(currentDay).(comp).semSnip_CBV_PFC   = std(allSnip_CBV_PFC)/sqrt(numel(allSnip_CBV_PFC));
            data.(currentGroup).(currentDay).(comp).AUCSnip_CBV_PFC = mean(data.(currentGroup).(currentDay).(comp).meanSnip_CBV_PFC);

            data.(currentGroup).(currentDay).(comp).allSnip_GFP_PFC  = allSnip_GFP_PFC;
            data.(currentGroup).(currentDay).(comp).meanSnip_GFP_PFC  = mean(allSnip_GFP_PFC);
            data.(currentGroup).(currentDay).(comp).semSnip_GFP_PFC   = std(allSnip_GFP_PFC)/sqrt(numel(allSnip_GFP_PFC));
            data.(currentGroup).(currentDay).(comp).AUCSnip_GFP_PFC = mean(data.(currentGroup).(currentDay).(comp).meanSnip_GFP_PFC);
        end
    end
end

% statistics - generalized linear mixed effects model. This is from Kevin.

% Data of interest
snipFields = {'DataSnip_CBV_S1BF_mouse', 'DataSnip_GFP_S1BF_mouse', 'DataSnip_CBV_PFC_mouse', 'DataSnip_GFP_PFC_mouse'};

% Initialize vectors
Mouse   = {};
Group   = {};
Day     = {};
Region  = {};
Signal  = {};
AUC     = [];

% Loop groups × days × animals × signals
for gg = 1:length(groups)
    currentGroup = groups{gg};

    for dd = 1:length(days)
        currentDay = days{dd};
        if ~isfield(data.(currentGroup), currentDay)
            continue;
        end

        FPanimalIDs = fieldnames(Results_Evoked.(currentGroup).(currentDay));
        if isempty(FPanimalIDs)
            continue;
        end

        for ff = 1:length(snipFields)
            fld = snipFields{ff};

            for aa = 1:length(FPanimalIDs)
                animalID = FPanimalIDs{aa};
                if ~isfield(data.(currentGroup).(currentDay), animalID)
                    continue;
                end

                % grab the single scalar AUC
                val = data.(currentGroup).(currentDay).(animalID).Contra.(fld);

                % parse out Region & Signal/ split the name so I can
                % allocate
                parts = split(fld, '_');
                % parts = {'DataSnip','CBV','S1BF','mouse'} or
                % {'DataSnip','CBV','S1BF','mouse'} This will be the result
                signal    = parts{2};
                region = parts{3};

                % Start prepping the table
                Mouse{end+1,1}  = animalID;
                Group{end+1,1}  = currentGroup;
                Day{end+1,1}    = currentDay;
                Region{end+1,1} = region;
                Signal{end+1,1} = signal;
                AUC(end+1,1)    = val;
            end
        end
    end
end

% Build the table for fitglme
T = table( categorical(Mouse), categorical(Group), categorical(Day), categorical(Region),categorical(Signal), AUC, ...
    'VariableNames',{'Mouse','Group','Day','Region','Signal','AUC'} );

% % Fit a mixed-effects model (from Kevin and all interactions)
% %    Here we include main effects of Group, Day, Region, Signal
% %    plus their interactions, with a random intercept per Mouse:
% formula = 'AUC ~ 1 + Group*Day*Region*Signal + (1|Mouse)';
% glme     = fitglme(T, formula);
% 
% % 6) View results
%disp(glme);

%% Within‐Group Day1 vs Day2–Day5 (GFP and CBV separted)
% if you do not have all the days it will crash and I do not know who to
% skip them

%% PFC WATER
%Water, CBV
T_W_CBV_PFC = T( T.Group == categorical("Water") & T.Region == categorical("PFC") & T.Signal == categorical("CBV"), : );
T_W_CBV_PFC.Day = categorical(T_W_CBV_PFC.Day, {'Day1','Day2','Day3','Day4','Day5'}, 'Ordinal', false);
glme_W_CBV_PFC = fitglme( T_W_CBV_PFC, 'AUC ~ 1 + Day + (1|Mouse)' );
disp('Water, CBV: DayK vs Day1')
disp(glme_W_CBV_PFC)

% Water, GFP
T_W_GFP_PFC = T( T.Group == categorical("Water") & T.Region == categorical("PFC") & T.Signal == categorical("GFP"), : );
T_W_GFP_PFC.Day = categorical(T_W_GFP_PFC.Day, {'Day1','Day2','Day3','Day4','Day5'}, 'Ordinal', false);
glme_W_GFP_PFC = fitglme( T_W_GFP_PFC, 'AUC ~ 1 + Day + (1|Mouse)' );
disp('Water, GFP: DayK vs Day1')
disp(glme_W_GFP_PFC)

%% S1BF WATER
%Water, CBV
T_W_CBV_S1BF = T( T.Group == categorical("Water") & T.Region == categorical("S1BF") & T.Signal == categorical("CBV"), : );
T_W_CBV_S1BF.Day = categorical(T_W_CBV_S1BF.Day, {'Day1','Day2','Day3','Day4','Day5'}, 'Ordinal', false);
glme_W_CBV_S1BF = fitglme( T_W_CBV_S1BF, 'AUC ~ 1 + Day + (1|Mouse)' );
disp('Water, CBV: DayK vs Day1')
disp(glme_W_CBV_S1BF)

% Water, GFP
T_W_GFP_S1BF = T( T.Group == categorical("Water") & T.Region == categorical("S1BF") & T.Signal == categorical("GFP"), : );
T_W_GFP_S1BF.Day = categorical(T_W_GFP_S1BF.Day, {'Day1','Day2','Day3','Day4','Day5'}, 'Ordinal', false);
glme_W_GFP_S1BF = fitglme( T_W_GFP_S1BF, 'AUC ~ 1 + Day + (1|Mouse)' );
disp('Water, GFP: DayK vs Day1')
disp(glme_W_GFP_S1BF)

%% PFC ALCOHOL
%Alcohol, CBV
T_A_CBV_PFC = T( T.Group == categorical("Alcohol") & T.Region == categorical("PFC") & T.Signal == categorical("CBV"), : );
T_A_CBV_PFC.Day = categorical(T_A_CBV_PFC.Day, {'Day1','Day2','Day3','Day4','Day5'}, 'Ordinal', false);
glme_A_CBV_PFC = fitglme( T_A_CBV_PFC, 'AUC ~ 1 + Day + (1|Mouse)' );
disp('Alcohol, CBV: DayK vs Day1')
disp(glme_A_CBV_PFC)

%Alcohol, GFP
T_A_GFP_PFC = T( T.Group == categorical("Alcohol") & T.Region == categorical("PFC") & T.Signal == categorical("GFP"), : );
T_A_GFP_PFC.Day = categorical(T_A_GFP_PFC.Day, {'Day1','Day2','Day3','Day4','Day5'}, 'Ordinal', false);
glme_A_GFP_PFC = fitglme( T_A_GFP_PFC, 'AUC ~ 1 + Day + (1|Mouse)' );
disp('Alcohol, GFP: DayK vs Day1')
disp(glme_A_GFP_PFC)

%% S1BF ALCOHOL
%Alcohol, CBV
T_A_CBV_S1BF = T( T.Group == categorical("Alcohol") & T.Region == categorical("S1BF") & T.Signal == categorical("CBV"), : );
T_A_CBV_S1BF.Day = categorical(T_A_CBV_S1BF.Day, {'Day1','Day2','Day3','Day4','Day5'}, 'Ordinal', false);
glme_A_CBV_S1BF = fitglme( T_A_CBV_S1BF, 'AUC ~ 1 + Day + (1|Mouse)' );
disp('Alcohol, CBV: DayK vs Day1')
disp(glme_A_CBV_S1BF)

%Alcohol, GFP
T_A_GFP_S1BF = T( T.Group == categorical("Alcohol") & T.Region == categorical("S1BF") & T.Signal == categorical("GFP"), : );
T_A_GFP_S1BF.Day = categorical(T_A_GFP_S1BF.Day, {'Day1','Day2','Day3','Day4','Day5'}, 'Ordinal', false);
glme_A_GFP_S1BF = fitglme( T_A_GFP_S1BF, 'AUC ~ 1 + Day + (1|Mouse)' );
disp('Alcohol, GFP: DayK vs Day1')
disp(glme_A_GFP_S1BF)

%% Between‐Group on Each Day (CBV vs GFP)

%Water vs Alcohol on Day1, CBV
T_D1_CBV = T( T.Day == categorical("Day1") & T.Signal == categorical("CBV"), : );
T_D1_CBV.Group = categorical(T_D1_CBV.Group, {'Alcohol','Water'}, 'Ordinal', false);
glme_D1_CBV = fitglme( T_D1_CBV, 'AUC ~ 1 + Group + (1|Mouse)' );
disp('--- Day1, CBV: Water vs Alcohol ---')
disp(glme_D1_CBV)

% Day2, CBV
T_D2_CBV = T( T.Day == categorical("Day2") & T.Signal == categorical("CBV"), : );
T_D2_CBV.Group = categorical(T_D2_CBV.Group, {'Alcohol','Water'}, 'Ordinal', false);
glme_D2_CBV = fitglme( T_D2_CBV, 'AUC ~ 1 + Group + (1|Mouse)' );
disp('--- Day2, CBV: Water vs Alcohol ---')
disp(glme_D2_CBV)

% Day3, CBV
T_D3_CBV = T( T.Day == categorical("Day3") & T.Signal == categorical("CBV"), : );
T_D3_CBV.Group = categorical(T_D3_CBV.Group, {'Alcohol','Water'}, 'Ordinal', false);
glme_D3_CBV = fitglme( T_D3_CBV, 'AUC ~ 1 + Group + (1|Mouse)' );
disp('--- Day3, CBV: Water vs Alcohol ---')
disp(glme_D3_CBV)

%Day4, CBV
T_D4_CBV = T( T.Day == categorical("Day4") & T.Signal == categorical("CBV"), : );
T_D4_CBV.Group = categorical(T_D4_CBV.Group, {'Alcohol','Water'}, 'Ordinal', false);
glme_D4_CBV = fitglme( T_D4_CBV, 'AUC ~ 1 + Group + (1|Mouse)' );
disp('--- Day4, CBV: Water vs Alcohol ---')
disp(glme_D4_CBV)

%Day5, CBV
T_D5_CBV = T( T.Day == categorical("Day5") & T.Signal == categorical("CBV"), : );
T_D5_CBV.Group = categorical(T_D5_CBV.Group, {'Alcohol','Water'}, 'Ordinal', false);
glme_D5_CBV = fitglme( T_D5_CBV, 'AUC ~ 1 + Group + (1|Mouse)' );
disp('--- Day5, CBV: Water vs Alcohol ---')
disp(glme_D5_CBV)

%%
%Day1, GFP
T_D1_GFP = T( T.Day == categorical("Day1") & T.Signal == categorical("GFP"), : );
T_D1_GFP.Group = categorical(T_D1_GFP.Group, {'Alcohol','Water'}, 'Ordinal', false);
glme_D1_GFP = fitglme( T_D1_GFP, 'AUC ~ 1 + Group + (1|Mouse)' );
disp('--- Day1, GFP: Water vs Alcohol ---')
disp(glme_D1_GFP)

%Day2, GFP
T_D2_GFP = T( T.Day == categorical("Day2") & T.Signal == categorical("GFP"), : );
T_D2_GFP.Group = categorical(T_D2_GFP.Group, {'Alcohol','Water'}, 'Ordinal', false);
glme_D2_GFP = fitglme( T_D2_GFP, 'AUC ~ 1 + Group + (1|Mouse)' );
disp('--- Day2, GFP: Water vs Alcohol ---')
disp(glme_D2_GFP)

%Day3, GFP
T_D3_GFP = T( T.Day == categorical("Day3") & T.Signal == categorical("GFP"), : );
T_D3_GFP.Group = categorical(T_D3_GFP.Group, {'Alcohol','Water'}, 'Ordinal', false);
glme_D3_GFP = fitglme( T_D3_GFP, 'AUC ~ 1 + Group + (1|Mouse)' );
disp('--- Day3, GFP: Water vs Alcohol ---')
disp(glme_D3_GFP)

%Day4, GFP
T_D4_GFP = T( T.Day == categorical("Day4") & T.Signal == categorical("GFP"), : );
T_D4_GFP.Group = categorical(T_D4_GFP.Group, {'Alcohol','Water'}, 'Ordinal', false);
glme_D4_GFP = fitglme( T_D4_GFP, 'AUC ~ 1 + Group + (1|Mouse)' );
disp('--- Day4, GFP: Water vs Alcohol ---')
disp(glme_D4_GFP)

%Day5, GFP
T_D5_GFP = T( T.Day == categorical("Day5") & T.Signal == categorical("GFP"), : );
T_D5_GFP.Group = categorical(T_D5_GFP.Group, {'Alcohol','Water'}, 'Ordinal', false);
glme_D5_GFP = fitglme( T_D5_GFP, 'AUC ~ 1 + Group + (1|Mouse)' );
disp('--- Day5, GFP: Water vs Alcohol ---')
disp(glme_D5_GFP)


end


