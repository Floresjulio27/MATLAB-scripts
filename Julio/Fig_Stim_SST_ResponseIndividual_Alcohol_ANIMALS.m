function [] = Fig_Stim_SST_ResponseIndividual_Alcohol_ANIMALS(rootFolder,saveState,delim)
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
            data.(currentGroup).(currentDay).GetTimeVector.(comp).mean_timeVector = ...
                data.(currentGroup).(currentDay).(firstID).(comp).timeVector;
        end
    end
end


% Compute Day-by-Day averages (and SEM) across animals

% ————————————————————————————————————————————————————————————————
% Compute Day-by-Day averages (and SEM) across animals. This is other way
% to organize it. 
% ————————————————————————————————————————————————————————————————
% Loop through groups  
% for ee = 1:length(groups)
%     currentGroup = groups{1,ee};
% 
%     % Loop through animals
%     for  xx = 1:length(days)
%         currentDay = days{1,xx};
%         % Skip if this day doesn't exist
%         if ~isfield(Results_Evoked.(currentGroup), currentDay)
%             continue;
%         end
%         FPanimalIDs = fieldnames(Results_Evoked.(currentGroup).(currentDay));
% 
%         if isempty(FPanimalIDs)
%             continue
%         end
% 
%         for ff = 1:length(compDataTypes)
%             comp = compDataTypes{ff};
% 
%             % gather all animals' raw traces for this day
%             allAChCBV = [];
%             allNECBV  = [];
%             allAChGFP = [];
%             allNEGFP  = [];
%             for aa = 1:length(FPanimalIDs)
%                 animalID = FPanimalIDs{aa};
%                 if ~isfield(data.(currentGroup).(currentDay), animalID)
%                     continue
%                 end
%                 S = data.(currentGroup).(currentDay).(animalID).(comp);
%                 allAChCBV = [allAChCBV; S.P_AChCBVRaw];   %#ok<AGROW>
%                 allNECBV  = [allNECBV;  S.P_NECBVRaw];    %#ok<AGROW>
%                 allAChGFP = [allAChGFP; S.P_AChGFPRaw];   %#ok<AGROW>
%                 allNEGFP  = [allNEGFP;  S.P_NEGFPRaw];    %#ok<AGROW>
%             end
% 
%             % compute and store stats
%             N = size(allAChCBV,1);
%             data.(currentGroup).(currentDay).(comp).P_AChCBV_N_Exp  = N;
%             data.(currentGroup).(currentDay).(comp).P_AChCBV_Mean   = mean(allAChCBV,   1);
%             data.(currentGroup).(currentDay).(comp).P_AChCBV_SEM    = std(allAChCBV,0,1) / sqrt(N);
% 
%             N = size(allNECBV,1);
%             data.(currentGroup).(currentDay).(comp).P_NECBV_N_Exp   = N;
%             data.(currentGroup).(currentDay).(comp).P_NECBV_Mean    = mean(allNECBV,    1);
%             data.(currentGroup).(currentDay).(comp).P_NECBV_SEM     = std(allNECBV,0,1)  / sqrt(N);
% 
%             N = size(allAChGFP,1);
%             data.(currentGroup).(currentDay).(comp).P_AChGFP_N_Exp  = N;
%             data.(currentGroup).(currentDay).(comp).P_AChGFP_Mean   = mean(allAChGFP,   1);
%             data.(currentGroup).(currentDay).(comp).P_AChGFP_SEM    = std(allAChGFP,0,1) / sqrt(N);
% 
%             N = size(allNEGFP,1);
%             data.(currentGroup).(currentDay).(comp).P_NEGFP_N_Exp   = N;
%             data.(currentGroup).(currentDay).(comp).P_NEGFP_Mean    = mean(allNEGFP,    1);
%             data.(currentGroup).(currentDay).(comp).P_NEGFP_SEM     = std(allNEGFP,0,1)  / sqrt(N);
%         end
%     end
% end

% --- Compute group stats from per-animal means ---
for ee = 1:length(groups)
    currentGroup = groups{ee};

    for xx = 1:length(days)
        currentDay = days{xx};

        % Skip missing days
        if ~isfield(data.(currentGroup), currentDay)
            continue;
        end

        FPanimalIDs = fieldnames(data.(currentGroup).(currentDay));
        if isempty(FPanimalIDs)
            continue;
        end

        for ff = 1:length(compDataTypes)
            comp = compDataTypes{ff};

            allAChCBV = [];
            allNECBV  = [];
            allAChGFP = [];
            allNEGFP  = [];

            for aa = 1:length(FPanimalIDs)
                animalID = FPanimalIDs{aa};
                if ~isfield(data.(currentGroup).(currentDay).(animalID), comp)
                    continue
                end

                S = data.(currentGroup).(currentDay).(animalID).(comp);

                if isfield(S, 'P_AChCBV'), allAChCBV = [allAChCBV; S.P_AChCBV]; end
                if isfield(S, 'P_NECBV'),  allNECBV  = [allNECBV;  S.P_NECBV];  end
                if isfield(S, 'P_AChGFP'), allAChGFP = [allAChGFP; S.P_AChGFP]; end
                if isfield(S, 'P_NEGFP'),  allNEGFP  = [allNEGFP;  S.P_NEGFP];  end
            end

            % store group stats (mean ± SEM across animals)
            N = size(allAChCBV,1);
            data.(currentGroup).(currentDay).(comp).P_AChCBV_N_Animals = N;
            data.(currentGroup).(currentDay).(comp).P_AChCBV_Mean      = mean(allAChCBV, 1);
            data.(currentGroup).(currentDay).(comp).P_AChCBV_SEM       = std(allAChCBV,0,1) / sqrt(N);

            N = size(allNECBV,1);
            data.(currentGroup).(currentDay).(comp).P_NECBV_N_Animals  = N;
            data.(currentGroup).(currentDay).(comp).P_NECBV_Mean       = mean(allNECBV, 1);
            data.(currentGroup).(currentDay).(comp).P_NECBV_SEM        = std(allNECBV,0,1) / sqrt(N);

            N = size(allAChGFP,1);
            data.(currentGroup).(currentDay).(comp).P_AChGFP_N_Animals = N;
            data.(currentGroup).(currentDay).(comp).P_AChGFP_Mean      = mean(allAChGFP, 1);
            data.(currentGroup).(currentDay).(comp).P_AChGFP_SEM       = std(allAChGFP,0,1) / sqrt(N);

            N = size(allNEGFP,1);
            data.(currentGroup).(currentDay).(comp).P_NEGFP_N_Animals  = N;
            data.(currentGroup).(currentDay).(comp).P_NEGFP_Mean       = mean(allNEGFP, 1);
            data.(currentGroup).(currentDay).(comp).P_NEGFP_SEM        = std(allNEGFP,0,1) / sqrt(N);
        end
    end
end




%% plot the comparison of GCaMP and TRITC with SST Ca2+ activity and TRITC


gfp_Day1_color = [0 0 0];
gfp_Day2_color = [0.38, 0.56, 0.74];
gfp_Day3_color = [0.867, 0.110, 0.467];
gfp_Day4_color = [0.70, 1.00, 1.0];
gfp_Day5_color = [0.55, 0.45, 0.70];
gfp_Day6_color = [0.70, 0.40, 0.10];

CBV_Day1_color = [0 0 0];
CBV_Day2_color = [0.38, 0.56, 0.74];
CBV_Day3_color = [0.867, 0.110, 0.467];
CBV_Day4_color = [0.70, 1.00, 1.0];
CBV_Day5_color = [0.55, 0.45, 0.70];
CBV_Day6_color = [0.70, 0.40, 0.10];

%% Contra stim PFC fiber GFP (Right hemisphere, left whisker stim) -
ax1 = subplot(2,2,1);
sgtitle('Stimulus Evoked Responses(Alcohol Group) ANIMALS')

% GFP day 1
x = data.Alcohol.Day1.GetTimeVector.Contra.mean_timeVector (:)';
y = data.Alcohol.Day1.Contra.P_NEGFP_Mean(:)';
e = data.Alcohol.Day1.Contra.P_NEGFP_SEM(:)';
fill([x fliplr(x)], [y+e fliplr(y-e)], gfp_Day1_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;
p1 = plot(x, y, 'Color', gfp_Day1_color, 'LineWidth', 2);

% GFP day 2 
x2 = data.Alcohol.Day2.GetTimeVector.Contra.mean_timeVector(:)';
y2 = data.Alcohol.Day2.Contra.P_NEGFP_Mean(:)';
e2 = data.Alcohol.Day2.Contra.P_NEGFP_SEM(:)';
fill([x2 fliplr(x2)], [y2+e2 fliplr(y2-e2)], gfp_Day2_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x2, y2, 'Color', gfp_Day2_color, 'LineWidth', 2);

% GFP day 3 
x3 = data.Alcohol.Day3.GetTimeVector.Contra.mean_timeVector(:)';
y3 = data.Alcohol.Day3.Contra.P_NEGFP_Mean(:)';
e3 = data.Alcohol.Day3.Contra.P_NEGFP_SEM(:)';
fill([x3 fliplr(x3)], [y3+e3 fliplr(y3-e3)], gfp_Day3_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p3 = plot(x3, y3, 'Color', gfp_Day3_color, 'LineWidth', 2);

% GFP day 4 
x4 = data.Alcohol.Day4.GetTimeVector.Contra.mean_timeVector(:)';
y4 = data.Alcohol.Day4.Contra.P_NEGFP_Mean(:)';
e4 = data.Alcohol.Day4.Contra.P_NEGFP_SEM(:)';
fill([x4 fliplr(x4)], [y4+e4 fliplr(y4-e4)], gfp_Day4_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p4 = plot(x4, y4, 'Color', gfp_Day4_color, 'LineWidth', 2);

% GFP day 5 
x5 = data.Alcohol.Day5.GetTimeVector.Contra.mean_timeVector(:)';
y5 = data.Alcohol.Day5.Contra.P_NEGFP_Mean(:)';
e5 = data.Alcohol.Day5.Contra.P_NEGFP_SEM(:)';
fill([x5 fliplr(x5)], [y5+e5 fliplr(y5-e5)], gfp_Day5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p5 = plot(x5, y5, 'Color', gfp_Day5_color, 'LineWidth', 2);

% GFP day 6
x6 = data.Alcohol.Day6.GetTimeVector.Contra.mean_timeVector(:)';
y6 = data.Alcohol.Day6.Contra.P_NEGFP_Mean(:)';
e6 = data.Alcohol.Day6.Contra.P_NEGFP_SEM(:)';
fill([x6 fliplr(x6)], [y6+e6 fliplr(y6-e6)], gfp_Day6_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p6 = plot(x6, y6, 'Color', gfp_Day6_color, 'LineWidth', 2);


title('Contra stim PFC fiber GFP')
ylabel('\DeltaF/F PFC SST Ca2+ activity (%)')
xlabel('Time (s)')
legend([p1 p2 p3 p4 p5 p6],'BL','DID1','DID2','DID3','DID4','Post-DID')
%legend([p1 p2 p3 p4 p5],'Week1','Week2','Week3','Week4','Week5')

axis square; 
axis tight;
xlim([-5 10]);

%% Contra stim PFC fiber CBV (Right hemisphere, left whisker stim) -
ax2 = subplot(2,2,2);

% CBV day 1
x = data.Alcohol.Day1.GetTimeVector.Contra.mean_timeVector (:)';
y = data.Alcohol.Day1.Contra.P_NECBV_Mean(:)';
e = data.Alcohol.Day1.Contra.P_NECBV_SEM(:)';
fill([x fliplr(x)], [y+e fliplr(y-e)], CBV_Day1_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;
p1 = plot(x, y, 'Color', CBV_Day1_color, 'LineWidth', 2);

% CBV day 2 
x2 = data.Alcohol.Day2.GetTimeVector.Contra.mean_timeVector(:)';
y2 = data.Alcohol.Day2.Contra.P_NECBV_Mean(:)';
e2 = data.Alcohol.Day2.Contra.P_NECBV_SEM(:)';
fill([x2 fliplr(x2)], [y2+e2 fliplr(y2-e2)], CBV_Day2_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x2, y2, 'Color', CBV_Day2_color, 'LineWidth', 2);

% CBV day 3 
x3 = data.Alcohol.Day3.GetTimeVector.Contra.mean_timeVector(:)';
y3 = data.Alcohol.Day3.Contra.P_NECBV_Mean(:)';
e3 = data.Alcohol.Day3.Contra.P_NECBV_SEM(:)';
fill([x3 fliplr(x3)], [y3+e3 fliplr(y3-e3)], CBV_Day3_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p3 = plot(x3, y3, 'Color', CBV_Day3_color, 'LineWidth', 2);

% CBV day 4 
x4 = data.Alcohol.Day4.GetTimeVector.Contra.mean_timeVector(:)';
y4 = data.Alcohol.Day4.Contra.P_NECBV_Mean(:)';
e4 = data.Alcohol.Day4.Contra.P_NECBV_SEM(:)';
fill([x4 fliplr(x4)], [y4+e4 fliplr(y4-e4)], CBV_Day4_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p4 = plot(x4, y4, 'Color', CBV_Day4_color, 'LineWidth', 2);

% CBV day 5 
x5 = data.Alcohol.Day5.GetTimeVector.Contra.mean_timeVector(:)';
y5 = data.Alcohol.Day5.Contra.P_NECBV_Mean(:)';
e5 = data.Alcohol.Day5.Contra.P_NECBV_SEM(:)';
fill([x5 fliplr(x5)], [y5+e5 fliplr(y5-e5)], CBV_Day5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p5 = plot(x5, y5, 'Color', CBV_Day5_color, 'LineWidth', 2);

% CBV day 6
x6 = data.Alcohol.Day6.GetTimeVector.Contra.mean_timeVector(:)';
y6 = data.Alcohol.Day6.Contra.P_NECBV_Mean(:)';
e6 = data.Alcohol.Day6.Contra.P_NECBV_SEM(:)';
fill([x6 fliplr(x6)], [y6+e6 fliplr(y6-e6)], CBV_Day6_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p6 = plot(x6, y6, 'Color', CBV_Day6_color, 'LineWidth', 2);


title('Contra stim PFC fiber CBV')
ylabel('\DeltaF/F PFC CBV (%)')
xlabel('Time (s)')
legend([p1 p2 p3 p4 p5 p6],'BL','DID1','DID2','DID3','DID4','Post-DID')
%legend([p1 p2 p3 p4 p5],'Week1','Week2','Week3','Week4','Week5')
axis square; 
axis tight;
xlim([-5 10]);

%% Contra stim S1BF fiber GFP (Left hemisphere, Right whisker stim) -
ax3 = subplot(2,2,3);

% GFP day 1
x = data.Alcohol.Day1.GetTimeVector.Contra.mean_timeVector (:)';
y = data.Alcohol.Day1.Contra.P_AChGFP_Mean(:)';
e = data.Alcohol.Day1.Contra.P_AChGFP_SEM(:)';
fill([x fliplr(x)], [y+e fliplr(y-e)], gfp_Day1_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;
p1 = plot(x, y, 'Color', gfp_Day1_color, 'LineWidth', 2);

% GFP day 2 
x2 = data.Alcohol.Day2.GetTimeVector.Contra.mean_timeVector(:)';
y2 = data.Alcohol.Day2.Contra.P_AChGFP_Mean(:)';
e2 = data.Alcohol.Day2.Contra.P_AChGFP_SEM(:)';
fill([x2 fliplr(x2)], [y2+e2 fliplr(y2-e2)], gfp_Day2_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x2, y2, 'Color', gfp_Day2_color, 'LineWidth', 2);

% GFP day 3 
x3 = data.Alcohol.Day3.GetTimeVector.Contra.mean_timeVector(:)';
y3 = data.Alcohol.Day3.Contra.P_AChGFP_Mean(:)';
e3 = data.Alcohol.Day3.Contra.P_AChGFP_SEM(:)';
fill([x3 fliplr(x3)], [y3+e3 fliplr(y3-e3)], gfp_Day3_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p3 = plot(x3, y3, 'Color', gfp_Day3_color, 'LineWidth', 2);

% GFP day 4 
x4 = data.Alcohol.Day4.GetTimeVector.Contra.mean_timeVector(:)';
y4 = data.Alcohol.Day4.Contra.P_AChGFP_Mean(:)';
e4 = data.Alcohol.Day4.Contra.P_AChGFP_SEM(:)';
fill([x4 fliplr(x4)], [y4+e4 fliplr(y4-e4)], gfp_Day4_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p4 = plot(x4, y4, 'Color', gfp_Day4_color, 'LineWidth', 2);

% GFP day 5 
x5 = data.Alcohol.Day5.GetTimeVector.Contra.mean_timeVector(:)';
y5 = data.Alcohol.Day5.Contra.P_AChGFP_Mean(:)';
e5 = data.Alcohol.Day5.Contra.P_AChGFP_SEM(:)';
fill([x5 fliplr(x5)], [y5+e5 fliplr(y5-e5)], gfp_Day5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p5 = plot(x5, y5, 'Color', gfp_Day5_color, 'LineWidth', 2);

% GFP day 6
x6 = data.Alcohol.Day6.GetTimeVector.Contra.mean_timeVector(:)';
y6 = data.Alcohol.Day6.Contra.P_AChGFP_Mean(:)';
e6 = data.Alcohol.Day6.Contra.P_AChGFP_SEM(:)';
fill([x6 fliplr(x6)], [y6+e6 fliplr(y6-e6)], gfp_Day6_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p6 = plot(x6, y6, 'Color', gfp_Day6_color, 'LineWidth', 2);


title('Contra stim S1BF fiber GFP')
ylabel('\DeltaF/F S1BF SST Ca2+ activity (%)')
xlabel('Time (s)')
legend([p1 p2 p3 p4 p5 p6],'BL','DID1','DID2','DID3','DID4','Post-DID')
%legend([p1 p2 p3 p4 p5],'Week1','Week2','Week3','Week4','Week5')
axis square; 
axis tight;
xlim([-5 10]);

%% Contra stim S1BF fiber CBV (Left hemisphere, right whisker stim) -
ax4 = subplot(2,2,4);

% CBV day 1
x = data.Alcohol.Day1.GetTimeVector.Contra.mean_timeVector (:)';
y = data.Alcohol.Day1.Contra.P_AChCBV_Mean(:)';
e = data.Alcohol.Day1.Contra.P_AChCBV_SEM(:)';
fill([x fliplr(x)], [y+e fliplr(y-e)], CBV_Day1_color,'FaceAlpha', 0.2, 'EdgeColor', 'none'); 
hold on;
p1 = plot(x, y, 'Color', CBV_Day1_color, 'LineWidth', 2);

% CBV day 2 
x2 = data.Alcohol.Day2.GetTimeVector.Contra.mean_timeVector(:)';
y2 = data.Alcohol.Day2.Contra.P_AChCBV_Mean(:)';
e2 = data.Alcohol.Day2.Contra.P_AChCBV_SEM(:)';
fill([x2 fliplr(x2)], [y2+e2 fliplr(y2-e2)], CBV_Day2_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x2, y2, 'Color', CBV_Day2_color, 'LineWidth', 2);

% CBV day 3 
x3 = data.Alcohol.Day3.GetTimeVector.Contra.mean_timeVector(:)';
y3 = data.Alcohol.Day3.Contra.P_AChCBV_Mean(:)';
e3 = data.Alcohol.Day3.Contra.P_AChCBV_SEM(:)';
fill([x3 fliplr(x3)], [y3+e3 fliplr(y3-e3)], CBV_Day3_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p3 = plot(x3, y3, 'Color', CBV_Day3_color, 'LineWidth', 2);

% CBV day 4 
x4 = data.Alcohol.Day4.GetTimeVector.Contra.mean_timeVector(:)';
y4 = data.Alcohol.Day4.Contra.P_AChCBV_Mean(:)';
e4 = data.Alcohol.Day4.Contra.P_AChCBV_SEM(:)';
fill([x4 fliplr(x4)], [y4+e4 fliplr(y4-e4)], CBV_Day4_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p4 = plot(x4, y4, 'Color', CBV_Day4_color, 'LineWidth', 2);

% CBV day 5 
x5 = data.Alcohol.Day5.GetTimeVector.Contra.mean_timeVector(:)';
y5 = data.Alcohol.Day5.Contra.P_AChCBV_Mean(:)';
e5 = data.Alcohol.Day5.Contra.P_AChCBV_SEM(:)';
fill([x5 fliplr(x5)], [y5+e5 fliplr(y5-e5)], CBV_Day5_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p5 = plot(x5, y5, 'Color', CBV_Day5_color, 'LineWidth', 2);

% CBV day 6
x6 = data.Alcohol.Day6.GetTimeVector.Contra.mean_timeVector(:)';
y6 = data.Alcohol.Day6.Contra.P_AChCBV_Mean(:)';
e6 = data.Alcohol.Day6.Contra.P_AChCBV_SEM(:)';
fill([x6 fliplr(x6)], [y6+e6 fliplr(y6-e6)], CBV_Day6_color,'FaceAlpha', 0.2, 'EdgeColor', 'none');
p6 = plot(x6, y6, 'Color', CBV_Day6_color, 'LineWidth', 2);


title('Contra stim PFC fiber CBV')
ylabel('\DeltaF/F S1BF CBV (%)')
xlabel('Time (s)')
legend([p1 p2 p3 p4 p5 p6],'BL','DID1','DID2','DID3','DID4','Post-DID')
%legend([p1 p2 p3 p4 p5],'Week1','Week2','Week3','Week4','Week5')
axis square; 
ylim([-1 1])
xlim([-5 10]);
%% Save figure
if saveState == true
    dirpath = [rootFolder delim 'MATLAB Figures Water vs Alcohol' delim 'PuffEvokedResponses'];
    if ~exist(dirpath,'dir')
        mkdir(dirpath);
    end
    savefig(fullfile(dirpath, 'Stimulus_Evoked_Responses_Alcohol_ANIMALS'));
end
