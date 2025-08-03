function [] = AnalyzeArousalTransitions_SST_project_handler(rootFolder,delim,runFromStart)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner

%Modified by: Julio Flores-Cuadra
%----------------------------------------------------------------------------------------------------------

% 1) Initialize as structs so isfield works:
if runFromStart == true
    Results_Transitions.Water   = struct();
    Results_Transitions.Alcohol = struct();
else
    cd(fullfile(rootFolder,'Results_SST_project'));
    if exist('Results_Transitions.mat','file')
        load('Results_Transitions.mat','Results_Transitions')
    else
        Results_Transitions.Water   = struct();
        Results_Transitions.Alcohol = struct();
    end
end

% 2) Define groups, set, days
dataPath = fullfile(rootFolder,'Data');
groups   = {'Water','Alcohol'};
setName  = 'Analysis';
days     = {'Day1','Day2','Day3','Day4','Day5','Day6'};

% 3) Compute total iterations for the waitbar
waitBarLength = 0;
for ia = 1:length(groups)
  for id = 1:length(days)
    dpath = fullfile(dataPath, groups{ia}, setName, days{id});
    if isfolder(dpath)
      dd = dir(dpath);
      dd = dd([dd.isdir] & ~startsWith({dd.name},'.'));
      waitBarLength = waitBarLength + length(dd);
    end
  end
end

% 4) Run analysis
cc = 1;
multiWaitbar('Analyzing Arousal State Transitions', 0, 'Color','P');
for ia = 1:length(groups)
  for id = 1:length(days)
    dayPath = fullfile(dataPath, groups{ia}, setName, days{id});
    if ~isfolder(dayPath)
      continue
    end
    dd = dir(dayPath);
    dd = dd([dd.isdir] & ~startsWith({dd.name},'.'));
    animalIDs = {dd.name};
    for ib = 1:length(animalIDs)
      grp      = groups{ia};
      dayName  = days{id};
      aID      = animalIDs{ib};

      % only run if this day/animal hasn’t been done yet
      if ~isfield(Results_Transitions.(grp), dayName) || ~isfield(Results_Transitions.(grp).(dayName), aID)

        Results_Transitions = AnalyzeArousalTransitions_SST_project(aID, grp, setName, rootFolder, delim, Results_Transitions, dayName);

        % ensure we’re back in Data folder
        cd(dataPath);
      end

      multiWaitbar('Analyzing Arousal State Transitions', cc/waitBarLength);
      pause(0.1);
      cc = cc + 1;
    end
  end
end
multiWaitbar('close all');