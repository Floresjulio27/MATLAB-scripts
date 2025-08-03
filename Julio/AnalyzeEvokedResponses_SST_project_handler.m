function [] = AnalyzeEvokedResponses_SST_project_handler(rootFolder,delim,runFromStart)
%----------------------------------------------------------------------------------------------------------
% Written by Kevin L. Turner
% Pennsylvania State University, Dept. of Biomedical Engineering
% https://github.com/KL-Turner

%Modified by: Julio Flores-Cuadra
%----------------------------------------------------------------------------------------------------------

% 1) Initialize as structs so isfield works:
if runFromStart == true
    Results_Evoked.Water   = struct();
    Results_Evoked.Alcohol = struct();
else
    cd(fullfile(rootFolder,'Results_SST_project'));
    if exist('Results_Evoked.mat','file')
        load('Results_Evoked.mat','Results_Evoked')
    else
        Results_Evoked.Water   = struct();
        Results_Evoked.Alcohol = struct();
    end
end

% 2) Define groups, set, days
dataPath = fullfile(rootFolder,'Data');
groups   = {'Water','Alcohol'};
setName  = 'Analysis'; %previously evoked 
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
%multiWaitbar('Analyzing evoked responses', 0, 'Color','P');
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
      if ~isfield(Results_Evoked.(grp), dayName) || ~isfield(Results_Evoked.(grp).(dayName), aID)

        Results_Evoked = AnalyzeEvokedResponses_SST_project(aID, grp, setName, rootFolder, delim, Results_Evoked, dayName);

        % ensure we’re back in Data folder
        cd(dataPath);
      end

      multiWaitbar('Analyzing evoked responses', cc/waitBarLength);
      pause(0.1);
      cc = cc + 1;
    end
  end
end
multiWaitbar('close all');
