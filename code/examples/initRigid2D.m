function [options,prams] = initRigid2D(options,prams)
% set the path and assign options and prams to default values if they
% have not been assigned

P = path; ii = find(pwd == filesep); ii = ii(end);
subPath = pwd; subPath = [subPath(1:ii) 'src'];
if isempty(strfind(P, subPath))
  addpath(subPath)
end

PramList = {'N','nb','T','number_steps'};
defaultPram.N = 64;
defaultPram.nb = 1;
defaultPram.T = 1;
defaultPram.rho = 1.0;
defaultPram.number_steps = 10;

for k = 1:length(PramList)
  if ~isfield(prams,PramList{k})
    eval(['prams.' PramList{k} '=defaultPram.' PramList{k} ';'])
    % Set any unassigned parameters to a default value
  end
end


OptionList = {'timeOrder','inear','farField','verbose','usePreco',...
      'gmresTol','saveData','usePlot','plotAxis'};
    
defaultOption.timeOrder = 1;
defaultOption.inear = true;
defaultOption.farField = 'shear';
defaultOption.verbose = true;
defaultOption.usePreco = false;
defaultOption.gmresTol = 1e-6;
defaultOption.saveData = true;
defaultOption.usePlot = true;
defaultOption.plotAxis = [-1 1 -1 1];

for k = 1:length(OptionList)
  if ~isfield(options,OptionList{k})
    eval(['options.' OptionList{k} '=defaultOption.' ...
        OptionList{k} ';'])
    % Set any unassigned parameters to a default value
  end
end


