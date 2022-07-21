function elecNames=readiELVisElecNames(fsSubj)

% pierre.megevand@gmail.com

global globalFsDir;

fileID=fopen(fullfile(globalFsDir,fsSubj,'elec_recon',[fsSubj '.electrodeNames']),'r');
elecNames=textscan(fileID,'%s','HeaderLines',2);
elecNames=elecNames{1};
elecNames=reshape(elecNames,3,numel(elecNames)/3)';

fclose(fileID);
