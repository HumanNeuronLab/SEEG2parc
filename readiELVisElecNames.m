function elecNames=readiELVisElecNames(datadir,patID)

% pierre.megevand@unige.ch

fileID=fopen(fullfile(datadir,patID,'elec_recon',[patID '.electrodeNames']),'r');
elecNames=textscan(fileID,'%s','HeaderLines',2);
elecNames=elecNames{1};
elecNames=reshape(elecNames,3,numel(elecNames)/3)';

fclose(fileID);
