function elecCoord=readiELVisElecCoord(fsSubj,elecCoordType)

% pierre.megevand@gmail.com

global globalFsDir;

fileID=fopen(fullfile(globalFsDir,fsSubj,'elec_recon',[fsSubj '.' upper(elecCoordType)]),'r');
elecCoord=textscan(fileID,'%f','HeaderLines',2);
elecCoord=cell2mat(elecCoord);
elecCoord=reshape(elecCoord',3,numel(elecCoord)/3)';

fclose(fileID);
