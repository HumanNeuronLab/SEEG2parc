function elecCoord=readiELVisElecCoord(datadir,patID,elecCoordType)

% pierre.megevand@unige.ch

fileID=fopen(fullfile(datadir,patID,'elec_recon',[patID '.' upper(elecCoordType)]),'r');
elecCoord=textscan(fileID,'%f','HeaderLines',2);
elecCoord=cell2mat(elecCoord);
elecCoord=reshape(elecCoord',3,numel(elecCoord)/3)';

fclose(fileID);
