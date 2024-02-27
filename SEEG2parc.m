function elecInfo=SEEG2parc(cfg)

% function elecInfo=SEEG2parc(cfg)
% Attribute cerebral parcellation labels to stereo-EEG electrodes.
%
% In its current version, the function is compatible both with iELVis and
% with the BIDS output of Voxeloc. In order to identify the folder where
% the relevant information is, it first opens a system UI box so that the
% user directs to the correct folder.
%
% Inputs:
%
%   cfg:            (structure) configuration strucure with the following
%                   fields:
%
%     dataType:     (string) either 'iELVis' or 'BIDS'
%
%     parcType:     (string) either 'DK' (default) for the Desikan-Killiany
%                   parcellation or 'Des' for the Destrieux parcellation
%
%     tissueType:   (string) 'cortex' if only cortical tissue
%                   hippocampus and amygdala) is to be considered,
%                   otherwise (default) all FreeSurfer parcellation labels
%                   are used
%
% Outputs:
%
%   elecInfo:       N-by-7 cell array, where N is the number of electrodes.
%
%     column 1:     electrode name
%
%     column 2:     electrode type ('D' or 'depth')
%
%     column 3:     hemisphere ('L' or 'R')
%
%     column 4:     electrode coordinates, in patient space, in ILA
%                   coordinates
%
%     column 5:     tissue label(s), deterministic; returns 'unknown' if no
%                   valid tissue label was found near the electrode
%
%     column 6:     tissue label(s), probabilistic; if there is only 1
%                   label, column 6 is equal to column 5; returns 'unknown'
%                   if no valid tissue label was found near the electrode
%
%     column 7:     weight of each tissue label in column 6 ([] if tissue
%                   label is 'unknown')
%
% The attribution of tissue labels is based on a 3x3x3-voxel cube centered
% around the electrode coordinate. Tissue labels are weighted like so:
%
%          +--------+--------+--------+
%         /    _   /    _   /    _   /|
%        /   \/3  /   \/2  /   \/3  / |
%       +--------+--------+--------+  |
%      /    _   /        /    _   /|  |
%     /   \/2  /    1   /   \/2  / |  +
%    +--------+--------+--------+  | /
%   /    _   /    _   /    _   /|  |/
%  /   \/3  /   \/2  /   \/3  / |  +
% +--------+--------+--------+  | /
% |        |        |        |  |/
% |        |        |        |  +
% |        |        |        | /
% |        |        |        |/
% +--------+--------+--------+
%
%          +--------+--------+--------+
%         /    _   /        /    _   /|
%        /   \/2  /    1   /   \/2  / |
%       +--------+--------+--------+  |
%      /        /        /        /|  |
%     /    1   /    1   /    1   / |  +
%    +--------+--------+--------+  | /
%   /    _   /        /    _   /|  |/
%  /   \/2  /    1   /   \/2  / |  +
% +--------+--------+--------+  | /
% |        |        |        |  |/
% |        |        |        |  +
% |        |        |        | /
% |        |        |        |/
% +--------+--------+--------+
%
%          +--------+--------+--------+
%         /    _   /    _   /    _   /|
%        /   \/3  /   \/2  /   \/3  / |
%       +--------+--------+--------+  |
%      /    _   /        /    _   /|  |
%     /   \/2  /    1   /   \/2  / |  +
%    +--------+--------+--------+  | /
%   /    _   /    _   /    _   /|  |/
%  /   \/3  /   \/2  /   \/3  / |  +
% +--------+--------+--------+  | /
% |        |        |        |  |/
% |        |        |        |  +
% |        |        |        | /
% |        |        |        |/
% +--------+--------+--------+
%
%
% https://github.com/HumanNeuronLab/SEEG2parc
% Pierre Mégevand, Human Neuron Lab, University of Geneva, Switzerland. 2022-2024.
% pierre.megevand@unige.ch; https://www.unige.ch/medecine/neucli/en/groupes-de-recherche/1034megevand/

if ~isfield(cfg,'dataType'), error('Please specify data type (iELVis or BIDS).'); end
if ~isfield(cfg,'parcType'), parcType='DK'; else, parcType=cfg.parcType; end
if ~isfield(cfg,'tissueType'), tissueType='all'; else, tissueType=cfg.tissueType; end

switch parcType
    case 'DK' % Desikan-Killiany parcellation: Desikan et al., Neuroimage 2006; https://doi.org/10.1016/j.neuroimage.2006.01.021
        myParcVol='aparc+aseg.mgz';
        onlyCtx=[17,18,53,54,1001:1035,2001:2035]'; % only hippocampus, amygdala and cortical labels from DK parcellation
    case 'Des' % Destrieux parcellation: Destrieux et al., Neuroimage 2010; https://doi.org/10.1016/j.neuroimage.2010.06.010
        myParcVol='aparc.a2009s+aseg.mgz';
        onlyCtx=[17,18,53,54,11101:11175,12101:12175]'; % only hippocampus, amygdala and cortical labels from Destrieux parcellation
    otherwise
        error('I do not recognize this parcellation: %s.',parcType);
end

filepath=uigetdir([],'Select patient folder');
[~,patID]=fileparts(filepath);


switch cfg.dataType
    
    case 'iELVis'
        
        % load electrode names
        elecNames=readiELVisElecNames(patID);
        elecType=elecNames(:,2);
        elecHem=elecNames(:,3);
        elecNames=elecNames(:,1);
        
        % load electrode coordinates (LIP)
        elecCoord=readiELVisElecCoord(patID,'LEPTOVOX');
        
        % load volume parcellations (ILA)
        wmparc=MRIread(fullfile(filepath,'mri','wmparc.mgz')); % load white matter parcellation
        parcVol=MRIread(fullfile(filepath,'mri',myParcVol)); % load volume parcellation
        parcVol.vol(parcVol.vol==2|parcVol.vol==41)=wmparc.vol(parcVol.vol==2|parcVol.vol==41); % replace white matter labels with DK white matter parcellation
        
    case 'BIDS'
        
        % load electrode table
        elecTable=readtable(fullfile(filepath,'ieeg',[patID '_electrodes.tsv']),'FileType','text','Delimiter','\t','ReadVariableNames',true,'ReadRowNames',false);
        elecNames=elecTable.name;
        elecType=elecTable.type;
        elecHem=elecTable.hemisphere;
        
        % voxeloc coords are LIA
        elecCoord=table2array(elecTable(:,2:4));
        elecCoord(:,3)=256-elecCoord(:,3); % make them LIP like iELVis
        
        % load volume parcellations (ILA)
        wmparc=MRIread(fullfile(filepath,'anat',[patID '_wmparc.mgz'])); % load white matter parcellation
        parcVol=MRIread(fullfile(filepath,'anat',[patID '_' myParcVol])); % load volume parcellation
        parcVol.vol(parcVol.vol==2|parcVol.vol==41)=wmparc.vol(parcVol.vol==2|parcVol.vol==41); % replace white matter labels with DK white matter parcellation
end

nElec=size(elecCoord,1);

% start populating elecInfo
elecInfo=cell(nElec,7);
elecInfo(:,1)=elecNames;
elecInfo(:,2)=elecType;
elecInfo(:,3)=elecHem;
elecInfo(:,4)=num2cell(elecCoord,2); % storing the original coordinates (LIP)

% match elecCoord dimension order and direction with scans
LIP2ILA=[0 1 0 0;1 0 0 0;0 0 -1 256;0 0 0 1];
elecCoord=(LIP2ILA*[elecCoord';ones(1,nElec)])';
elecCoord=elecCoord(:,1:3);

% load FreeSurfer lookup table
fid=fopen('FreeSurferColorLUTnoFormat.txt','r'); % file provided with iELVis
fsLut=textscan(fid,'%d %s %d %d %d %d');
fsLut=cat(2,num2cell(fsLut{1}),fsLut{2},num2cell([fsLut{3:5}],2)); % ugly code to rearrange cell array
fclose(fid);clear('fid');
fsLutIdx=[fsLut{:,1}]'; % get FreeSurfer indices for parcellations

% define weights for neighboring voxels
voxelWeights=1/cat(3, ...
    [sqrt(3) sqrt(2) sqrt(3); sqrt(2) 1 sqrt(2); sqrt(3) sqrt(2) sqrt(3)], ...
    [sqrt(2) 1 sqrt(2); 1 1 1; sqrt(2) 1 sqrt(2)], ...
    [sqrt(3) sqrt(2) sqrt(3); sqrt(2) 1 sqrt(2); sqrt(3) sqrt(2) sqrt(3)]);

for ctElec=1:nElec
    if strncmpi(elecType{ctElec},'D',1) % double-check that this is a depth electrode
        
        voxelCubeCoord=[ ...
            round(elecCoord(ctElec,1))-1, round(elecCoord(ctElec,2))-1, round(elecCoord(ctElec,3))-1; ...
            round(elecCoord(ctElec,1))+1, round(elecCoord(ctElec,2))+1, round(elecCoord(ctElec,3))+1];
        
        % REMmed by Pierre 20220721
        %{
        % deal with case where voxel cube would reach beyond the limits of the MRI (unlikely in practice)
        voxelCubeCoord(voxelCubeCoord<1)=1;
        voxelCubeCoord(voxelCubeCoord>256)=256;
        %}
        
        voxelCube=parcVol.vol(voxelCubeCoord(1,1):voxelCubeCoord(2,1),voxelCubeCoord(1,2):voxelCubeCoord(2,2),voxelCubeCoord(1,3):voxelCubeCoord(2,3));
        
        % REMmed by Pierre 20220721
        %{
        % deal with case where voxel cube would reach to the other cerebral hemisphere
        switch elecNames{ctElec,3}
            case 'L' % left hemisphere
                voxelCube(~ismember(voxelCube,myLOI_L))=0;
            case 'R' % right hemisphere
                voxelCube(~ismember(voxelCube,myLOI_R))=0;
        end
        %}
        
        % if only cortical tissue is required by user
        if strcmp(tissueType,'cortex')
            voxelCube(~ismember(voxelCube,onlyCtx))=NaN;
            % in case there is no cortical tissue
            if all(isnan(voxelCube),'all')
                elecInfo{ctElec,5}='unknown';
                elecInfo{ctElec,6}='unknown';
                elecInfo{ctElec,7}=[];
                continue
            end
        end
        
        % weighted score for each label's occurrence
        voxelCubeUnique=unique(voxelCube);
        voxelCubeUnique(isnan(voxelCubeUnique))=[];
        voxelCubeUniqueWeight=zeros(size(voxelCubeUnique));
        for ctUniq=1:numel(voxelCubeUniqueWeight)
            voxelCubeUniqueWeight(ctUniq)=sum(sum(sum((voxelCube==voxelCubeUnique(ctUniq)).*voxelWeights)));
        end
        
        % sort weights and labels in descending order
        [voxelCubeUniqueWeightSorted,sortIdx]=sort(voxelCubeUniqueWeight,'descend');
        voxelCubeUniqueSorted=voxelCubeUnique(sortIdx);
        
        % express weights in fractions of 1
        voxelCubeUniqueWeightSorted=voxelCubeUniqueWeightSorted./sum(voxelCubeUniqueWeightSorted(:));
        
        % get FreeSurfer labels
        nLabels=numel(voxelCubeUniqueSorted);
        voxelCubeLabels=cell(nLabels,1);
        for ctLbl=1:nLabels
            voxelCubeLabels{ctLbl}=fsLut{voxelCubeUniqueSorted(ctLbl)==fsLutIdx,2};
        end
        
        elecInfo{ctElec,5}=voxelCubeLabels{1};
        elecInfo{ctElec,6}=voxelCubeLabels;
        elecInfo{ctElec,7}=voxelCubeUniqueWeightSorted;
        
    else % this electrode is not a depth
        elecInfo{ctElec,5}='unknown';
        elecInfo{ctElec,6}='unknown';
        elecInfo{ctElec,7}=[];
    end
    
end

%% plotting (useful for debug)

%{
switch cfg.dataType
    case 'iELVis'
        CT=MRIread(fullfile(filepath,'elec_recon','postInPre.nii'));  % load post-implant CT
    case 'BIDS'
        CT=MRIread(fullfile(filepath,'anat',[patID '_postInPre.nii'])); % load post-implant CT
end
elecCoordRound=round(elecCoord);
figure;
for ctElec=1:size(elecCoord,1)
    colormap bone;
    subplot(1,3,1);
    imagesc(squeeze(CT.vol(elecCoordRound(ctElec,1),:,:)));
    hold all;
    scatter(elecCoordRound(ctElec,3),elecCoordRound(ctElec,2),'MarkerEdgeColor',[1 1 0],'MarkerFaceColor',[1 1 0],'MarkerFaceAlpha',0.3);
    hold off;
    axis image;
    xlabel('P to A');
    ylabel('R to L');
    subplot(1,3,2);
    imagesc(squeeze(CT.vol(:,elecCoordRound(ctElec,2),:)));
    hold all;
    scatter(elecCoordRound(ctElec,3),elecCoordRound(ctElec,1),'MarkerEdgeColor',[1 1 0],'MarkerFaceColor',[1 1 0],'MarkerFaceAlpha',0.3);
    hold off;
    axis image;
    xlabel('P to A');
    ylabel('S to I');
    title(elecNames{ctElec});
    subplot(1,3,3);
    imagesc(squeeze(CT.vol(:,:,elecCoordRound(ctElec,3))));
    hold all;
    scatter(elecCoordRound(ctElec,2),elecCoordRound(ctElec,1),'MarkerEdgeColor',[1 1 0],'MarkerFaceColor',[1 1 0],'MarkerFaceAlpha',0.3);
    hold off;
    axis image;
    xlabel('R to L');
    ylabel('S to I');
    pause;
end
%}
