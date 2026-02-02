function [elecTable,tissueLabels,tissueWeights]=SEEG2parc(cfg)

% function elecTable=SEEG2parc(cfg)
% Attribute cerebral parcellation labels to intracranial EEG electrodes.
%
% In its current version, the function is compatible both with iELVis and
% with the BIDS-iEEG output of Voxeloc. However, the attribution of
% parcellation labels to subdural electrodes (ECOG) only works with the
% iELVis data type.
%
% Inputs:
%
%   cfg:            (structure) configuration strucure with the following
%                   fields:
%
%     filePath:     (string) full path to folder containing patient's iEEG 
%                   anatomical information (no trailing slash). If missing,
%                   a system UI box opens so that the user directs to the
%                   correct folder.
%
%     dataType:     (string) either 'iELVis' or 'BIDS'
%
%     parcType:     (string) either 'DK' (default) for the Desikan-Killiany
%                   parcellation or 'Des' for the Destrieux parcellation
%
% Outputs:
%
%   elecTable:      N-by-K table, where N is the number of electrodes. The
%                   table follows the BIDS-iEEG convention. Therefore, the 
%                   number of variables K is flexible. However, the first 5
%                   columns are REQUIRED and must be provided in this exact
%                   order.
%
%     column 1:     electrode name
%
%     column 2:     electrode coordinate x
%
%     column 3:     electrode coordinate y
%
%     column 4:     electrode coordinate z
%
%     column 5:     electrode size (surface in mm^2)
%                   if not provided (iELVis): takes the value NaN
%
%     column ki:    hemisphere (L or R)
%
%     column ki:    electrode type (e.g. D or depth)
%
%     column ki:    electrode material
%                   if not provided (iELVis): the variable is not added to 
%                   the table
%
%     column ki:    electrode manufacturer
%                   if not provided (iELVis): the variable is not added to 
%                   the table
%
%     column ki:    electrode group (derived from electrode names)
%                   if not provided (iELVis): the variable is not added to 
%                   the table
%
%     column ki:    dimension of electrode group (e.g. [1x12])
%                   if not provided (iELVis): the variable is not added to 
%                   the table
%
%     column ki:    tissue label(s)
%                   returns 'unknown' if no valid tissue label was found 
%                   anywhere near the electrode
%
%     column ki:    tissue weights(s)
%
%   tissueLabels:   provided as a N-by-1 cell array of strings
%
%   tissueWeights:  provided as a N-by-1 numeric cell array
%
%
% For stereo-EEG electrodes, the attribution of tissue labels is based on a
% 3x3x3-voxel cube centered around the electrode coordinate. Tissue labels 
% are weighted like so:
%
%          +--------+--------+--------+
%         /    _   /    _   /    _   /|
%        /   \/3  /   \/2  /   \/3  / |
%       +--------+--------+--------+  |
%      /    _   /        /    _   /|  +
%     /   \/2  /    1   /   \/2  / | /
%    +--------+--------+--------+  |/
%   /    _   /    _   /    _   /|  +
%  /   \/3  /   \/2  /   \/3  / | /
% +--------+--------+--------+  |/
% |        |        |        |  +
% |        |        |        | /
% |        |        |        |/
% +--------+--------+--------+
%
%          +--------+--------+--------+
%         /    _   /        /    _   /|
%        /   \/2  /    1   /   \/2  / |
%       +--------+--------+--------+  |
%      /        /        /        /|  +
%     /    1   /    1   /    1   / | /
%    +--------+--------+--------+  |/
%   /    _   /        /    _   /|  +
%  /   \/2  /    1   /   \/2  / | /
% +--------+--------+--------+  |/
% |        |        |        |  +
% |        |        |        | /
% |        |        |        |/
% +--------+--------+--------+
%
%          +--------+--------+--------+
%         /    _   /    _   /    _   /|
%        /   \/3  /   \/2  /   \/3  / |
%       +--------+--------+--------+  |
%      /    _   /        /    _   /|  +
%     /   \/2  /    1   /   \/2  / | /
%    +--------+--------+--------+  |/
%   /    _   /    _   /    _   /|  +
%  /   \/3  /   \/2  /   \/3  / | /
% +--------+--------+--------+  |/
% |        |        |        |  +
% |        |        |        | /
% |        |        |        |/
% +--------+--------+--------+
%
% For subdural electrodes, the attribution of tissue labels is based on a
% 4-mm sphere surrounding the surface vertex closest to the electrode
% coordinates. Tissue labels within 2 mm have a weight of 1, those between
% 2 and 3 mm have a weight of 1/sqrt(2), those further than 3 mm have a
% weight of 1/sqrt(3).
%
% https://github.com/HumanNeuronLab/SEEG2parc
% Pierre Megevand, Human Neuron Lab, University of Geneva, Switzerland. 2022-2025.
% pierre.megevand@unige.ch; https://www.unige.ch/medecine/neucli/en/groupes-de-recherche/1034megevand/

if ~isfield(cfg,'dataType'), error('Please specify data type (iELVis or BIDS).'); end
if ~isfield(cfg,'parcType'), parcType='DK'; else, parcType=cfg.parcType; end
if ~isfield(cfg,'filePath'), filePath=uigetdir([],'Select patient folder'); else, filePath=cfg.filePath; end
if ~exist(filePath,'dir'), error('Patient folder not found.'); end

[datadir,patID]=fileparts(filePath);

switch parcType
    case 'DK' % Desikan-Killiany parcellation: Desikan et al., Neuroimage 2006; https://doi.org/10.1016/j.neuroimage.2006.01.021
        myParcVol='aparc+aseg';
        myParcSurf='h.aparc.annot';
        %onlyCtx=[17,18,53,54,1001:1035,2001:2035]'; % only hippocampus, amygdala and cortical labels from DK parcellation
    case 'Des' % Destrieux parcellation: Destrieux et al., Neuroimage 2010; https://doi.org/10.1016/j.neuroimage.2010.06.010
        myParcVol='aparc.a2009s+aseg';
        myParcSurf='h.aparc.a2009s.annot';
        %onlyCtx=[17,18,53,54,11101:11175,12101:12175]'; % only hippocampus, amygdala and cortical labels from Destrieux parcellation
    otherwise
        error('I do not recognize this parcellation: %s.',parcType);
end

switch cfg.dataType
    
    case 'iELVis'
        
        fileExt='.mgz';
        
        % load electrode names
        elecNames=readiELVisElecNames(datadir,patID);
        elecType=elecNames(:,2);
        elecHem=elecNames(:,3);
        elecNames=elecNames(:,1);
        
        % load volume electrode coordinates (LIP)
        elecCoordVol=readiELVisElecCoord(datadir,patID,'LEPTOVOX');
        
        % load surface electrode coordinates (RAS?)
        elecCoordSurf=readiELVisElecCoord(datadir,patID,'LEPTO');
        
        % load volume parcellations (ILA)
        wmparc=MRIread(fullfile(filePath,'mri',['wmparc' fileExt])); % load white matter parcellation
        parcVol=MRIread(fullfile(filePath,'mri',[myParcVol fileExt])); % load volume parcellation
        parcVol.vol(parcVol.vol==2|parcVol.vol==41)=wmparc.vol(parcVol.vol==2|parcVol.vol==41); % replace white matter labels with DK white matter parcellation
        
        % load surface parcellations
        [vertCoordL,facesCoordL]=read_surf(fullfile(filePath,'surf','lh.pial'));
        [vertCoordR,facesCoordR]=read_surf(fullfile(filePath,'surf','rh.pial'));
        [avertsL,alblL,actblL]=read_annotation(fullfile(filePath,'label',['l' myParcSurf]));
        avertsL=avertsL+1; % FreeSurfer indexing is 0-based
        [avertsR,alblR,actblR]=read_annotation(fullfile(filePath,'label',['r' myParcSurf]));
        avertsR=avertsR+1; % FreeSurfer indexing is 0-based
        
    case 'BIDS'
        
        fileExt='.nii';
        
        % load electrode table
        elecTable=readtable(fullfile(filePath,'ieeg',[patID '_electrodes.tsv']),'FileType','text','Delimiter','\t','ReadVariableNames',true,'ReadRowNames',false);
        elecNames=elecTable.name;
        elecType=elecTable.type;
        elecHem=elecTable.hemisphere;
        
        % voxeloc coords are LIA
        elecCoordVol=table2array(elecTable(:,2:4));
        elecCoordVol(:,3)=256-elecCoordVol(:,3); % make them LIP like iELVis
        
        % load volume parcellations (ILA)
        wmparc=MRIread(fullfile(filePath,'anat',[patID '_wmparc' fileExt])); % load white matter parcellation
        parcVol=MRIread(fullfile(filePath,'anat',[patID '_' myParcVol fileExt])); % load volume parcellation
        parcVol.vol(parcVol.vol==2|parcVol.vol==41)=wmparc.vol(parcVol.vol==2|parcVol.vol==41); % replace white matter labels with DK white matter parcellation
end

nElec=size(elecCoordVol,1);

% match elecCoordVol dimension order and direction with scans
LIP2ILA=[0 1 0 0;1 0 0 0;0 0 -1 256;0 0 0 1];
elecCoordVol=(LIP2ILA*[elecCoordVol';ones(1,nElec)])';
elecCoordVol=elecCoordVol(:,1:3);

% load FreeSurfer label lookup table
% provided by FreeSurfer here: https://surfer.nmr.mgh.harvard.edu/fswiki/BIDS
% file bundled with SEEG2parc, must be in the path
fsLut=readtable('labels.tsv','FileType','text','ReadVariableNames',true,'Delimiter','\t');

% define weights for neighboring voxels
voxelWeights=1/cat(3, ...
    [sqrt(3) sqrt(2) sqrt(3); sqrt(2) 1 sqrt(2); sqrt(3) sqrt(2) sqrt(3)], ...
    [sqrt(2) 1 sqrt(2); 1 1 1; sqrt(2) 1 sqrt(2)], ...
    [sqrt(3) sqrt(2) sqrt(3); sqrt(2) 1 sqrt(2); sqrt(3) sqrt(2) sqrt(3)]);

% prepare containers for tissue labels
tissueLabels=cell(nElec,1);
tissueWeights=cell(nElec,1);

% main loop over electrodes
for ctElec=1:nElec
    if strncmpi(elecType{ctElec},'D',1) % double-check that this is a depth electrode
        
        voxelCubeCoord=[ ...
            round(elecCoordVol(ctElec,1))-1, round(elecCoordVol(ctElec,2))-1, round(elecCoordVol(ctElec,3))-1; ...
            round(elecCoordVol(ctElec,1))+1, round(elecCoordVol(ctElec,2))+1, round(elecCoordVol(ctElec,3))+1];
        
        % REMmed by Pierre 20220721
        % deal with case where voxel cube would reach beyond the limits of the MRI (unlikely in practice)
        %{
        voxelCubeCoord(voxelCubeCoord<1)=1;
        voxelCubeCoord(voxelCubeCoord>256)=256;
        %}
        
        voxelCube=parcVol.vol(voxelCubeCoord(1,1):voxelCubeCoord(2,1),voxelCubeCoord(1,2):voxelCubeCoord(2,2),voxelCubeCoord(1,3):voxelCubeCoord(2,3));
        
        % REMmed by Pierre 20220721
        % deal with case where voxel cube would reach to the other cerebral hemisphere
        %{
        switch elecNames{ctElec,3}
            case 'L' % left hemisphere
                voxelCube(~ismember(voxelCube,myLOI_L))=0;
            case 'R' % right hemisphere
                voxelCube(~ismember(voxelCube,myLOI_R))=0;
        end
        %}
        
        % REMmed by Pierre 20240227
        % if only cortical tissue is required by user
        %{
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
        %}
        
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
            voxelCubeLabels{ctLbl}=fsLut.name{voxelCubeUniqueSorted(ctLbl)==fsLut.index};
        end
        
        tissueLabels{ctElec}=voxelCubeLabels;
        tissueWeights{ctElec}=voxelCubeUniqueWeightSorted;
        
    elseif any(strncmpi(elecType{ctElec},{'S','G'},1)) % double-check that this is a subdural electrode
        
        % determine whether working on L or R hemisphere
        switch lower(elecHem{ctElec})
            case 'l'
                vertCoord=vertCoordL;
                facesCoord=facesCoordL;
                actbl=actblL;
                averts=avertsL;
                albl=alblL;
            case 'r'
                vertCoord=vertCoordR;
                facesCoord=facesCoordR;
                actbl=actblR;
                averts=avertsR;
                albl=alblR;
        end
        
        % compute distance between electrode and all vertices
        vertDist=pdist2(vertCoord,elecCoordSurf(ctElec,:),'euclidean');
        
        % sort these distances in ascending order and keep track of the indices to the corresponding vertices
        [vertDistSort,vertDistIdx]=sort(vertDist);
        
        % sort labels and retain only those within 4 mm of the one closest to the electrode
        alblSort=albl(vertDistIdx);
        alblSortSelect=alblSort(1:find(vertDistSort<vertDistSort(1)+4,1,'last'));
        
        % determine label weights (2 cutoffs at 2 and 3 mm)
        vertDistSortSelect=vertDistSort(1:find(vertDistSort<vertDistSort(1)+4,1,'last'));
        pialSurfWeights=ones(size(vertDistSortSelect));
        pialSurfWeights(vertDistSortSelect>vertDistSortSelect(1)+2)=1/sqrt(2);
        pialSurfWeights(vertDistSortSelect>vertDistSortSelect(1)+3)=1/sqrt(3);
        
        % determine unique labels and weights
        [alblUniq,~,alblUniqIdx]=unique(alblSortSelect,'stable');
        nLabels=numel(alblUniq);
        pialSurfLabels=cell(nLabels,1);
        pialSurfWeightsUniq=zeros(size(pialSurfLabels));
        for ctUniq=1:nLabels
            pialSurfLabels(ctUniq)=actbl.struct_names(alblUniq(ctUniq)==actbl.table(:,5));
            pialSurfWeightsUniq(ctUniq)=sum(pialSurfWeights(alblUniqIdx==ctUniq));
        end
        
        % express weights in fractions of 1
        pialSurfWeightsUniq=pialSurfWeightsUniq./sum(pialSurfWeightsUniq(:));
        
        tissueLabels{ctElec}=pialSurfLabels;
        tissueWeights{ctElec}=pialSurfWeightsUniq;
                
    else
        error('I did not recognize this type of electrode: %s, %s.',elecNames{ctElec},elecType{ctElec});
    end
    
end

% prepare cell array with JSON notation for tissue labels and weights
tissueLabelsWeights=cell(nElec,2);
for ctElec=1:nElec
    tmpLabel='';
    tmpWeight='';
    for ctLabel=1:numel(tissueLabels{ctElec})
        %tmpLabel=[tmpLabel tissueLabels{ctElec}{ctLabel} '=' num2str(tissueWeights{ctElec}(ctLabel),'%.3f') ';'];
        tmpLabel=[tmpLabel '"' tissueLabels{ctElec}{ctLabel} '",'];
        tmpWeight=[tmpWeight num2str(tissueWeights{ctElec}(ctLabel),'%.3f') ','];
    end
    tissueLabelsWeights{ctElec,1}=['{' tmpLabel(1:end-1) '}'];
    tissueLabelsWeights{ctElec,2}=['[' tmpWeight(1:end-1) ']'];
end

% prepare elecTable
if exist('elecTable','var')==1
    
    % update electrode coordinates
    elecTable.x=elecCoordVol(:,1);
    elecTable.y=elecCoordVol(:,2);
    elecTable.z=elecCoordVol(:,3);
    
    % add tissue labels/weights variable
    elecTable.tissueLabels=tissueLabelsWeights(:,1);
    elecTable.tissueWeights=tissueLabelsWeights(:,2);
    
else
    
    % create table from scratch
    elecTable=table(elecNames,elecCoordVol(:,1),elecCoordVol(:,2),elecCoordVol(:,3),NaN(nElec,1),elecHem,elecType,tissueLabelsWeights(:,1),tissueLabelsWeights(:,2), ...
        'VariableNames',{'name','x','y','z','size','hemisphere','type','tissueLabels','tissueWeights'});
    
end

% to convert tissueWeights back into a numeric cell array:
% tissueWeights=cellfun(@str2num,elecTable.tissueWeights,'UniformOutput',false)

% to write elecTable as a .tsv file:
% writetable(elecTable,'electrodes.tsv','FileType','text','Delimiter','\t');

%% plotting volume (useful for debug)

%{
switch cfg.dataType
    case 'iELVis'
        CT=MRIread(fullfile(filePath,'elec_recon','postInPre.nii'));  % load post-implant CT
        if isempty(CT)
            CT=MRIread(fullfile(filePath,'elec_recon','ctINt1.nii'));  % load post-implant CT
        end
    case 'BIDS'
        CT=MRIread(fullfile(filePath,'anat',[patID '_postInPre.nii'])); % load post-implant CT
end
elecCoordRound=round(elecCoord);
figure;
for ctElec=1:nElec
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
