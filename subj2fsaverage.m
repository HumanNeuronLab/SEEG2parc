function elecCoordAvg=subj2fsaverage(datadir,subj)

% function elecCoordAvg=subj2fsaverage(datadir,subj)
% Get iEEG electrode coordinates in the FreeSurfer average brain (MRI305).
%
% This is an adaptation from the sub2AvgBrain function of iELVis, by David
% Groppe and myself. https://github.com/iELVis/iELVis/tree/iELVis_pm
%
% Inputs:
%   datadir         full path to patient's iEEG anatomical information
%   subj            patient code name
%
% Outputs:
%   elecCoordAvg    N-by-3 matrix with electrode coordinates in the MNI305
%                   (fsaverage) space
%
% Requirements:
%   FreeSurfer (for the MATLAB functions): https://surfer.nmr.mgh.harvard.edu/
%
% https://github.com/HumanNeuronLab/SEEG2parc
% Pierre Mégevand, Human Neuron Lab, University of Geneva, Switzerland. 2026.
% pierre.megevand@unige.ch; https://www.unige.ch/medecine/neucli/en/groupes-de-recherche/1034megevand/


% get electrode names and coordinates
elecNames=readiELVisElecNames(subj);
elecCoord=readiELVisElecCoord(subj,'LEPTO');
isDepth=strcmp(elecNames(:,2),'D');
isLeft=strcmp(elecNames(:,3),'L');

elecCoordAvg=zeros(size(elecCoord));


%% depth electrodes

if any(isDepth)
    elecCoordDepth=elecCoord(isDepth,:);
    
    % get MRI transform matrices
    MRIhdr=MRIread(fullfile(datadir,subj,'mri','brainmask.mgz'),true);
    Norig=MRIhdr.vox2ras;
    Torig=MRIhdr.tkrvox2ras;
    TalTransform=freesurfer_read_talxfm(fullfile(datadir,subj,'mri','transforms','talairach.xfm'));
    
    % Pierre took this code from http://surfer.nmr.mgh.harvard.edu/fswiki/CoordinateSystems
    % Section #2 of "Transforms within a subject's anatomical space"
    % "2. I have an RAS point on the surface (tkrR tkrA tkrS) ("Vertex RAS"
    % from tksurfer) and want to compute the MNI305 RAS that corresponds to this point:"
    elecCoordDepthAvg=(TalTransform*Norig*(Torig\[elecCoordDepth'; ones(1, sum(isDepth))]))';
    elecCoordDepthAvg=elecCoordDepthAvg(:,1:3);

    elecCoordAvg(isDepth,:)=elecCoordDepthAvg;
end


%% subdural electrodes

if any(~isDepth)
    elecCoordSubdural=elecCoord(~isDepth,:);
    elecCoordSubduralAvg=zeros(size(elecCoordSubdural));
    vertSubj=zeros(sum(~isDepth),1);
    vertAvg=zeros(sum(~isDepth),1);
    isLeftSubdural=isLeft(~isDepth);
    
    % load pial surfaces
    [pialSubjL.vert,pialSubjL.tri]=read_surf(fullfile(datadir,subj,'surf','lh.pial'));
    [pialSubjR.vert,pialSubjR.tri]=read_surf(fullfile(datadir,subj,'surf','rh.pial'));
    [sphereSubjL.vert,sphereSubjL.tri]=read_surf(fullfile(datadir,subj,'surf','lh.sphere.reg'));
    [sphereSubjR.vert,sphereSubjR.tri]=read_surf(fullfile(datadir,subj,'surf','rh.sphere.reg'));
    [pialAvgL.vert,pialAvgL.tri]=read_surf(fullfile(datadir,'fsaverage','surf','lh.pial'));
    [pialAvgR.vert,pialAvgR.tri]=read_surf(fullfile(datadir,'fsaverage','surf','rh.pial'));
    [sphereAvgL.vert,sphereAvgL.tri]=read_surf(fullfile(datadir,'fsaverage','surf','lh.sphere.reg'));
    [sphereAvgR.vert,sphereAvgR.tri]=read_surf(fullfile(datadir,'fsaverage','surf','rh.sphere.reg'));
    
    euclDist=@(x1,x2) sqrt(sum((x1-x2).^2,2));
    
    for ctSubdural=1:sum(~isDepth)
        
        % find nearest vertex on subject's pial surface
        if isLeftSubdural(ctSubdural)
            dst=euclDist(pialSubjL.vert,elecCoordSubdural(ctSubdural,:));
        else
            dst=euclDist(pialSubjR.vert,elecCoordSubdural(ctSubdural,:));
        end
        [~,vertSubj(ctSubdural)]=min(dst);
        
        % find nearest vertex between subject's and fsaverage's sphere surfaces
        if isLeftSubdural(ctSubdural)
            dst=euclDist(sphereAvgL.vert,sphereSubjL.vert(vertSubj(ctSubdural),:));
        else
            dst=euclDist(sphereAvgR.vert,sphereSubjR.vert(vertSubj(ctSubdural),:));
        end
        [~,vertAvg(ctSubdural)]=min(dst);
        
        % get coordinates of fsaverage's pial surface vertex
        if isLeftSubdural(ctSubdural)
            elecCoordSubduralAvg(ctSubdural,:)=pialAvgL.vert(vertAvg(ctSubdural),:);
        else
            elecCoordSubduralAvg(ctSubdural,:)=pialAvgR.vert(vertAvg(ctSubdural),:);
        end
        
    end
    
    elecCoordAvg(~isDepth,:)=elecCoordSubduralAvg;
    
end
