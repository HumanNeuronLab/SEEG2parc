function [coord_LEPTOVOX,coord_LEPTOSURF]=voxeloc_iELVis_coordConvert(coord);

% function [coord_LEPTOVOX,coord_LEPTOSURF]=voxeloc_iELVis_coordConvert(coord)
% Converts stereo-EEG electrode coordinates from the Voxeloc coordinate
% system to the iELVis ones.

% coordinate systems:
% 
% Voxeloc and elecTable are VOX and
% x: dorsal to ventral (superior to inferior)
% y: right to left
% z: caudal to rostral (posterior to anterior)
% 
% .LEPTOVOX is VOX and
% x: right to left
% y: dorsal to ventral (superior to inferior)
% z: rostral to caudal (anterior to posterior)
% 
% .LEPTO is SURFACE and
% x: left to right
% y: caudal to rostral (posterior to anterior)
% z: ventral to dorsal (inferior to superior)

% https://github.com/HumanNeuronLab/xxxxxxx
% Pierre Megevand, Human Neuron Lab, University of Geneva, Switzerland. 2026.
% pierre.megevand@unige.ch; https://www.unige.ch/medecine/neucli/en/groupes-de-recherche/1034megevand/


% convert from Voxeloc to iELVis .LEPTOVOX
VOXELOC2LEPTOVOX=[0 1 0 0; 1 0 0 0; 0 0 -1 256; 0 0 0 1];
coord_LEPTOVOX=(VOXELOC2LEPTOVOX*[coord';ones(1,size(coord,1))])';
coord_LEPTOVOX(:,4)=[];

% convert from Voxeloc to iELVis .LEPTO
VOXELOC2LEPTOSURF=[0 -1 0 128; 0 0 1 -128; -1 0 0 128; 0 0 0 1];
coord_LEPTOSURF=(VOXELOC2LEPTOSURF*[coord';ones(1,size(coord,1))])';
coord_LEPTOSURF(:,4)=[];
