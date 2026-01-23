% This script plots the cortica lsurface of both hemispheres, together with
% iEEG electrodes.
% Requirements:
%   - iELVis: https://github.com/iELVis/iELVis/tree/iELVis_pm
%   - SEEG2parc: 


% provide access path to "FreeSurfer" folder with patient's anatomy
global globalFsDir;
globalFsDir='your\path';

% select patient to plot
fsSubj='PAT_6619';
[str,rem]=strtok(fsSubj,'_');

% set basic figure plotting parameters
iEEG_fig=figure;
cfgPlot=[];
cfgPlot.title=[str ' ' rem(2:end)];
cfgPlot.backgroundColor=[1 1 1];

% set cortical surface plotting parameters
cfgPlot.opaqueness=0.3; % set < 1 to see through brain surface and show depth electrodes
% cfgPlot.overlayParcellation='DK';

% set electrode plotting parameters
cfgPlot.ignoreDepthElec='n'; % show depth electrodes
cfgPlot.showLabels='y';
cfgPlot.elecShape='marker';
cfgPlot.pullOut=0;

% get electrode names and coordinates
cfgPlot.elecNames=readiELVisElecNames(fsSubj);
cfgPlot.elecCoord=readiELVisElecCoord(fsSubj,'LEPTO');
cfgPlot.elecCoord(:,4)=strcmp('L',cfgPlot.elecNames(:,3));
cfgPlot.elecNames=cfgPlot.elecNames(:,1);

% do the plotting
cfgPlot.figId=iEEG_fig;
cfgPlot.clearFig='n';
cfgPlot.view='lf';
cfgOut=plotPialSurf(fsSubj,cfgPlot);
cfgPlot.view='rf';
cfgOut=plotPialSurf(fsSubj,cfgPlot);
