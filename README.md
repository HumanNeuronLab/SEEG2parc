# SEEG2parc

#### Attribute cerebral parcellation labels to intracranial EEG electrodes

In its current version, the function is compatible both with iELVis and with the BIDS-iEEG output of Voxeloc. However, the attribution of parcellation labels to subdural electrodes (ECOG) only works with the iELVis data type.

Requisites:
- iELVis, https://github.com/iELVis/iELVis/tree/iELVis_pm
- Voxeloc, https://github.com/HumanNeuronLab/voxeloc
- FreeSurfer, https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferWiki

This function was inspired by Mercier et al.'s Proximal Tissue Density (PTD) index, which is implemented in iELVis.

Pierre Mégevand, Human Neuron Lab, University of Geneva, Switzerland. 2022-2025.  
pierre.megevand@unige.ch  
https://www.unige.ch/medecine/neucli/en/groupes-de-recherche/1034megevand/


#### References:
- Groppe DM, Bickel S, Dykstra AR, Wang X, Mégevand P, Mercier MR, Lado FA, Mehta AD, Honey CJ. iELVis: An open source MATLAB toolbox for localizing and visualizing human intracranial electrode data. Journal of Neuroscience Methods. 2017. https://doi.org/10.1016/j.jneumeth.2017.01.022
- Monney J, Dallaire SE, Stoutah L, Fanda L, Mégevand P. Voxeloc: A time-saving graphical user interface for localizing and visualizing stereo-EEG electrodes. Journal of Neuroscience Methods. 2024. https://doi.org/10.1016/j.jneumeth.2024.110154
- Fischl B. FreeSurfer. NeuroImage. 2012. https://doi.org/10.1016/j.neuroimage.2012.01.021
- Desikan RS, Ségonne F, Fischl B, Quinn BT, Dickerson BC, Blacker D, Buckner RL, Dale AM, Maguire RP, Hyman BT, Albert MS, Killiany RJ. An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest. NeuroImage. 2006. https://doi.org/10.1016/j.neuroimage.2006.01.021
- Destrieux C, Fischl B, Dale A, Halgren E. Automatic parcellation of human cortical gyri and sulci using standard anatomical nomenclature. NeuroImage. 2010. https://doi.org/10.1016/j.neuroimage.2010.06.010
- Mercier MR, Bickel S, Megevand P, Groppe DM, Schroeder CE, Mehta AD, Lado FA. Evaluation of cortical local field potential diffusion in stereotactic electro-encephalography recordings: A glimpse on white matter signal. NeuroImage. 2017. https://doi.org/10.1016/j.neuroimage.2016.08.037
