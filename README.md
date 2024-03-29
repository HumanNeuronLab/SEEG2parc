# SEEG2parc
Attribute cerebral parcellation labels to stereo-EEG electrodes

This function uses stereo-EEG electrode coordinates, obtained from the iELVis (1) pipeline, and the cerebral parcellations generated by FreeSurfer (2), to attribute tissue labels to each electrode.

There are 2 main options:
- parcellation type: can use either the Desikan-Killiany (3) or Destrieux (4) atlas
- tissue type: can be restricted to cortex (+ hippocampus and amygdala) or unrestricted, taking into account all FreeSurfer parcellation labels

Requisites: iELVis, https://github.com/iELVis/iELVis/tree/iELVis_pm.

The tissue labels are taken from a 3x3x3-voxel cube surrounding the electrode coordinate and are weighted like so:

<pre>         +--------+--------+--------+
        /    _   /    _   /    _   /|
       /   \/3  /   \/2  /   \/3  / |
      +--------+--------+--------+  |
     /    _   /        /    _   /|  |
    /   \/2  /    1   /   \/2  / |  +
   +--------+--------+--------+  | /
  /    _   /    _   /    _   /|  |/
 /   \/3  /   \/2  /   \/3  / |  +
+--------+--------+--------+  | /
|        |        |        |  |/
|        |        |        |  +
|        |        |        | /
|        |        |        |/
+--------+--------+--------+

         +--------+--------+--------+
        /    _   /        /    _   /|
       /   \/2  /    1   /   \/2  / |
      +--------+--------+--------+  |
     /        /        /        /|  |
    /    1   /    1   /    1   / |  +
   +--------+--------+--------+  | /
  /    _   /        /    _   /|  |/
 /   \/2  /    1   /   \/2  / |  +
+--------+--------+--------+  | /
|        |        |        |  |/
|        |        |        |  +
|        |        |        | /
|        |        |        |/
+--------+--------+--------+

         +--------+--------+--------+
        /    _   /    _   /    _   /|
       /   \/3  /   \/2  /   \/3  / |
      +--------+--------+--------+  |
     /    _   /        /    _   /|  |
    /   \/2  /    1   /   \/2  / |  +
   +--------+--------+--------+  | /
  /    _   /    _   /    _   /|  |/
 /   \/3  /   \/2  /   \/3  / |  +
+--------+--------+--------+  | /
|        |        |        |  |/
|        |        |        |  +
|        |        |        | /
|        |        |        |/
+--------+--------+--------+</pre>

The tissue labels are returned both in "deterministic" (whichever tissue label is the most common in the weighted cube) and "probabilistic" fashion (all tissue labels, as well as their weights, are returned).


This function was inspired by Mercier et al.'s Proximal Tissue Density (PTD) index (5), which is implemented in iELVis.

Pierre Mégevand, Human Neuron Lab, University of Geneva, Switzerland. 2022.
pierre.megevand@unige.ch; https://www.unige.ch/medecine/neucli/en/groupes-de-recherche/1034megevand/


#### References:
1. Groppe DM, Bickel S, Dykstra AR, Wang X, Mégevand P, Mercier MR, Lado FA, Mehta AD, Honey CJ. iELVis: An open source MATLAB toolbox for localizing and visualizing human intracranial electrode data. Journal of Neuroscience Methods. 2017. https://doi.org/10.1016/j.jneumeth.2017.01.022
2. Fischl B. FreeSurfer. NeuroImage. 2012. https://doi.org/10.1016/j.neuroimage.2012.01.021
3. Desikan RS, Ségonne F, Fischl B, Quinn BT, Dickerson BC, Blacker D, Buckner RL, Dale AM, Maguire RP, Hyman BT, Albert MS, Killiany RJ. An automated labeling system for subdividing the human cerebral cortex on MRI scans into gyral based regions of interest. NeuroImage. 2006. https://doi.org/10.1016/j.neuroimage.2006.01.021
4. Destrieux C, Fischl B, Dale A, Halgren E. Automatic parcellation of human cortical gyri and sulci using standard anatomical nomenclature. NeuroImage. 2010. https://doi.org/10.1016/j.neuroimage.2010.06.010
5. Mercier MR, Bickel S, Megevand P, Groppe DM, Schroeder CE, Mehta AD, Lado FA. Evaluation of cortical local field potential diffusion in stereotactic electro-encephalography recordings: A glimpse on white matter signal. NeuroImage. 2017. https://doi.org/10.1016/j.neuroimage.2016.08.037
