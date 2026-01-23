% a quick script to help convert SEEG2parc output and build iEEG-BIDS compatible electrodes.tsv
%
% https://github.com/HumanNeuronLab/SEEG2parc
% Pierre Mégevand, Human Neuron Lab, University of Geneva, Switzerland. 2026.
% pierre.megevand@unige.ch; https://www.unige.ch/medecine/neucli/en/groupes-de-recherche/1034megevand/


% first, run SEEG2parc while collecting optional outputs:
% [elecTable,tissueLabels,tissueWeights]=SEEG2parc(cfg);


% for tissueWeights
N=numel(tissueWeights);
n=zeros(N,1);
for ct1=1:N
    n(ct1)=numel(tissueWeights{ct1});
end

w=zeros(N,max(n));
for ct1=1:N
    for ct2=1:n(ct1)
        w(ct1,ct2)=tissueWeights{ct1}(ct2);
    end
end

% w can then be pasted into Excel and cells containing only zeros erased


% for tissueLabels
l=strings(N,1);
for ct1=1:N
    if n(ct1)==1
        l(ct1)=tissueLabels{ct1}{1};
    else
        tmp=tissueLabels{ct1}{1};
        for ct2=2:n(ct1)
            tmp=[tmp ' ' tissueLabels{ct1}{ct2}];
        end
        l(ct1)=tmp;
    end
end

% l can then be pasted into Excel as is



            