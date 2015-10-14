function res = SpectralParcellation(subjectlist,hemisphere,weight,Nclus,pathToData,wbCommand)

% % Groupwise parcellation using multiscale spectral clustering


%%%%%% INPUTS %%%%%%
% % subjectlist: cell array of strings listing the ID of the subjects to be parcellated
% % hemisphere: 'L' or 'R', whether parcellating the left or right hemisphere
% % weight: multiplicative weight for the inter-subjects links
% % pathToData: base path where the data is stored
% % wbCommand: workbench command, typically
% 'path/to/workbench/installation/wb_command'
% % Nclus: Number of desired parcels

%%%%%% OUTPUTS %%%%%%
% % res: structure containing the obtained parcellation and merged
% connectivity matrices for all subjects
% res structure fields:
% -parcellation: number of vertices x number of subjects, describes for each subject
% the parcel assignment of each vertex
% -Mat: number of parcels x number of parcels x number of subjects, describes the low
% dimensionality connectivity networl obtained after parcellation for each subject.
% -Inlist: cell array, number of subjects x 1, each cell Inlist{i}{j} lists which
% vertices belong to parcel j for subject i. Another way of describing the parcel assignments.

% Copyright (C) Sarah Parisot, Imperial College London, 2015
%
%This program is free software: you can redistribute it and/or modify
%it under the terms of the GNU General Public License as published by
%the Free Software Foundation, either version 3 of the License, or
%(at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program.  If not, see <http://www.gnu.org/licenses/>.


%% Preprocessing: compute the correlation matrix and load structural data

% path to tractography data, TO MODFIFY based on data organisation
pathToCorMat = '/diffusion/preprocessed/T1w/probtrack/';

% path where temporary files are saved, has to be subject specific
for i=1:numel(subjectlist)
    savePath(i) = [pathToData subjectlist{i} pathToCorMat hemisphere];
end

for i=1:numel(subjectlist)

    subjectID = subjectlist{i};

    %%% load useful structural data
    data = loadSubjectData(subjectID,hemisphere,pathToData);

    if ~exist([savePath(i) '/SCM.mat'])
    %%% Get the correlation matrix
    CorMat = ComputeCorrelationMatrix(subjectID,pathToData,pathToCorMat,wbCommand);

    %%% Create supervertices
    SCmatrix = MultiAffinity(data.Sphere,data.BrainSurf,data.roi,CorMat,data.Mat,1,data.LV);

    %%% saving results (better for memory than loading all data and SCmatrices at once)
    cd(savePath(i))
    save('SCM.mat','SCmatrix','data')
    clear SCmatrix
    clear CorMat
    end

    %% Construction of inter-subjects links

    R=cell(numel(subjectlist),numel(subjectlist));
    match=cell(numel(subjectlist),numel(subjectlist));

    for j=1:numel(subjectlist)

        Subject2=subjectlist{j};
        disp(Subject2)

        if strcmp(subjectID,Subject2)
            continue
        end

        cd(savePath(i))

        load SCM.mat

        SCmatrix1=SCmatrix;

        cd(savePath(j))

        load SCM.mat

        SCmatrix2=SCmatrix;

        clear SCMatrix
        data.Mat=[];


        [R{i,j},match{i,j}] = ComputePWlinks(SCmatrix1,SCmatrix2,data);
    end

end

%% Do the groupwise parcellation


%%% load all subject data and set up some parameters
for i=1:numel(subjectlist)
    cd(savePath(i))
    load('SCM.mat','SCmatrix');
    SCmatrices{i}=SCmatrix;
    clear SCmatrix
end

Nlayers=length(SCmatrix.baseSeg.Seeds);
layer =1; %% supervertex resolution on which the parcellation is visualised


for i=1:Nlayers
    nC(i)=numel(SCmatrix.baseSeg.Seeds{i});
end

tmp=[SCmatrices{:}];
C={tmp(:).C};
tmp = [tmp(:).baseSeg];

%%% Parcellation

[Wtot,C,d]=BuildMatrices(C,[tmp(:).Cor],R,nC,weight);
[NcutDiscrete,NcutEigenvectors,Eigenvalues]=Ncut_multi(Wtot,C,d,Nclus);


clear tmp

%%% Project parcellation on the brain surface

for j=1:numel(subjectlist)
    cd(savePath(j))
    load('SCM.mat','data');

    sN=j-1;

    if sN==0
        Map=zeros(length(data.BrainSurf.vertices),numel(subjectlist));
    end

    %%%  offset to find the subject's parcellation in the global parcellation
    %%% matrix
    offset=sN*sum(nC(1:Nlayers))+sum(nC(1:layer-1));
    for i=offset+1:nC(layer)+offset
        Map(SCmatrices{sN+1}.baseSeg.Maps{layer}==i-offset,j)=find(NcutDiscrete(i,:)==1);
    end

    %%% Compute the merged connectivity profiles and normalise them
    [Corr(:,:,j),Inlist{j}] = MergeConnectivityMatrix2(1:Nclus,Map(data.Inc,j),data.Mat(data.Inc,data.Inc),data.Inc);

    for i=1:Nclus
        if numel(Inlist{j}{i})==0
            dg(i)=1;
        else
            dg(i)=numel(Inlist{j}{i});
        end
    end

    Corr(:,:,j)=bsxfun(@rdivide,Corr(:,:,j),dg');

end

res.parcellation=Map;
res.Mat=Corr;
res.Inlist=Inlist;

end
