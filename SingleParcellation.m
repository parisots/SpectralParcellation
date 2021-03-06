function res = SingleParcellation(subjectID,hemisphere,Nclus,pathToData,wbCommand)

% % Single subject parcellation using multiscale spectral clustering
% % Assumes an organisation of the data similar to the Human Connectome
% % Project data

%%%%%% INPUTS %%%%%%
% % subjectID: identification of the subjects to be parcellated
% % hemisphere: 'L' or 'R', whether parcellating the left or right hemisphere
% % pathToData: base path where the data is stored
% % wbCommand: workbench command, typically
% 'path/to/workbench/installation/wb_command'
% % Nclus: number of desired parcels

%%%%%% OUTPUTS %%%%%%
% % res: structure containing the obtained parcellation and merged
% connectivity matrix
% res structure fields: 
% -parcellation: number of vertices x 1, describes the parcel assignment of each vertex
% -Mat: number of parcels x number of parcels, describes the low dimensionality 
% 	  connectivity networl obtained after parcellation. 
% -Inlist: cell array number of parcels x 1, each cell Inlist{i} lists which 
% vertices belong to parcel i. Another way of describing the parcel assignments.

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

%%% path where temporary files are saved, has to be subject specific
savePath= [pathToData subjectID pathToCorMat hemisphere];

%%% load useful structural data
data = loadSubjectData(subjectID,hemisphere,pathToData);

%%% Get the correlation matrix
if ~exist([savePath '/SCM.mat']) %% if matrices are saved, only compute them once
CorMat = ComputeCorrelationMatrix(subjectID,pathToData,pathToCorMat,wbCommand);

%%% Create supervertices

    SCmatrix = MultiAffinity(data.Sphere,data.BrainSurf,data.roi,CorMat,data.Mat,1,data.LV);


%%% OPTIONAL: saving results
cd(savePath)
save('SCM.mat','SCmatrix','data')
clear CorMat

end

%% Do the parcellation

% Parameters
Nlayers=length(SCmatrix.baseSeg.Seeds);
layer =1; %% supervertex resolution on which the parcellation is visualised

for i=1:Nlayers
    nC(i)=numel(SCmatrix.baseSeg.Seeds{i});
end

[NcutDiscrete,NcutEigenvectors,Eigenvalues]=Ncut_multi(...
    SCmatrix.W(1:sum(nC),1:sum(nC)),SCmatrix.C(1:sum(nC(2:end)),1:sum(nC)),SCmatrix.d(1:sum(nC)),Nclus);

%%% Project parcellation on the brain surface%%% load useful structural data

cd([pathToData subjectID pathToCorMat hemisphere])

Map=zeros(length(data.BrainSurf.vertices),1);


%%%  Offset to find the subject's parcellation in the multiscale parcellation
%%% matrix
offset=sum(nC(1:layer-1));
for i=offset+1:nC(layer)+offset
    Map(SCmatrix.baseSeg.Maps{layer}==i-offset)=find(NcutDiscrete(i,:)==1);
end

%%% Compute the merged connectivity profiles and normalise them
[M,Inlist] = MergeConnectivityMatrix2(1:Nclus,Map(data.Inc),data.Mat(data.Inc,data.Inc),data.Inc);

for i=1:Nclus
    if numel(Inlist{i})==0
        dg(i)=1;
    else
        dg(i)=numel(Inlist{i});
    end
end

M=bsxfun(@rdivide,M,dg');


res.parcellation=Map;
res.Mat=M;
res.Inlist=Inlist;

end



