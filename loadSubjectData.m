function data = loadSubjectData(subject,hemisphere,pathToData)

% % load useful structural data for further parcellation tasks
% % assumes HCP style organisation of the data

%%%%%% INPUTS %%%%%%
% % subject: subject ID, string
% % hemisphere: 'L' or 'R', left or right hemisphere
% % surface: string, options for HCP data are: midthickness, pial, white
% % pathToData: base path where the data is stored

%%%%%% OUTPUT %%%%%%
% % data: structure regrouping all loaded information

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


%%


surface = 'midthickness';
folder = 'MNINonLinear';

%%% Identify medial wall vertices
cd([pathToData subject '/structural/' folder '/fsaverage_LR32k'])
roi=gifti([ subject '.' hemisphere '.atlasroi.32k_fs_LR.shape.gii']);
Nvertices=size(roi.cdata,1);
MWlist = find(roi.cdata==0); % list of medial wall vertices
data.Inc=setdiff(1:Nvertices,MWlist)'; 
data.roi=roi;
% roi list the vertices that are excluded the mesh (medial wall vertices),
% it should be all zeros if no vertices are excluded

%%% load surface data
cd([pathToData subject '/structural/' folder '/fsaverage_LR32k'])
%actual cortical surface
data.BrainSurf=gifti([subject '.' hemisphere '.' surface '.32k_fs_LR.surf.gii']);
%spherical inflation of the surface (necessary for supervertex parcellation)
data.Sphere=gifti([subject '.' hemisphere '.sphere.32k_fs_LR.surf.gii']);
%inflated mesh (for visualisation purposes)
data.Inflated=gifti([subject '.' hemisphere '.inflated.32k_fs_LR.surf.gii']);

%%% describes the neighboring vertices of each surface vertex. 
list_vertices=cell(Nvertices,1);
for S=1:Nvertices
    [idx,~]=find(data.BrainSurf.faces==S);
    list_vertices{S}=setdiff(unique(data.BrainSurf.faces(idx,:)),S)';
end

data.LV = list_vertices;

%%% connectivity data, output of probtrackX
cd([pathToData subject '/diffusion/preprocessed/T1w/probtrack/' hemisphere])
data.Mat = spconvert(load('fdt_matrix1.dot'));

