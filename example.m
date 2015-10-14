
% Example usage of the SingleParcellation and SpectralParcellation functions

% Copyright (C) Sarah Parisot, Imperial College London, 2015

%% Set the parameters

% list all subjects identifications here. 
subjectlist={'100307','103414','105115'};

% depends on the number of parcels and number of subjects in the group
% example weight values for groups of 50 subjects are given in the 
% 'param.mat' matrix, with the corresponding number of parcels (Nclus).
% weights too strong create isolated supervertices, too weak prevent the
% coherence between subjects. 
weight = 0.5;

% number of desired parcels
Nclus = 100; 

% Left or right hemisphere to parcellate. Could be whole brain if data
% dependent function (loadSubjectData.m and ComputeCorrelationMatrix.m are 
% adapted.)
hemisphere = 'L';

% Setting the paths, TO BE MODIFIED
pathToData='/homes/user/data/';
wb_command='/homes/user/workbench/bin_linux64/wb_command';

%% Do the parcellation

%%%%%%% IMPORTANT: PATHS SETTING %%%%%%%
% Both parcellation functions MUST be run from the main installation folder for the 
% ComputeCorrelationMatrix function to work properly. 
% Alternatively, the paths in ComputeCorrelationMatrix.m have to be edited accordingly.

% res_single and res_groups are structures containing the obtained parcellations and merged
% connectivity matrices.
% res_single structure fields: 
% -parcellation: number of vertices x 1, describes the parcel assignment of each vertex
% -Mat: number of parcels x number of parcels, describes the low dimensionality 
% 	  connectivity networl obtained after parcellation. 
% -Inlist: cell array number of parcels x 1, each cell Inlist{i} lists which 
% vertices belong to parcel i. Another way of describing the parcel assignments.

% res_group structure fields: 
% -parcellation: number of vertices x number of subjects, describes for each subject 
% the parcel assignment of each vertex
% -Mat: number of parcels x number of parcels x number of subjects, describes the low 
% dimensionality connectivity network obtained after parcellation for each subject. 
% -Inlist: cell array, number of subjects x 1, each cell Inlist{i}{j} lists which 
% vertices belong to parcel j for subject i. Another way of describing the parcel assignments.

%%%% Single subject parcellation
res_single = SingleParcellation(subjectlist{1},hemisphere,Nclus,pathToData,wbCommand);

%%%% Group-wise parcellation
res_group = SpectralParcellation(subjectlist,hemisphere,weight,Nclus,pathToData,wbCommand);


%%  Compute the group average parcellation and network

% Average parcellation, computed from all subjects' parcellations through majority voting. 
AvMap = mode(res_group.parcellation,2);

% Average connectivity matrix, average of all subjects' connectivity matrices
AvMat = mean(res_group.Mat,3);

%% Visualise

% Access the connectivity matrix of subject 100307
Mat_single = res_single.Mat;
Mat_group = res_group.Mat(:,:,1);

% Access the parcel assignment of subject 100307
Parc_single = res_single.parcellation;
Parc_group = res_group.parcellation(:,1);

% Visualise the parcellation using the gifti toolbox
filename='BrainSurf.gii'; %name of the gifti file representing the cortical surface
BrainSurf = gifti([pathToData subjectlist{1} '/' filename ]);
surf.cdata = Parc_group; % or Parc_single

figure;plot(BrainSurf,surf);



