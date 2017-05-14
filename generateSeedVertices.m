function [ seeds, idx ] = generateSeedVertices( vSphere, fSphere, vGray, labels, ratio )
%GENERATESEEDVERTICES Summary of this function goes here
%   Generate seed points using mesh resampling and knn search. results are
%   then masked to exlude vertices whose timeseries data are not available.
%   To run the functions iso2mesh library has to be added to the path.
%   vSphere: sphere vertices (obtained from 103414.L.sphere.32k_fs_LR.surf)
%   fSphere: sphere faces  (obtained from 103414.L.sphere.32k_fs_LR.surf)
%   vGray: real surface vertices  (obtained from 103414.L.midthickness.32k_fs_LR.surf) (for finding the coordinates of the seed points)
%   labels: functional vertices on the surface (1s cortical vertices, 0s non-cortical medial wall vertices)
%   ratio: downsample ratio

    [v,f]=meshresample(vSphere,fSphere,ratio);
    idx = knnsearch(vSphere,v);
    elim = find(labels == 0);
    neighs = arrayfun(@(x) find(idx == x, 1,'first'), elim, 'UniformOutput', false);
    A = cell2mat(neighs);
    idx(A) = [];
    seeds = vGray(idx,:);
end




% ratio = 0.01783; n = 500;
% ratio =  0.022; n = 600;
% ratio = 0.0255; n = 750;
% ratio = 0.0269; n = 800;
% ratio = 0.0336; n = 1000;
% ratio = 0.06725; n = 2000;
% ratio = 0.10097; n = 3000;