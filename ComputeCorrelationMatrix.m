function CorMat = ComputeCorrelationMatrix(subjectID,pathToData,pathToCorMat,wbCommand)

% % Compute the correlation matrix, assuming FSL's probtrackX and HCP data
% are used (or that the data has the same organisation)

% % This function must be run from the main installation folder, or the path to the SeedSpaceMetric.func.gii 
% must be modified accordingly

%%%%%% INPUTS %%%%%%
% % subjectID: identification of the subjects to be parcellated
% % pathToData: base path where the data is stored
% % wbCommand: workbench command, typically
% 'path/to/workbench/installation/wb_command'
% % pathToCorMat: path where the tractography matrix is stored

%%%%%% OUTPUTS %%%%%%
% % CorMat: NverticesxNvertices Correlation matrix

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
    
    savePath=[pathToData subjectID pathToCorMat hemisphere];
    
    if ~exist('fdt_matrix1log.dconn.nii')
        
        %%% convert the output of probtrackX
        %%% SeedSpaceMetric.func.gii (provided) is a metric file with ones on all
        %%% vertices
        unix([wbCommand '-probtrackx-dot-convert' savePath '/fdt_matrix1.dot' savePath 'fdt_matrix1.dconn.nii -row-surface ./Dependencies/SeedSpaceMetric.func.gii -col-surface ./Dependencies/SeedSpaceMetric.func.gii -transpose']);
        
        %%% compute the log transform of the matrix
        cd savePath
        cii = ciftiopen('fdt_matrix1.dconn.nii',wbCommand);
        cii.cdata=log(cii.cdata);
        cii.cdata(cii.cdata==-Inf)=0; %% minimal non zero value of the matrix is 1
        
        
        ciftisave(cii,'fdt_matrix1log.dconn.nii',wbCommand)
        
        delete fdt_matrix1.dconn.nii
        
        
    end
    
    cd savePath
    %%% compute the correlation matrix: correlate the rows of the matrix
    if ~exist('Corlog.dconn.nii')
        % can be replaced by Matlab's corrcoef function. 
        unix([wbCommand '-cifti-correlation fdt_matrix1log.dconn.nii Corlog.dconn.nii']); 
    end
    
    
    %%% load correlation matrix
    cii = ciftiopen('Corlog.dconn.nii',wbCommand);
    CorMat=cii.cdata(data.Inc,data.Inc);
    clear cii
    CorMat(isnan(CorMat))=min(CorMat(:));
    