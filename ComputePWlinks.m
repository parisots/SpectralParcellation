function [R,match] = ComputePWlinks(SCmatrix1,SCmatrix2,data,Adj)

% % Compute inter-subject connections based on supervertex similarity in
% terms of connectivity

%%%%%% INPUTS %%%%%%
% % SCmatrix1: structure containing all supervertex parcellation
% information (relevant here: supervertex assignments for all vertices and 
% merged connectivity matrices, supervertex centres), subject 1. Output of 
% MultiAffinity.m
% % SCmatrix2: same as SCmatrix2, subject 2

% % data: structure listing all relevant structural information on the
% surface mesh (relevant here: mesh neighbours and surface information for
% visualisation purposes). Output of loadSubjectData.m
% % Adj : OPTIONAL. NsuperverticesxNsupervertices, binary matrix. Has a 1
% if two supervertices are neighbours, 0 otherwise. 

%%%%%% OUTPUTS %%%%%%
% % W: Global affinity matrix
% % Ctot: Global constraint matrix
% % d: Degree matrix of W

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


%% pairwise matching SVlevel

Inc=data.Inc;

%%% only connect the coarsest resolution
Map1=SCmatrix1.baseSeg.Maps{end};
Map2=SCmatrix2.baseSeg.Maps{end};

Cor1=corrcoef(SCmatrix1.baseSeg.Mat{end}');
Cor2=corrcoef(SCmatrix2.baseSeg.Mat{end}');


%%% Construct adjacency matrix, Adj(i,j)=1 specifies that supervertices i
%%% and j are neighbours
if nargin<4
    Adj = zeros(size(Cor2));
    for j=1:numel(Inc)
        k=Inc(j);
        neigh = intersect(Inc,data.LV{k});
        Cnn=Map2(neigh);
        if numel(unique([Cnn;Map2(k)]))>1
            Adj(Map2(k),unique(Cnn))=1;
        end
    end
end



Nseeds=size(Adj,1);

% for a given SV in subject1, find the one in subject2 that has the highest
% amount of overlap (number of vertices)
Connections = zeros(Nseeds,1);
for i=1:Nseeds
    
    
    tmp=zeros(32492,1);
    tmp(Map1==i)=Nseeds;
    a=find((tmp+Map2)>Nseeds);
    clus=unique(Map2(a));
    n1(i)=0;
    for j=1:numel(clus)
        nb=numel(find(Map2(a)==clus(j)));
        if nb>n1(i)
            n1(i)=nb;
            Connections(i)=clus(j);
        end
    end
end

%%% find which SV (closest SV or one in the immediate neighbourhood) has
%%% the largest correlation with subject1's SV
R=zeros(Nseeds,Nseeds);
for i=1:Nseeds
    newMap=zeros(numel(Inc),1);
    
    list=find(Adj(Connections(i),:));
    %%% Uncomment to look at a larger neighbourhood
    %     [a,b]=find(Adj2(list,:)); 
    %     list=[list,b'];
    %     list=unique(list);
    
    newMap2=zeros(numel(Inc),numel(list));
    
    for l=1:Nseeds
        in1=find(Map1(Inc)==l);
        in2=find(Map2(Inc)==l);
        newMap(in1)=Cor1(i,l);
        for j=1:numel(list)
            newMap2(in2,j)=Cor2(list(j),l);
        end
    end
    
    
    
    
    Ccoef= corr(newMap,newMap2);
    
    
    [~,idx]=nanmax(Ccoef);
    match(i,1)=list(idx);
    match(i,2)=nanmax(Ccoef);
    
    R(i,match(i,1))=match(i,2);
    
%%% Visualisation    
%     surf.cdata=zeros(32492,1);
%     for j=1:numel(list)
%         surf.cdata(Map2==list(j))=Ccoef(j);
%     end
%     surf.cdata(SCmatrix1.baseSeg.Maps{3}==i)=surf.cdata(SCmatrix1.baseSeg.Maps{3}==i)+1;
%     figure;plot(data.Inflated,surf)
%     view(-94,2)
%     caxis([min(Ccoef)-0.01 max(Ccoef)]);
%     hold on
%     scatter3(data.Inflated.vertices((SCmatrix1.baseSeg.Seeds{end}(i)),1),data.Inflated.vertices((SCmatrix1.baseSeg.Seeds{end}(i)),2),data.Inflated.vertices((SCmatrix1.baseSeg.Seeds{end}(i)),3));hold on;
%     scatter3(data.Inflated.vertices((SCmatrix2.baseSeg.Seeds{end}(match(i,1))),1),data.Inflated.vertices((SCmatrix2.baseSeg.Seeds{end}(match(i,1))),2),data.Inflated.vertices((SCmatrix2.baseSeg.Seeds{end}(match(i,1))),3))
%     
%     %

    
    
    
end
