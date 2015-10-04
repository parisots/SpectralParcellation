function SCmatrices = MultiAffinity(Sphere,BrainSurf,roi,CorMat,Mat,sp,list_vertices)

% Constructs the single subject supervertex parcellations (multiple 
% resolutions) and corresponding affinity and constraint matrices

%%%%%% INPUTS %%%%%%
% % Sphere: as loaded from a gifti file, sphere representation of the brain
% % surface
% % BrainSurf: as loaded from a gifti file, the brain surface
% % roi: metric file, zero valued on medial wall vertices, one on the rest
% % CorMat: NverticesxNvertices, where Nvertices is the number of non medial wall vertices.
% % Correlation matrix
% % sp: 0/1 specifies whether we are constructing spatially constrained connectivity matrices
% % list_vertices: cell array Nverticesx1, each cell lists the neighbours 
%  of a given vertex   

%%%%%% OUTPUTS %%%%%%
% % SCmatrices: structure. Contains the base parcellation (assignments, 
% centres, merged connectivity matrices), and the multiresolution affinity and 
% constraint matrices

% Copyright (C) Sarah Parisot, Imperial College London, 2015

%% Parameters
nlayers=3;
Inc = find(roi.cdata==1);
MWlist=find(roi.cdata==0);



%% Create multilevel supervertices
%%%% 3 levels: 3000 2000 1000

%%% Parameters to sample the desired number of uniform seeds on the surface

% Right Hemisphere
%%%% 3000 seeds:  0.10097
%%%% 2000 seeds:  0.06718
%%%% 1000 seeds:  0.03365


% Left Hemisphere
%%%% 3000 seeds:  0.10097
%%%% 2000 seeds:  0.06725
%%%% 1000 seeds:  0.03352

if length(Inc)==29716 % right hemisphere
    ratio = [0.10097,0.06718,0.03365];
else % left hemisphere
    ratio = [0.10097,0.06725,0.03352];
end

for i=1:nlayers
    %%% sample uniform seeds on the surface, use the sphere for this
    [ Sd, idx ] = generateSeedVertices( Sphere.vertices, Sphere.faces, BrainSurf.vertices, roi.cdata, ratio(i) );
    %%% construct the supervertex parcellations
    [Map{i}, Seeds{i}, Mn{i}] = Iterative_SV(CorMat,double(idx),list_vertices,Inc,MWlist, double(BrainSurf.vertices), double(BrainSurf.faces),3);
    %%% merge the tractography matrix to obtain a NseedxNseed connectivity
    %%% matrix (one connectivity profile per supervertex)
    [MatSV{i},Inlist{i},CM{i}] = MergeConnectivityMatrix(Seeds{i},Map{i},Mat(Inc,Inc),Inc);
end


% number of supervertices in each layer
for i=1:nlayers
    nC(i)=numel(Seeds{i});
end


%% Spatially constrained lower layer
%%% only keep connection between supervertices that are neighbours to
%%% obtain spatially contiguous parcels
if sp
    
    

for i=1:nlayers
SurMat = sparse(size(CM{i}));
for j=1:numel(Inc)
    k=Inc(j);
    neigh = intersect(Inc,list_vertices{k});
    Cnn=Map{i}(neigh);
    if numel(unique([Cnn;Map{i}(k)]))>1
        SurMat(Map{i}(k),unique(Cnn))=1;
    end
end

for j=1:length(CM{i})
    lv{j}=find(SurMat(j,:)==1);
end


order =1; %% only keep immediate neighbours
CoM=single(zeros(size(CM{i})));
for k=1:size(CM{i},1)
    LV{k,1}=lv{k};
    LVtot=LV{k,1};
    for l=2:order        
        LV{k,l} = setdiff([lv{[LV{k,:}]}],LVtot);
        LVtot = [LV{k,:}];
    end
     CoM(k,LVtot)=CM{i}(k,LVtot);
end
CoM = max(CoM,0);
 CM{i}=double(CoM);
end

end



%% Inter layer connections

%%% Connect to supervertices at different resolution if they share vertices
%%% on the orginal surface. Strength of the connection is the amount of
%%% overlap

Conn = cell(1,nlayers-1);
for k=1:nlayers-1
    l=k+1;
    Nclus=numel(unique(Map{k}))-1;
    Nclow=numel(unique(Map{l}))-1;
    Conn{k}=zeros(nC(l),nC(k));
    for i=1:Nclow
        tmp=zeros(size(Map{l}));
        tmp(Map{l}==i)=Nclus;
        a=find((tmp+Map{k})>Nclus);
        clus=unique(Map{k}(a));
        for j=1:numel(clus)
        Conn{k}(i,clus(j))=numel(find(Map{k}(a)==clus(j)))/numel(find(Map{l}==i));%max(numel(find(Map{k}==i)),numel(find(Map{l}==clus(j))));
        end
    end
end

%% Constraint matrix


Ctot = 1*blkdiag(Conn{:});
offsetx=0;
offsety=0;
for i=1:nlayers-1
    
    Ctot(offsetx+1:offsetx+nC(i+1),offsety+nC(i)+1:offsety+nC(i+1)+nC(i))=-speye(nC(i+1));
    offsetx=offsetx+nC(i+1);
    offsety=offsety+nC(i);
end

Ctot=sparse(Ctot);

%% Block weight matrix

Wtot = blkdiag(CM{:});
Wtot=sparse(Wtot);
Wtot(Wtot<0)=0;

n = size(Wtot,1);
offset=0.02;
d = sum((Wtot),2);
offset=offset*mean(d);
d = d + offset * 2;
Wtot = Wtot + spdiags(offset*ones(n,1),0,n,n);

SCmatrices=struct;
SCmatrices.W = Wtot;
SCmatrices.d = d;
SCmatrices.C = Ctot;
SCmatrices.Neigh = Conn;
SCmatrices.baseSeg.Seeds = Seeds;
SCmatrices.baseSeg.Maps = Map;
SCmatrices.baseSeg.Mn = Mn;
SCmatrices.baseSeg.Mat = MatSV;
SCmatrices.baseSeg.Cor = CM;
SCmatrices.baseSeg.Inlist = Inlist;

    
