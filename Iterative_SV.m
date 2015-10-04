function [Map, Seeds, Mn] = Iterative_SV(CoM,Sd,list_vertices,Inc,MWlist,vertex,faces,S)

% constructs a supervertex parcellation of a triangular surface mesh from a given affinity matrix

%%%%%% INPUTS %%%%%%
% % CoM: Affinity matrix NverticesxNvertices
% % Sd: List of initial seeds (must correspond to vertices on the mesh)
% % list_vertices: cell array Nverticesx1, each cell lists the neighbours 
%  of a given vertex
% % Inc: list of vertices included in the parcellation scheme, useful to 
%  avoid the medial wall on a cortical surface. Lists all vertices if no
% exclusions
% % MWlist: list of excluded vertices (complementary to Inc). Emtpy if no
% exclusions
% % vertex: Nverticesx3, 3D coordinates of all vertices on the mesh
% % faces: Ntrianglesx3, each row lists the vertices corresponding to a mesh triangle
% % S: weighting parameter for the fast marching speed function

%%%%%% OUTPUTS %%%%%%
% % Map: Nverticesx1, supervertex assignment of all vertices
% % Seeds: final supervertices centres (updated version of the Sd input)
% % Mn: mean correlation within a supervertex

% Copyright (C) Sarah Parisot, Imperial College London, 2015

%%

Nseeds=numel(Sd); %% number of seeds 
Nvertices=numel(list_vertices);

%%% lists the new centres and whether they differ from the previous iteration
NCenter = zeros(Nseeds,2); 

for it=1:20 %% maximum number of iterations
    if it==1
        Seeds=Sd;
    else
        Seeds = NCenter(:,1);
    end

    %%% set fast marching options
    options.constraint_map = Inf.*ones(Nvertices,1);
    options.constraint_map(setdiff(1:Nvertices,Inc)) = -Inf;
    options.dmax=20; %% limit maximum geodesic distance computation
    options.W=ones(Nvertices,1);

    %%% compute the geodesic distance from new centres
    New=find(NCenter(:,2)==0);
    for i=1:numel(New)
        options.W(Inc) = double(exp(-S.*CoM(Inc==Seeds(New(i)),:)));
        Dn(:,New(i)) = perform_fast_marching_mesh(vertex, faces, Seeds(New(i)),options);
    end
    
    
    %%% Do the parcellation by minimising the geodesic distance from all
    % centres
    [~,Map]=min(Dn(:,1:Nseeds),[],2);
    Map(MWlist)=0;
    [~,tmp]=sort(Map(Inc));
    
    %%% Compute the average correlation within a supervertex and update the
    %%% center
    offset = 0;
    for i=1:Nseeds
        LC = find(Map(Inc)==i);
        SzSV = numel(LC);
        Cor = CoM(tmp(offset+1:offset+SzSV),tmp(offset+1:offset+SzSV));
        Mx= find(sum(Cor)==nanmax(sum(Cor)));
        NCenter(i,1) = Inc(tmp(offset+Mx(1)));
        NCenter(i,2) = isequal(NCenter(i,1),Seeds(i));
        offset=offset+SzSV;
        Mn{i,it}=Cor;
        tmpval(i,it)=median(median(Cor));
    end
    Ns(it)=sum(NCenter(:,2));
    disp(sum(NCenter(:,2)))
    disp(mean(tmpval(:,1:it)))
    
    %%% Break condition: no centre updates
    if sum(NCenter(:,2))==Nseeds
        %%% optional: check whether supervertices are fully continuous
        for i=1:numel(Seeds)
            inclus = find(Map==i);
            LV=list_vertices(inclus);
            LV=unique([LV{:}]);
            
            iscont(i)=isempty(setdiff(inclus,intersect(LV,inclus)));
        end
        
        disp(sum(iscont)) 
        break
    end
    
    
    %%% optional: check whether supervertices are fully continuous
    if it==20 %% maximum number of iterations
        
        for i=1:numel(Seeds)
            inclus = find(Map==i);
            LV=list_vertices(inclus);
            LV=unique([LV{:}]);
            
            iscont(i)=isempty(setdiff(inclus,intersect(LV,inclus)));
        end
        
        disp(sum(iscont))
        
    end
    
end



