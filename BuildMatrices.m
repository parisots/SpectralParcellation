function [W,Ctot,d]=BuildMatrices(C,Aff,R,nC,weight)

% % Construct global affinity and constraint matrices from individual matrices

%%%%%% INPUTS %%%%%%
% % C: cell matrix of all subjects' constraint matrices 1xNsubjects
% % Aff: cell matrix of all subjects' affinity matrices. 1xNlayers*Nsubjects
% % R: inter-subjects links, cell matrix NsubjectsxNsubjects
% % weight: multiplicative weights for inter-subjects links

%%%%%% OUTPUTS %%%%%%
% % W: Global affinity matrix
% % Ctot: Global constraint matrix
% % d: Degree matrix of W

% Copyright (C) Sarah Parisot, Imperial College London, 2015

%%

Nim= numel(C);
nlayers=numel(nC);

if nargin==4
    weight = 1;
end

%% Constraint matrix

Ctot=blkdiag(C{:});

%% Affinity matrix with inter image contraints: block matrix of all subjects affinities

  tic
 sizes=repmat(nC,1,Nim);
 idx=cumsum(sizes);
 allsz=(3000*3000+1000*1000+2000*2000)*50;
 offset=0;
 offset2=0;
 nx=zeros(1,allsz);
 ny=zeros(1,allsz);
 val=zeros(1,allsz);
 for i=1:nlayers*Nim
     [X,Y]=meshgrid(1:sizes(i),1:sizes(i));
     nx(offset2+1:offset2+length(Y(:)'))=[Y(:)'+offset];
     ny(offset2+1:offset2+length(Y(:)'))=[X(:)'+offset];
     val(offset2+1:offset2+length(Y(:)'))=sparse(Aff{i}(:));
     offset=idx(i);
     offset2=offset2+length(Y(:)');
 end
  neg=val<=0;
 nx(neg)=[];ny(neg)=[];val(neg)=[];

 tic
 [X,Y]=meshgrid(1:nC(end),1:nC(end));
 Npairs=Nim*(Nim-1)/2;
 off=length(val);
 nx=[nx,zeros(1,2*Npairs*nC(end)*nC(end))];ny=[ny,zeros(1,2*Npairs*nC(end)*nC(end))];val=[val,zeros(1,2*Npairs*nC(end)*nC(end))];
 for i=1:Nim-1
     for j=i+1:Nim
         
         SM=sum(nC(1:nlayers));
         Ri1 = R{j,i};
         Ri2 = R{i,j};
         
         offset=(j-1)*(SM)+sum(nC(1:2));
         offset2=(i-1)*(SM)+sum(nC(1:2));
         nx(off+1:off+length(Y(:)'))=[Y(:)'+offset];
         ny(off+1:off+length(Y(:)'))=[X(:)'+offset2];
         val(off+1:off+length(Y(:)'))=weight*Ri1(:);
         
         off=off+length(Y(:));
         
         nx(off+1:off+length(Y(:)'))=[Y(:)'+offset2];
         ny(off+1:off+length(Y(:)'))=[X(:)'+offset];
         val(off+1:off+length(Y(:)'))=weight*Ri2(:);
         off=off+length(Y(:));

     end
 end
 toc
 tic
  neg=val<=0;
 nx(neg)=[];ny(neg)=[];val(neg)=[];
 W=sparse(nx,ny,val,idx(end),idx(end));
 toc
  clear nx ny val

n = size(W,1);
offset=0.02;
d = sum((W),2);
offset=offset*mean(d);
d = d + offset * 2;
W = W + spdiags(offset*ones(n,1),0,n,n);
