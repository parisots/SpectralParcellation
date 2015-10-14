function [NcutDiscrete,NcutEigenvectors,Eigenvalues]=Ncut_multi(Wtot,C,d,Nclus)

% Eigenvector decomposition of the affinity matrix Wtot subject to the
% constraint matrix C following the multiscale normalised cuts method 
% (Cour et al. CVPR 2005)

%%%%%% INPUTS %%%%%%
% % Wtot: Affinity matrix
% % C: Constraint matrix
% % d: degree matrix of Wtot
% % Nclus: number of desired parcels

%%%%%% OUTPUTS %%%%%%
% % NcutDiscrete: Discrete parcel assignments
% % NcutEigenvectors: Eigenvectors of the 
% % CorMat: Correlation Matrix of MatM

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

n=size(Wtot,1);
Dinvsqrt = 1./sqrt(d+eps);
P = spmtimesd(Wtot,Dinvsqrt,Dinvsqrt);

Ctot = spmtimesd(C,[],Dinvsqrt);
Prod = Ctot*Ctot';

opts.type = 'ict';
opts.droptol = 1e-12;
opts.shape='upper';
Rc=ichol(Prod,opts);

options=getDefaultOptionsEigs(n,Nclus);
[result,convergence] = eigs_compatible_with_eigs_optimized(@mex_projection_QR_symmetric,[],Nclus,options,tril(P),Ctot',Rc);

Eigenvalues = result.lambda;
Eigenvectors = spdiags(Dinvsqrt,0,n,n) * result.X;



%%% normalise eigenvectors
for  i=1:size(Eigenvectors,2)
    Eigenvectors(:,i) = (Eigenvectors(:,i) / norm(Eigenvectors(:,i))  )*norm(ones(n,1));
    if Eigenvectors(1,i)~=0
        Eigenvectors(:,i) = - Eigenvectors(:,i) * sign(Eigenvectors(1,i));
    end
end

%%% discretisation of the solution to the relaxed problem
[NcutDiscrete,NcutEigenvectors] =discretisation(Eigenvectors);