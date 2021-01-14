%%%%%%%%% PBC data fitting
load pbc
pbc(:,1)=log(pbc(:,1)/365.25); % First column: log(T)
pbc(:,5:7)=log(pbc(:,5:7)); % log-transformation for some covariates
pbc(:,2)=[pbc(:,2)==2]; % Indicate censoring status
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pbc_data=[pbc(:,1),pbc(:,[3,4,5,7,6]),pbc(:,2)];
% Fit PBC data
[beta,support,weight,cons]=aftmix(pbc_data);

% beta: parameter estimates
% support: support for mixing distribution
% weight: weight for each support in the mixing distribution
% cons: the lower bound of the support for the mixing distirbution
