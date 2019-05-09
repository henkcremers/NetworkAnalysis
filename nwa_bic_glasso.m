function [conn bicdat] = bic_glasso(tsMat)
%==========================================================================
% USE
%  IN
%  OUT
%  EXAMPLE
% =========================================================================
plotdat = 0;
lambda_range = [0.01:0.01:1];
nts  = size(tsMat,1);
nvar = size(tsMat,2);
for l = 1:length(lambda_range)
    lambda = lambda_range(l);
    [w, theta] = GraphicalLasso(tsMat,lambda);
    [bic,ll] = calcbic(theta,w,lambda,nts);
    BIC(l) = bic;
    LL(l) = ll;
    theta_temp(:,:,l) = theta;
    
    minbic = find(BIC==min(BIC));
    if l-minbic>10
        theta_best = theta_temp(:,:,minbic);
        lambda_best = lambda_range(minbic);
        clear theta_temp
        break
    else
        theta_best = theta_temp(:,:,minbic);
        lambda_best = lambda_range(minbic);
    end
end

% Convert to partial correlations (see e.g. Cribben, 2013)
% ------------------------------------------
D = diag(theta_best);
N = size(D,1);
PrD1 = repmat(D,1,N);
PrD2 = repmat(D',N,1);
conn = (-1*theta_best)./sqrt(PrD1.*PrD2);
conn(logical(eye(N,N))) = 1;

%% save the data
bicdat.minbic = minbic;
bicdat.lambdabest = lambda_best;
bicdat.lambdarange = lambda_range;
bicdat.bic = BIC;
bicdat.LL = LL;

% plot trajectory
if plotdat == 1;
    figure
    plot(lambda_range(1:l),BIC)
    hold on
    plot(lambda_range(minbic),BIC(minbic),'x')
end

end

function [bic,LL] = calcbic(theta,w,lambda,nts)
nvar = size(theta,1);
E    = nwa_reshape(theta,'mat2vec');
crit = -log(det(theta))+trace(w*theta)+sum(sum(abs(lambda.*theta)));
LL   = -(nvar/2)*crit;
% then BIC
% based on: https://www4.stat.ncsu.edu/~reich/BigData/code/glasso.html
% Foygel, "Extended Bayesian Information Criteria for Gaussian Graphical
% Models"
bic = -2*LL+sum(E~=0)*log(nts);
end




