function [conn bicdat] = nwa_bic_glasso(tsMat,varargin)
% gLasso network estimation with the optimal lambda assessed with the
% bayesian information criterion (bic).
% USE: [conn bicdat] = nwa_bic_glasso(tsMat,varargin)
% =========================================================================
% IN: (time-series) data, time x nodes.
%   optional: 
%    'plot' - plot the trajectory
%    'lambda' - range of lambda values
%
% OUT: 
%   conn: connectivity matrix, converted to partial correlations
%   bicdat: data on the lamdbda estimation.
% =========================================================================

% defaults
plotdat = false;
lambda_range = [0.01:0.01:1];
pendiag = true;

% input
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'plot', plotdat = true;
            case 'lambda', lambda_range = varargin{i+1};
            case 'pendiag',pendiag = varargin{i+1};    
        end
    end
end

% start loop
nts  = size(tsMat,1);
nvar = size(tsMat,2);
for l = 1:length(lambda_range)
    lambda = lambda_range(l);
    [w, theta] = GraphicalLasso(tsMat,lambda);
    c = cov(tsMat);
    % glasso in r uses the covariance for LL calculations. 
    [bic,ll] = calcbic(theta,c,lambda,nts,pendiag);
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

% In case of a tie, pick the first
if length(size(theta_best))>2
   theta_best = theta_best(:,:,1);
end

% Convert to partial correlations (e.g. Cribben, 2013)
D = diag(theta_best);
N = size(D,1);
PrD1 = repmat(D,1,N);
PrD2 = repmat(D',N,1);
conn = (-1*theta_best)./sqrt(PrD1.*PrD2);
conn(logical(eye(N,N))) = 1;

% gather the output
bicdat.minbic = minbic;
bicdat.lambdabest = lambda_best;
bicdat.lambdarange = lambda_range;
bicdat.bic = BIC;
bicdat.LL = LL;

% plot trajectory
if plotdat 
    figure
    plot(lambda_range(1:l),BIC)
    hold on
    plot(lambda_range(minbic),BIC(minbic),'x')
end

end

function [bic,LL] = calcbic(theta,w,lambda,nts,diag0)
% https://www4.stat.ncsu.edu/~reich/BigData/code/glasso.html
nvar = size(theta,1);
E    = nwa_reshape(theta,'mat2vec');
% optional: set diagonal to zero
if diag0
    temp = nwa_proc_conn(theta,'diag0');
else
    temp = theta;
end
crit = -log(det(theta))+trace(w*theta)+sum(sum(abs(lambda.*temp)));
LL   = -(nvar/2)*crit;
bic = -2*LL+sum(E~=0)*log(nts);
end





