function [BIC,LL] = nwa_bic(theta,w,labda,nts)
nvar = size(theta,1);
E = nwa_reshape(theta,'mat2vec');
%% first calc the log likelihood
% from the glasso (R) code: 
%     critfun = function(Sigmahati, s, rho, penalize.diagonal = TRUE) {
%         d = det(Sigmahati)
%         temp = Sigmahati
%         if (!penalize.diagonal) {
%             diag(temp) = 0
%         }
%         val = -log(d) + sum(diag(s %*% Sigmahati)) + sum(abs(rho * 
%             temp))
%         return(val)
%     }
%     crit = critfun(xx, s, rho, penalize.diagonal)
%     loglik = -(n/2) * crit

% theta_d0 = nwa_procmat(theta,'diag0');
crit = -log(det(theta))+trace(w*theta)+sum(sum(abs(labda.*theta)));
LL   = -(nvar/2)*crit;

%% then BIC
% based on: https://www4.stat.ncsu.edu/~reich/BigData/code/glasso.html
% Foygel, "Extended Bayesian Information Criteria for Gaussian Graphical
% Models"
BIC = -2*LL+sum(E~=0)*log(nts);

end