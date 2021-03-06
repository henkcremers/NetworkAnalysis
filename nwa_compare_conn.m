function [eCorr cCorr compdat] = nwa_compare_conn(M1,M2,varargin);
% Compare two connectivity matrices M1 and M2
% USE: [eCorr cCorr compdat] = nwa_compare_conn(M1,M2,varargin);
% =========================================================================
% IN: M1 - matrix 1 
%     M2 - Matrix 2 
% Optional:
%     'matpart' - part of the tmatir 'l' or 'u' (lower/upper triangel)
%     'plot' - plot the overalap 
%     'central' -  also compare the node centrality
% OUT: eCorr - element-wise spearman correlation 
%      cCorr -  column wise spearman correlation
%      compdat - some additional data. 
%
% =========================================================================


% Defaults
% -----------------------
matpart = {'ll'}; % lower (l) or upper (u) triangle of matrix.
plotfig = false;  % plot figure(s)
central = false;  % also test the centrality of the connectivity

% User input
%-------------------------
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'matpart', matpart = varargin{i+1};   % part of the matrix (lower/upper triangle)
            case 'plot', plotfig = true;               % plot output
            case 'central', central = true;            % analyse centrality
                %case 'M2', M2 = varargin{i+1};
        end
    end
end

%% check the input.
nn = size(M1);
L = tril(ones(nn(1),nn(1)),-1);
U = triu(ones(nn(1),nn(1)),1);
if isequal(M1,M2); % isempty(M2)
    disp('Same matrix - Compare upper and lower triangle')
    
    v1dat = M1(logical(L));
    M2    = rot90(fliplr(M1),1);
    v2dat = M2(logical(L));
    
    if plotfig
        ddat = v1dat - v2dat;
        figure
        hist(ddat)
    end
    
    if isequal(v1dat,v2dat);
        error('this is a symmetric matrix, nothing to compare')
        
    else
        Mmean = mean([v1dat,v2dat],2);
        compdat.Mmean = nwa_reshape(Mmean,'vec2mat');
    end
    
else
    
    if isequal('l',matpart{1}(1))
        v1dat = M1(logical(L));
    else
        %v1dat = M1(logical(U));
        M1flip = rot90(fliplr(M1),1);
        v1dat  = M1flip(logical(L));
    end
    
    if isequal('l',matpart{1}(2))
        v2dat = M2(logical(L));
    else
        %v2dat = M2(logical(U));
        M2flip = rot90(fliplr(M2),1);
        v2dat  = M2flip(logical(L));
    end
end

%% Main comparison
for j = 1:nn(1);
    a = M1(:,j); a(j) = [];
    b = M2(:,j); b(j) = [];
    cdat(j) = corr(a,b);
end

% gather output
eCorr = corr(v1dat,v2dat);
cCorr = mean(cdat);
d1 = v1dat~=0;dens1 = sum(d1)/length(v1dat); compdat.dens1 = dens1;
d2 = v2dat~=0;dens2 = sum(d2)/length(v2dat); compdat.dens2 = dens2;

% network dice coef. 
cd = [d1 d2];
dc = (sum(sum(cd,2)>1)*2)/(sum(d1)+sum(d2));
compdat.dc = dc;

if plotfig | plotfig== 1
    b = [-1:0.05:1];
    [h1 x1] = histc(v1dat,b);
    [h2 x2] = histc(v2dat,b);
    
    figure('Name','Matrix Comparison');
    subplot(3,1,1);
    p1 = plot(b,h1,'k');
    hold on
    p2 = plot(b,h2,'r');
    legend(['Matrix 1: dens' num2str(dens1)],['Matrix 2: dens' num2str(dens2)])
    subplot(3,1,2);
    scatter(v1dat,v2dat);
    xlabel('values matrix 1');
    ylabel('values matrix 2');
    title(['element-wise correlation: ' num2str(eCorr)])
    subplot(3,1,3)
    hist(cdat);
    title(['colum-wise correlation: ' num2str(cCorr)])
end

%% centrality
if central | central==1
    
    % compare node strength
    strength1 = sum(abs(M1))';
    strength2 = sum(abs(M2))';
    sCorr = corr(strength1,strength2); compdat.sCorr = sCorr;
    
    % eigen vector
    eig1 = eigenvector_centrality_und(abs(M1));
    eig2 = eigenvector_centrality_und(abs(M2));
    eigCorr = corr(eig1,eig2); compdat.eigCorr = eigCorr;
    
    % efficiency
    eff1 = efficiency_wei(abs(M1),2);
    eff2 = efficiency_wei(abs(M2),2);
    effCorr = corr(eff1,eff2); compdat.effCorr = effCorr;
    
    % clustering
    clus1 = clustering_coef_wu(abs(M1));
    clus2 = clustering_coef_wu(abs(M2));
    clusCorr = corr(clus1,clus2); compdat.clusCorr = clusCorr;
    
    % modularity and particpation coef.
    [Ci Q1] = modularity_und(M1);
    part1 = participation_coef(M1,Ci);
    [Ci Q2] = modularity_und(M2);
    part2 = participation_coef(M2,Ci);
    partCorr = corr(part1,part2); compdat.partCorr = partCorr;
end

if (plotfig | plotfig== 1) & (central | central==1)
    
    figure('Name','Centrality Measures');
    subplot(5,1,1)
    scatter(strength1,strength2);
    title(['Node-strength correlation: ' num2str(sCorr)])
    
    subplot(5,1,2)
    scatter(eig1,eig2);
    title(['Eigenvector correlation: ' num2str(eigCorr)])
    
    subplot(5,1,3)
    scatter(eff1,eff2);
    title(['Efficiency correlation: ' num2str(effCorr)])
    
    subplot(5,1,4)
    scatter(clus1,clus2);
    title(['Clustering correlation: ' num2str(clusCorr)])
    
    subplot(5,1,5)
    scatter(part1,part2);
    title(['Participation correlation: ' num2str(partCorr)])
end

end

