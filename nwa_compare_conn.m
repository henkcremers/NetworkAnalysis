function [eCorr cCorr Mmean] = nwa_comp_connec(M1,varargin);

% defaults
% -----------------------
matpart = {'ll'};
plotfig = 0;
Mmean = [];
M2 = [];

% user input
%-------------------------
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'matpart', matpart = varargin{i+1};   % part of the matrix (lower/upper triangle)
            case 'plot', plotfig = 1;   % part of the matrix (lower/upper triangle)
            case 'M2', M2 = varargin{i+1};
        end
    end
end

nn = size(M1);
L = tril(ones(nn(1),nn(1)),-1);
U = triu(ones(nn(1),nn(1)),1);

if isempty(M2)
    
    v1dat = M1(logical(L));
    M2    = rot90(fliplr(M1),1);
    v2dat = M2(logical(L));
    
    if plotfig==1
    ddat = v1dat - v2dat;
    figure
    hist(ddat)
    end
    
    if isequal(v1dat,v2dat);
        error('this is a symmetric matrix, nothing to compare')
        
    else
        Mmean = mean([v1dat,v2dat],2);
        Mmean = nwa_reshape(Mmean,'vec2mat');    
    end

else
    % M2 = varargin{1};
    
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

% Column wise correlations
for j = 1:nn(1);
    a = M1(:,j); a(j) = [];
    b = M2(:,j); b(j) = [];
    cdat(j) = corr(a,b);
end

dens1 = sum(v1dat~=0)/length(v1dat);
dens2 = sum(v2dat~=0)/length(v2dat);

b = [-1:0.05:1];
[h1 x1] = histc(v1dat,b);
[h2 x2] = histc(v2dat,b);

% main output
eCorr = corr(v1dat,v2dat);
cCorr = mean(cdat);

% compare node strength 
strength1 = sum(abs(M1))';
strength2 = sum(abs(M2))';
sCorr = corr(strength1,strength2);

% eig1 = eigenvector_centrality_und(M1);
% eig2 = eigenvector_centrality_und(M2);
% eigCorr = corr(eig1,eig2)
% 
% eff1 = efficiency_wei(M1,2);
% eff2 = efficiency_wei(M2,2);
% effCorr = corr(eff1,eff2)

if plotfig == 1
    figure;
    subplot(4,1,1);
    p1 = plot(b,h1,'k');
    hold on
    p2 = plot(b,h2,'r');
    legend(['Matrix 1: dens' num2str(dens1)],['Matrix 2: dens' num2str(dens2)])
    subplot(4,1,2);
    scatter(v1dat,v2dat);
    xlabel('values matrix 1');
    ylabel('values matrix 2');
    title(['element-wise correlation: ' num2str(eCorr)])
    subplot(4,1,3)
    hist(cdat);
    title(['colum-wise correlation: ' num2str(cCorr)])
    subplot(4,1,4)
    scatter(strength1,strength2);
    title(['node-strength correlation: ' num2str(sCorr)])
end
end

