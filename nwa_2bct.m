function [NWA] = nwa_bct(NWA,varargin)

%% data format
if  iscell(NWA)
    nC = size(NWA,2);
elseif isstruct(NWA)
    C = NWA.data;
    nC = size(C,2);
else
    nC = 1;
    Cres = C;
    C = {};
    C{1} = Cres;
end

%% check and filter data
for j = 1:nC
    % NaN
    datacheck.nan(j,:) = squeeze(sum(sum(isnan(C{j}))));
    % zeros
    zerocheck = squeeze(sum(sum(C{j})))==0;
    datacheck.zeros(j,:) = double(zerocheck);
    % nMat(j) = size(C{j},3);
end
filter = [datacheck.nan; datacheck.zeros];
filter = sum(filter)>0;
for j = 1:nC; C{j} = C{j}(:,:,~filter); end

% add warining
if sum(filter)>0;
    fs = find(filter==1);
    for s = 1:length(fs)
        if isfield(NWA,'subjects');
            sub = NWA.subjects{fs(s)}; if isempty(sub); sub = 'jane doe'; end
        else
            sub = num2str(fs(s));
        end
        disp(['WARNING: subject ' sub ' is removed from the analyses'])
    end
end

NWA.filter = filter;
nMat = nC*sum(filter==0);

% defaults
% ---------
lmet = {'edges' 'strength' 'deg' 'bet' 'clus' 'eig' 'mod' 'eff'};
gmet = {'clus'};
% nrun  = size(PHYS,3);

% input
%-----------
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'lmet', lmet = varargin{i+1}; % local  measures
            case 'gmet', gmet = varargin{i+1}; % global measures
        end
    end
end

%% start analyses
progressbar_new('running the network analyses:')
count = 0;
for j = 1:nC
    
    cdat = C{j};
    nmat = size(cdat,3);
    
    for n = 1:nmat
        
        conn = cdat(:,:,n);
        nNode = size(conn,1);
        count = count + 1;
        progressbar_new(count/nMat)
        
        % Global measures
        % ====================
        for nC = 1:length(gmet);
            met = gmet{nC};
            switch met
                case 'clus'
                    NWA.global.clus{j}(n,nC) = mean(clustering_coef_bu(conn));
            end
        end
        
        % Local Measures
        % ======================
        for nL = 1:length(lmet);
            met = lmet{nL};
            switch met
                
                case 'edges'
                    edges = nwa_reshape(conn,'mat2vec');
                    NWA.local.edges{j}(n,:) = edges;
                    
                case 'deg'
                    % binarize values !! i.e all non-zero values set to 1
                    deg = degrees_und(conn);
                    NWA.local.deg{j}(n,:) = deg;
                    degdist = histc(deg,[1:length(deg)]);
                    NWA.local.degdist{j}(n,:) = degdist;
                    
                case 'bet'
                    % convert to connection-length matrix ! - normalize?
                    ind = conn~=0;
                    L = conn;
                    L(ind) = 1./conn(ind);
                    bet = betweenness_wei(L);
                    bet = bet;
                    NWA.local.bet{j}(n,:) = bet;
                    
                case 'strength'
                    c = nwa_proc_conn(conn,'abs','diag0');
                    s = sum(c);
                    NWA.local.strength{j}(n,:) = s;
                    
                case 'eig'
                    eig = eigenvector_centrality_und(conn);
                    NWA.local.eig{j}(n,:) = eig;
                    
                case 'eff'
                    % same as closeness centrality?
                    Eloc = efficiency_wei(conn,2);
                    NWA.local.eff{j}(n,:) = Eloc;
                    
                case 'clus'
                    clus = clustering_coef_wu(conn);
                    NWA.local.clus{j}(n,:) = clus;
                    
                case 'mod'
                    % note that the participation coeff seems less stable
                    % over different thresholds..
                    [Ci Q] = modularity_und(conn);
                    NWA.global.mod{j}(n,nC) = Q;
                    NWA.local.mod{j}(n,:) = Ci;
                    NWA.local.part{j}(n,:) = participation_coef(conn,Ci);
                
                otherwise
                warning(['unknown input: ' met])
            end
        end
    end
end
