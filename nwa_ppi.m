function [ppinw] = nwa_ppi(PHYS,PSY,varargin)
% =========================================================================
% ppi network analyses
% =========================================================================
% USE: [ppinw] = nwa_ppi(PHYS,PSY,varargin)
% IN: PHYS - timeseries data
%     PSY -  task regressors 
%   optional 
%   'contrast', contrast vector or matrix
%   'confound', - confounding regressors 
%   'regression'
%    correction with a AR(1) model (default). or select other: 
%    'armodel', armodel =  varargin{i+1};
%   'physadjust', correction method for the physiological data
%   'psyadjust',  correction method for the psychological data
%   'taskthr',  theshold to define the task 
%   'taskbin', binarize task regressor
%   'connest', connectivity estimation
% 
% OUT: ppinw - data strucuture
% =========================================================================

[nt nreg nrun] = size(PHYS);
[nc ndv nruncheck]  = size(PSY);

if nt ~= nc
    error('number of time points for PSY and PHYS need to be the same')
end

if nrun ~= nruncheck
    error('number of runs for PSY and PHYS need to be the same')
end

%% defaults
% --------------------
Xconf = []; % confound regressors (motion etc.)
ncon  = 1;
contrast   = [];
nrun  = size(PHYS,3);

% adjustment of the regressors: 
% demean ('dm'); decenter ('dc'); or zscore ('z'), or 'none'
physadjust = 'dm';
psyadjust  = 'dc';
taskthr = 0.1;   % theshold to define the task 
taskbin = false; % binarize the task regressor

% regression method
% 'ar' or 'standard) and (in case of ar), specify the model order  (default is 1).
regression = 'ar';  % or 'standard'
armodel = 1;

% connectivity estimation of the resdiudal time-series
connest    = 'glasso';

% high pass filter for times-series
hpf = 128;

%% user input
%--------------------
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'contrast', contrast = varargin{i+1}; ncon = size(contrast,1);   % contrast vector or matrix
            case 'confound', Xconf = varargin{i+1};                % high pass filter ??
            case 'method', method = varargin{i+1};
            case 'physadjust', physadjust =  varargin{i+1};
            case 'psyadjust',  psyadjust =  varargin{i+1};
            case 'taskthr',  taskthr = varargin{i+1}; % theshold to define the task 
            case 'taskbin', taskbin = true;
            case 'regression', regression =  varargin{i+1};
            case 'armodel', armodel =  varargin{i+1};
            case 'connest', connest =  varargin{i+1};
        end
    end
end

%% add settings to structure
ppinw.settings.contrast = contrast;
ppinw.settings.adjust.phys = physadjust;
ppinw.settings.adjust.psy = psyadjust;
ppinw.settings.adjust.taskbin = taskbin;
ppinw.settings.adjust.taskthr = taskthr;
ppinw.settings.regression.method = regression;
ppinw.settings.regression.ar= armodel;
ppinw.settings.connest = connest;
ppinw.settings.hpf = hpf;

npsy  = size(PSY,2);
disp('running the analyses: ')

for run = 1:nrun
    disp([' run:' num2str(run)])
    
    % get the data for the run
    PHYSrun  = PHYS(:,:,run);
    PSYrun   = PSY(:,:,run);
    Xconfrun = Xconf(:,:,run);
    
    %% Rank Intensity
    disp('  ..intensity rank')
    m = mean(PHYSrun);
    [rank] = tiedrank(m);
    irank = (rank./nreg)*100;
    ppinw.run(run).irank = irank;
    
    % high-pas filter of PHYS
    
    TR = 2;
    K.row=1:nt;
    K.RT=TR;
    K.HParam=hpf;
    % removed filtering of the task regressor
    % PSYrun  = spm_filter(K,PSYrun);
    PHYSrun = spm_filter(K,PHYSrun);
    
    % time-series adjustment
    switch physadjust
        case 'dm'
            mPHYS = mean(PHYSrun);
            mPHYS = repmat(mPHYS,nt,1);
            PHYSrun = PHYSrun - mPHYS;
        case 'dc'
            maxval          = max(PHYSrun);
            minval          = min(PHYSrun);
            subval          = mean([maxval;minval]);
            PHYSrun   = PHYSrun - repmat(subval,nt,1);
        case 'z'
            PHYSrun = zscore(PHYSrun);
        case 'none'
    end
    
    % task  adjustment
    
    % first define the task/rest
    taskloc = logical(sum(PSYrun>taskthr,2));
    
    if taskbin
    disp('NOTE: Binarizing the task regressors')    
    PSYrun = double(PSYrun>taskthr);
    end
    
    switch psyadjust
        case 'dm'
            mPSY = mean(PSYrun);
            mPSY = repmat(mPSY,nt,1);
            PSYrun = PSYrun - mPHYS;
        case 'dc'
            maxval = max(PSYrun);
            minval = min(PSYrun);
            subval = mean([maxval;minval]);
            PSYrun = PSYrun - repmat(subval,nt,1);
        case 'z'
            PSYrun = zscore(PSYrun);
        case 'none'
    end
    
    
    % set up:  residual time-series, PPI terms and check task/rest distribution
    for n = 1:nreg;
        phys = PHYSrun(:,n);
        [b,bint,PHYSr] = regress(phys,[PSYrun Xconfrun]);
        PHYSres(:,n) = PHYSr;
        for np = 1:npsy;
            ppi = phys.*PSYrun(:,np);
            PPI(:,n,np)      = ppi;
        end
        
        % compare the distribution
        physt = phys(taskloc); physr = phys(~taskloc);
        mt = mean(physt); mr = mean(physr);
        vt = var(physt); vr = var(physr);
        st = skewness(physt); sr = skewness(physr);
        kt = kurtosis(physt); kr = kurtosis(physr);
        [h,p,ks2stat] = kstest2(physt,physr);
        ppinw.run(run).phystest.kstest(n) = p;
        ppinw.run(run).phystest.mean(1,n) = mt; ppinw.run(run).phystest.mean(2,n) = mr;
        ppinw.run(run).phystest.var(1,n)  = vt; ppinw.run(run).phystest.var(2,n)  = vr;
        ppinw.run(run).phystest.skew(1,n) = st; ppinw.run(run).phystest.skew(2,n) = sr;
        ppinw.run(run).phystest.kurt(1,n) = kt; ppinw.run(run).phystest.kurt(2,n) = kr;
    end
    
    % add to structure
    ppinw.run(run).data.PHYSres = PHYSres;
    ppinw.run(run).data.PHYS    = PHYSrun;
    ppinw.run(run).data.PPI     = PPI;
    ppinw.run(run).data.PSY     = PSYrun;
    
    %% task effects
    disp('  ..task regression')
    for r = 1:nreg
        X = [ones(nt,1) PSYrun Xconfrun];
        y = PHYSrun(:,r);
        con_psy = [zeros(ncon,1) contrast zeros(ncon,size(Xconf,2))];
        
        switch regression
            case 'ar'
                b = fit_gls(y,X,con_psy',armodel);
            case 'standard'
                [b,bint,e] = regress(y,X);
        end
        psycope =  con_psy*b;
        ppinw.run(run).task(:,r) = psycope;
    end
    
    %% times-series connectivity
    disp(['  ..task-unrelated connectivity method: ' connest])
    switch connest
        case 'partial'
            RHO = partialcorr(PHYSres);
            ppinw.run(run).conn.tspar = RHO;
        case 'corr'
            Cr = corr(PHYSres);
            ppinw.run(run).conn.tscorr = Cr;
        case 'glasso'
            connMat = nwa_bic_glasso(PHYSres); 
            ppinw.run(run).conn.tsglasso = connMat;
    end
    
    %% task-dependent connectivity.
    disp('  ..PPI regression')
    progressbar_new('running the PPI regression:')
    nconnec = (nreg*(nreg-1))/2;
    count = 0;
    
    for r1 = 1:(nreg-1)
        for r2 = (r1+1):nreg
            count = count + 1;
            progressbar_new(count/nconnec)
            
            % gPPI approach
            X1 = [ones(nt,1) PHYSrun(:,r1) PSYrun squeeze(PPI(:,r1,:)) Xconfrun];
            y1 = PHYSrun(:,r2);
            X2 = [ones(nt,1) PHYSrun(:,r2) PSYrun squeeze(PPI(:,r2,:)) Xconfrun];
            y2 = PHYSrun(:,r1);
            con_ppi = [zeros(ncon,2) zeros(ncon,npsy) contrast zeros(ncon,size(Xconf,2))];
            
            switch regression
                case 'ar'
                    b1 = fit_gls(y1,X1,con_ppi',armodel,pinv(X1));
                    b2 = fit_gls(y2,X2,con_ppi',armodel,pinv(X2));
                case 'standard'
                    [b1,bint1,e1] = regress(y1,X1);
                    [b2,bint2,e2] = regress(y2,X2);
            end
            
            ppicope1 =  con_ppi*b1;
            ppicope2 =  con_ppi*b2;
            ppinw.run(run).conn.gppi(r1,r2,:) = ppicope1;
            ppinw.run(run).conn.gppi(r2,r1,:) = ppicope2;
            
        end
    end
end

if nrun >1
[ppinw] = nwa_ppi_summary(ppinw,connest);
end
%return


    function [ppinw] = nwa_ppi_summary(ppinw,connest)
        
        connmethod = connest;
        
        nrun = size(ppinw.run,2);
        ncon = size(ppinw.run(1).conn.gppi,3);
        
        % restructure data
        for r = 1:nrun
            task(:,:,r) = ppinw.run(r).task;
            irank(r,:) = ppinw.run(r).irank;
            % ts connectivity
            switch connmethod
                case 'partial'
                    conn(:,:,r) = ppinw.run(run).conn.tspar
                case 'corr'
                    conn(:,:,r) = ppinw.run(run).conn.tscorr;
                case 'glasso'
                    conn(:,:,r) = ppinw.run(run).conn.tsglasso;
            end
            % task ppi
            for c = 1:ncon;
                p = ppinw.run(r).conn.gppi(:,:,c);
                [eCorr] = nwa_comp_connec(p);
                % ppicon(c).dat(:,:,r) = pM;
                % ppinw.diag.gppisym(c,r) = eCorr;
                [pM,d] = nwa_procmat(p,'avesym');
                ppicon(c).dat(:,:,r) = pM;
                ppinw.diag.gppisym(c,r) = d(1);
                ppinw.diag.gppicorr(c,r) = eCorr;
                
                eCorr = nwa_comp_connec(pM,'M2',conn(:,:,r));
                ppinw.diag.task2glasso(c,r) = eCorr;
            end
        end
        % summarize over runs
        ppinw.mean.conn = mean(conn,3);
        ppinw.mean.task = mean(task,3);
        ppinw.mean.irank = mean(irank,1);
        % ppi
        for c = 1:ncon;
            ppinw.mean.gppi(:,:,c) = mean(ppicon(c).dat,3);
        end
    end

end





