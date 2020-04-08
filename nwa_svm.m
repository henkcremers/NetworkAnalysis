function [svmstats] = nwa_svm(NWA,varargin)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Defaults
% ===============================================
niter = 500;
nperm = 1000; %00;
features = {'strength'};
compare = {'BPD','NPC'};
con = [1];
kfold = 10;
conname = 'subject_contrast';
subsam = false;
Conf = [];

% the deafult svm function and input arguments
% --------------------------------------
svmfunc = @(X,Y)fitcsvm(X,Y,...
    'Standardize',true,...
    'KernelFunction','linear',...
    'KernelScale','auto',...
    'prior','uniform',...
    'OutlierFraction',0.05);

%     'BoxConstraint',1,...

%% Input
% ===============================================
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'features', features = varargin{i+1};
            case 'niter', niter = varargin{i+1};
            case 'nperm', nperm = varargin{i+1};
            case 'compare', compare = varargin{i+1};
            case 'con', con = varargin{i+1};
            case 'conname', conname = varargin{i+1};
            case 'svmfunc', svmfunc = varargin{i+1};
            case 'kfold', kfold = varargin{i+1};
            case 'subsam', subsam = true;
            case 'Conf',Conf = varargin{i+1};
                % otherwise error('dont recognize the input');
        end
    end
end
%svmstats.feature 
% features

%% eval input: NWA structure or X/Y for testing
if isfield(NWA,'data')
    %% select the data
    selectdat = nwa_selectdata(NWA,'groups',compare,'contrast',con,'features',features);
    svmstats.ftlabels = selectdat.ftlabels;
    svmstats.ftnum = selectdat.ftnum;
    X = selectdat.X;
    Y = selectdat.Y;
    Conf = Conf(selectdat.index',:);
    
else
    X = NWA;
    Y = varargin{1};
end

nn = size(X,1);
nf = size(X,2);

%% adjust X for potential confounders
% xConf = NWA.site.num;
if ~isempty(Conf);
    if size(Conf,1)~=length(Y); Conf = Conf'; end
%     % xConf1 = xConf(g1loc,:); xConf2 = xConf(g2loc,:); xConf = [xConf1; xConf2];
%     try
%     
%     catch
%     end
    xC = [];
    count = 0;
    
    % recode into dummy variables
    for c = 1:size(Conf,2);
        creg = Conf(:,c);
        U = unique(creg);
        if length(U)<5;
            warning('recoding into dummy variable');
            for nU = 1:(length(U)-1);
                count = count +1;
                dum = zeros(nn,1);
                dum(creg==U(nU),1) = 1;
                xC(:,count) = dum;
            end
        else
            l = size(xC,2);
            xC(:,l+1) = creg;
        end
    end
    svmstats.Conf.xC = xC;
    svmstats.Conf.Xunadj = X;
    % this can be done more easy...!! lukas/steven.
    for j = 1:nf
        s = regstats(X(:,j),xC);
        H(1,j) = s.fstat;
        Xadj(:,j) = s.r;
    end
    svmstats.Conf.effect = H;
    X = Xadj;
end
svmstats.X = X;
svmstats.Y = Y;


%% Run the SVM
% ======================================================
progressbar_new(['running the SVM'])
for i = 1:niter;
    progressbar_new(i/niter)
    [stats svmMDl boostpred] = svm(X,Y,kfold,svmfunc);
    
    evalm = fields(stats);
    for e = 1:length(evalm);
        edat(e,i) = getfield(stats,evalm{e});
    end
    
    % boosting 
    bdat(i,:) = boostpred;
    
    try
        ensmodel.bias(i)    = svmMdl.bias;
        ensmodel.scale(i)   = svmMdl.scale;
        ensmodel.beta(i,:)  = svmMdl.beta;
    catch
    end
end

%% add the data to structure.
evalm = fields(stats);
for e = 1:size(edat,1);
    d = edat(e,:);
    svmstats = setfield(svmstats,evalm{e},'dat',d);
    svmstats = setfield(svmstats,evalm{e},'mean',nanmean(d));
    svmstats = setfield(svmstats,evalm{e},'CI95',prctile(d,[2.5 97.5]));
end

% boosting
% prob = sum(bdat==1)./niter;
% majority = double((prob>0.5))';
% boostacc = sum(majority==svmstats.Y)/nn;

% svmstats.boost.dat  = bdat;
% svmstats.boost.prob = prob';
% svmstats.boost.majority = majority;
% svmstats.boost.acc = boostacc;


%% Permutation "null model

% only run when its worth checking (e.g. svmstats.bacc.mean>0.6)
if nperm > 0 & svmstats.bacc.mean>0.6;
    progressbar_new(['running the SVM permutation'])
    for p = 1:nperm;
        progressbar_new(p/nperm)
        
        % permutation test
        r = randperm(length(Y));
        Yperm = Y(r);
        [stats svmMDl] = svm(X,Yperm,kfold,svmfunc);
        
        evalm = fields(stats);
        for e = 1:length(evalm);
            pdat(e,p) = getfield(stats,evalm{e});
        end
    end
    
    % add the permutation data
    for e = 1:size(pdat,1);
        pd = pdat(e,:);
        svmstats = setfield(svmstats,evalm{e},'permdat',pd);
        m = getfield(svmstats,evalm{e},'mean');
        p = sum(pd>m)/nperm;
        if p == 0; p = 1*10^-(log10(nperm)); end
        svmstats = setfield(svmstats,evalm{e},'pval',p);
    end
    
    %% effect size ??
    %     m1 = svmstats.acc.mean;
    %     m2 = mean(svmstats.acc.permdat);
    %     sd1 = std(svmstats.acc.dat);
    %     sd2 = std(svmstats.acc.permdat);
    %     svmstats.acc.ES = (m1-m2)/sqrt((sd1^2+sd2^2)/2);
else  
    % add the permutation data
    evalm = fields(stats);
    for e = 1:length(evalm);
        p = 1;
        svmstats = setfield(svmstats,evalm{e},'pval',p);
    end
end

%% run the subsamples
if subsam & svmstats.bacc.pval<0.05
    subsam_n = 10;
    subsam_frac = [0.5:0.1:1];
    count = 0;
    progressbar_new(['running the SVM subsampling'])
    for sl = 1:length(subsam_frac);
        for s = 1:subsam_n
            count = count +1;
            progressbar_new(count/(length(subsam_frac)*subsam_n))
            y1 = find(Y==1);
            y0 = find(Y==0);
            
            y1ss = randsample(y1,round(subsam_frac(sl)*length(y1)));
            y0ss = randsample(y0,round(subsam_frac(sl)*length(y0)));
            locss = [y1ss; y0ss];
            Yss = Y(locss);
            Xss = X(locss,:);
            
            [stats svmMDl] = svm(Xss,Yss,kfold,svmfunc);
            subsambacc(s,sl) = stats.bacc;
            
            evalm = fields(stats);
            for e = 1:length(evalm);
                sdat(e).d(s,sl) = getfield(stats,evalm{e});
            end
            
        end
    end
    
    % add the data.
    evalm = fields(stats);
    svmstats.subsam.frac = subsam_frac;
    for e = 1:length(sdat);
        d = sdat(e).d;
        svmstats = setfield(svmstats,'subsam',evalm{e},'mean',mean(d));
        svmstats = setfield(svmstats,'subsam',evalm{e},'CI95',prctile(d,[2.5 97.5]));
    end
end

%% ensemble model
try
    svmstats.ensmodel.bias  = mean(ensmodel.bias);
    svmstats.ensmodel.scale = mean(ensmodel.scale);
    svmstats.ensmodel.beta  = mean(ensmodel.beta);
    
    % sanity check - fitting the ensemble model to the entire data. 
    Xz = zscore(X);
    Y_pred_raw = (Xz/svmstats.ensmodel.scale)*svmstats.ensmodel.beta'+svmstats.ensmodel.bias;
    pos = Y_pred_raw>0;
    Y_test = zeros(length(Y_pred_raw),1);
    Y_test(pos) = 1;
    svmstats.ensmodel.acc = 1-mean(abs(Y-Y_test));
catch
end

%% print basic stats
for e = 1:size(edat,1);
    d = getfield(svmstats,evalm{e});
    prM   = num2str(round(d.mean,2));
    prCIL = num2str(round(d.CI95(1),2));
    prCIH = num2str(round(d.CI95(2),2));
    if isfield(d,'pval');
        prP   = num2str(round(d.pval,3));
    else
        prP = 1;
    end
    disp([evalm{e} ' mean: ' prM ' (CI95: ' prCIL  ' - ' prCIH ') pval: ' prP ]);
end

% main svm function
    function [stats svmMdl boostpred] = svm(X,Y,kfold,svmfunc);
        svmMdl = [];
        
        % kfold cross validation
        ny  = length(Y);
        CVP = cvpartition(Y,'kfold',kfold);
        boostpred = zeros(ny,1);
        for k = 1:kfold;
            kloc = CVP.test(k);
            
            % split the data
            Xtrain = X(~kloc,:);
            Xtest  = X(kloc,:);
            Ytrain = Y(~kloc,:);
            Ytest  = Y(kloc,:);
            
            % train the model
            Mdl = svmfunc(Xtrain,Ytrain);
            
            % save model parameters (only available for linear model)
            try
                svmMdl.beta(k,:)= Mdl.Beta;
                svmMdl.bias(k)  = Mdl.Bias;
                svmMdl.scale(k) = Mdl.KernelParameters.Scale;
            catch
            end
            
            % prediction
            Ypred = predict(Mdl,Xtest);
            boostpred(kloc)=Ypred;
            %acc1(k) = mean(Ytest==Ypred);
            
            % test evaluation features
            % --------------------------
            TP  = sum(Ypred==1 & Ytest==1);
            FP  = sum(Ypred==1 & Ytest==0);
            TN  = sum(Ypred==0 & Ytest==0);
            FN  = sum(Ypred==0 & Ytest==1);
            tot = TP + FP + TN + FN;
            
            sens(k) = TP/(TP + FN);
            spec(k) = TN/(TN + FP);
            prec(k) = TP/(TP+FP);
            F1(k)   = 2*TP/(2*TP+FP+FN);
            acc(k)  = (TP+TN)/tot;
            bacc(k) = 0.5*(sens(k)+spec(k));
        end
        stats.acc  = mean(acc);
        stats.spec = mean(spec);
        stats.sens = mean(sens);
        stats.prec = mean(prec);
        stats.F1   = mean(F1);
        stats.bacc = mean(bacc);
        % stats.bpred = boostpred;
        
        % model summary
        try
            svmMdl.beta  = mean(svmMdl.beta);
            svmMdl.bias  = mean(svmMdl.bias);
            svmMdl.scale = mean(svmMdl.scale);
        catch
        end
    end

end