function [svmstats] = nwa_classify(NWA,varargin)
% Classification through cross-validation of NWA data.
% =========================================================================
% USE: stats = nwa_classify(NWA,varargin)
% IN:
%
% OUT:
%
% Check/todo:
% - data partioning
% - bootstrapping - overestimation ?
% - save all the model parameters vs. fit all data.. ?
% - X single/multiple models.. 
% =========================================================================

% Defaults
% ----------------------
rep = 500;
nperm = 1000; %00;
%features = {'strength'};
%compare = {'BPD','NPC'};
con = 1;
kfold = 10;
conname = 'subject_contrast';
subsam = false;
Conf = [];
confreg = 1; % 1 = overall regression, 2 = per fold.
repres = 'reg';

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
            case 'rep', rep = varargin{i+1};
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

% eval input: NWA structure or X/Y for testing
% --------------------------------------------
if isfield(NWA,'features')
    % select the data
    selectdat = nwa_selectdata(NWA,'groups',compare,'contrast',con,'features',features);
    svmstats.ftlabels = selectdat.ftlabels;
    svmstats.ftnum = selectdat.ftnum;
    Conf = Conf(selectdat.index',:);
    X = selectdat.X;
    Y = selectdat.Y;
else
    X = NWA;
    Y = varargin{1};
end

% set-up the confound regressors
% ------------------------------
if ~isempty(Conf);
    if size(Conf,1)~=length(Y); Conf = Conf'; end
    xC = [];
    count = 0;
    nn = size(Conf,1);
    % recode into dummy variables
    for c = 1:size(Conf,2);
        creg = Conf(:,c);
        U = unique(creg);
        if length(U)<5;
            warning('recoding into dummy variable');
            D = dummyvar(creg);
            xC(:,1:(length(U)-1)) = D(:,1:(length(U)-1));
        else
            l = size(xC,2);
            xC(:,l+1) = creg;
        end
    end
    Conf = xC;
    svmstats.Conf.xC = xC;
end

if confreg == 1 & ~isempty(Conf);
    % option 1 - regress out all together
    if iscell(X)
        nx = length(X);
        svmstats.Conf.Xunadj = X;
        for n = 1:nx
            Xadj = Correct(X{n},Conf);
            X{n} = Xadj;
        end
    else
        svmstats.Conf.Xunadj = X;
        Xadj = Correct(X,Conf);
        X = Xadj;
    end
end
svmstats.X = X;
svmstats.Y = Y;

%% Run the SVM
% ======================================================
progressbar_new(['running the classification'])
boostdat.pred  = ones(length(Y),rep)*2;
boostdat.model = zeros(length(Y),rep);
for i = 1:rep;
    progressbar_new(i/rep)
    
    switch repres
        
        % regular
        case 'reg'
            Xi = X;
            Yi = Y;
            Conf_i = Conf;
            
            %[stats svmMDl bin bout Ypr_out Ypr_model] = svm(X,Y,kfold,svmfunc); %,Conf);
            %             [stats svmMDl bin bout Ypr_out Ypr_model] = svm(X,Y,kfold,svmfunc); %,Conf);
            %             % bdat(:,i) = pred;
            %             boostdat.pred(:,i) = Ypr_out;
            %             boostdat.model(:,i) = Ypr_model';
            
            % bootstrap
        case 'boot'
            [Xdump,idx] = datasample(X{1},size(X{1},1),'Replace',true);
            for j = 1:length(X)
                Xi{j} = X{j}(idx,:);
            end
            Yi = Y(idx);
            Conf_i = Conf(idx,:);
            
            % boosting - subsampling
        case 'boost'
            f = 0.9;
            bsLoc = randsample(1:length(Y),round(f*length(Y)));
            Yi = Y(bsLoc);
            
            if iscell(X)
                nx = length(X);
                for n = 1:nx
                    Xi{n} = X{n}(bsLoc,:);
                end
            else
                Xbs = X(bsLoc,:);
            end
            Conf_i = Conf(bsloc,:);
            % bdat(:,i) = ones(nn,1)*2;
    end
    
    
    % nested cv for multiple models.
    if confreg == 2 & ~isempty(Conf)
        [stats baccin svmparm Ypr_out Ypr_model modemodel] = svm(Xi,Yi,kfold,svmfunc,Confi);
    else
        [stats baccin svmparm Ypr_out Ypr_model modemodel] = svm(Xi,Yi,kfold,svmfunc);
    end
    
    % save some data
    svmstats.model.BaccIn(:,:,i)  = baccin;
    svmstats.model.predict(:,i) = Ypr_out;
    svmstats.model.models(:,i) = Ypr_model';
    try
    svmstats.model.win(:,i) = modemodel;
    catch
    end
    
    evalm = fields(stats);
    for e = 1:length(evalm);
        edat(e,i) = getfield(stats,evalm{e});
    end
    
    
    try
        for h = 1:length(X)
            avemodel(h).bias(i)    = svmparm(h).bias;
            avemodel(h).scale(i)   = svmparm(h).scale;
            avemodel(h).beta(i,:)  = svmparm(h).beta;
        end
    catch
    end
end

% add the data to structure.
% --------------------------
evalm = fields(stats);
for e = 1:size(edat,1);
    d = edat(e,:);
    svmstats = setfield(svmstats,evalm{e},'dat',d);
    svmstats = setfield(svmstats,evalm{e},'mean',nanmean(d));
    svmstats = setfield(svmstats,evalm{e},'CI95',prctile(d,[2.5 97.5]));
end


%% Permutation null model
% --------------------------------
% only run when its worth checking (e.g. svmstats.bacc.mean>0.6)
ndat = length(Y);
if nperm > 0 & svmstats.bacc.mean>0.4;
    progressbar_new(['running the SVM permutation'])
    for p = 1:nperm;
        progressbar_new(p/nperm)
        
        % permutation test
        r = randperm(ndat);
        Yperm = Y(r);
        [stats svmMDl] = svm(X,Yperm,kfold,svmfunc); % this needs to be updated for nested cv
        
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
    
else
    % add the permutation data
    evalm = fields(stats);
    for e = 1:length(evalm);
        p = 1;
        svmstats = setfield(svmstats,evalm{e},'pval',p);
    end
end

%% run the subsamples
% ---------------------------------
if subsam %& svmstats.bacc.pval<0.05
    subsam_n = 2;
    subsam_frac = [0.5:0.1:0.6];
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
            
            %Xss = X(locss,:);
            
            for j = 1:length(X)
                Xss{j} = X{j}(locss,:);
            end
            
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
% ---------------
try
    for h = 1:length(X)
        svmstats.avemodel(h).bias  = mean(avemodel(h).bias);
        svmstats.avemodel(h).scale = mean(avemodel(h).scale);
        svmstats.avemodel(h).beta  = mean(avemodel(h).beta);
    end
    %     % sanity check - fitting the ensemble model to the entire data.
    %     Xz = zscore(X);
    %     Y_pred_raw = (Xz/svmstats.avemodel.scale)*svmstats.avemodel.beta'+svmstats.avemodel.bias;
    %     pos = Y_pred_raw>0;
    %     Y_test = zeros(length(Y_pred_raw),1);
    %     Y_test(pos) = 1;
    %     svmstats.avemodel.acc = 1-mean(abs(Y-Y_test));
catch
end

% compare against model fit of all data 

% % try
%     for h = 1:length(X)
%         
%         % train the model
%         Mdlfit = svmfunc(X{h},Y);
%         
% %         % save model parameters (only available for linear models)
% %         %                         try
% %         svmpar(nx).beta(kin,:)= Mdlfit.Beta;
% %         svmpar(nx).bias(kin)  = Mdl_in.Bias;
% %         svmpar(nx).scale(kin) = Mdl_in.KernelParameters.Scale;
%         
%         svmstats.modelfit(h).bias  = Mdlfit.Bias;
%         svmstats.modelfit(h).scale = Mdlfit.KernelParameters.Scale;
%         svmstats.modelfit(h).beta  = Mdlfit.Beta;
%         fit = predict(Mdlfit,X{h});
%         M = metrics(Y,fit);
%         svmstats.modelfit(h).bacc = M.bacc;
%     end
    
% catch
% end




%% print basic stats
% ---------------------
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

% main classify function
    function [stats baccin svmparm Ypr_out Ypr_model modemodel] = svm(X,Y,kfold,svmfunc,varargin);
        svmparm = [];
        baccin = [];
        modemodel = [];
        
        if ~isempty(varargin); Ccv = varargin{1}; end
        
        % kfold cross validation
        % ----------------------
        
        % Outer loop
        % ------------------------------
        ny  = length(Y);
        CVPol = cvpartition(Y,'kfold',kfold);
        Ypr_out = zeros(ny,1);
        for kout = 1:kfold;
            kout_loc = CVPol.test(kout);
            
            % split the data
            Ytrain_out = Y(~kout_loc,:);
            Ytest_out  = Y(kout_loc,:);
            
            if ~isempty(varargin);
                Ctrain_out = Ccv(~kout_loc);
                Ctest_out  = Ccv(kout_loc);
            end
            
            if iscell(X)
                
                % Inner loop - loop over possible feature models.
                % -----------------------------------------------
                CVPil = cvpartition(Ytrain_out,'kfold',5);
                for nx = 1:length(X)
                    Xn = X{nx}; Xn = Xn(~kout_loc,:);
                    
                    Ypr_in = zeros(length(Ytrain_out),1);
                    for kin = 1:5; % max 5 for computational reasons?
                        kin_loc = CVPil.test(kin);
                        
                        Xtrain_in    = Xn(~kin_loc,:);
                        Xtest_in     = Xn(kin_loc,:);
                        if ~isempty(varargin);
                            Ctrain_in  = Ctrain_out(~kin_loc);
                            Ctest_in  = Ctrain_out(kin_loc);
                            Xtrain_in = Correct(Xtrain_in,Ctrain_in);
                            Xtest_in  = Correct(Xtest_in,Ctest_in);
                        end
                        
                        Ytrain_in    = Ytrain_out(~kin_loc,:);
                        Ytest_in     = Ytrain_out(kin_loc,:);
                        
                        % train the model
                        Mdl_in = svmfunc(Xtrain_in,Ytrain_in);
                        
                        % save model parameters (only available for linear models)
                        %                         try
                        svmpar(nx).beta(kin,:)= Mdl_in.Beta;
                        svmpar(nx).bias(kin)  = Mdl_in.Bias;
                        svmpar(nx).scale(kin) = Mdl_in.KernelParameters.Scale;
                        %                         catch
                        %                         end
                        
                        % predict
                        pred_in = predict(Mdl_in,Xtest_in);
                        Ypr_in(kin_loc) = pred_in;
                        
                    end
                    M_in = metrics(Ytrain_out,Ypr_in);
                    model_bacc_in(nx) = M_in.bacc;
                end
                
                % pick the winning model.
                [v wmodel] = max(model_bacc_in); % winning model
                Xtrain_out = X{wmodel}(~kout_loc,:);
                Xtest_out  = X{wmodel}(kout_loc,:);
                
                % save the data for all models.
                baccin(kout,:) = model_bacc_in;
                model(kout) = wmodel;
            else
                
                Xtrain_out = X(~kout_loc,:);
                Xtest_out  = X(kout_loc,:);
                
                % correct
                if ~isempty(varargin);
                    Xtrain_out = Correct(Xtrain_out,Ctrain_out);
                    Xtest_out  = Correct(Xtest_out,Ctest_out);
                end
                
            end

            % train the winning model
            Mdl_out = svmfunc(Xtrain_out,Ytrain_out);
            
            % prediction
            pred_out = predict(Mdl_out,Xtest_out);
            Ypr_out(kout_loc)=pred_out;
            
            try
            Ypr_model(kout_loc)=wmodel;
            catch
            Ypr_model(kout_loc)=1;
            end
            
        end
        stats = metrics(Y,Ypr_out);
        % bacc_out = M_out.bacc;
        if iscell(X); modemodel = mode(model); end
        
        % model summary
        try
            for j = 1:length(X)
                svmparm(j).beta  = mean(svmpar(j).beta);
                svmparm(j).bias  = mean(svmpar(j).bias);
                svmparm(j).scale = mean(svmpar(j).scale);
            end
            
        catch
        end
    end

    function M = metrics(Y,Yp)
        % test evaluation features
        TP  = sum(Yp==1 & Y==1);
        FP  = sum(Yp==1 & Y==0);
        TN  = sum(Yp==0 & Y==0);
        FN  = sum(Yp==0 & Y==1);
        tot = TP+FP+TN+FN;
        
        M.sens = TP/(TP+FN);
        M.spec = TN/(TN+FP);
        M.prec = TP/(TP+FP);
        M.F1   = 2*TP/(2*TP+FP+FN);
        M.acc  = (TP+TN)/tot;
        M.bacc = 0.5*(M.sens+M.spec);
    end

    function Xc = Correct(X,Conf)
        %  regress out confounds
        Conf1 = [ones(size(Conf,1),1) Conf];
        bcf = pinv(Conf1)*X;
        Xc = X - Conf1*bcf;
    end

end