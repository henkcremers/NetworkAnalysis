function [stats] = nwa_ttest(NWA,varargin)


%% Defaults
% ===============================================
features = {'strength'};
compare = {1,2};
con = [1];
conname = 'subject_contrast';
Conf = [];
alpha = 0.05;
nroi =0;
makeim=0;

%% Input
% ===============================================
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'features', features = varargin{i+1};
            case 'compare', compare = varargin{i+1};
            case 'con', con = varargin{i+1};
            case 'conname', conname = varargin{i+1};
            case 'Conf',Conf = varargin{i+1};
            case 'alpha',alpha = varargin{i+1};
                % otherwise error('dont recognize the input');
        end
    end
end
%tstats.feature
features
%% select the data
selectdat = nwa_selectdata(NWA,'groups',compare,'contrast',con,'features',features);
tstats.ftlabels = selectdat.ftlabels;
tstats.ftnum = selectdat.ftnum;
X = selectdat.X;
Y = selectdat.Y;
nn = size(X,1);
nf = size(X,2);

%% adjust X for potential confounders
% xConf = NWA.site.num;
if ~isempty(Conf);
    if size(Conf,1)~=size(NWA.local.edges{1},1); Conf = Conf'; end
    % xConf1 = xConf(g1loc,:); xConf2 = xConf(g2loc,:); xConf = [xConf1; xConf2];
    Conf = Conf(selectdat.index',:);
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
    tstats.Conf.xC = xC;
    tstats.Conf.Xunadj = X;
    for j = 1:nf
        s = regstats(X(:,j),xC);
        H(1,j) = s.fstat;
        Xadj(:,j) = s.r;
    end
    tstats.Conf.effect = H;
    X = Xadj;
end
tstats.X = X;
tstats.Y = Y;

%% T-test

g1=X(Y==1,:);
g2=X(Y==0,:);

%% ttest
[H P CI STATS] = ttest2(g1,g2,'alpha',alpha);
T = STATS.tstat;
stats.T = T;
%     [p_fdr, p_masked] = fdr(P, alpha);

plotdat(:,1) = mean(g1);
plotdat(:,2) = mean(g2);

% else
%     g = find(gcon~=0);
%     gdat = [];
%     for ng = 1:length(g);
%         gdat = [gdat;cdat(group==g(ng),:)];
%     end
%     [H P CI STATS] = ttest(gdat);
%     T = STATS.tstat;
%     stats.T = T;
%     H = P<alpha;
%     plotdat(:,1) = mean(gdat);
%
% end




%% plot
figure;
% plotdat(:,1) = mean(g1);
% plotdat(:,2) = mean(g2);
p = plot(plotdat);
% legend({'NPC','BPD'})
%set(gca,'XTick',1:nreg)
%set(gca,'XTickLabel',NWA.atlas.RegionList(:,1))
%xticklabel_rotate([],90,[],'Fontsize',8)
hold on;
sdat(1,:) = find(H>0);
sdat(2,:) = ones(1,sum(H))*0.01;
scatter(sdat(1,:),sdat(2,:),'r')
% % ROI analyses
rcolors = {'k','g','y','b','k','g','y','b'};
if nroi>0
    for r = 1:nroi
        rdat(1,:) = roiloc{r};
        rdat(2,:) = ones(1,length(roiloc{r}))*(0.1*r);
        scatter(rdat(1,:),rdat(2,:),rcolors{r})
        clear rdat
        %         txt = '\leftarrow sin(\pi) = 0';
        %         text(pi,sin(pi),txt)
    end
end

%% print significant results
if sum(H)>0;
    loc = find(H>0);
    for h = 1:sum(H);
        t = STATS.tstat(loc(h));
        p = P(loc(h));
        reg = num2str(loc(h));
        if isfield(NWA,'atlas')
            reg = ['Region: ' reg ' : '  NWA.atlas.RegionList{loc(h)}];
        end
        %        con = num2str(nc);
        disp(['region: ' reg ' T= ' num2str(t) ' p= ' num2str(p)])
        G1m = mean(g1(:,loc(h)));
        G2m = mean(g2(:,loc(h)));
        disptext = [' .. G1: ' num2str(G1m) ' G2: ' num2str(G2m) ];
        disp(disptext)
    end
end

%% create image of data
if makeim == 1;
    [MaskInfo, dat] = iimg_read_img(atlas, 2);
    Data = iimg_get_data(MaskInfo, atlas);
    imdata = zeros(1,length(Data));
    Udata = unique(Data);
    Udata = Udata(Udata~=0);
    for j = 1:nreg
        locclus = (Data==Udata(j));
        imdata(1,locclus) = T(j);
    end
    imfile = [met '.nii'];
    iimg_reconstruct_vols(imdata', MaskInfo, 'outname',imfile);
end


%% ROI analyses
if nroi>0
    figure
    for r = 1:nroi
        %  group differences:
        g1rdat = mean(cdat(group==g1,roiloc{r}),2);
        g2rdat = mean(cdat(group==g2,roiloc{r}),2);
        [H P CI STATS] = ttest2(g1rdat,g2rdat,'alpha',alpha);
        T = STATS.tstat;
        for g = 1:ng
            violintdat{g} = mean(cdat(group==g,roiloc{r}),2);
        end
        subplot(nroi,1,r)
        violin(violintdat,'xlabel',{'NPC','CLC','BPD'});
        roititle = [roiname{r} ' T=' num2str(STATS.tstat) ' p=' num2str(P)];
        title(roititle)
    end
end

return