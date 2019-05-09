function [selectdat] = nwa_selectdata(NWA,varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% defaults 
gcompare = {'BPD','NPC'};
features = {'strength'};
Xt = [];
Y = [];
con = [1];

for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch arg
            case 'groups', gcompare = varargin{i+1};
            case 'features', features = varargin{i+1};
            case 'contrast', con = varargin{i+1};
            case 'conname', conname = varargin{i+1};
%            otherwise error('dont recognize the input');
        end
    end
end


%% get all the variables/feature
ftlabels = {};
for m = 1:length(features); 
    try
        dat = getfield(NWA.local,features{m});
        
        % this is ugly, but it works: apply the "low"-level contrast
        % number of dimensions
        nd = size(dat,2);
        
        % zero pad contrast if needed
        nc = length(con);  if nd>nc; con = [con zeros(1,nd-nc)]; end
        
        % apply the contrast
        for n = 1:nd; tempdat(:,:,n) = dat{n}.*con(n); end
        dat = sum(tempdat,3); clear tempdat
       
        % update the feature data
        Xt = [Xt dat];
        
        % variable names and numbers
        sloc = length(ftlabels);
        for j = 1:length(dat);
        ftlabels{sloc+j} = [features{m} '_' num2str(j)];
        ftnum(sloc+j)  = j;
        end
    catch
        error([features{m} ' is not a field '])
    end
end
selectdat.ftlabels = ftlabels;
selectdat.ftnum   = ftnum;

%% select the subjects
g1 = gcompare{1};
if isstr(g1)
    gname = NWA.group.name;
    g1loc = find(strcmp(gname,g1));
    Xg1 = Xt(g1loc,:);

    g2loc = [];
    for j = 2:length(gcompare);
        g2loc = [g2loc find(strcmp(gname, gcompare{j}))];
    end
    Xg2 = Xt(g2loc,:);
else
    gnum = NWA.group.num;
    g1loc = find(gnum==g1);
    Xg1 = svmdat(g1loc,:);

    g2loc = [];
    for j = 2:length(gcompare);
        g2loc = [g2loc find(gnum==compare{j})];
    end
    Xg2 = Xt(g2loc,:);
end

%% combine the data
X = [Xg1; Xg2];
Y = ones(size(X,1),1);
Y(size(Xg1,1)+1:end) = 0;
selectdat.X = X;
selectdat.Y = Y;
selectdat.index = [g1loc g2loc];

end

