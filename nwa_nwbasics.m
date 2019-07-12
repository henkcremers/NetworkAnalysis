function [nwtable] = nwa_nwbasics(NWA)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ngroup = length(unique(NWA.group.num));

%% basic network metrics
clus  = NWA.global.clus{1};
mod   = NWA.global.mod{1};
Meff  = mean(NWA.local.eff{1},2);
Meig  = mean(NWA.local.eig{1},2);
Mst   = mean(NWA.local.strength{1},2);

nwdat = [clus mod Meff Meig Mst];
for j = 1:ngroup
    nwtable(j,:) = mean(nwdat(NWA.group.num==j,:));
end

% nwtable = nwtable';
try
nwtable = array2table(nwtable,...
    'VariableNames',{'ClusCoef' 'Modularity' 'glEff' 'mEig' 'mStrength'});
catch
end

% 'RowNames',
% nwtable = table(nwtable,'VariableNames',{'NPC','CLC','BPD'})

%% distribution of strength
nwa_colors;

%figure;
st = NWA.local.strength{1};
%range = linspace(min(min(st)),max(max(st)),50);
range = linspace(0,max(max(st)),20);
disp(max(max(st)))
for i = 1:size(st,1);
sthist(i,:) = histcounts(st(i,:),range);
end
xvals = range(2:end);
plot(xvals,sthist','Color',colors{6,3});

% overlay averages per group
hold on;
for j = 1:ngroup;
    d= mean(sthist(NWA.group.num==j,:));
    p(j) = plot(xvals,d,'Color',colors{j,1});
    g = NWA.group.name(NWA.group.num==j);
    label = g{1};
    if isnumeric(label); label = num2str(label); end
    l{j} = label;
end
legend(p,l)
return
% xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
% title('Degree distribution')
