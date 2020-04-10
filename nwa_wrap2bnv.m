function [] = nwa_wrap2bnv(C,varargin)
% Wrapper function get a network matric into the format of the Brain
% Network Viewer (BNV).
% USE: nwa_wrap2bnv(C,varargin)
% =========================================================================
% =========================================================================

nn = size(C,1);

% defaults:
%-------------------
dispmode = 'plotnw';
M = ones(nn,1); 
S = ones(nn,1); 
atinfo = load('/Volumes/WD2T/BPD/BPD_info/networkQC/20180824/NCUT121_BPD.mat');
%atinfo = load('/Volumes/WD2T/HCP_PTN820/NWA/hcpica300.mat');
nwfile = 'nwgraph';
nodefile = 'node';
edgefile = 'edge';
thr = 0; % threshold 
% template = '/Volumes/WD2T/BCS_UofC/SET_MFX/NWA/ncutTcorr22_glassoLabda0.01_20150630/networktemplate3.mat';
template = '/Volumes/WD2T/BPD/BPD_info/Plots/bnv_settings.mat';
% template = '/Volumes/WD2T/Tools/CCP/NWA/Checks/template_bnv_color.mat';
% surftemplate 
surf = '/Volumes/WD2T/Tools/BrainNetViewer_20150807/Data/SurfTemplate/BrainMesh_Ch2.nv';


% get the user imput 
%--------------------
 for i = 1:length(varargin)
  arg = varargin{i};
  if ischar(arg)
      switch (arg)
         case 'dispmode', dispmode = varargin{i+1};
         case 'M', M = varargin{i+1};
         case 'S', S = varargin{i+1};  
         case 'nwfile', nwfile = varargin{i+1};  
         case 'node', nodefile = varargin{i+1}; 
         case 'edge', edgefile = varargin{i+1}; 
         case 'atlas', atinfo = varargin{i+1};   
         case 'surf',  surf =   varargin{i+1};   
         case 'thr',  thr =   varargin{i+1}; 
         case 'template',  template =   varargin{i+1};    
       end
   end
 end


% node list. 
% -----------------------
coord = num2cell(round(atinfo.DataList(:,1:3)),3); 
inf = num2cell([M S],3);
labels = atinfo.RegionList(:,2);
nodedata = [coord inf labels];
nodefile = [nodefile '.node'];
dlmcell(nodefile,nodedata,' ') 


% edge list 
% ----------------------
C(abs(C)<thr) = 0;
spar = sum(sum(abs(C)>thr))*100/(nn^2);
disp(['The network has ' num2str(spar) ' % edges'])
C(logical(eye(nn,nn))) = - 1;
edgename = [edgefile '.edge'];
dlmwrite(edgename,C,'delimiter','\t');

% switch to display mode
% ----------------------
switch dispmode
    case 'plotnw'       
    BrainNet_MapCfg(surf,nodefile,edgename,template)
    %BrainNet_MapCfg(surf,nodefile,edgename)
    
    case 'setupnw'
    copyfile(surf,'.')
    BrainNet

    case 'savenw'
    BrainNet_MapCfg(surf,nodefile,edgename,[nwfile '.eps'],template)   
end

return




