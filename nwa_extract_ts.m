function [tsMat ddata] = nwa_extract_ts(image,roi,mask,varargin)
% Extract time-series data from an (preprocessed) image. 
% USE  [tsMat] = nwa_extract_ts(image,roi,mask,varargin)
% works with 3d and 4d atlas and ROI data
% =========================================================================
% IN:   
%  image - a preprocessed nifto image. 
%  roi - a brain region (binarized) or atlas structure (3D or 4D)
%  mask - a mask image to read in the image data (this can be the same as
%  the image)
% OUT: 
%  tsdat - time series data per region in the roi
%  ddata - some diagnostic data
% =========================================================================

diagnose = 0;
gmdata = [];
ddata = [];

%% read the images. 
disp( '.. reading the data')

[maskInfo, dat] = iimg_read_img(mask, 2);

% load the roi/atlas data 
roidata = iimg_get_data(maskInfo,roi);

for i = 1:length(varargin)
  arg = varargin{i};
  if ischar(arg)
      switch arg
         case 'diagnose', diagnose = 1;  
         case 'gm', gm = varargin{i+1}; gmdata  = iimg_get_data(maskInfo,gm);  
       end
   end
 end


% determine the type of image 3d/4d
ndim = size(roidata,1);
if ndim > 1;  
   nreg = ndim;
   is4d = logical(1);
else
   udata = unique(roidata); 
   udata = udata(udata~=0); % filter out zeros
   nreg = length(udata);
   is4d = logical(0);
end

% load the image data. 
imdata  = iimg_get_data(maskInfo,image);
nt = size(imdata,1); 

% check the size
if size(roidata,2) ~= size(imdata,2)
    error('the number of voxels in the roi and image should match')
end

% loop over all regions in the roi;
tsMat = zeros(nt,nreg);   
data = zeros(nreg,2);

progressbar_new('time-series extration:')
count = 0;
for n = 1:nreg 
    progressbar_new(n/nreg)
    %disp(n)
    if is4d
    loc = roidata(n,:)>0;
    regData = imdata(:,loc);      
    else     
    loc = roidata==udata(n);
    regData = imdata(:,loc);
    end
    
    if ~isempty(gmdata)
        
        regGM   = gmdata(1,loc);
        
        % clean up the data
        
        nanvec = sum(isnan(regData),1)>0;
        noGm   = regGM==0;
        datuse = (nanvec~=1 & noGm~=1); % usable data
        
        if sum(datuse)==0;            
            warning('no valid data points')
        else
            count = count+1;
            regData = regData(:,datuse);
            regGM   = regGM(1,datuse);
            
            % weight the ts data according to the Gm value
            w = regGM;
            nt = size(regData,1);
            wmat = repmat(w,nt,1);
            regData = regData .*wmat;
            
            % mean
            tsm = mean(regData,2);   
        end  
    else
        tsm = nanmean(regData,2);
    end
    tsMat(:,n) = tsm;
    
    if diagnose == 1;
        % PCA
        [COEFF, SCORE, LATENT] = pca(regData);
        tsPC1 = SCORE(:,1);
        expVar = LATENT(1)/(sum(LATENT));
        
        % compare the mean and PC1;
        c = corr(tsm,tsPC1);
        ddata(n,1) = expVar;
        ddata(n,2) = c;
%        ddata(n,3) = sum(noGm)/sum(loc);
    end
end

progressbar_new(1)
% for non-time serie data return all voxel values
if nreg == 1; tsMat = regData; end

end

