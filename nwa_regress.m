function [] = nwa_regress(image,mask,varargin)
%USE  [matlabbatch] = nwa_regress(image,mask,varargin)
% =========================================================================
% IN:   Detailed explanation goes here
% OUT: 
% ERAMPLE: 

%% read the images. 

method = 'spm';
hp = 128;
fileid = 'res4D';
disp( '.. regression')

% get the user input 
%------------------------------------------------------------
 for i = 1:length(varargin)
  arg = varargin{i};
  if ischar(arg)
      switch arg
         case 'method', method = varargin{i+1};   
         case 'R', R = varargin{i+1};    % can be a matrix or path to *.mat file for spm analysis
         case 'hp', hp = varargin{i+1};   
         case 'fileid', fileid = varargin{i+1};    
       end
   end
 end

disp( '.. reading the data')
[maskInfo, dat] = iimg_read_img(mask, 2);
% load the image data. 
imdata  = iimg_get_data(maskInfo,image);
nt   = size(imdata,1);
nvox = size(imdata,2); 

[path file ext] = fileparts(image);

respath = [path '/statsResiduals'];
if ~isdir(respath)
    mkdir(respath)
end

cd(respath);

switch method
    case 'basic'
        
        % create new matrix with residuals;
        regDat = zeros(nt,nvox);
        
        % add filtering to the model? 
        
        for n = 1:nvox
            %     loc = roidata(n,:)>0;
            %
            %     if sum(loc)==0
            %         warning('some of the frames of the roi seem empty..')
            %     end
            y = imdata(:,n);
            [b,r] = y_regress_ss(y,R);
            regDat(:,n) = r;
        end
        
        % add constant to the residual image (see FSL mailinglist on this).
        addCon = 1000;
        regDat = regDat + addCon;

        % write new residual image.

        resim = ([respath '/' imname fileid]);
        disp('..write residual image');
        iimg_reconstruct_vols(regDat', maskInfo, 'outname',resim);

        
    case 'spm'

        %extend the file name for spm.
        files = {};
        for j = 1:nt;
            files{j,1} = [image ',' num2str(j)];
        end
        newfile = [fileid '_hp' num2str(hp) '.nii'];
        
        matlabbatch{1}.spm.stats.fmri_spec.dir = {respath};
        matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;
        matlabbatch{1}.spm.stats.fmri_spec.sess.scans = files;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
        matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {R};
        matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = hp;
        matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
        matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
        matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
        matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
        matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
        matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
        matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 1;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        matlabbatch{3}.spm.util.cat.vols(1) = cfg_dep('Model estimation: Residual Images', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','res'));
        matlabbatch{3}.spm.util.cat.name = newfile;
        matlabbatch{3}.spm.util.cat.dtype = 4;
        
        %run 
        save([newfile '_hp' num2str(hp) '_batch.mat'],'matlabbatch');
        spm('defaults', 'FMRI');
        spm_jobman('initcfg'); %this is recomendent, gives warning otherwise
        spm_jobman('run',matlabbatch)
        clear matlabbatch
        
        % add 1000 to the image, in voxel only
        imdata  = iimg_get_data(maskInfo,newfile);
        imdata = imdata + 1000;
        iimg_reconstruct_vols(imdata', maskInfo, 'outname',newfile);
        
        % remove Residuals
        disp('cleanup')
        delete('Res_*')
        delete('beta*')
        delete('ResMS*');
        delete('RPV*');
        
        % move SPM.mat file 
        newSPMdir = [newfile '_hp' num2str(hp) '_SPMdir'];
        mkdir(newSPMdir);
        movefile('SPM.mat',newSPMdir);
 
end

end

