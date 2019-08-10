function [fem,mesh,segm] = segm_6C(pathout,path_t1,path_t2,path_toolbox)

% segm_6C() segments MRI data into a 6 compartment isotropic head model.
%   Data from T1 and T2 MRI-images is used to create a 6C hexahedral FEM 
%   head model (white matter (WM), grey matter (GM), skin, cerebrospinal 
%   fluid (CSF), bone compacta and bone spongiosa). This script provides
%   geometry-adapted hexahedral meshes in vista-format that can be used for
%   forward computation of current flow (using SimBio) in the framework of 
%   transcranial electric stimulation (tES) and associated inverse
%   optimization of stimulation montages.
%
% Cite: Radecke et al., "Simulating individually targeted transcranial 
%   electric stimulation for experimental application" and refer to the 
%   used MATLAB, FieldTrip and SPM versions. 
% Dependencies
%   - FielTrip 20170326 (http://www.fieldtriptoolbox.org)
%   - SPM 12 (www.fil.ion.ucl.ac.uk/spm/)
%   Make sure that both toolboxes are downloaded and stored in path_toolbox
% The SimBio toolbox (https://www.mrt.uni-jena.de/simbio) can be utilized
%   for forward modeling of tES-induced electric fields using the
%   vista mesh created in the current script.
% SciRun (http://www.sci.utah.edu/cibc-software/scirun.html) can be used
%   for high-resolution visualization of the estimated electric fields
%   using the scirun mat-file created in the current script.
% 
% INPUT     
%   pathout         directory; path to store resulting files
%   path_t1         directory; path to t1 nii-file (NIfTI)
%   path_t2         directory; path to t2 nii-file
%   path_toolbox    directory; path to FieldTrip and SPM toolboxes
%
% OUTPUT
%   fem             structure; FieldTrip structure holding the
%                   finite-element head model (FEM), resulting from the 
%                   segmentation. This structure is only computed e.g. for 
%                   the convenient further processing of EEG data using 
%                   FieldTrip using the same head model as used for tES 
%                   forward modeling in e.g. the SimBio toolbox.
%   mesh            structure; FieldTrip structure holding the
%                   geometry-adapted hexahedral mesh. The mesh is
%                   additionally stored as vista v-file and scirun mat-file
%                   for further processing
%   segm            structure; FieldTrip structure holding the segmented
%                   head volume integrating information of T1 and T2 MRI
%                   data. 
% 
% EXAMPLE 
%   pathout = '/FILEPATH/';
%   addpath(genpath(pathout));      % segm_6C.m, segm_interp.m and checkSeg.m are stored in pathout/scripts/
%   path_t1 = [pathout,'t1.nii'];   % nii files were stored in pathout
%   path_t2 = [pathout,'t2.nii'];
%   path_toolbox = '/Applications/MATLAB/toolboxes/';
%   [fem,mesh,segm] = segm_6C(pathout,path_t1,path_t2,path_toolbox);
%
%                                               by Jan-Ole Radecke, 05/2017
%                                                                    v1.0.0
% v1.0.0 (26.07.2019)
% Reduced code for distribution. Cite: Radecke et al., "Simulating 
%   individually targeted transcranial electric stimulation for 
%   experimental application".
% Please note that segm_6C was optimized for the MRI data aquisition and 
%   the specific stimulation targets in the above-mentioned study. 
%   Therefore, it cannot be adapted to any other application without
%   extensive validation of the segmentation results.
% The function was run and tested using MATLAB 9.1.0.441655 (2016b) for MAC
%   on a macOS 10.12.6 (Sierra) machine.
% Error-handling
%   Due to the different versions of FieldTrip and SPM (and 
%   co-dependencies) that were available to the time of  implementation, 
%   some typical errors may occure that are related to version confusions. 
%   For example, FielTrip 20170326 uses legacy SPM versions for the 
%   segmentation which is done using SPM8 or even SPM2. Newer versions of 
%   FieldTrip and SPM-MATLAB scripts might shortcut some of the operations
%   below, but were not tested. If an error occurs (typical in 
%   spm_vol_hdr.m or ft_preamble.m) try the following workarounds:
%   A) restart MATLAB (reset paths).
%   B) check for fixed paths to SPM or FieldTrip versions in the MATLAB 
%       preferences. Delete these and use dynamic definition of paths and
%       toolbox versions in the function below.
%   C) check the respective versions of FieldTrip and SPM, after adding
%       them in the respective lines below.
%   D) If the segmentation volumes are cut off, replace l.130 in the 
%       FieldTrip function ft_volumereslice: 
%           range = [-127.5 127.5] * cfg.resolution;
%       by e.g.
%           range = [-150 150] * cfg.resolution; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters 
Tws = {'t1','t2'};
compartsSPM = {'grey','white','csf','bone','skin'}; % SPM compartments
pMinSpongiosa   = 0.5;  % integer; p-value threshold for bone spongiosa
pSegThresh      = 0.4;  % integer; p-value threshold for SPM probability maps
interpSize      = 2;    % integer [>= 1]; number of neighboring voxels that are included in the iterative interpolation of tissue labels
hexaShift       = 0.3;  % integer; shifting parameter for geometry-adapted hexhedral mesh 

compartsFT = {'grey','white','csf','skin','spongiosa','compacta'}; % final compartments
outputTissueLabel = compartsFT;
tissueConductivity = [0.33 0.14 1.79 0.43 0.025 0.007]; % S/m  % Wagner et al., 2016, Siam J Appl Math

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% coregister t1 and t2 images
addpath(genpath([path_toolbox,'spm12/']));
hdr = [];
path_t = {path_t1,path_t2};
for iT = 1:2
    hdrName{iT} = path_t{iT};
    cd(pathout);
    hdr{iT} = spm_vol(hdrName{iT}); % get nii header
end % for iT

% coregister t1 on t2
invhdr = hdr;
invhdr(1) = hdr(2); invhdr(2) = hdr(1);
x = spm_coreg(invhdr{1},invhdr{2});

% reslice t1 image to the size of the t2 image
% OUTPUT: r*.nii files for t1 and t2 each
cfg.mean = 0;
invhdrNames = {invhdr{1}.fname,invhdr{2}.fname};
spm_reslice(invhdrNames,cfg);


%%  segment data using SPM 12 batch-manager
% OUTPUT: five c*.nii files for t1 and t2 each, stored in pathout: c1 = grey matter, c2 = white matter, c3 = CSF, c4 = bone, c5 = skin/soft tissue
addpath(genpath([pathout,'/segment_batch/']));
for iT = 1:2
    matlabbatch = segment_batch_job(path_toolbox); % create batch-file for current data
    cd(pathout);
    disp(' '); disp(['Please select preprocessed data (',Tws{iT},') manually: '])
    % run batch-file and check for missing input serially. Automatic filename assignment ran into error
    spm_jobman('serial',matlabbatch); % select nii-file that was realigned and resliced from the pop-up batch GUI (r*.nii) 
end % for iT


%% import *.nii-files to FieldTrip format  
addpath([path_toolbox,'fieldtrip-20170326/']);
ft_defaults;
addpath(genpath([path_toolbox,'spm12/']));

for iT = 1:2
    filename = dir(pathout); filename = {filename.name}; filename = filename(3:end);
    % list filenames and find *.nii-files
    filenameOrig = filename(cellfun(@(x)strcmp('r',x(1)),filename))';  % get filenames starting with 'rsPRI' (resliced and realigned T1/T2 weighted images)
    filenameOrig = filenameOrig(cellfun(@(x)strcmp('.nii',x(end-3:end)),filenameOrig))'; % extract *.nii (NIfTI)-files
    filenameCmp = filename(cellfun(@(x)strcmp('c',x(1)),filename))';  % get filenames starting with 'c' (segmented datasets produced by SPM)
    filenameCmp = filenameCmp(cellfun(@(x)strcmp(Tws{iT},x(end-5:end-4)),filenameCmp));
    data.(Tws{iT}) = ft_read_mri([pathout,filenameOrig{1}]);
    
    defineCoordsys = 0;
    if  defineCoordsys % determine coordinate system for current dataset
        data.(Tws{iT}) = ft_determine_coordsys(data.(Tws{iT})); close(gcf);
    else data.(Tws{iT}).coordsys = 'spm'; % let FieldTrip know how to interpret the coordinate system
    end
    
    for iCmp = 1:size(compartsSPM,2)
        disp(filenameCmp(iCmp,:));                
        dataCmp = ft_read_mri([pathout,filenameCmp{iCmp}]);
        data.(Tws{iT}).(compartsSPM{iCmp}) = dataCmp.anatomy;
    end % for iCmp
% checkSeg(data.(Tws{1}),data.(Tws{2}));
end % for iT


%% build 6C model
addpath(genpath([path_toolbox,'spm8/']));
addpath(genpath([path_toolbox,'fieldtrip-20170326/'])); % error in ft_preamble when using ft_defaults

for iT = 1:length(Tws)
    cfg = [];
    cfg.output = {'scalp'};
    headMask = ft_volumesegment(cfg,data.(Tws{iT})); % binary scalp mask
    data.(Tws{iT}).scalp = headMask.scalp;
end % for iTw
% checkSeg(data.t1,data.t2.anatomy);                       

% headmask
headmask = data.(Tws{1}).scalp;                                      
headmask(data.(Tws{2}).scalp) = 1;

s = strel_bol(3);                                                  
headmask = imdilate(headmask,s);
headmask = imclose(headmask,s);
headmask = imfill(headmask,'holes');
headmask = imerode(headmask,s); 

% temporary datamatrix dataCmp (X x Y x Z x tissue) with: comparts = [T1grey,T1white,T2csf,mean(T1/T2skin),mean(T1/T2bone)]
dataCmp = [];                                                                           % element in dim 4:
dataCmp = data.(Tws{1}).grey;                                                           % 1
dataCmp = cat(4,dataCmp,data.(Tws{1}).white);                                           % 2
dataCmp = cat(4,dataCmp,data.(Tws{2}).csf);                                             % 3
dataCmp = cat(4,dataCmp,squeeze(mean(cat(4,data.(Tws{1}).skin,data.(Tws{2}).skin),4))); % 4
dataCmp = cat(4,dataCmp,zeros(size(dataCmp(:,:,:,1))));                                 % 5 zeros, for p.spongiosa
dataCmp = cat(4,dataCmp,zeros(size(dataCmp(:,:,:,1))));                                 % 6 zeros, for p.compacta
% dataCmp = cat(4,dataCmp,data.(Tws{1}).bone);                                          % 7 (bone = 7 based on T1 images)
% dataCmp = cat(4,dataCmp,data.(Tws{2}).bone);                                          % 7 (bone = 7 based on T2 images)
dataCmp = cat(4,dataCmp,squeeze(mean(cat(4,data.(Tws{1}).bone,data.(Tws{2}).bone),4))); % 7 (bone = 7)

% threshold p-values
disp(['Threshold probability maps to p = ',num2str(pSegThresh),'...']);
dataCmp(dataCmp < pSegThresh) = 0;
% define background
dataCmp(repmat(headmask,1,1,1,size(dataCmp,4)) == 0) = 0; 

% find voxel-wise maximum p across tissues
segm = data.(Tws{1});
segm = rmfield(segm,compartsSPM); segm.hdr.segm_hdr = []; segm = rmfield(segm,'scalp');
segm.headmask_h = headmask;
% checkSeg(segm,headmask);                                          
segm.segment_p.tissue_p = dataCmp;
segm.segment_p.tissue_label = {'T1w_grey','T1w_white','T2w_csf','mean_T1wT2w_skin',[],[],'mean_T1wT2w_skull'};
segm.hdr.segm_hdr.pSegThresh = pSegThresh;
segm.hdr.segm_hdr.pMinSpongiosa = pMinSpongiosa;
% redefine relevant tissues from boolean to 1-7   
dataMax = repmat(max(dataCmp,[],4),[1,1,1,size(dataCmp,4)]);
dataCmp(dataCmp ~= dataMax) = 0;
for iDim = 1:size(dataCmp,4)                                       % relabel compartments (1-4 = grey, white, csf, skin; 7 = bone; 5-6 = later p.compacta/spongiosa)
    dataCmpTemp2 = []; dataCmpTemp2 = dataCmp(:,:,:,iDim);
    dataCmpTemp2(dataCmpTemp2 ~= 0) = iDim;
    dataCmp(:,:,:,iDim) = dataCmpTemp2;
end % for iDim

% brain mask
brainmask = (dataCmp(:,:,:,1) > 0|dataCmp(:,:,:,2) > 0|dataCmp(:,:,:,3) > 0); % GM | WM | CSF  
s = strel_bol(2);
brainmask = imdilate(brainmask,s);
brainmask = imclose(brainmask,s);
brainmask = imfill(brainmask,'holes');
brainmask = imerode(brainmask,s);                   
segm.brainmask_h = brainmask;
% checkSeg(segm,brainmask);                 
% checkSeg(segm,headmask);                 

% separate bone compacta and spongiosa
bonemask = squeeze(dataCmp(:,:,:,7));
bonemask(bonemask ~= 0) = 1; % boolean labels
bonemask = imclose(bonemask,strel_bol(5));
bonemask = imfill(bonemask, 'holes');
% checkSeg(data,bonemask);                                          
brainmaskTemp = imerode(brainmask,strel_bol(2));
bonemask(brainmaskTemp == 1) = 0;
% limit bone to being max 19mm from the brain surface
brainmaskTemp2 = imdilate(brainmask,strel_bol(19));
bonemask(brainmaskTemp2 == 0) = 0;
segm.bonemask_h = bonemask;
% checkSeg(segm,bonemask);                                          
bonemaskErode  = imerode(bonemask,strel_bol(3)); % erode bonemask to avoid leakage artifacts
% checkSeg(segm,bonemaskErode);                                          

% bone spongiosa defined based on the eroded bone mask thresholded by p-values from t2 images
boneROI = data.(Tws{2}).bone;
boneROI(bonemaskErode == 0) = 0;
boneROI(boneROI >= pMinSpongiosa & bonemaskErode == 1) = 5;
tempBone = boneROI; tempBone(tempBone ~= 5) = 0;
tempBone = imclose(tempBone,strel_bol(3));
tempBone = imfill(tempBone,'holes');
dataCmp(:,:,:,5) = tempBone;
% checkSeg(segm,dataCmp(:,:,:,5));                                          
tempBone = zeros(size(dataCmp(:,:,:,6)));
tempBone(dataCmp(:,:,:,5) ~= 5 & bonemask == 1) = 6; % define bone compacta
dataCmp(:,:,:,6) = tempBone;
% checkSeg(segm,dataCmp(:,:,:,6));                                          

% skin
bonebrainmask = bonemask > 0 | brainmask > 0;
s = strel_bol(2);
bonebrainmask = imdilate(bonebrainmask,s);
bonebrainmask = imclose(bonebrainmask,s);
bonebrainmask = imfill(bonebrainmask,'holes');
bonebrainmask = imerode(bonebrainmask,s);                
% checkSeg(segm,bonebrainmask);                                          
tempSkin = zeros(size(dataCmp(:,:,:,4)));
tempSkin(headmask == 1 & bonebrainmask == 0) = 4;
tempSkin = imclose(tempSkin,strel_bol(2));
dataCmp(:,:,:,4) = tempSkin;
% checkSeg(segm,tempSkin);                                          

disp('Finished computing masks...');


%% integration of tissue compartments (interpolation)
disp('Interpolate ambigue and missing tissue labels...');
dataCmp = dataCmp(:,:,:,1:6);
segm.tissue = segm_interp(dataCmp,interpSize,segm);
segm.tissuelabel = compartsFT;

% checkSeg(segm,segm.tissue);
% for iCmp = 1:size(dataCmp,4)                             
%     checkSeg(segm,dataCmp(:,:,:,iCmp));
% end % for iCmp


%% prepare 3C segmentation, prepare trihedral or hexahedral mesh and compute FEM and/or BEM headmodel
% addpath(genpath([path_toolbox,'fieldtrip-20170326/']));
addpath([path_toolbox,'fieldtrip-20170326/']); ft_defaults;

% transform datastructure to something fieldtrip understands
segmTmp = segm;
segm = rmfield(segm,{'tissuelabel','bonemask_h','headmask_h','brainmask_h','segment_p','tissue'});
if  sum(ismember(fieldnames(segm),'inside')) > 0 segm = rmfield(segm,'inside'); end
for iCmp = 1:length(segmTmp.tissuelabel)
    tempTissue = segmTmp.tissue;
    tempTissue(tempTissue ~= iCmp) = 0;
    tempTissue(tempTissue > 0) = 1;
    segm.(segmTmp.tissuelabel{iCmp}) = logical(tempTissue); % fieldTrip extracts logical matrices only
end % for iCmp

% transform to CTF
cfg         = [];
cfg.method  = 'interactive';
cfg.coordsys = 'ctf';            
segm = ft_volumerealign(cfg,segm);
% tempfid = segm.cfg.fiducial; 
% reslice data
cfg     = [];
segm    = ft_volumereslice(cfg,segm);

% cut the neck to minimize the mesh volume (in CTF space) 
maxCut = floor(abs(min(segm.cfg.zrange)-min(segm.cfg.zrange)/3));
segm.anatomy(:,:,1:maxCut) = 0;
segm.inside(:,:,1:maxCut) = 0;
for iCmp = 1:size(compartsFT,2)
    segm.(compartsFT{iCmp})(:,:,1:maxCut) = 0;
end
segm.fiducial.vox = segm.cfg.previous.fiducial; % store fiducials in voxel space      

% OUTPUT (RAS) using data
% export fieldTrip 6C segmentation datastructure as single nii volume before reslicing 
xSegm           = segmTmp;
xSegm.anatomy   = segmTmp.tissue;
cfg.parameter   = 'anatomy';
cfg.filename    = [pathout,'6C_segm_ras'];
cfg.filetype    = 'nifti';
ft_volumewrite(cfg,xSegm);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Define hexahedral mesh to prepare head model (Shift: ',num2str(hexaShift),')...']);
cfg        = [];
cfg.tissue = outputTissueLabel;
cfg.shift  = hexaShift;
cfg.method = 'hexahedral';
mesh = ft_prepare_mesh(cfg,segm);

% OUTPUT: saves vista mesh for forward modeling using SimBio (and adjacent inverse optimization of tES montages)
vistaName = ['6C_mesh_vista.v'];
write_vista_mesh([pathout,vistaName],mesh.pos,mesh.hex,mesh.tissue); % vista-format is needed to compute the forward model with SimBio
% OUTPUT: saves SciRun mesh for visualization
scirunName      = ['6C_mesh_scirun'];
scirun.node     = mesh.pos;
scirun.cell     = mesh.hex;
scirun.field    = mesh.tissue;
save([pathout,scirunName],'scirun','-v7'); % -v7 needed for SciRun

% compute FielTrip FEM for further processing of EEG data etc.
cfg        = [];
cfg.method = 'simbio';
cfg.conductivity = tissueConductivity;
fem = ft_prepare_headmodel(cfg,mesh);

disp('Finished 6C segmentation, mesh and fem computation.')

end % function
