clear; close all;

%% set up paths
procRoot='/proc5';
addpath('/usr/local/spm/spm12');
addpath([procRoot '/software/relaxometry/HIFI_3a']);
addpath([procRoot '/software/UTILITIES']);

spm('defaults', 'FMRI'); spm_jobman('initcfg');

%% set parameters
opts.scanner='SIEMENS'; %works for Siemens and GE
opts.series=[8:13]; %series numbers (dicom subdirectories should be named either (e.g.) 8 or 8_spgr
opts.isIR=[1 1 0 0 0 0]; %indicate which acquisitions are inverision-recovery
opts.ESP=[5e-3 5e-3 5e-3 5e-3 5e-3 5e-3]; %echo spacing for IR-SPGR (Siemens only - obtain this from protocol)
opts.scaleFactor=[1 1 1 1 1 1]; %scale images before fitting (depends on sequences used)
opts.fit= [ 1 1 1 1 1 1 ]; %indicate which sequences to include in fitting
opts.slices={'all' 'all' 'all'}; %indicate which voxels to fit (useful for testing)
%opts.NSlices=nan; %specify number of k-space partition lines - if not specified, code will assume this is equal to NSlices (which may be incorrect if oversampling used)
opts.threshold=50; % the unscaled maximum signal for a voxel (across all acquisitions used for fitting) must be >= this threshold, otherwise it will not be fitted
opts.NTry=1; %number of fitting attempts with different starting values
opts.NCores=1; %number of cores to use when fitting data
opts.dicomExamDir=[pwd '/../../HIFI_test_data/dicom']; %dicom exam dir
opts.niftiDir='./nifti'; %output dir for raw nifti images
opts.niftiRegDir='./nifti_reg'; %output dir for co-registered nifti images
opts.mapDir='./maps_T1'; %output for parameter maps

%% run pipeline steps
pipeline_R1_convert(opts); %convert dicoms
pipeline_R1_reg(opts); %co-register
pipeline_R1_create_map(opts); %fit

save('./options','opts'); %save settings
