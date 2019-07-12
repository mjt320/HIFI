function pipeline_R1_create_map(opts)
% fit HIFI/VFA T1 mapping data

if ~isfield(opts,'fitOptions')
    fitOptions=struct('tolFun',1e-6,'tolX',1e-6);
else
    fitOptions=opts.fitOptions;
end

load([opts.niftiDir filesep 'acqPars'],'acqPars'); %load acquisition parameters
acqPars=acqPars; %necessary for par loop to work

mkdir(opts.mapDir); delete([opts.mapDir filesep '*.*']); %create output dir/delete contents


%% load 4D magnitude data
[signal,xyz]=spm_read_vols(spm_vol([opts.niftiRegDir '/r4D.nii']));
if any(opts.scaleFactor~=1)
    for n=1:size(signal,4); signal(:,:,:,n)=opts.scaleFactor(n)*signal(:,:,:,n); end %scale signal
end

if length(opts.isIR) ~= length(opts.fit) % check opts consistency
    error('opts.isIR and opts.fit are of different size')
else
    isIR=logical(opts.isIR) & logical(opts.fit); %images that are IR and should be fitted
    isFit=logical(opts.fit); %images that should be fitted
end

%% initialise output arrays

volTemplate =spm_vol([opts.niftiDir filesep 'series' num2str(opts.series(1),'%02d') '.nii']); %use this header as template for 3D output files
T1     = nan(volTemplate.dim);
S0     = nan(volTemplate.dim); 
k      = nan(volTemplate.dim); 
RSq    = nan(volTemplate.dim); 
model  = nan([volTemplate.dim sum(isFit)]); %initialise output arrays
R1_LCI = nan(volTemplate.dim); 
R1_UCI = nan(volTemplate.dim);

%% do the fitting

if isfield(opts,'NCores'); NCores=opts.NCores; else NCores=1; end %check number of cores

if NCores<=1 % non-parallel version (allows fitting selected voxels)
    % --------------------------------------------------------------
    for iDim=1:3
        if strcmp(opts.slices{iDim},'all'); slices{iDim}=1:size(signal,iDim); else slices{iDim}=opts.slices{iDim}; end %determine which indices to fit
    end
    for i1=slices{1}; for i2=slices{2}; for i3=slices{3}; % loop through voxels (only loop through indices to be fitted)
                
                if max(signal(i1,i2,i3,isFit)./opts.scaleFactor(1))<opts.threshold; continue; end %skip voxels that don't pass threshold criterion
                
                % run the fitting kernel
                [T1(i1,i2,i3),S0(i1,i2,i3),k(i1,i2,i3),model(i1,i2,i3,:),R1_LCI(i1,i2,i3),R1_UCI(i1,i2,i3),RSq_temp(i1,i2,i3),exitFlag]=...
                    fit_R1_2(squeeze(signal(i1,i2,i3,:)).',isIR,isFit,acqPars.TR,acqPars.FA,acqPars.TI,acqPars.PECentre,acqPars.NReadout,opts.NTry,fitOptions);
                
            end
        end
        disp([num2str(i1) filesep num2str(size(signal,1))]); %display progress
    end
    
else % parallel version - fit all voxels
    % ----------------------------------------
    delete(gcp('nocreate')); poolobj = parpool('local',NCores);
    slices={1:size(signal,1) 1:size(signal,2) 1:size(signal,3)};
    parfor i1=slices{1} %parallel loop for first dimension
        if ~strcmp(opts.slices{1},'all') && isempty(find(opts.slices{1}==i1)); continue; end
        
        %%create temporary variables for parfor look
        T1_temp     = nan(volTemplate.dim([2 3]));
        S0_temp     = nan(size(T1_temp));
        k_temp      = nan(size(T1_temp));
        R1_LCI_temp = nan(size(T1_temp));
        R1_UCI_temp = nan(size(T1_temp));
        RSq_temp    = nan(size(T1_temp));
        model_temp  = nan([volTemplate.dim([2 3]) sum(isFit)]);
        
        %loop through plane of voxels
        for i2=slices{2};
            if ~strcmp(opts.slices{2},'all') && isempty(find(opts.slices{2}==i2)); continue; end %skip some voxels
            for i3=slices{3};
                if ~strcmp(opts.slices{3},'all') && isempty(find(opts.slices{3}==i3)); continue; end %skip some voxels
                
                if max(signal(i1,i2,i3,isFit)./opts.scaleFactor(1))<opts.threshold; continue; end %skip voxels that don't pass threshold criterion
                
                % run the fitting kernel
                [T1_temp(i2,i3),S0_temp(i2,i3),k_temp(i2,i3),model_temp(i2,i3,:),R1_LCI_temp(i2,i3),R1_UCI_temp(i2,i3),RSq_temp(i2,i3),exitFlag]=...
                    fit_R1(squeeze(signal(i1,i2,i3,:)).',isIR,isFit,acqPars.TR,acqPars.FA,acqPars.TI,acqPars.PECentre,acqPars.NReadout,opts.NTry);
                
            end
        end
        
        % assign sliced variables for this i1 from temporary variables
        T1(i1,:,:)      = T1_temp;
        S0(i1,:,:)      = S0_temp;
        k(i1,:,:)       = k_temp;
        R1_LCI(i1,:,:)  = R1_LCI_temp;
        R1_UCI(i1,:,:)  = R1_UCI_temp;
        RSq(i1,:,:)     = RSq_temp;
        model(i1,:,:,:) = model_temp;
        disp([num2str(i1) filesep num2str(size(slices{1},2))]); %display progress
    end  
end


%% write output images
iEchoFit=0;
for iEcho=1:size(opts.series,2)
    volModel=volTemplate; volModel.dt=[16 0]; volModel.fname=[opts.mapDir filesep 'model_series_' num2str(iEcho,'%02d') '.nii'];
    volSignal=volTemplate; volSignal.dt=[16 0]; volSignal.fname=[opts.mapDir filesep 'signal_series_' num2str(iEcho,'%02d') '.nii'];
    if isFit(iEcho); %if this echo was used in fitting then write values to files
        iEchoFit=iEchoFit+1;
        spm_write_vol(volModel,model(:,:,:,iEchoFit)); %NB model array only includes echoes used for fitting
        spm_write_vol(volSignal,signal(:,:,:,iEcho));
    else %if echo not used in fitting then write nans to files
        spm_write_vol(volModel,nan(volTemplate.dim));
        spm_write_vol(volSignal,nan(volTemplate.dim));
    end
end
spm_file_merge(sort(getMultipleFilePaths([opts.mapDir filesep 'model_series_*.nii'])),[opts.mapDir filesep 'model.nii'],0);
spm_file_merge(sort(getMultipleFilePaths([opts.mapDir filesep 'signal_series_*.nii'])),[opts.mapDir filesep 'signal.nii'],0);
delete([opts.mapDir filesep 'model_series_*.nii']);
delete([opts.mapDir filesep 'signal_series_*.nii']);

paramNames={'T1' 'S0' 'k' 'RSq' 'R1' 'R1_LCI' 'R1_UCI'};
outputs={T1 S0 k RSq 1./T1 R1_LCI R1_UCI};

for iOutput=1:size(outputs,2)
    volOutput=volTemplate;
    volOutput.fname=[opts.mapDir '/' paramNames{iOutput} '.nii'];
    volOutput.dt=[16 0];
    spm_write_vol(volOutput,outputs{iOutput});
end

spm_file_merge({[opts.mapDir filesep 'R1.nii'] [opts.mapDir filesep 'R1_LCI.nii'] [opts.mapDir filesep 'R1_UCI.nii']},[opts.mapDir filesep 'R1_CI'],0);

end
