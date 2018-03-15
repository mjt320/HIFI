function pipeline_R1_create_map(opts)
%fit HIFI/VFA T1 mapping data

load([opts.niftiDir '/acqPars'],'acqPars'); %load acquisition parameters

mkdir(opts.mapDir); delete([opts.mapDir '/*.*']); %create output dir/delete contents

%% load 4D magnitude data
[signal,xyz]=spm_read_vols(spm_vol([opts.niftiRegDir '/r4D.nii']));
for n=1:size(signal,4); signal(:,:,:,n)=opts.scaleFactor(n)*signal(:,:,:,n); end %scale signal

isIR=logical(opts.isIR) & logical(opts.fit); %images that are IR and should be fitted
isFit=logical(opts.fit); %images that should be fitted

%% initialise output arrays
volTemplate=spm_vol([opts.niftiDir '/series' num2str(opts.series(1),'%02d') '.nii']); %use this header as template for 3D output files
T1=nan(volTemplate.dim); S0=nan(volTemplate.dim); k=nan(volTemplate.dim); RSq=nan(volTemplate.dim); model=nan([volTemplate.dim sum(isFit)]); %initialise output arrays
R1_LCI=nan(volTemplate.dim); R1_UCI=nan(volTemplate.dim);

%% do the fitting
for iDim=1:3
    if strcmp(opts.slices{iDim},'all'); slices{iDim}=1:size(signal,iDim); else slices{iDim}=opts.slices{iDim}; end %determine which indices to fit
end
for i3=slices{3}; for i1=slices{1}; for i2=slices{2}; % loop through voxels (only loop through indices to be fitted)
                
            if max(signal(i1,i2,i3,isFit)./opts.scaleFactor(1))<opts.threshold; continue; end %skip voxels that don't pass threshold criterion
            
            %% run the fitting kernel
            [T1(i1,i2,i3),S0(i1,i2,i3),k(i1,i2,i3),model(i1,i2,i3,:),R1_LCI(i1,i2,i3),R1_UCI(i1,i2,i3),RSq(i1,i2,i3),exitFlag]=...
                fit_R1(squeeze(signal(i1,i2,i3,:)).',isIR,isFit,acqPars.TR,acqPars.FA,acqPars.TI,acqPars.PECentre,acqPars.NReadout,opts.NTry);
            
        end;
    end;
    disp([num2str(i3) '/' num2str(size(signal,3))]); %display progress
end;

%% write output images
iEchoFit=0;
for iEcho=1:size(opts.series,2)
    volModel=volTemplate; volModel.dt=[16 0]; volModel.fname=[opts.mapDir '/model_series_' num2str(iEcho,'%02d') '.nii'];
    volSignal=volTemplate; volSignal.dt=[16 0]; volSignal.fname=[opts.mapDir '/signal_series_' num2str(iEcho,'%02d') '.nii'];
    if isFit(iEcho); %if this echo was used in fitting then write values to files
        iEchoFit=iEchoFit+1;
        spm_write_vol(volModel,model(:,:,:,iEchoFit)); %NB model array only includes echoes used for fitting
        spm_write_vol(volSignal,signal(:,:,:,iEcho));
    else %if echo not used in fitting then write nans to files
        spm_write_vol(volModel,nan(volTemplate.dim));
        spm_write_vol(volSignal,nan(volTemplate.dim));
    end
end
spm_file_merge(sort(getMultipleFilePaths([opts.mapDir '/model_series_*.nii'])),[opts.mapDir '/model.nii'],0);
spm_file_merge(sort(getMultipleFilePaths([opts.mapDir '/signal_series_*.nii'])),[opts.mapDir '/signal.nii'],0);
delete([opts.mapDir '/model_series_*.nii']);
delete([opts.mapDir '/signal_series_*.nii']);

paramNames={'T1' 'S0' 'k' 'RSq' 'R1' 'R1_LCI' 'R1_UCI'};
outputs={T1 S0 k RSq 1./T1 R1_LCI R1_UCI};

for iOutput=1:size(outputs,2)
    volOutput=volTemplate;
    volOutput.fname=[opts.mapDir '/' paramNames{iOutput} '.nii'];
    volOutput.dt=[16 0];
    spm_write_vol(volOutput,outputs{iOutput});
end

spm_file_merge({[opts.mapDir '/R1.nii'] [opts.mapDir '/R1_LCI.nii'] [opts.mapDir '/R1_UCI.nii']},[opts.mapDir '/R1_CI'],0);

end
