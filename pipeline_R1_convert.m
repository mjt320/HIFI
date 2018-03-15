function pipeline_R1_convert(opts)
%convert dicoms to NIFTI and get acquisition parameters

tempDir=[opts.niftiDir '/temp'];
mkdir(opts.niftiDir); delete([opts.niftiDir '/*.*']);
mkdir(tempDir); delete([tempDir '/*.*']);

%% initialise variables
acqPars.NSeries=size(opts.series,2);
acqPars.TR=nan(acqPars.NSeries,1); %time between excitation pulses in acquisition block
acqPars.TE=nan(acqPars.NSeries,1);
acqPars.FA=nan(acqPars.NSeries,1);
acqPars.FADeg=nan(acqPars.NSeries,1);
acqPars.TI=nan(acqPars.NSeries,1); %delay between inversion and start of excitation block
acqPars.NReadout=nan(acqPars.NSeries,1);
acqPars.PECentre=nan(acqPars.NSeries,1);
acqPars.NSlices=nan;

%% for each series convert dicoms and record acquisition parameters
for iSeries=1:size(opts.series,2)
    
    temp2=dir([opts.dicomExamDir '/']);
    temp3=~cellfun(@isempty,regexp({temp2.name},['^' num2str(opts.series(iSeries)) '_'])) | strcmp({temp2.name},num2str(opts.series(iSeries))); %look for directories names 'iSeries' or beginning with 'iSeries_'
    if sum(temp3)~=1; error('Cannot find single unique dicom directory for this series.'); end
    dicomDir=[opts.dicomExamDir '/' temp2(temp3).name];
    
    dicoms=dir([dicomDir '/*.dcm']); %check for dcm files, failing that ima files
    if isempty(dicoms); dicoms=dir([dicomDir '/*.IMA']); end;
    if isempty(dicoms); error(['No dicoms found in ' dicomDir]); end;
    
    %% get acquisition parameters
    if iSeries==1;
        acqPars.NSlices=size(dicoms,1); %get number of slices, assuming this is number of dicom files found
    elseif size(dicoms,1)~= acqPars.NSlices; error('Number of slices is not equal across acquisitions.');
    end
    temp=dicominfo([dicoms(1).folder '/' dicoms(1).name]); %use first dicom header to get TE and FA
    acqPars.TE(iSeries)=0.001*temp.EchoTime;
    acqPars.FADeg(iSeries)=temp.FlipAngle;
    acqPars.FA(iSeries)=((2*pi)/360)*temp.FlipAngle;

    if opts.isIR(iSeries); %is this one of the inversion-recovery SPGR acquisitions?
        switch opts.scanner
            case 'GE'
                acqPars.TI(iSeries)=0.001*temp.InversionTime; %GE "Inversion time" is TI (time between inversion and start of acquisition train)
                acqPars.TR(iSeries)=0.001*temp.RepetitionTime; %GE "Repetition time" is TR (time between pulses in acquisition train)
                acqPars.NReadout(iSeries)=acqPars.NSlices/2; %GE scans half of the slice encodes per inversion pulse;
                acqPars.PECentre(iSeries)=0; %i.e. centric slice encoding
            case 'SIEMENS'
                %Siemens "Inversion time" is time to centre of acquisition train old formula: 0.001*(2/acqPars.NSlices)*(temp.RepetitionTime-temp.InversionTime)
                %Siemens "Repetition time" is time between inversion pulses
                acqPars.TR(iSeries)=opts.ESP(iSeries); %get TR for pulse train
                acqPars.TI(iSeries)=0.001*temp.InversionTime - 0.5*acqPars.NSlices*opts.ESP(iSeries); %calculate TI (delay before pulse train)
                acqPars.NReadout(iSeries)=acqPars.NSlices; %Siemens scans all of the slice encode lines following a single inversion pulse;
                acqPars.PECentre(iSeries)=0.5; %Siemens uses linear slice encoding
            otherwise
                error('opts.scanner unrecognised.');
        end
    elseif strcmp(opts.scanner,'SIEMENS') && ~isnan(opts.ESP(iSeries)) %SPGR acquisitions obtained using MPRAGE with IR switched off
        acqPars.TR(iSeries)=opts.ESP(iSeries); %this is to allow specification of TR when MPRAGE sequence without inversion is used
    else %SPGR sequences
        acqPars.TR(iSeries)=0.001*temp.RepetitionTime; %get TR
    end
    
    %% convert dicoms to 3D niftis
    system(['dcm2niix -f series' num2str(opts.series(iSeries),'%02d') ' -o ' opts.niftiDir ' ' dicomDir]);
    
end
rmdir(tempDir);

%% check whether opts.NSlices is defined - if so, overwrite acqPars.NSlices since it may be necessary to specify slice oversampling
if isfield(opts,'NSlices'); acqPars.NSlices=opts.NSlices; disp('Using user-specified number of slices.'); else disp('Warning: assuming no slice oversampling!'); end

%% display acquisition parameters
disp(['TR: ' num2str(acqPars.TR.')]); disp(['TE: ' num2str(acqPars.TE.')]); disp(['FA (deg): ' num2str(acqPars.FADeg.')]); disp(['TI: ' num2str(acqPars.TI.')]); disp(['Slices: ' num2str(acqPars.NSlices)]);

%% make a 4D nifti containing all series
outFileList=[];
for iSeries=1:size(opts.series,2); outFileList=[outFileList ' ' opts.niftiDir '/series' num2str(opts.series(iSeries),'%02d')]; end
system(['fslmerge -t ' opts.niftiDir '/4D' outFileList]);

save([opts.niftiDir '/acqPars'],'acqPars');
