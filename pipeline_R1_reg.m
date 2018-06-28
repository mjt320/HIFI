function pipeline_R1_reg(opts)
% co-register nifti images prior to fitting

mkdir(opts.niftiRegDir); delete([opts.niftiRegDir '/*.*']); %make directory/delete contents

regRefIdx=1; %use first series listed in "series" vector as registration target
outputFiles={};
outputFiles{regRefIdx}=[opts.niftiRegDir '/rSeries' num2str(opts.series(regRefIdx),'%02d') '.nii']; %names of output files
copyfile([opts.niftiDir '/series' num2str(opts.series(regRefIdx),'%02d') '.nii'],outputFiles{regRefIdx}); %registration target doesn't change so is just copied

%% co-reg each image to target
for n=1:size(opts.series,2)
    if n==regRefIdx; continue; end %skip reference image
    
    refFile=[opts.niftiDir '/series' num2str(opts.series(regRefIdx),'%02d') '.nii'];
    inputFile=[opts.niftiDir '/series' num2str(opts.series(n),'%02d') '.nii'];
    outputFiles{n}=[opts.niftiRegDir '/rSeries' num2str(opts.series(n),'%02d') '.nii'];
    matFile=[opts.niftiRegDir '/series' num2str(opts.series(n),'%02d') '.mat'];
    
    system([ 'flirt -in ' inputFile ' -ref ' refFile ' -out ' outputFiles{n} ' -omat ' matFile ' -dof 6 -cost normmi']); %FLIRT
    system([ 'fslchfiletype NIFTI ' outputFiles{n}]); %change file type
end


%% make a 4D nifti containing all series
outFileList=[];
for iSeries=1:size(opts.series,2); outFileList=[outFileList ' ' outputFiles{iSeries} ]; end
system(['fslmerge -t ' opts.niftiRegDir '/r4D' outFileList]);
system(['fslchfiletype NIFTI ' opts.niftiRegDir '/r4D']);

end