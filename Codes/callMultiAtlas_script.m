% set the target filename
target_fname = ['/tmp/MR-CT-x-260-x-260_MR-x-T1-x-Multi_Atlas/260_1.nii.gz'];
% set the directory where all of the output information will be stored
mdir = ['/tmp/MR-CT-x-260-x-260_MR-x-T1-x-Multi_Atlas/'];
% set the input directory where everything is installed
in_dir = '/data/mcr/full-multi-atlas/';
% add all of the important directories to the path
addpath([in_dir, 'masi-fusion/src/shared/']);
addpath([in_dir, 'masi-fusion/src/ext/']);
addpath([in_dir, 'masi-fusion/src/simulations/']);
addpath([in_dir, 'masi-fusion/src/multi-atlas/']);
% set the important directories for running each component
in_opts.niftyreg_loc = [in_dir, 'niftyreg/bin/'];
in_opts.ants_loc = [in_dir, 'ANTs-bin/'];
in_opts.atlas_loc = [in_dir, 'atlas-processing/'];
in_opts.mni_loc = [in_dir, 'MNI/'];
in_opts.mipav_loc = [in_dir, 'mipav/'];
% run the pipeline
multi_atlas_pipeline(target_fname, mdir, in_opts);
