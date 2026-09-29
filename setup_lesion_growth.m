function setup_lesion_growth()
%SETUP_LESION_GROWTH Put this project's folders on the MATLAB path.
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root, 'data'), fullfile(root, 'analysis'), fullfile(root, 'tests'));
if isempty(which('spm_vol'))
    warning('lesion_growth:spm', 'SPM is not on the path; it is needed to read lesion images.');
end
end
