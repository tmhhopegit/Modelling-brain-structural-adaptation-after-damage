function [lesions, ok, voxels, dims, template] = read_lesion_images(files, brain_mask, dims)
%READ_LESION_IMAGES Binary lesion images as a patients x voxels matrix.
%   [lesions, ok, voxels, dims, template] = read_lesion_images(files, brain_mask)
%   [...] = read_lesion_images(files, brain_mask, dims)   require this image size
%
%   files:      cell array of image paths ('' = no image for that patient)
%   brain_mask: image whose non-zero voxels are kept ('' = every voxel)
%   lesions:    N x V logical (voxel > 0; NaN counts as 0)
%   ok:         N x 1, false where the file is missing or unreadable. Those rows
%               are all false and are listed in a warning. (PLORAS_Raw turned
%               an unreadable image into an empty lesion silently, i.e. a
%               patient with volume 0.)
%   voxels:     linear indices of the kept voxels in a volume of size dims
%   template:   SPM header of the first readable image
%   Every image must be in the same space (e.g. 2 mm MNI, 91 x 109 x 91). Needs SPM.
if nargin < 2, brain_mask = ''; end
if nargin < 3, dims = []; end
n = numel(files);
ok = false(n, 1);
headers = cell(n, 1);
for i = 1:n
    if isempty(files{i}) || ~exist(files{i}, 'file'), continue; end
    try
        headers{i} = spm_vol(files{i});
        ok(i) = true;
    catch
    end
end
if ~any(ok)
    error('lesion_growth:images', 'None of the %d lesion images could be opened.', n);
end
template = headers{find(ok, 1)};
if isempty(dims), dims = template.dim; end
if ~isempty(brain_mask)
    M = spm_read_vols(spm_vol(brain_mask));
    if ~isequal(size(M), dims)
        error('lesion_growth:space', 'The brain mask is %s; the lesions are %s.', mat2str(size(M)), mat2str(dims));
    end
    voxels = find(M > 0);
else
    voxels = (1:prod(dims))';
end
lesions = false(n, numel(voxels));
for i = find(ok)'
    if ~isequal(headers{i}.dim, dims)
        error('lesion_growth:space', '%s is %s, not %s like the other images.', files{i}, ...
            mat2str(headers{i}.dim), mat2str(dims));
    end
    try
        X = spm_read_vols(headers{i});
        X(isnan(X)) = 0;
        lesions(i, :) = X(voxels)' > 0;
    catch
        ok(i) = false;
    end
end
if ~all(ok)
    bad = files(~ok);
    bad(cellfun(@isempty, bad)) = {'(no file given)'};
    warning('lesion_growth:images', '%d of %d lesion images missing or unreadable:\n  %s', ...
        nnz(~ok), n, strjoin(bad(1:min(end, 20)), '\n  '));
end
end
