function [volume, left, right] = lesion_volumes(lesions, voxels, dims, midline_x)
%LESION_VOLUMES Lesion size in voxels, in total and per hemisphere.
%   [volume, left, right] = lesion_volumes(lesions, voxels, dims, midline_x)
%   Voxels with x index < midline_x count as right hemisphere (in the
%   91 x 109 x 91 MNI template x index 1 is the right side), the rest as left,
%   as PLORAS_Raw did (right = x 1-45, left = x 46-91).
[x, ~, ~] = ind2sub(dims, voxels);
is_right = x(:)' < midline_x;
volume = double(sum(lesions, 2));
right = double(sum(lesions(:, is_right), 2));
left = volume - right;
end
