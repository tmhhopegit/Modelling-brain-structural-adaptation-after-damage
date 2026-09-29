function x = code_demographic(values, kind)
%CODE_DEMOGRAPHIC Turn text demographics into 0/1 predictors, the same way for both datasets.
%   x = code_demographic(values, 'right_handed')  1 if the text starts with 'right'
%       ('Right', 'Right now left'); 0 otherwise ('Left', 'Ambi', 'Left now right').
%       Pre-stroke handedness, in other words.
%   x = code_demographic(values, 'female')        1 for 'F' or 'Female'
%   x = code_demographic(values, 'english')       1 if the text contains 'english'
%   Case is ignored. Empty or missing text gives NaN.
%
%   The originals coded handedness in OPPOSITE directions: PLORAS_Raw set 1
%   for left/ambidextrous/'...now Right', create_chile_mat set 1 for right.
%   So the Chilean patients' handedness was reversed relative to the model
%   trained on PLORAS.
values = cellstr(string(values(:)));
x = nan(numel(values), 1);
for i = 1:numel(values)
    v = strtrim(lower(values{i}));
    if isempty(v) || strcmp(v, 'nan') || strcmp(v, '<missing>'), continue; end
    switch kind
        case 'right_handed'
            x(i) = double(startsWith(v, 'right'));
        case 'female'
            x(i) = double(any(strcmp(v, {'f', 'female'})));
        case 'english'
            x(i) = double(contains(v, 'english'));
        otherwise
            error('lesion_growth:code', 'Unknown demographic ''%s''.', kind);
    end
end
end
