function scores = lesion_pca_project(model, lesions)
%LESION_PCA_PROJECT Scores of lesions on components from lesion_pca_fit.
%   scores = lesion_pca_project(model, lesions)  ->  (lesions - mean) * coeff
%   (StepwiseLM_Nested.m projected test lesions with Les'\XL, a least-squares
%   solve that is not a projection, without removing the mean first.)
scores = double(single(lesions) * single(model.coeff)) - model.mean * model.coeff;
end
