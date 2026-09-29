function model = lesion_pca_fit(lesions, min_explained)
%LESION_PCA_FIT Principal components of binary lesion maps.
%   model = lesion_pca_fit(lesions, min_explained)
%   lesions: N x V (logical or numeric); min_explained: keep components that
%   explain at least this percentage of the variance (the original used 1).
%   model.mean (1 x V), model.coeff (V x K), model.explained (all components, %),
%   model.kept (K), model.scores (N x K, the training rows' scores).
%
%   Same components as MATLAB's pca() (up to sign), computed from the N x N
%   Gram matrix, so the V x V covariance is never formed and memory stays
%   around N x V single precision: pca() on 500 whole-brain images needs
%   several GB. Use lesion_pca_project for new lesions.
X = single(lesions);
n = size(X, 1);
if n < 2, error('lesion_growth:pca', 'Need at least two lesions.'); end
mu = mean(X, 1);
G = double(X * X');
H = eye(n) - ones(n) / n;
G = H * G * H;                      % Gram matrix of the centred data
G = (G + G') / 2;
[U, L] = eig(G, 'vector');
[L, order] = sort(L, 'descend');
U = U(:, order);
L(L < 0) = 0;
explained = 100 * L / sum(L);
keep = find(explained >= min_explained & L > 1e-10 * max(L));
U = U(:, keep);
sv = sqrt(L(keep))';                % singular values of the centred data
coeff = double(X' * single(U)) ./ sv;   % V x K, unit length (U is orthogonal to the mean, so no centring needed)
model.mean = double(mu);
model.coeff = coeff;
model.explained = explained;
model.kept = numel(keep);
model.scores = U .* sv;             % = (X - mean) * coeff for the training rows
end
