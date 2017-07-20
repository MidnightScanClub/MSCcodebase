function [r p] = paircorr_mod(a,b)
%PAIRCORR Computes pairwise Pearson's linear correlation coefficient with
% optional significance. Returns r, a p1-by-p2 matrix containing the
% pairwise correlation coefficient between each pair of columns in the
% n-by-p1 and n-by-p2 matrices a and b. r is calculated as the dot
% product between two vectors divided by the product of their magnitudes.
% If a second output argument is provided, like so:
% [r p] = paircorr(a,b)
% then p is the two-tailed significance.
% TOL 03/01/11.
% Added single input functionality TOL, 04/01/12.

if nargin<2
    b = a;
end

a = bsxfun(@minus, a, mean(a));
b = bsxfun(@minus, b, mean(b));

mag_a = sqrt(sum(a.^2, 1));
mag_b = sqrt(sum(b.^2, 1));

r = (a' * b) ./ (mag_a' * mag_b);

if nargout > 1
    [n p1] = size(a);
    
    % calculate t-statistic
    t = r ./ sqrt((1 - r.^2)/(n - 2));
    % calculate significance, two-tailed
    p = 2 * tcdf(-abs(t), n - 2);
end

