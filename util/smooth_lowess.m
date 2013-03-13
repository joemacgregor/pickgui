function ys                 = smooth_lowess(y, span)
% SMOOTH_LOWESS Smooth data using Lowess method.
%   YS = SMOOTH(Y,SPAN) smooths vector Y using across SPAN using the Lowess
%   method. Adapted from SMOOTH in the Curve Fitting Toolbox.

% make column vector
y                           = y(:);

% For problems where x is uniform, there's a faster way
% For uniform data, an even span actually covers an odd number of
% points.  For example, the four closest points to 5 in the
% sequence 1:10 are {3,4,5,6}, but 7 is as close as 3.
% Therfore force an odd span.
span                        = (2 * floor(span / 2)) + 1;

% Omit points at the extremes, which have zero weight
halfw                       = (span - 1) / 2; % halfwidth of entire span
d                           = abs(((1 - halfw):(halfw - 1))); % distances to pts with nonzero weight

% set up weighted Vandermonde matrix using equally spaced X values
x1                          = (2:(span - 1)) - (halfw + 1);
weight                      = (1 - ((d / halfw) .^ 3)) .^ 1.5; % tri-cubic weight
v                           = [ones(length(x1), 1) x1(:)];
V                           = v .* repmat(weight', 1, size(v, 2));

% QR decomposition
[Q, ~]                      = qr(V, 0);

% The projection matrix is Q*Q'. We want to project onto the middle point,
% so we can take just one row of the first factor.
alpha                       = Q(halfw, :) * Q';

% This alpha defines the linear combination of the weighted y values that
% yields the desired smooth values.  Incorporate the weights into the
% coefficients of the linear combination, then apply filter.
alpha                       = alpha .* weight;
ys                          = filter(alpha, 1, y);

% We need to slide the values into the center of the array.
ys((halfw + 1):(end - halfw)) = ys((span - 1):(end - 1));

% Now we have taken care of everything except the end effects.  Loop over
% the points where we don't have a complete span.  Now the Vandermonde
% matrix has span-1 points, because only 1 has zero weight.
x1                          = 1:(span - 1);
v                           = [ones(length(x1), 1) x1(:)];

for ii = 1:halfw;
    % Compute weights based on deviations from the jth point,
    % then compute weights and apply them as above.
    weight                  = (1 - ((abs((1:(span - 1)) - ii) / (span - ii)) .^ 3)) .^ 1.5;
    V                       = v .* repmat(weight(:), 1, size(v, 2));
    [Q, ~]                  = qr(V, 0);
    alpha                   = Q(ii, :) * Q';
    alpha                   = alpha .* weight;
    ys(ii)                  = alpha * y(1:(span - 1));
    % These coefficients can be applied to the other end as well
    ys(end + 1 - ii)        = alpha * y(end:-1:(end - span + 2));
end;