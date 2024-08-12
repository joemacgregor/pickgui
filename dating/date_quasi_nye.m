function varargout          = date_quasi_nye(depth_bound, age_bound, strain_rate, depth, var_tol, iter_max, return_type)
% DATE_QUASI_NYE Date an isochrone within a bounded layer using a self-consistent intra-layer uniform vertical strain rate.
%   
%   [AGE/DEPTH,STRAIN_RATE] = DATE_QUASI_NYE(DEPTH_BOUND,AGE_BOUND,STRAIN_RATE,DEPTH/AGE,TOL,ITER_MAX,RETURN_TYPE)
%   determines the age AGE (yr) or DEPTH (m) of the isochrone at depth
%   DEPTH (m) or AGE (yr), respectively, within the layer sandwiched
%   between two dated isochrones with depths DEPTH_BOUND (m) and AGE_BOUND
%   (yr) by assuming that the vertical strain rate is uniform within that
%   layer. This approach is termed "quasi-Nye" dating by MacGregor et al.
%   [2015]. It can also return the best-fit strain rate. STRAIN_RATE (1/yr)
%   is the initial guess at the uniform vertical strain rate within the
%   layer, VAR_TOL (m or yr) is the tolerance on the depth/age residual and
%   ITER_MAX is the maximum number of iterations. RETURN_TYPE is the string
%   that determines which return is desired ('age' or 'depth' only).
%   
%   This approach is described in two studies:
%   
%   MacGregor, J.A. et al., 2012, Spatial variation of englacial radar
%   attenuation: Modeling approach and application to the Vostok flowline
%   (2012), J. Geophys. Res., 117, F03022.
%   
%   MacGregor, J.A. et al, 2015, Radiostratigraphy and age structure
%   of the Greenland Ice Sheet, J. Geophys. Res. Earth Surf., 120
%   
% Joe MacGregor (NASA), Ed Waddington (UW)
% Last updated: 13 May 2024

if (nargin ~= 7)
    error('date_strain:nargin', ['Number of arguments (' num2str(nargin) ') is not equal to 7.'])
end
if (~isnumeric(depth_bound) || ~isvector(depth_bound) || (numel(depth_bound) ~= 2) || ~issorted(depth_bound))
    error('date_strain:depthbound', 'DEPTH_BOUND is not a numeric 2-element sorted vector.')
end
if (~isnumeric(age_bound) || ~isvector(age_bound) || (numel(age_bound) ~= 2) || ~issorted(age_bound))
    error('date_strain:agebound', 'AGE_BOUND is not a numeric 2-element sorted vector.')
end
if (~isnumeric(depth) || ~isscalar(depth) || (depth <= 0))
    error('date_strain:depth', 'DEPTH is not a numeric positive scalar.')
end
if (~isnumeric(var_tol) || ~isscalar(var_tol) || (var_tol <= 0))
    error('date_strain:tol', 'VAR_TOL is not a numeric positive scalar.')
end
if (~isnumeric(strain_rate) || ~isscalar(strain_rate))
    error('date_strain:strainrate', 'STRAIN_RATE is not a numeric scalar.')
end
if (~isnumeric(iter_max) || ~isscalar(iter_max) || mod(iter_max, 1) || (iter_max < 1))
    error('date_strain:itermax', 'ITER_MAX is not a numeric scalar positive integer.')
end
if (~ischar(return_type) || ~any(strcmp(return_type, {'age' 'depth'})))
    error('date_strain:returntype', 'RETURN_TYPE is not a string that is either ''age'' or ''depth''.')
end
if (nargout > 2)
    error('date_strain:nargout', ['Number of outputs (' num2str(nargout) ') is greater than 2.'])
end

% grid search for local minimum in residual rate of change, staying well above trivial residual, i.e., e = 0)
strain_rate_test            = logspace((log10(strain_rate) - 2), (log10(strain_rate) + 2));
[~, ind_test_min]           = min(gradient(abs(R(depth_bound, age_bound, strain_rate_test)), strain_rate_test));
strain_rate                 = strain_rate_test(ind_test_min);

% initialize loop counter
ii                          = 0;

% iterate for strain rate; e = e - r/(dr/de); dr/de = (A_b * t) - (A_t * b)
while ((abs(R(depth_bound, age_bound, strain_rate)) > var_tol) && (ii < iter_max))
    strain_rate             = strain_rate - ((R(depth_bound, age_bound, strain_rate) / ((age_bound(2) * T(depth_bound(1), age_bound(2), strain_rate)) - (age_bound(1) * B(depth_bound(2), age_bound(1), strain_rate)))));
    ii                      = ii + 1; % increment counter
end

% age/depth of isochrone if strain rate converged; a = (-1/e) * log(1 - z/H_eff); H_eff = z_t / (1 - exp(-A_t * e)); Equations B5/7/9 of MacGregor et al. [2012, JGR]
varargout{1}                = NaN;
if (ii < iter_max)
    switch return_type
        case 'age'
            if (depth < (depth_bound(1) / (1 - exp(-age_bound(1) * strain_rate))))
                varargout{1}= (-1 / strain_rate) * log(1 - (depth / (depth_bound(1) / (1 - exp(-age_bound(1) * strain_rate)))));
            end
        case 'depth'
            age             = depth;
            varargout{1}    = (depth_bound(1) / (1 - exp(-age_bound(1) * strain_rate))) * (1 - exp(-age * strain_rate));
    end
end

if (nargout == 2)
    varargout{2}            = strain_rate;
end

% top depth exponential; t = z_t * exp(-A_b * e)
function t                  = T(depth_top, age_bot, strain_rate)
t                           = depth_top .* exp(-age_bot .* strain_rate);

% bottom depth exponential; b = z_b * exp(-A_t * e)
function b                  = B(depth_bot, age_top, strain_rate)
b                           = depth_bot .* exp(-age_top .* strain_rate);

% residual; r = (z_t * (1 - exp(-A_b * e'))) - (z_b * (1 - exp(-A_t * e))) = -t + b + z_t - z_b; Equation B8 of MacGregor et al. [2012, JGR]
function r                  = R(depth_bound, age_bound, strain_rate)
r                           = -T(depth_bound(1), age_bound(2), strain_rate) + B(depth_bound(2), age_bound(1), strain_rate) + depth_bound(1) - depth_bound(2);