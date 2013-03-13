% INTERSECTI Linearly interpolated intersections for two vectors.
% 
%   [XINT,YINT,IND1,IND2] = INTERSECTI(X1,Y1,X2,Y2) finds the linearly
%   interpolated intersections XINT and YINT between vectors (X1,Y1) and
%   (X2,Y2). IND1 and IND2 are the indices of the nearest points in (X1,Y1)
%   and (X2,Y2) to the those intersections, respectively.
% 
% Joe MacGregor (UTIG)
% Last updated: 02/01/13

function [x_int, y_int, ind1, ind2] ...
                        = intersecti(x1, y1, x2, y2)

if (length(x1) ~= length(y1))
    error('intersecti:xy1size', 'X1 and Y1 are not the same size.')
end
if (length(x2) ~= length(y2))
    error('intersecti:xy2size', 'X2 and Y2 are not the same size.')
end

% want both (x1,y1) and (x2,y2) to be row vectors
if ~isrow(x1)
    x1                  = x1';
end
if ~isrow(y1)
    y1                  = y1';
end
if ~isrow(x2)
    x2                  = x2';
end
if ~isrow(y2)
    y2                  = y2';
end

% lengths needed later
n1                      = length(x1) - 1; 
n2                      = length(x2) - 1;

% slopes
m1                      = repmat((diff(y1) ./ diff(x1))', 1, n2); 
m2                      = repmat((diff(y2) ./ diff(x2)), n1, 1);

% replicated position vectors offset appropriately
x1a                     = repmat(x1(1:(end - 1))', 1, n2);
x1b                     = repmat(x1(2:end)', 1, n2);
y1a                     = repmat(y1(1:(end - 1))', 1, n2);
y1b                     = repmat(y1(2:end)', 1, n2);
x2a                     = repmat(x2(1:(end - 1)), n1, 1);
x2b                     = repmat(x2(2:end), n1, 1);
y2a                     = repmat(y2(1:(end - 1)), n1, 1);
y2b                     = repmat(y2(2:end), n1, 1);

% linearly interpolated positions for each segment
x                       = (y2a - y1a + (m1 .* x1a) - (m2 .* x2a)) ./ (m1 - m2); 
y                       = (m1 .* (x - x1a)) + y1a;

% fix slopes that are vertical/infinite/x-unchanged by calculating x(y) instead of y(x), so that slope is zero
if any(isinf(m1(:)))
    ind_fix             = find(isinf(m1))';
    y(ind_fix)          = (m2(ind_fix) .* (x1a(ind_fix) - x2a(ind_fix))) + y2a(ind_fix);
    x(ind_fix)          = x1a(ind_fix);
end
if any(isinf(m2(:)))
    ind_fix             = find(isinf(m2))';
    y(ind_fix)          = (m1(ind_fix) .* (x2a(ind_fix) - x1a(ind_fix))) + y1a(ind_fix);
    x(ind_fix)          = x2a(ind_fix);
end

% NaN out any interpolations that are outside the bounds of the segments' x/y ranges
x((x < min(x1)) | (x > max(x1)) | (y < min(y1)) | (y > max(y1)) | (x < min(x2)) | (x > max(x2)) | (y < min(y2)) | (y > max(y2))) ...
                        = NaN;
y(isnan(x))             = NaN;

% find intersections within x/y bounds, separated by quadrant
test_all                = cell(4, 2);
test_all{1, 1}          = find((x >= x1a) & (x <= x1b) & (y >= y1a) & (y <= y1b)); % assumes both x and y increasing, vector 1
test_all{2, 1}          = find((x <= x1a) & (x >= x1b) & (y >= y1a) & (y <= y1b)); % assumes x decreasing, y increasing, vector 1
test_all{3, 1}          = find((x <= x1a) & (x >= x1b) & (y <= y1a) & (y >= y1b)); % assumes both x and y decreasing, vector 1
test_all{4, 1}          = find((x >= x1a) & (x <= x1b) & (y <= y1a) & (y >= y1b)); % assumes x increasing, y decreasing, vector 1
test_all{1, 2}          = find((x >= x2a) & (x <= x2b) & (y >= y2a) & (y <= y2b)); % repeat for vector 2
test_all{2, 2}          = find((x <= x2a) & (x >= x2b) & (y >= y2a) & (y <= y2b));
test_all{3, 2}          = find((x <= x2a) & (x >= x2b) & (y <= y2a) & (y >= y2b));
test_all{4, 2}          = find((x >= x2a) & (x <= x2b) & (y <= y2a) & (y >= y2b));

% find intersections that pass for both segments
ind_int                 = [];
for ii = 1:4
    for jj = 1:4
        ind_curr        = intersect(test_all{ii, 1}, test_all{jj, 2}); % intersections between index sets that pass
        if ~isempty(ind_curr) % only keep if not empty, can't simplify or a deprecation warning results
            ind_int     = [ind_int; ind_curr]; %#ok<AGROW>
        end
    end
end
ind_int                 = unique(ind_int); % keep only unique intersections

% return x/y positions at intersections and the nearest indices for those vectors
x_int                   = x(ind_int);
y_int                   = y(ind_int);
[ind1, ind2]            = ind2sub([n1 n2], ind_int);

% correct ind1/2 by finding which index is closest to the linearly interpolated intersection
for ii = 1:length(x_int)
    ind_test1           = (ind1(ii) - 1):(ind1(ii) + 1); % range about ind1 to test
    ind_test1           = ind_test1(ind_test1 > 0); % only indices that make sense
    [~, ind_test2]      = min(sqrt(((x_int(ii) - x1(ind_test1)) .^ 2) + ((y_int(ii) - y1(ind_test1)) .^ 2))); % find which of the 2 or 3 points are closest to the intersection
    ind1(ii)            = ind1(ii) + ind_test2 - (length(ind_test1) - 1); % correct ind1 to the closest index
    ind_test1           = (ind2(ii) - 1):(ind2(ii) + 1);
    ind_test1           = ind_test1(ind_test1 > 0);
    [~, ind_test2]      = min(sqrt(((x_int(ii) - x2(ind_test1)) .^ 2) + ((y_int(ii) - y2(ind_test1)) .^ 2)));
    ind2(ii)            = ind2(ii) + ind_test2 - (length(ind_test1) - 1);
end