% Computes the empirical cumulative distribution function (ECDF), accounting for
% the presence of ties. Useful for bootstrap statistics, which often contain ties.
%
% -- Function File: [x, F] = bootcdf (y)
% -- Function File: [x, F] = bootcdf (y, trim)
% -- Function File: [x, F] = bootcdf (y, trim, m)
% -- Function File: [x, F] = bootcdf (y, trim, m, tol)
% -- Function File: [x, F, P] = bootcdf (...)
%
%     '[x, F] = bootcdf (y)' computes the empirical cumulative distribution
%     function (ECDF) of the vector y of length N. This funtction accounts for
%     the presence of ties and so is suitable for computing the ECDF of
%     bootstrap statistics.
%
%     '[x, F] = bootcdf (y, trim)' removes redundant rows of the ECDF when trim
%     is true. When trim is false, x and F are are the same length as y. The
%     default is true.
%
%     '[x, F] = bootcdf (y, trim, m)' specifies the denominator in the
%     calculation of F as (N + m). Accepted values of m are 0 or 1, with the
%     default being 0. When m is 1, quantiles formed from x and F are akin to
%     qtype 6 in the R quantile function.
%
%     '[x, F] = bootcdf (y, trim, m, tol)' applies a tolerance for the absolute
%     difference in y values that constitutes a tie. The default tolerance
%     is 1e-12 for double precision, or 1e-6 for single precision.
%
%     '[x, F, P] = bootcdf (...)' also returns the distribution of P values.
%
%  bootcdf (version 2023.07.05)
%  Author: Andrew Charles Penn
%  https://www.researchgate.net/profile/Andrew_Penn/
%
%  Copyright 2019 Andrew Charles Penn
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see http://www.gnu.org/licenses/


function [x, F, P] = bootcdf (y, trim, m, Tol)

  % Computes empirical cumulative distribution function and p-value distribution
  % in the presence of ties
  % brainder.org/2012/11/28/competition-ranking-and-empirical-distributions/

  % Check input arguments
  if (nargin > 4)
    error ('bootcdf: too many input arguments provided.');
  end
  if (~ isa (y, 'numeric'))
    error ('bootcdf: y must be numeric.');
  end
  if (all (size (y) > 1))
    error ('bootcdf: y must be a vector.');
  end
  if (size (y, 2) > 1)
    y = y.';
  end
  if (nargin < 2)
    trim = true;
  end
  if ( (~ islogical (trim)) && (~ ismember (trim, [0, 1])) )
    error ('bootcdf: m must be scalar.');
  end
  if (nargin < 3)
    % Denominator in calculation of F is (N + m)
    % When m is 1, quantiles formed from x and F are akin to qtype 6
    % www.rdocumentation.org/packages/stats/versions/3.6.2/topics/quantile
    % Hyndman and Fan (1996) Am Stat. 50(4):361-365
    m = 0;
  end
  if (~ isscalar (m))
    error ('bootcdf: m must be scalar.');
  end
  if (~ ismember (m, [0, 1]))
    error ('bootcdf: m must be either 0 or 1');
  end
  if ( (nargin < 4) || isempty (Tol) )
    if isa (y, 'double')
      Tol = 1e-12;
    elseif isa (y, 'single')
      Tol = 1e-6;
    else
      Tol = 0;
    end
  else
    if (~ any (isa (Tol, {'double','single'})))
      error ('bootcdf: tol must be single or double precision.');
    end
  end

  % Check output arguments
  if (nargout > 3)
    error ('bootcdf: too many output arguments requested.');
  end

  % Discard NaN values
  ridx = isnan (y);
  y(ridx) = [];

  % Set Inf values to the maximum finite value
  infidx = isinf (y);
  y(infidx) = max (y (~ infidx));

  % Get size of y
  N = numel (y);

  % Sort values and apply tolerance for ties
  x = sort (y);
  belowTol = find ((x(2:end) - x(1:end-1)) < Tol); % abs() would be redundant
  for i = 1:numel(belowTol)
    x(belowTol(i) + 1) = x(belowTol(i));
  end

  % Create empirical CDF accounting for ties by competition ranking
  [jnk, IA, IC] = unique (x, 'first');
  R = cat (1, IA(2:end) - 1, N);
  F = arrayfun (@(i) R(IC(i)), (1 : N)') / (N + m);

  % Create p-value distribution accounting for ties by competition ranking
  P = 1 - arrayfun (@(i) IA(IC(i)) - 1, (1 : N)') / N;

  % Remove redundancy
  if trim
    M = unique ([x, F, P], 'rows', 'last');
    x = M(:,1); F = M(:,2); P = M(:,3);
  end

end


%!test
%! 
%! % Example 1 from: 
%! % brainder.org/2012/11/28/competition-ranking-and-empirical-distributions/
%!
%! y = [75; 76; 79; 80; 84; 85; 86; 88; 90; 94];
%! [x, F, P] = bootcdf (y, false, 0);
%! F_ref = [0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1];
%! assert (max (abs (F - F_ref)), 0, 1e-10);
%! P_ref = [1; 0.9; 0.8; 0.7; 0.6; 0.5; 0.4; 0.3; 0.2; 0.1];
%! assert (max (abs (P - P_ref)), 0, 1e-10);

%!test
%! 
%! % Example 2 from: 
%! % brainder.org/2012/11/28/competition-ranking-and-empirical-distributions/
%!
%! y = [81; 81; 82; 83; 83; 83; 84; 85; 85; 85];
%! [x, F, P] = bootcdf (y, false, 0);
%! F_ref = [0.2; 0.2; 0.3; 0.6; 0.6; 0.6; 0.7; 1; 1; 1];
%! assert (max (abs (F - F_ref)), 0, 1e-10);
%! P_ref = [1; 1; 0.8; 0.7; 0.7; 0.7; 0.4; 0.3; 0.3; 0.3];
%! assert (max (abs (P - P_ref)), 0, 1e-10);
