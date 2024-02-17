% Calculates a smoothed version of the median.
%
% -- Function File: M = smoothmedian (X)
% -- Function File: M = smoothmedian (X, DIM)
% -- Function File: M = smoothmedian (X, DIM, TOL)
%
%     If X is a vector, find the univariate smoothed median (M) of X. If X is a
%     matrix, compute the univariate smoothed median value for each column and
%     return them in a row vector.  If the optional argument DIM is given,
%     operate along this dimension. Arrays of more than two dimensions are not
%     NaN values currently supported. The MEX file versions of this function
%     ignore (omit) whereas the m-file includes NaN in it's calculations. Use
%     the 'which' command to establish which version of the function is being
%     used.
%
%     The smoothed median is a slightly smoothed version of the ordinary median
%     and is an M-estimator that is both robust and efficient:
%
%     | Asymptotic                            | Mean |    Median  |    Median  |
%     | properties                            |      | (smoothed) | (ordinary) |
%     |---------------------------------------|------|------------|------------|
%     | Breakdown point                       | 0.00 |      0.341 |      0.500 |
%     | Pitman efficacy                       | 1.00 |      0.865 |      0.637 |
%
%     Smoothing the median is achieved by minimizing the objective function:
%
%           S (M) = sum (((X(i) - M).^2 + (X(j) - M).^2).^ 0.5)
%                  i < j
% 
%     where i and j refers to the indices of the Cartesian product of each
%     column of X with itself. 
%
%     With the ordinary median as the initial value of M, this function
%     minimizes the above objective function by finding the root of the first
%     derivative using a fast, but reliable, Newton-Bisection hybrid algorithm.
%     The tolerance (TOL) is the maximum value of the step size that is
%     acceptable to break from optimization. By default, TOL = range * 1e-04.
%
%     The smoothing works by slightly reducing the breakdown point of the median.
%     Bootstrap confidence intervals using the smoothed median have good
%     coverage for the ordinary median of the population distribution and can be
%     used to obtain second order accurate intervals with Studentized bootstrap
%     and calibrated percentile bootstrap methods [1]. When the population
%     distribution is thought to be strongly skewed, coverage errors can be
%     reduced by improving symmetry through appropriate data transformation.
%     Unlike kernel-based smoothing approaches, bootstrapping smoothmedian does
%     not require explicit choice of a smoothing parameter or a probability
%     density function.
%
%  Bibliography:
%  [1] Brown, Hall and Young (2001) The smoothed median and the
%       bootstrap. Biometrika 88(2):519-534
%
%  smoothmedian (version 2023.05.02)
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


function M = smoothmedian (x, dim, Tol)

  % Evaluate input arguments
  if (nargin < 1) || (nargin > 3)
    error ('smoothmedian: Invalid number of input arguments')
  end

  if (nargout > 1)
    error ('smoothmedian: Invalid number of output arguments')
  end

  if (numel (size (x)) > 2)
    error ('smoothmedian: X cannot have more than 2 dimensions')
  end

  % Check data dimensions
  if (nargin < 2)
    if (size (x, 2) == 1)
      dim = 1;
    elseif (size (x, 1) == 1)
      dim = 2;
    else
      dim = 1;
    end
  end
  if (~ ismember (dim, [1, 2]))
    error ('smoothmedian: DIM must be a valid dimension');
  end

  % If applicable, switch dimension
  if (dim > 1)
    x = x.';
  end

  % Check input data type
  if (~ isa (x, 'double'))
    error ('smoothmedian: X must be double precision')
  end

  % Obtain data dimensions
  s = size (x);
  m = s(1);
  n = s(2);
  l = m * (m - 1) / 2;
  
  % Sort the data and calculate the median for each column of the data
  x = sort(x, 1);
  mid = 0.5 * m;
  M = x(fix (mid + 1), 1 : n); % Median when m is odd
  if ( mid == fix (mid) ) 
      % Median when m is even
      M = M + x(fix (mid), 1 : n);
      M = M * 0.5;
  end

  % Set initial bracket bounds and calculate range along each column
  a = x(1,:);           % Minimum
  b = x(m,:);           % Maximum
  range = (b - a) / 2;  % Range
  
  % Check/set tolerance
  if ((nargin < 3) || isempty(Tol))
    Tol = range * 1e-04;
  else 
    Tol = Tol * ones (1, n);
  end

  % Obtain m(m-1)/2 pairs from the Cartesian product of each column of
  % x with itself by enforcing the restriction i < j on xi and xj
  q = logical (triu (ones (m, m), 1));
  i = uint32 ((1:m)' * ones (1, m));
  xi = x(i(q), :);
  j = uint32 (ones (m,1) * (1:m));
  xj = x(j(q), :);
  idx = 1:n;

  % Nonlinear root finding by Newton-Bisection hybrid algorithm
  % Set starting value as the median
  p = M;
  
  % Calculate commonly used operations and assign them to new variables
  z = (xi - xj).^2;
  y = xi + xj;
  
  % Minimize objective function (vectorized)
  MaxIter = 20;
  for Iter = 1:MaxIter
  
    % Compute derivatives
    temp = ones (l, 1) * p;
    D = (xi - temp).^2 + (xj - temp).^2;
    D (D == 0) = 1; % Ensures that no NaN values occur when the
                    % objective function is not differentiable
    R = sqrt(D); 
    % Objective function (S)
    %S = sum(R);
    % First derivative (T)
    T = sum ((2 * temp - y) ./R, 1); 
    % Second derivative (U)
    % Equivalent to (but much faster to compute than):
    % U = sum ( (xi-xj).^2 .* ((xi-temp).^2 + (xj-temp).^2).^(-3/2) )   
    U = sum(z .* R ./ D.^2, 1);
    
    % Reduce memory usage
    temp = []; %#ok<NASGU> Faster than using clear.
    D = [];    %#ok<NASGU> Faster than using clear.
    R = [];    %#ok<NASGU> Faster than using clear.
    
    % Compute Newton step (fast quadratic convergence but unreliable)
    step = T./U; 
    
    % Evaluate convergence
    cvg = ((abs(step) <= Tol) | (range <= Tol));
    if any(cvg)
      % Export converged parameters
      M(idx(cvg)) = p(cvg);
      % Avoid excess computations in following iterations
      idx(cvg) = [];
      xi(:, cvg) = [];
      xj(:, cvg) = [];
      z(:, cvg) = [];
      y(:, cvg) = [];
      a(cvg) = [];
      b(cvg) = [];
      p(cvg) = [];
      step(cvg) = [];
      Tol(cvg) = [];  
    end
    
    % Break from loop when all optimizations have converged
    if (all(cvg))
      break
    end
    
    % Update bracket bounds
    a(step < 0) = p(step < 0) + Tol(step < 0);
    b(step > 0) = p(step > 0) - Tol(step > 0);
                
    % Update the range with the distance between bracket bounds
    range = b - a;
    
    % Preview new value of the smoothed median
    nwt = p - step;
    
    % Prefer Newton step if it is within brackets
    I = (nwt > a) & (nwt < b);
    p(I) = nwt(I);
    
    % Otherwise, compute Bisection step (slow linear convergence but very safe)
    p(~I) = 0.5 * (a(~I) + b(~I));

    % Tidy up
    nwt = [];  %#ok<NASGU> Faster than using clear.
    I = [];    %#ok<NASGU> Faster than using clear.
    T = [];    %#ok<NASGU> Faster than using clear.
    U = [];    %#ok<NASGU> Faster than using clear.
    
  end

  % Set the smoothmedian to NaN where columns/rows contain NaN
  M(any (isnan (x))) = NaN;

  if (Iter == MaxIter)
    fprintf('Warning: Root finding failed to reach the specified tolerance.\n');
  end

  % If applicable, switch dimension
  if (dim > 1)
    M  = M.';
  end
