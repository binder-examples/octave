% Performs a permutation or randomization test to compare the distributions of 
% two independent or paired data samples. 
%
% -- Function File: PVAL = randtest2 (X, Y)
% -- Function File: PVAL = randtest2 (X, Y, PAIRED)
% -- Function File: PVAL = randtest2 (X, Y, PAIRED, NREPS)
% -- Function File: PVAL = randtest2 (X, Y, PAIRED, NREPS)
% -- Function File: PVAL = randtest2 (X, Y, PAIRED, NREPS, FUNC)
% -- Function File: PVAL = randtest2 (X, Y, PAIRED, NREPS, FUNC, SEED)
% -- Function File: PVAL = randtest2 ([X, GX], [Y, GY], ...)
% -- Function File: [PVAL, STAT] = randtest (...)
% -- Function File: [PVAL, STAT, FPR] = randtest (...)
% -- Function File: [PVAL, STAT, FPR, PERMSTAT] = randtest (...)
%
%     'PVAL = randtest2 (X, Y)' performs a randomization (a.k.a. permutation)
%     test to ascertain whether data samples X and Y come from populations with
%     the same distribution. Distributions are compared using the Wasserstein
%     metric [1,2], which is the area of the difference between the empirical
%     cumulative distribution functions of X and Y. The data in X and Y should
%     be column vectors that represent measurements of the same variable. The
%     value returned is a 2-tailed p-value against the null hypothesis computed
%     using the absolute values of the test statistics.
%
%     'PVAL = randtest2 (X, Y, PAIRED)' specifies whether X and Y should be
%     treated as independent (unpaired) or paired samples. PAIRED accepts a
%     logical scalar:
%        o false (default): As above.
%        o true: Performs a randomization or permutation test to ascertain
%                whether paired or matched data samples X and Y come from
%                populations with the same distribution. The vectors X and Y
%                must contain the same number of sampling units.
%
%     'PVAL = randtest2 (X, Y, PAIRED, NREPS)' specifies the number of resamples
%     without replacement to take in the randomization test. By default, NREPS
%     is 5000. If the number of possible permutations is smaller than NREPS, the
%     test becomes exact. For example, if the number of sampling units across
%     two independent samples is 6, then the number of possible permutations is
%     factorial (6) = 720, so NREPS will be truncated at 720 and sampling will
%     systematically evaluate all possible permutations. If the number of
%     sampling units in each paired sample is 12, then the number of possible
%     permutations is 2^12 = 4096, so NREPS will be truncated at 4096 and
%     sampling will systematically evaluate all possible permutations. 
%
%     'PVAL = randtest2 (X, Y, PAIRED, NREPS, FUNC)' also specifies a custom
%     function calculated on the original samples, and the permuted or
%     randomized resamples. Note that FUNC must compute a difference statistic
%     between samples X and Y, and should either be a:
%        o function handle or anonymous function,
%        o string of function name, or
%        o a cell array where the first cell is one of the above function
%          definitions and the remaining cells are (additional) input arguments 
%          to that function (other than the data arguments).
%        See the built-in demos for example usage with the mean or vaiance.
%
%     'PVAL = randtest2 (X, Y, PAIRED, NREPS, FUNC, SEED)' initialises the
%     Mersenne Twister random number generator using an integer SEED value so
%     that that the results of 'randtest2' results are reproducible when the
%     test is approximate (i.e. when using randomization if not all permutations
%     can be evaluated systematically).
%
%     'PVAL = randtest2 ([X, GX], [Y, GY], ...)' also specifies the sampling
%     units (i.e. clusters) using consecutive positive integers in GX and GY
%     for X and Y respectively. Defining the sampling units has applications
%     for clustered resampling, for example in the cases of nested experimental 
%     designs. If PAIRED is false, numeric identifiers in GX and GY must be
%     unique (e.g. 1,2,3 in GX, 4,5,6 in GY). If PAIRED is true, numeric
%     identifiers in GX and GY must by identical (e.g. 1,2,3 in GX, 1,2,3 in
%     GY). Note that when sampling units contain different numbers of values,
%     function evaluations after sampling cannot be vectorized. If the parallel
%     computing toolbox (Matlab) or package (Octave) is installed and loaded,
%     then the function evaluations will be automatically accelerated by
%     parallel processing on platforms with multiple processors.
%
%     '[PVAL, STAT] = randtest2 (...)' also returns the test statistic.
%
%     '[PVAL, STAT, FPR] = randtest2 (...)' also returns the minimum false
%     positive risk (FPR) calculated for the p-value, computed using the
%     Sellke-Berger approach.
%
%     '[PVAL, STAT, FPR, PERMSTAT] = randtest2 (...)' also returns the
%     statistics of the permutation distribution.
%
%  Bibliography:
%  [1] Dowd (2020) A New ECDF Two-Sample Test Statistic. arXiv.
%       https://doi.org/10.48550/arXiv.2007.01360
%  [2] https://en.wikipedia.org/wiki/Wasserstein_metric
%
%  randtest2 (version 2023.09.16)
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

function [pval, stat, fpr, STATS] = randtest2 (x, y, paired, nreps, func, seed)

  % Check if we are running Octave or Matlab
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Store local function wass_stat in a stucture (for parallel processes)
  localfunc = struct ('wass_stat', @wass_stat);

  % Check if we have parallel processing capabilities
  PARALLEL = false; % Default
  if (ISOCTAVE)
    software = pkg ('list');
    names = cellfun (@(S) S.name, software, 'UniformOutput', false);
    status = cellfun (@(S) S.loaded, software, 'UniformOutput', false);
    index = find (~ cellfun (@isempty, regexpi (names, '^parallel')));
    if ( (~ isempty (index)) && (logical (status{index})) )
      PARALLEL = true;
    end
  else
    try 
      pool = gcp ('nocreate'); 
      PARALLEL = ~ isempty (pool);
    catch
      % Do nothing
    end
  end

  % Check the number of function arguments
  if (nargin < 2)
    error ('randtest2: X and Y must be provided');
  end
  if (nargin > 6)
    error ('randtest2: Too many input arguments')
  end
  if (nargout > 4)
    error ('randtest2: Too many output arguments')
  end
  
  % Set defaults
  if ( (nargin < 3) || isempty (paired) )
    paired = false;
  end
  if ( (nargin < 4) || isempty (nreps) )
    nreps = 5000;
  end
  if ( (nargin > 4) && (~ isempty (func)) )
    if (iscell (func))
      if (ischar (func{1}))
        % Convert character string of a function name to a function handle
        func = str2func (func{1});
      else
        func = func{1};
      end
      args = func(2:end);
      func = @(varargin) func (varargin{:}, args{:});
    elseif (ischar (func))
      % Convert character string of a function name to a function handle
      func = str2func (func);
    elseif (isa (func, 'function_handle'))
      % Do nothing
    else
      error ('randtest2: FUNC must be a function name or function handle')
    end
  else
    func = localfunc.wass_stat;
  end
  if ( (nargin > 5) && (~ isempty (seed)) )
    % Set random seed
    rand ('seed', seed);
  end

  % Remove NaN data values
  x(isnan (x(:,1)), :) = [];
  y(isnan (y(:,1)), :) = [];
  if ( any (isinf (x(:,1))) || any (isinf (y(:,1))) )
    error ('randtest2: x and y cannot not contain inf values')
  end

  % Get size of the data
  szx = size (x);
  szy = size (y);

  % Evaluate definition of the sampling units (e.g. clusters) of x and y
  if (szx(2) > 1)
    gx = x(:, 2);
    [ux, ix] = unique_stable (gx); % ix are indices where ux appear last in gx
    nx = numel (ux);
  else
    gx = (1 : szx(1))';
    ux = gx;
    ix = gx;
    nx = szx(1);
  end
  if (szy(2) > 1)
    gy = y(:, 2);
    [uy, iy] = unique_stable (gy);  % iy are indices where uy appear last in gy
    ny = numel (uy);
  else
    switch paired
      case false
        gy = (nx + 1 : nx + szy(1))';
        iy = gy - nx;
      case true
        gy = gx;
        iy = gy;
    end
    uy = gy;
    ny = szy(1);
  end

  switch paired

    case false

      % Error checking
      if ( any (ismember (ux, uy)) )
        error (cat (2, 'randtest2: sampling units defined in GX and GY', ...
                       ' must be unique when PAIRED is false'))
      end
      nz = nx + ny;
      if ( ~ all (ismember ((1 + nx) : nz, uy)) )
        error (cat (2, 'randtest2: sampling units defined in GY must', ...
                       ' continue numbering from GX when PAIRED is false'))
      end
      if ( (any ((ux ~= (1 : nx)'))) || ...
           (any ((uy ~= (1 : ny)' + nx))) )
        error (cat (2, 'randtest2: sampling units must be defined as', ...
                       ' consecutive positive integers (1, 2, 3, etc.)'))
      end
      if ( (any (ix ~= cumsum (accumarray (gx, 1)))) || ...
           (any (iy ~= cumsum (accumarray (gy - nx, 1)))) )
        error ('randtest2: clustered observations must be grouped together')
      end

      % Compute test statistic on the original data
      stat = func (x(:, 1), y(:, 1));

      % Create cell array of x and y samples
      z = cat (1, mat2cell (x(:, 1), accumarray (gx, 1)),...
                  mat2cell (y(:, 1), accumarray (gy - nx, 1)));
      gz = cat (1, gx, gy);

      % Create permutations or randomized samples
      if (factorial (nz) < nreps)
        I = perms (1:nz).';                 % For exact permutation test
        nreps = factorial (nz);
      else 
        [jnk, I] = sort (rand (nz, nreps)); % For approximate randomization test
      end
      X = arrayfun (@(i) z(I(i, :)), 1 : nx, 'UniformOutput', false);
      X = [X{:}]';
      Y = arrayfun (@(i) z(I(i, :)), (1 + nx) : nz, 'UniformOutput', false);
      Y = [Y{:}]';

      % Check if we can vectorize function evaluations
      VECTORIZED = ~ any (diff (accumarray (gz, 1)));

    case true

      % Error checking
      if (nx ~= ny)
        error (cat (2, 'randtest2: X and Y must have the same number of', ...
                       ' sampling units when PAIRED is true'))
      end
      if (any (ux ~= uy))
        error (cat (2, 'randtest2: GX and GY must use the same IDs for', ...
                       ' sampling units when PAIRED is true'))
      end
      if ( (~ all (ismember (1 : nx, ux))) || ...
           (~ all (ismember (1 : ny, uy))) )
        error (cat (2, 'randtest2: sampling units must be defined as', ...
                       ' consecutive positive integers (1, 2, 3, etc.)'))
      end
      nz = nx;
      if ( (any (ix ~= cumsum (accumarray (gx, 1)))) || ...
           (any (iy ~= cumsum (accumarray (gy, 1)))) )
        error ('randtest2: clustered observations must be grouped together')
      end
      if (any (gx ~= gy))
        error (cat (2, 'randtest2: cluster definitions must be identical', ...
                       ' for x and y when PAIRED is true'))
      end

      % Compute test statistic on the original data
      stat = func (x(:, 1), y(:, 1));

      % Create cell array of x and y samples
      z = cat (1, mat2cell (x(:, 1), accumarray (gx, 1))', ...
                  mat2cell (y(:, 1), accumarray (gy, 1))');

      % Create permutations or perform randomization
      if (2^nz < nreps)
        I = (dec2bin (0 : 2^nz - 1).' - '0') + 1; % For exact permutation test
        nreps = 2^nz;
      else 
        I = (rand (nz, nreps) > 0.5) + 1; % For approximate randomization test
      end
      X = arrayfun (@(i) z(I(i, :), i), 1 : nz, 'UniformOutput', false);
      X = [X{:}]';
      Y = arrayfun (@(i) z(3 - (I(i, :)), i), 1 : nz, 'UniformOutput', false);
      Y = [Y{:}]';

      % Check if we can vectorize function evaluations
      VECTORIZED = and (~ any (diff (accumarray (gx, 1))), ...
                        ~ any (diff (accumarray (gy, 1))));

  end

  % Perform function evaluations
  if VECTORIZED
    X = reshape (vertcat (X{:}), [], nreps);
    Y = reshape (vertcat (Y{:}), [], nreps);
    STATS = func (X, Y);
  else
    if (PARALLEL)
      if (ISOCTAVE)
        STATS = pararrayfun (inf, @(b) func (vertcat (X{:,b}), ...
                                             vertcat (Y{:,b})), 1:nreps);
      else
        STATS = zeros (1, nreps);
        parfor b = 1:nreps
          STATS(b) = func (vertcat (X{:,b}), vertcat (Y{:,b}))
        end
      end
    else
      STATS = arrayfun (@(b) func (vertcat (X{:,b}), ...
                                   vertcat (Y{:,b})), 1:nreps);
    end
  end

  % Calculate two-tailed p-value(s) by linear interpolation
  [x, jnk, P] = bootcdf (abs (STATS), true);
  res_lim = 1 / nreps;
  if (abs (stat) < x(1))
    pval = interp1 (x, P, abs (stat), 'linear', 1);
  else
    pval = interp1 (x, P, abs (stat), 'linear', res_lim);
  end
  if (nargout > 2)
    % Compute minimum false positive risk
    fpr = pval2fpr (pval);
  end

end

%--------------------------------------------------------------------------

function W = wass_stat (x, y)

  % Vectorized function to compute the Wasserstein metric

  % Get sample sizes
  % The number of columns in x and y are assumed to be the same
  [nx, ncols] = size (x);
  ny = size (y, 1);

  % Compute the difference in the areas under the empirical cumulative
  % distribution functions of x and y
  z = cat (1, x, y);
  [zs, I] = sort (z);
  D = cat (1, zs(2:end,:) - zs(1:end-1,:), zeros (1, ncols));
  E = cumsum (I <= nx) / nx;
  F = cumsum (I > nx) / ny;
  W = sum (D .* abs (E - F));

end

%--------------------------------------------------------------------------

% FUNCTION THAT RETURNS UNIQUE VALUES IN THE ORDER THAT THEY FIRST APPEAR

function [U, IA, IC] = unique_stable (A, varargin)

  % Subfunction used for backwards compatibility

  % Error checking
  if any (ismember (varargin, {'first', 'last', 'sorted', 'stable'}))
    error ('unique_stable: the only option available is ''rows''')
  end
  if (iscell (A) && ismember ('rows', varargin))
    error ('unique_stable: ''rows'' option not supported for cell arrays')
  end

  % Flatten A to a column vector if 'rows' option is not specified
  if (~ ismember ('rows', varargin))
    A = A(:);
  end

  % Obtain sorted unique values
  [u, ia, ic] = unique (A, 'last', varargin{:});

  % Sort index of last occurence of unique values as they first appear
  IA = sort (ia);

  % Get unique values in the order of appearace (a.k.a. 'stable')
  U = A(IA,:);

  % Create vector of numeric identifiers for unique values in A
  n = numel (IA);
  if iscell (A)
    IC = sum (cell2mat (arrayfun (@(i) i * ismember (A, U(i,:)), ...
                        (1:n), 'UniformOutput', false)), 2);
  elseif isnumeric (A)
    IC = sum (cell2mat (arrayfun (@(i) i * (all (bsxfun (@eq, A, U(i,:)), ...
                        2)), (1:n), 'UniformOutput', false)), 2);
  end

end

%--------------------------------------------------------------------------

% FUNCTION TO COMPUTE MINIMUM FALSE POSITIVE RISK (FPR)

function fpr = pval2fpr (p)

  % Subfunction to compute minimum false positive risk. These are calculated
  % from a Bayes factor based on the sampling distributions of the p-value and
  % that H0 and H1 have equal prior probabilities. This is called the Sellke-
  % Berger approach.
  % 
  % References:
  %  Held and Ott (2018) On p-Values and Bayes Factors. 
  %    Annu. Rev. of Stat. Appl. 5:393-419
  %  David Colquhoun (2019) The False Positive Risk: A Proposal 
  %    Concerning What to Do About p-Values, The American Statistician, 
  %    73:sup1, 192-201, DOI: 10.1080/00031305.2018.1529622 

  % Calculate minimum Bayes Factor (P(H0) / P(H1)) by the Sellke-Berger method 
  logp = min (log (p), -1);
  minBF = exp (1 + logp + log (-logp));

  % Calculate the false-positive risk from the minumum Bayes Factor
  L10 = 1 ./ minBF;      % Convert to Maximum Likelihood ratio L10 (P(H1)/P(H0))
  fpr = max (0, 1 ./ (1 + L10));  % Calculate minimum false positive risk 
  fpr(isnan(p)) = NaN; 

end

%--------------------------------------------------------------------------

%!demo
%!
%! % Mouse data from Table 2 (page 11) of Efron and Tibshirani (1993)
%! treatment = [94 197 16 38 99 141 23]';
%! control = [52 104 146 10 51 30 40 27 46]';
%!
%! % Randomization test comparing the distributions of observations from two
%! % independent samples (assuming i.i.d and exchangeability) using the
%! % Wasserstein metric
%! pval = randtest2 (control, treatment, false, 5000)
%!
%! % Randomization test comparing the difference in means between two
%! % independent samples (assuming i.i.d and exchangeability) 
%! pval = randtest2 (control, treatment, false, 5000, ...
%!                           @(x, y) mean (x) - mean (y))
%!
%! % Randomization test comparing the ratio of variances between two
%! % independent samples (assuming i.i.d and exchangeability). (Note that
%! % the log transformation is necessary to make the p-value two-tailed)
%! pval = randtest2 (control, treatment, false, 5000, ...
%!                           @(x, y) log (var (y) ./ var (x)))
%!

%!demo
%!
%! % Example data from: 
%! % https://www.biostat.wisc.edu/~kbroman/teaching/labstat/third/notes18.pdf
%! A = [117.3 100.1 94.5 135.5 92.9 118.9 144.8 103.9 103.8 153.6 163.1]';
%! B = [145.9 94.8 108 122.6 130.2 143.9 149.9 138.5 91.7 162.6 202.5]';
%!
%! % Randomization test comparing the distributions of observations from two
%! % paired or matching samples (assuming i.i.d and exchangeability) using the
%! % Wasserstein metric
%! pval = randtest2 (A, B, true, 5000)
%!
%! % Randomization test comparing the difference in means between two
%! % independent samples (assuming i.i.d and exchangeability) 
%! pval = randtest2 (A, B, true, 5000, @(A, B) mean (A) - mean (B))
%!
%! % Randomization test comparing the ratio of variances between two
%! % paired or matching samples (assuming i.i.d and exchangeability). (Note
%! % that the log transformation is necessary to make the p-value two-tailed)
%! pval = randtest2 (A, B,  true, 5000, @(A, B) log (var (A) ./ var (B)))
%!                           
%!

%!demo
%!
%! X = [21,26,33,22,18,25,26,24,21,25,35,28,32,36,38]';
%! GX = [1,1,1,1,2,2,2,2,2,2,3,3,3,3,3]';
%! Y = [26,34,27,38,44,34,45,38,31,41,34,35,38,46]';
%! GY = [4,4,4,5,5,5,5,5,6,6,6,6,6,6]';
%!
%! % Randomization test comparing the distributions of observations from two
%! % independent samples (assuming i.i.d) using the Wasserstein metric
%! pval = randtest2 (X, Y, false, 5000)
%!
%! % Randomization test comparing the distributions of clustered observations
%! % from two independent samples using the Wasserstein metric
%! pval = randtest2 ([X GX], [Y GY], false, 5000)
%!

%!demo
%!
%! X = [21,26,33,22,18,25,26,24,21,25,35,28,32,36,38]';
%! GX = [1,1,1,1,2,2,2,2,2,2,3,3,3,3,3]';
%! Y = [26,34,27,38,44,34,45,38,31,41,34,35,38,46,36]';
%! GY = [1,1,1,1,2,2,2,2,2,2,3,3,3,3,3]';
%!
%! % Randomization test comparing the distributions of observations from two
%! % paired or matched samples (assuming i.i.d) using the Wasserstein metric
%! pval = randtest2 (X, Y, true, 5000)
%!
%! % Randomization test comparing the distributions of clustered observations
%! % from two paired or matched using the Wasserstein metric
%! pval = randtest2 ([X GX], [Y GY], true, 5000)
%!

%!demo
%!
%! % Load example data from CSV file
%! data = csvread ('demo_data.csv');
%! trt = data(:,1); % Predictor: 0 = no treatment; 1 = treatment
%! grp = data(:,2); % Cluster IDs
%! val = data(:,3); % Values measured of the outcome
%! X = val(trt==0); GX = grp(trt==0);
%! Y = val(trt==1); GY = grp(trt==1);
%!
%! % Randomization test comparing the distributions of clustered observations
%! % from two independent samples using the Wasserstein metric
%! pval = randtest2([X, GX], [Y, GY], false)
%!


%!test
%!
%! % Test various capabilities of randtest2
%! X = randn (3,1);
%! Y = randn (3,1);
%! pval1 = randtest2 (X, Y);
%! pval2 = randtest2 (X, Y, false);
%! pval3 = randtest2 (X, Y, []);
%! randtest2 (X, Y, true);
%! randtest2 (X, Y, [], 500);
%! randtest2 (X, Y, [], []);
%! X = randn (9,1);
%! Y = randn (9,1);
%! pval5 = randtest2 (X, Y, false, 5000);
%! pval5 = randtest2 (X, Y, false, [], [], 1);
%! pval6 = randtest2 (X, Y, false, [], @(X, Y) mean (X) - mean (Y), 1);
%! pval7 = randtest2 (X, Y, false, [], @(A, B) log (var (A) ./ var (B)), 1);
