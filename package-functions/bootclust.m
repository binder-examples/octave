% Performs balanced bootstrap (or bootknife) resampling of clustered data and 
% calculates bootstrap bias, standard errors and confidence intervals.
%
% -- Function File: bootclust (DATA)
% -- Function File: bootclust (DATA, NBOOT)
% -- Function File: bootclust (DATA, NBOOT, BOOTFUN)
% -- Function File: bootclust ({D1, D2, ...}, NBOOT, BOOTFUN)
% -- Function File: bootclust (DATA, NBOOT, {BOOTFUN, ...})
% -- Function File: bootclust (DATA, NBOOT, BOOTFUN, ALPHA)
% -- Function File: bootclust (DATA, NBOOT, BOOTFUN, ALPHA, CLUSTID)
% -- Function File: bootclust (DATA, NBOOT, BOOTFUN, ALPHA, CLUSTSZ)
% -- Function File: bootclust (DATA, NBOOT, BOOTFUN, ALPHA, CLUSTID, LOO)
% -- Function File: bootclust (DATA, NBOOT, BOOTFUN, ALPHA, CLUSTID, LOO, SEED)
% -- Function File: STATS = bootclust (...)
% -- Function File: [STATS, BOOTSTAT] = bootclust (...)
%
%     'bootclust (DATA)' uses nonparametric balanced bootstrap resampling
%     to generate 1999 resamples from clusters of rows of the DATA (column
%     vector or matrix). By default, each rows is it's own cluster (i.e. no
%     clustering). The means of the resamples are then computed and the
%     following statistics are displayed:
%        - original: the original estimate(s) calculated by BOOTFUN and the DATA
%        - bias: bootstrap bias of the estimate(s)
%        - std_error: bootstrap standard error of the estimate(s)
%        - CI_lower: lower bound(s) of the 95% bootstrap confidence interval
%        - CI_upper: upper bound(s) of the 95% bootstrap confidence interval
%
%     'bootclust (DATA, NBOOT)' specifies the number of bootstrap resamples,
%     where NBOOT is a scalar, positive integer corresponding to the number
%     of bootstrap resamples. THe default value of NBOOT is the scalar: 1999.
%
%     'bootclust (DATA, NBOOT, BOOTFUN)' also specifies BOOTFUN: the function
%     calculated on the original sample and the bootstrap resamples. BOOTFUN
%     must be either a:
%       <> function handle or anonymous function,
%       <> string of function name, or
%       <> a cell array where the first cell is one of the above function
%          definitions and the remaining cells are (additional) input arguments 
%          to that function (other than the data arguments).
%        In all cases BOOTFUN must take DATA for the initial input argument(s).
%        BOOTFUN can return a scalar or any multidimensional numeric variable,
%        but the output will be reshaped as a column vector. BOOTFUN must
%        calculate a statistic representative of the finite data sample; it
%        should NOT be an estimate of a population parameter (unless they are
%        one of the same). If BOOTFUN is @mean or 'mean', narrowness bias of
%        the confidence intervals for single bootstrap are reduced by expanding
%        the probabilities of the percentiles using Student's t-distribution
%        [1]. By default, BOOTFUN is @mean.
%
%     'bootclust ({D1, D2, ...}, NBOOT, BOOTFUN)' resamples from the clusters
%     of rows of the data vectors D1, D2 etc and the resamples are passed onto
%     BOOTFUN as multiple data input arguments. All data vectors and matrices
%     (D1, D2 etc) must have the same number of rows.
%
%     'bootclust (DATA, NBOOT, BOOTFUN, ALPHA)', where ALPHA is numeric
%     and sets the lower and upper bounds of the confidence interval(s). The
%     value(s) of ALPHA must be between 0 and 1. ALPHA can either be:
%       <> scalar: To set the (nominal) central coverage of equal-tailed
%                  percentile confidence intervals to 100*(1-ALPHA)%.
%       <> vector: A pair of probabilities defining the (nominal) lower and
%                  upper percentiles of the confidence interval(s) as
%                  100*(ALPHA(1))% and 100*(ALPHA(2))% respectively. The
%                  percentiles are bias-corrected and accelerated (BCa) [2].
%        The default value of ALPHA is the vector: [.025, .975], for a 95%
%        BCa confidence interval.
%
%     'bootclust (DATA, NBOOT, BOOTFUN, ALPHA, CLUSTID)' also sets CLUSTID,
%     which are identifiers that define the grouping of the DATA rows for
%     cluster bootstrap case resampling. CLUSTID should be a column vector or
%     cell array with the same number of rows as the DATA. Rows in DATA with
%     the same CLUSTID value are treated as clusters of observations that are
%     resampled together.
%
%     'bootclust (DATA, NBOOT, BOOTFUN, ALPHA, CLUSTSZ)' groups consecutive
%     DATA rows into clusters of length CLUSTSZ. This is equivalent to block
%     bootstrap resampling. By default, CLUSTSZ is 1.
%
%     'bootclust (DATA, NBOOT, BOOTFUN, ALPHA, CLUSTID, LOO)' sets the
%     resampling method. If LOO is false, the resampling method used is
%     balanced bootstrap resampling. If LOO is true, the resampling method used
%     is balanced bootknife resampling [3]. Where N is the number of clusters,
%     bootknife cluster resampling involves creating leave-one-out jackknife
%     samples of size N - 1, and then drawing resamples of size N with
%     replacement from the jackknife samples, thereby incorporating Bessel's
%     correction into the resampling procedure. LOO must be a scalar logical
%     value. The default value of LOO is false.
%
%     'bootclust (DATA, NBOOT, BOOTFUN, ALPHA, CLUSTID, LOO, SEED)' initialises
%     the Mersenne Twister random number generator using an integer SEED value
%     so that bootclust results are reproducible.
%
%     'STATS = bootclust (...)' returns a structure with the following fields
%     (defined above): original, bias, std_error, CI_lower, CI_upper.
%
%     '[STATS, BOOTSTAT] = bootclust (...)' returns BOOTSTAT, a vector or matrix
%     of bootstrap statistics calculated over the bootstrap resamples.
%
%  REQUIREMENTS:
%    The function file boot.m (or better boot.mex) and bootcdf, which are
%    distributed with the statistics-resampling package.
%
%  BIBLIOGRAPHY:
%  [1] Hesterberg, Tim (2014), What Teachers Should Know about the 
%        Bootstrap: Resampling in the Undergraduate Statistics Curriculum, 
%        http://arxiv.org/abs/1411.5279
%  [2] Efron, and Tibshirani (1993) An Introduction to the Bootstrap. 
%        New York, NY: Chapman & Hall
%  [3] Hesterberg T.C. (2004) Unbiasing the Bootstrap—Bootknife Sampling 
%        vs. Smoothing; Proceedings of the Section on Statistics & the 
%        Environment. Alexandria, VA: American Statistical Association.
%
%  bootclust (version 2023.09.20)
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

function [stats, bootstat] = bootclust (x, nboot, bootfun, alpha, clustid, ...
                                        loo, seed)

  % Check if we are running Octave or Matlab
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Check the number of function arguments
  if (nargin < 1)
    error ('bootclust: DATA must be provided');
  end
  if (nargin > 7)
    error ('bootclust: Too many input arguments')
  end
  if (nargout > 2)
    error ('bootclust: Too many output arguments')
  end

  % NBOOT input argument
  if ((nargin < 2) || isempty (nboot))
    nboot = 1999;
  else
    if (~ isa (nboot, 'numeric'))
      error ('bootclust: NBOOT must be numeric');
    end
    if (numel (nboot) > 1)
      error ('bootclust: NBOOT cannot contain more than 1 value');
    end
    if (nboot ~= abs (fix (nboot)))
      error ('bootclust: NBOOT must contain positive integers');
    end    
  end
  if (~ all (size (nboot) == [1, 1]))
    error ('bootclust: NBOOT must be a scalar value')
  end

  % BOOTFUN input argument
  if ((nargin < 3) || isempty (bootfun))
    bootfun = @mean;
    bootfun_str = 'mean';
  else
    if (iscell (bootfun))
      if (ischar (bootfun{1}))
        % Convert character string of a function name to a function handle
        bootfun_str = bootfun{1};
        func = str2func (bootfun{1});
      else
        bootfun_str = func2str (bootfun{1});
        func = bootfun{1};
      end
      args = bootfun(2:end);
      bootfun = @(varargin) func (varargin{:}, args{:});
    elseif (ischar (bootfun))
      % Convert character string of a function name to a function handle
      bootfun_str = bootfun;
      bootfun = str2func (bootfun);
    elseif (isa (bootfun, 'function_handle'))
      bootfun_str = func2str (bootfun);
    else
      error ('bootclust: BOOTFUN must be a function name or function handle')
    end
  end

  % ALPHA input argument
  if ( (nargin < 4) || isempty (alpha) )
    alpha = [.025, .975];
  end
  nalpha = numel (alpha);
  if (~ isa (alpha, 'numeric') || (nalpha > 2))
    error (cat (2, 'bootclust: ALPHA must be a scalar (two-tailed', ...
                   'probability) or a vector (pair of probabilities)'))
  end
  if (size (alpha, 1) > 1)
    alpha = alpha.';
  end
  if (any (isnan (alpha)))
    error ('bootclust: ALPHA cannot contain NaN values');
  end
  if (any ((alpha < 0) | (alpha > 1)))
    error ('bootclust: Value(s) in ALPHA must be between 0 and 1');
  end
  if (nalpha > 1)
    % alpha is a pair of probabilities
    % Make sure probabilities are in the correct order
    if (alpha(1) > alpha(2) )
      error (cat (2, 'bootclust: The pair of probabilities must be', ...
                     ' in ascending numeric order'))
    end
    probs = alpha;
    alpha = 1 - probs(2) + probs(1);
  else
    probs = [alpha / 2 , 1 - alpha / 2];
  end

  % LOO input argument
  if ((nargin > 5) && ~ isempty (loo))
    if (~ islogical (loo))
      error ('bootclust: LOO must be a logical scalar value')
    end
  else
    loo = false;
  end
  
  % Initialise the random number generator with the SEED (if provided)
  if ( (nargin > 6) && (~ isempty (seed)) )
    boot (1, 1, true, seed);
  end

  % If DATA is a cell array of equal size colunmn vectors, convert the cell
  % array to a matrix and redefine bootfun to parse multiple input arguments
  if (iscell (x))
    szx = cellfun (@(x) size (x, 2), x);
    x = [x{:}];
    bootfun = @(x) col2args (bootfun, x, szx);
  else
    szx = size (x, 2);
  end

  % Determine properties of the DATA (x)
  [n, nvar] = size (x);
  if (n < 2)
    error ('bootclust: DATA must be numeric and contain > 1 row')
  end


  % Sort rows of CLUSTID and the DATA accordingly
  if ((nargin < 5) || isempty (clustid))
    clustid = (1 : n)';
  else
    if isscalar (clustid)
      % Group consecutive DATA rows into clusters of >= CLUSTID rows
      clustsz = clustid;
      if ( (~ isnumeric (clustsz)) || (clustsz ~= abs (clustsz)) || ...
                 (clustsz >= n) || (clustsz ~= fix (clustsz)) )
        error (cat (2, 'bootclust: CLUSTSZ must be a positive', ...
                       ' integer less than the number of DATA rows'))
      end
      nx = fix (n / clustsz);
      clustid = (nx + 1) * ones (n, 1);
      clustid(1:clustsz * nx, :) = reshape (ones (clustsz, 1) * (1:nx), [], 1);
      nx = clustid(end);
    end
    if ( any (size (clustid) ~= [n, 1]) )
      error (cat (2, 'bootclust: CLUSTID must be a column vector with', ...
                     ' the same number of rows as DATA'))
    end
    [clustid, idx] = sort (clustid);
    x = x(idx,:);
  end

  % Evaluate definition of the sampling units (e.g. clusters) of x 
  [ux, jnk, ic] = unique (clustid);
  nx = numel (ux);

  % Calculate the number of elements in the return value of bootfun and check
  % whether function evaluations can be vectorized
  T0 = bootfun (x);
  m = numel (T0);
  if (nvar > 1)
    M = cell2mat (cellfun (@(i) repmat (x(:, i), 1, 2), ... 
                 num2cell (1 : nvar), 'UniformOutput', false));
  else 
    M = repmat (x, 1, 2); 
  end
  if (any (szx > 1))
    VECTORIZED = false;
  else
    try
      chk = bootfun (M);
      if (all (size (chk) == [size(T0, 1), 2]) && all (chk == bootfun (x)))
        VECTORIZED = true;
      else
        VECTORIZED = false;
      end
    catch
      VECTORIZED = false;
    end
  end
  if (m > 1)
    % Vectorized along the dimension of the return values of bootfun so
    % reshape the output to be a column vector before proceeding with bootstrap
    if (size (T0, 2) > 1)
      bootfun = @(x) reshape (bootfun (x), [], 1);
      T0 = reshape (T0, [], 1);
      VECTORIZED = false;
    end
  end
  % Check if we can vectorize function evaluations
  if (any (diff (accumarray (ic, 1))))
    VECTORIZED = false;
  end

  % Convert x to a cell array of clusters
  x = mat2cell (x, accumarray (ic, 1));

  % Perform resampling of clusters
  bootsam = boot (nx, nboot, loo);
  X = arrayfun (@(i) x(bootsam(i, :)), 1 : nx, 'UniformOutput', false);
  X = [X{:}]';

  % Perform the function evaluations
  if VECTORIZED
    X = reshape (vertcat (X{:}), n, nboot * nvar);
    bootstat = bootfun (X);
  else
    bootstat = cell2mat (arrayfun (@(b) bootfun (vertcat (X{:,b})), ...
                                        1 : nboot, 'UniformOutput', false));
  end

  % Remove bootstrap statistics that contain NaN or inf
  ridx = any (or (isnan (bootstat), isinf (bootstat)) , 1);
  bootstat(:, ridx) = [];
  if (isempty (bootstat))
    error ('bootclust: BOOTFUN returned NaN or inf for all bootstrap resamples')
  end
  nboot = nboot - sum (ridx);

  % Bootstrap bias estimation
  bias = mean (bootstat, 2) - T0;

  % Bootstrap standard error
  se = std (bootstat, 0, 2);

  % Make corrections to the probabilities for the lower and upper bounds of the
  % confidence intervals.
  % First, if bootfun is the arithmetic meam, expand the probabilities of the 
  % percentiles for the confidence intervals using Student's t-distribution
  if (strcmpi (bootfun_str, 'mean'))
    probs = ExpandProbs (probs, nx - 1, loo);
  end
  % If requested, perform adjustments to the probabilities to correct for bias
  % and skewness
  switch (nalpha)
    case 1
      % No adjustements made
      probs = repmat (probs, m, 1);
    case 2
      % Create distribution functions
      stdnormcdf = @(x) 0.5 * (1 + erf (x / sqrt (2)));
      stdnorminv = @(p) sqrt (2) * erfinv (2 * p - 1);
      % Try using Jackknife resampling to calculate the acceleration constant (a)
      state = warning;
      if (ISOCTAVE)
        warning ('on', 'quiet');
      else
        warning ('off', 'all');
      end
      try
        jackfun = @(i) bootfun (vertcat (x{1 : nx ~= i, :}));
        % Evaluate bootfun on each jackknife resample
        T = cell2mat (arrayfun (jackfun, 1 : nx, 'UniformOutput', false));
        % Calculate empirical influence function
        U = (nx - 1) * bsxfun (@minus, mean (T, 2), T);
        a = sum (U.^3, 2) ./ (6 * sum (U.^2, 2) .^ 1.5);
      catch
        % Revert to bias-corrected (BC) bootstrap confidence intervals
        warning ('bootclust:jackfail', cat (2, 'BOOTFUN failed during', ... 
              ' jackknife calculations; acceleration constant set to 0.\n'))
        a = zeros (m, 1);
      end
      % Calculate the median bias correction constant (z0)
      z0 = stdnorminv (sum (bsxfun (@lt, bootstat, T0), 2) / nboot);
      if (~ all (isfinite (z0)))
        % Revert to percentile bootstrap confidence intervals
        warning ('bootclust:biasfail', ...
                 cat (2, 'Unable to calculate the bias correction', ...
                         ' constant; reverting to percentile intervals.\n'))
        z0 = zeros (m, 1);
        a = zeros (m, 1); 
      end
      % Calculate BCa or BC percentiles
      z = stdnorminv (probs);
      probs = stdnormcdf (bsxfun (@plus, z0, bsxfun (@plus, z0, z) ./ ...
                          (1 - (bsxfun (@times, a, bsxfun (@plus, z0, z))))));
  end

  % Intervals constructed from kernel density estimate of the bootstrap
  % statistics (with shrinkage correction)
  ci = nan (m, 2);
  for j = 1 : m
    try
      ci(j, :) = kdeinv (probs(j, :), bootstat(j, :), ...
                         se(j) * sqrt (1 / (nx - 1)), 1 - 1 / (nx - 1));
    catch
      % Linear interpolation (legacy)
      fprintf (strcat ('Note: Falling back to linear interpolation to', ...
                       ' calculate percentiles for interval pair %u\n'), j);
      [t1, cdf] = bootcdf (bootstat(j, :), true, 1);
      ci(j, 1) = interp1 (cdf, t1, probs(1), 'linear', min (t1));
      ci(j, 2) = interp1 (cdf, t1, probs(2), 'linear', max (t1));
    end
  end

  % Create STATS output structure
  stats = struct;
  stats.original = T0;
  stats.bias = bias;          % Bootstrap bias estimation
  stats.std_error = se;       % Bootstrap standard error
  stats.CI_lower = ci(:, 1);  % Lower percentile
  stats.CI_upper = ci(:, 2);  % Upper percentile

  % Print output if no output arguments are requested
  if (nargout == 0) 
    print_output (stats, nboot, nalpha, alpha, probs, m, bootfun_str, loo);
  else
    if (isempty (bootsam))
      [warnmsg, warnID] = lastwarn;
      if (ismember (warnID, {'bootclust:biasfail','bootclust:jackfail'}))
        warning ('bootclust:lastwarn', warnmsg);
      end
      lastwarn ('', '');
    end
  end

end

%--------------------------------------------------------------------------

function retval = col2args (func, x, szx)

  % Usage: retval = col2args (func, x, nvar)
  % col2args evaluates func on the columns of x. When nvar > 1, each of the
  % blocks of x are passed to func as a separate arguments. 

  % Extract columns of the matrix into a cell array
  [n, ncols] = size (x);
  xcell = mat2cell (x, n, ncols / sum (szx) * szx);

  % Evaluate column vectors as independent of arguments to bootfun
  retval = func (xcell{:});

end

%--------------------------------------------------------------------------

function X = kdeinv (P, Y, BW, CF)

  % Inverse of the cumulative density function (CDF) of a kernel density 
  % estimate (KDE)
  % 
  % The function returns X, the inverse CDF of the KDE of Y for the bandwidth
  % BW evaluated at the values in P. CF is a shrinkage factor for the variance
  % of the data in Y

  % Set defaults for optional input arguments
  if (nargin < 4)
    CF = 1;
  end

  % Create Normal CDF function
  pnorm = @(X, MU, SD) (0.5 * (1 + erf ((X - MU) / (SD * sqrt (2)))));

  % Calculate statistics of the data
  N = numel (Y);
  MU = mean (Y);

  % Apply shrinkage correction
  Y = ((Y - MU) * sqrt (CF)) + MU;

  % Set initial values of X0
  YS = sort (Y, 2);
  X0 = YS(fix ((N - 1) * P) + 1);

  % Perform root finding to get quantiles of the KDE at values of P
  findroot = @(X0, P) fzero (@(X) sum (pnorm (X - Y, 0, BW)) / N - P, X0);
  X = [-Inf, +Inf];
  for i = 1 : numel(P)
    if (~ ismember (P(i), [0, 1]))
      X(i) = findroot (X0(i), P(i));
    end
  end

end

%--------------------------------------------------------------------------

function PX = ExpandProbs (P, DF, LOO)

  % Modify ALPHA to adjust tail probabilities assuming that the kurtosis
  % of the sampling distribution scales with degrees of freedom like the
  % t-distribution. This is related in concept to ExpandProbs in the
  % R package 'resample':
  % www.rdocumentation.org/packages/resample/versions/0.6/topics/ExpandProbs

  % Get size of P
  sz = size (P);

  % Create required distribution functions
  stdnormcdf = @(X) 0.5 * (1 + erf (X / sqrt (2)));
  stdnorminv = @(P) sqrt (2) * erfinv (2 * P - 1);
  if ((exist ('betaincinv', 'builtin')) || (exist ('betaincinv', 'file')))
    studinv = @(P, DF) sign (P - 0.5) * ...
                sqrt ( DF ./ betaincinv (2 * min (P, 1 - P), DF / 2, 0.5) - DF);
  else
    % Earlier versions of Matlab do not have betaincinv
    % Instead, use betainv from the Statistics and Machine Learning Toolbox
    try 
      studinv = @(P, DF) sign (P - 0.5) * ...
                  sqrt ( DF ./ betainv (2 * min (P, 1 - P), DF / 2, 0.5) - DF);
    catch
      % Use the Normal distribution (i.e. do not expand probabilities) if
      % either betaincinv or betainv are not available
      studinv = @(P, DF) stdnorminv (P);
      warning ('bootclust:ExpandProbs', ...
          'Could not create studinv function; intervals will not be expanded.');
    end
  end
 
  % Calculate expanded probabilities
  if LOO
    PX = stdnormcdf (arrayfun (studinv, P, repmat (DF, sz)));
  else
    n = DF + 1;
    PX = stdnormcdf (sqrt (n / (n - 1)) * ...
                     arrayfun (studinv, P, repmat (DF, sz)));
  end

end

%--------------------------------------------------------------------------

function print_output (stats, nboot, nalpha, alpha, probs, m, bootfun_str, loo)

    fprintf (cat (2, '\nSummary of nonparametric cluster bootstrap', ...
                     ' estimates of bias and precision\n', ...
                     '*************************************************', ...
                     '*****************************\n\n'));
    fprintf ('Bootstrap settings: \n');
    fprintf (' Function: %s\n', bootfun_str);
    if loo
      fprintf (' Resampling method: Balanced, bootknife cluster resampling \n');
    else
      fprintf (' Resampling method: Balanced, bootstrap cluster resampling \n');
    end
    fprintf (' Number of resamples: %u \n', nboot(1));
    if (nalpha > 1)
      [jnk, warnID] = lastwarn;
      switch warnID
        case 'bootclust:biasfail'
          if (strcmpi (bootfun_str, 'mean'))
            fprintf (cat (2, ' Confidence interval (CI) type:', ...
                             ' Expanded percentile\n'));
          else
            fprintf (' Confidence interval (CI) type: Percentile\n');
          end
        case 'bootclust:jackfail'
          if (strcmpi (bootfun_str, 'mean'))
            fprintf (cat (2, ' Confidence interval (CI) type:', ...
                             ' Expanded bias-corrected (BC) \n'));
          else
            fprintf (cat (2, ' Confidence interval (CI) type:', ...
                             ' Bias-corrected (BC) \n'));
          end
        otherwise
          if (strcmpi (bootfun_str, 'mean'))
            fprintf (cat (2, ' Confidence interval (CI) type: Expanded', ...
                             ' bias-corrected and accelerated (BCa) \n'));
          else
            fprintf (cat (2, ' Confidence interval (CI) type: Bias-', ...
                             'corrected and accelerated (BCa) \n'));
          end
      end
    else
      if (strcmpi (bootfun_str, 'mean'))
        fprintf (cat (2, ' Confidence interval (CI) type: Expanded', ...
                         ' percentile (equal-tailed)\n'));
      else
        fprintf (cat (2, ' Confidence interval (CI) type: Percentile', ...
                         ' (equal-tailed)\n'));
      end
    end
    coverage = 100 * (1 - alpha);
    if (all (bsxfun (@eq, probs, probs(1, :))))
      fprintf (cat (2, ' Nominal coverage (and the percentiles used):', ...
                       ' %.3g%% (%.1f%%, %.1f%%)\n\n'), ...
                       coverage, 100 * probs(1,:));
    else
      fprintf (' Nominal coverage: %.3g%%\n\n', coverage);
    end
    fprintf ('Bootstrap Statistics: \n');
    fprintf (cat (2, ' original     bias         std_error    CI_lower', ...
                     '     CI_upper  \n'));
    for i = 1 : m
      fprintf (cat (2, ' %#-+10.4g   %#-+10.4g   %#-+10.4g   %#-+10.4g', ...
                     '   %#-+10.4g \n'), [stats.original(i), stats.bias(i), ...
                     stats.std_error(i), stats.CI_lower(i), stats.CI_upper(i)]);
    end
    fprintf ('\n');

end

%--------------------------------------------------------------------------

%!demo
%!
%! % Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! % 95% expanded BCa bootstrap confidence intervals for the mean
%! bootclust (data, 1999, @mean);

%!demo
%!
%! % Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%! clustid = {'a';'a';'b';'b';'a';'c';'c';'d';'e';'e';'e';'f';'f'; ...
%!            'g';'g';'g';'h';'h';'i';'i';'j';'j';'k';'l';'m';'m'};
%!
%! % 95% expanded BCa bootstrap confidence intervals for the mean with
%! % cluster resampling
%! bootclust (data, 1999, @mean, [0.025,0.975], clustid);

%!demo
%!
%! % Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! % 90% equal-tailed percentile bootstrap confidence intervals for
%! % the variance
%! bootclust (data, 1999, {@var, 1}, 0.1);

%!demo
%!
%! % Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%! clustid = {'a';'a';'b';'b';'a';'c';'c';'d';'e';'e';'e';'f';'f'; ...
%!            'g';'g';'g';'h';'h';'i';'i';'j';'j';'k';'l';'m';'m'};
%!
%! % 90% equal-tailed percentile bootstrap confidence intervals for
%! % the variance
%! bootclust (data, 1999, {@var, 1}, 0.1, clustid);

%!demo
%!
%! % Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%!
%! % 90% BCa bootstrap confidence intervals for the variance
%! bootclust (data, 1999, {@var, 1}, [0.05 0.95]);

%!demo
%!
%! % Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41].';
%! clustid = {'a';'a';'b';'b';'a';'c';'c';'d';'e';'e';'e';'f';'f'; ...
%!            'g';'g';'g';'h';'h';'i';'i';'j';'j';'k';'l';'m';'m'};
%!
%! % 90% BCa bootstrap confidence intervals for the variance
%! bootclust (data, 1999, {@var, 1}, [0.05 0.95], clustid);

%!demo
%!
%! % Input dataset
%! y = randn (20,1); x = randn (20,1); X = [ones(20,1), x];
%!
%! % 90% BCa confidence interval for regression coefficients 
%! bootclust ({y,X}, 1999, @(y,X) X\y, [0.05 0.95]); % Could also use @regress

%!demo
%!
%! % Input dataset
%! y = randn (20,1); x = randn (20,1); X = [ones(20,1), x];
%! clustid = [1;1;1;1;2;2;2;3;3;3;3;4;4;4;4;4;5;5;5;6];
%!
%! % 90% BCa confidence interval for regression coefficients 
%! bootclust ({y,X}, 1999, @(y,X) X\y, [0.05 0.95], clustid);

%!demo
%!
%! % Input bivariate dataset
%! x = [576 635 558 578 666 580 555 661 651 605 653 575 545 572 594].';
%! y = [3.39 3.3 2.81 3.03 3.44 3.07 3 3.43 ...
%!      3.36 3.13 3.12 2.74 2.76 2.88 2.96].';
%! clustid = [1;1;3;1;1;2;2;2;2;3;1;3;3;3;2];
%!
%! % 95% BCa bootstrap confidence intervals for the correlation coefficient
%! bootclust ({x, y}, 1999, @cor, [], clustid);
%!
%! % Please be patient, the calculations will be completed soon...

%!test
%! % Test for errors when using different functionalities of bootclust
%! y = randn (20,1); 
%! clustid = [1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;3;3;3;3;3];
%! stats = bootclust (y, 1999, @mean);
%! stats = bootclust (y, 1999, 'mean');
%! stats = bootclust (y, 1999, {@var,1});
%! stats = bootclust (y, 1999, {'var',1});
%! stats = bootclust (y, 1999, @mean, [], 4);
%! stats = bootclust (y, 1999, @mean, [], clustid);
%! stats = bootclust (y, 1999, {'var',1}, [], clustid);
%! stats = bootclust (y, 1999, {'var',1}, [], clustid, true);
%! stats = bootclust (y, 1999, {@var,1}, [], clustid, true, 1);
%! stats = bootclust (y, 1999, @mean, .1, clustid, true);
%! stats = bootclust (y, 1999, @mean, .1, clustid, true, 1);
%! stats = bootclust (y, 1999, @mean, [.05,.95], clustid, true);
%! stats = bootclust (y, 1999, @mean, [.05,.95], clustid, true, 1);
%! stats = bootclust (y(1:5), 1999, @mean, .1);
%! stats = bootclust (y(1:5), 1999, @mean, [.05,.95]);
%! Y = randn (20); 
%! clustid = [1;1;1;1;1;1;1;1;1;1;2;2;2;2;2;3;3;3;3;3];
%! stats = bootclust (Y, 1999, @mean);
%! stats = bootclust (Y, 1999, 'mean');
%! stats = bootclust (Y, 1999, {@var, 1});
%! stats = bootclust (Y, 1999, {'var',1});
%! stats = bootclust (Y, 1999, @mean, [], clustid);
%! stats = bootclust (Y, 1999, {'var',1}, [], clustid);
%! stats = bootclust (Y, 1999, {@var,1}, [], clustid, true);
%! stats = bootclust (Y, 1999, {@var,1}, [], clustid, true, 1);
%! stats = bootclust (Y, 1999, @mean, .1, clustid, true);
%! stats = bootclust (Y, 1999, @mean, .1, clustid, true, 1);
%! stats = bootclust (Y, 1999, @mean, [.05,.95], clustid, true);
%! stats = bootclust (Y, 1999, @mean, [.05,.95], clustid, true, 1);
%! stats = bootclust (Y(1:5,:), 1999, @mean, .1);
%! stats = bootclust (Y(1:5,:), 1999, @mean, [.05,.95]);
%! y = randn (20,1); x = randn (20,1); X = [ones(20,1), x];
%! stats = bootclust ({x,y}, 1999, @cor);
%! stats = bootclust ({x,y}, 1999, @cor, [], clustid);
%! stats = bootclust ({y,x}, 1999, @(y,x) pinv(x)*y); % Could use @regress
%! stats = bootclust ({y,X}, 1999, @(y,X) pinv(X)*y);
%! stats = bootclust ({y,X}, 1999, @(y,X) pinv(X)*y, [], clustid);
%! stats = bootclust ({y,X}, 1999, @(y,X) pinv(X)*y, [], clustid, true);
%! stats = bootclust ({y,X}, 1999, @(y,X) pinv(X)*y, [], clustid, true, 1);
%! stats = bootclust ({y,X}, 1999, @(y,X) pinv(X)*y, [.05,.95], clustid);

%!test
%! % Air conditioning failure times in Table 1.2 of Davison A.C. and
%! % Hinkley D.V (1997) Bootstrap Methods And Their Application. (Cambridge
%! % University Press)
%! x = [3, 5, 7, 18, 43, 85, 91, 98, 100, 130, 230, 487]';
%!
%! % Nonparametric 95% expanded percentile confidence intervals (equal-tailed)
%! % Balanced bootknife resampling
%! % Example 5.4 percentile intervals are 43.9 - 192.1
%! % Note that the intervals calculated below are wider because the narrowness
%! % bias was removed by expanding the probabilities of the percentiles using
%! % Student's t-distribution
%! stats = bootclust(x,1999,@mean,0.05,[],true,1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (stats.original, 108.0833333333333, 1e-08);
%!   assert (stats.bias, -2.842170943040401e-14, 1e-08);
%!   assert (stats.std_error, 38.21311346451331, 1e-08);
%!   assert (stats.CI_lower, 37.63178001821489, 1e-08);
%!   assert (stats.CI_upper, 200.846030428085, 1e-08);
%! end
%!
%! % Nonparametric 95% expanded BCa confidence intervals
%! % Balanced bootknife resampling
%! % Example 5.8 BCa intervals are 55.33 - 243.5
%! % Note that the intervals calculated below are wider because the narrowness
%! % bias was removed by expanding the probabilities of the percentiles using
%! % Student's t-distribution
%! stats = bootclust(x,1999,@mean,[0.025,0.975],[],true,1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (stats.original, 108.0833333333333, 1e-08);
%!   assert (stats.bias, -2.842170943040401e-14, 1e-08);
%!   assert (stats.std_error, 38.21311346451331, 1e-08);
%!   assert (stats.CI_lower, 49.68734232019933, 1e-08);
%!   assert (stats.CI_upper, 232.2854912518225, 1e-08);
%! end
%!
%! % Exact intervals based on an exponential model are 65.9 - 209.2
%! % (Example 2.11)

%!test
%! % Spatial test data from Table 14.1 of Efron and Tibshirani (1993)
%! % An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! % Probability 57 (Springer)
%! A = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!      0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! % Nonparametric 90% equal-tailed percentile confidence intervals
%! % Balanced bootknife resampling
%! % Table 14.2 percentile intervals are 100.8 - 233.9
%! stats = bootclust(A,1999,{@var,1},0.1,[],true,1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (stats.original, 171.534023668639, 1e-08);
%!   assert (stats.bias, -7.305657266503118, 1e-08);
%!   assert (stats.std_error, 43.17157379039285, 1e-08);
%!   assert (stats.CI_lower, 95.33589724383222, 1e-08);
%!   assert (stats.CI_upper, 237.1866652820803, 1e-08);
%! end
%!
%! % Nonparametric 90% BCa confidence intervals
%! % Balanced bootknife resampling
%! % Table 14.2 BCa intervals are 115.8 - 259.6
%! stats = bootclust(A,1999,{@var,1},[0.05 0.95],[],true,1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (stats.original, 171.534023668639, 1e-08);
%!   assert (stats.bias, -7.305657266503118, 1e-08);
%!   assert (stats.std_error, 43.17157379039285, 1e-08);
%!   assert (stats.CI_lower, 113.2858617321027, 1e-08);
%!   assert (stats.CI_upper, 264.0328613673329, 1e-08);
%! end
%!
%! % Exact intervals based on normal theory are 118.4 - 305.2 (Table 14.2)
%! % Note that all of the bootknife intervals are slightly wider than the
%! % nonparametric intervals in Table 14.2 because the bootknife (rather than
%! % standard bootstrap) resampling used here reduces small sample bias

%!test
%! % Law school data from Table 3.1 of Efron and Tibshirani (1993)
%! % An Introduction to the Bootstrap in Monographs on Statistics and Applied 
%! % Probability 57 (Springer)
%! LSAT = [576 635 558 578 666 580 555 661 651 605 653 575 545 572 594]';
%! GPA = [3.39 3.3 2.81 3.03 3.44 3.07 3 3.43 ...
%!        3.36 3.13 3.12 2.74 2.76 2.88 2.96]';
%!
%! % Nonparametric 90% equal-tailed percentile confidence intervals
%! % Balanced bootstrap resampling
%! % Percentile intervals on page 266 are 0.524 - 0.928
%! stats = bootclust({LSAT,GPA},1999,@cor,0.1,[],false,1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (stats.original, 0.7763744912894071, 1e-08);
%!   assert (stats.bias, -0.007614791775856333, 1e-08);
%!   assert (stats.std_error, 0.1355245146889644, 1e-08);
%!   assert (stats.CI_lower, 0.5146251204001586, 1e-08);
%!   assert (stats.CI_upper, 0.9531054945982934, 1e-08);
%! end
%!
%! % Nonparametric 90% BCa confidence intervals
%! % Balanced bootstrap resampling
%! % BCa intervals on page 266 are 0.410 - 0.923
%! stats = bootclust({LSAT,GPA},1999,@cor,[0.05 0.95],[],false,1);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (stats.original, 0.7763744912894071, 1e-08);
%!   assert (stats.bias, -0.007614791775856333, 1e-08);
%!   assert (stats.std_error, 0.1355245146889644, 1e-08);
%!   assert (stats.CI_lower, 0.4177828971982108, 1e-08);
%!   assert (stats.CI_upper, 0.9238952404759969, 1e-08);
%! end
