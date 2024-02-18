% Performs wild bootstrap and calculates bootstrap-t confidence intervals and 
% p-values for the mean, or the regression coefficients from a linear model.
%
% -- Function File: bootwild (y)
% -- Function File: bootwild (y, X)
% -- Function File: bootwild (y, X, CLUSTID)
% -- Function File: bootwild (y, X, BLOCKSZ)
% -- Function File: bootwild (y, X, ..., NBOOT)
% -- Function File: bootwild (y, X, ..., NBOOT, ALPHA)
% -- Function File: bootwild (y, X, ..., NBOOT, ALPHA, SEED)
% -- Function File: bootwild (y, X, ..., NBOOT, ALPHA, SEED, L)
% -- Function File: STATS = bootwild (y, ...)
% -- Function File: [STATS, BOOTSTAT] = bootwild (y, ...)
% -- Function File: [STATS, BOOTSTAT, BOOTSSE] = bootwild (y, ...)
% -- Function File: [STATS, BOOTSTAT, BOOTSSE, BOOTFIT] = bootwild (y, ...)
%
%     'bootwild (y)' performs a null hypothesis significance test for the
%     mean of y being equal to 0. This function implements wild bootstrap-t
%     resampling of Webb's 6-point distribution of the residuals, and computes
%     p-values after imposing the null hypothesis (H0) [1-4]. The following
%     statistics are printed to the standard output:
%        - original: the mean of the data vector y
%        - std_err: heteroscedasticity-consistent standard error(s)
%        - CI_lower: lower bound(s) of the 95% bootstrap-t confidence interval
%        - CI_upper: upper bound(s) of the 95% bootstrap-t confidence interval
%        - tstat: Student's t-statistic
%        - pval: two-tailed p-value(s) for the parameter(s) being equal to 0
%        - fpr: minimum false positive risk for the corresponding p-value
%          By default, the confidence intervals are symmetric, two-sided
%          bootstrap-t confidence intervals. The p-values are computed
%          following both of the guidelines by Hall and Wilson [5]. The minimum
%          false positive risk (FPR) is computed according to the Sellke-Berger
%          approach as described in [6,7].
%
%     'bootwild (y, X)' also specifies the design matrix (X) for least squares
%     regression of y on X. X should be a column vector or matrix the same
%     number of rows as y. If the X input argument is empty, the default for X
%     is a column of ones (i.e. intercept only) and thus the statistic computed
%     reduces to the mean (as above). The statistics calculated and returned in
%     the output then relate to the coefficients from the regression of y on X.
%
%     'bootwild (y, X, CLUSTID)' specifies a vector or cell array of numbers
%     or strings respectively to be used as cluster labels or identifiers.
%     Rows in y (and X) with the same CLUSTID value are treated as clusters
%     with dependent errors. Rows of y (and X) assigned to a particular
%     cluster will have identical resampling during wild bootstrap. If empty
%     (default), no clustered resampling is performed and all errors are
%     treated as independent. The standard errors computed are cluster robust.
%
%     'bootwild (y, X, BLOCKSZ)' specifies a scalar, which sets the block size
%     for bootstrapping when the residuals have serial dependence. Identical
%     resampling occurs within each (consecutive) block of length BLOCKSZ
%     during wild bootstrap. Rows of y (and X) within the same block are
%     treated as having dependent errors. If empty (default), no block
%     resampling is performed and all errors are treated as independent.
%     The standard errors computed are cluster robust.
%
%     'bootwild (y, X, ..., NBOOT)' specifies the number of bootstrap resamples,
%     where NBOOT must be a positive integer. If empty, the default value of
%     NBOOT is 1999.
%
%     'bootwild (y, X, ..., NBOOT, ALPHA)' is numeric and sets the lower and
%     upper bounds of the confidence interval(s). The value(s) of ALPHA must
%     be between 0 and 1. ALPHA can either be:
%        o scalar: To set the (nominal) central coverage of SYMMETRIC
%                  bootstrap-t confidence interval(s) to 100*(1-ALPHA)%.
%                  For example, 0.05 for a 95% confidence interval.
%        o vector: A pair of probabilities defining the (nominal) lower and
%                  upper bounds of ASYMMETRIC bootstrap-t confidence interval(s)
%                  as 100*(ALPHA(1))% and 100*(ALPHA(2))% respectively. For
%                  example, [.025, .975] for a 95% confidence interval.
%        The default value of ALPHA is the scalar: 0.05, for a symmetric 95%
%        bootstrap-t confidence interval.
%
%     'bootwild (y, X, ..., NBOOT, ALPHA, SEED)' initialises the Mersenne
%     Twister random number generator using an integer SEED value so that
%     'bootwild' results are reproducible.
%
%     'bootwild (y, X, ..., NBOOT, ALPHA, SEED, L)' multiplies the regression
%     coefficients by the hypothesis matrix L. If L is not provided or is empty,
%     it will assume the default value of 1 (i.e. no change to the design). 
%
%     'STATS = bootwild (...) returns a structure with the following fields:
%     original, std_err, CI_lower, CI_upper, tstat, pval, fpr and the sum-of-
%     squared error (sse).
%
%     '[STATS, BOOTSTAT] = bootwild (...)  also returns a vector (or matrix) of
%     bootstrap statistics (BOOTSTAT) calculated over the bootstrap resamples
%     (before studentization).
%
%     '[STATS, BOOTSTAT, BOOTSSE] = bootwild (...)  also returns a vector
%     containing the sum-of-squared error for the fit on each bootstrap 
%     resample.
%
%     '[STATS, BOOTSTAT, BOOTSSE, BOOTFIT] = bootwild (...)  also returns an
%     N-by-NBOOT matrix containing the N fitted values for each of the NBOOT
%     bootstrap resamples.
%
%  Bibliography:
%  [1] Wu (1986). Jackknife, bootstrap and other resampling methods in
%        regression analysis (with discussions). Ann Stat.. 14: 1261–1350. 
%  [2] Cameron, Gelbach and Miller (2008) Bootstrap-based Improvements for
%        Inference with Clustered Errors. Rev Econ Stat. 90(3), 414-427
%  [3] Webb (2023) Reworking wild bootstrap-based inference for clustered
%        errors. Can J Econ. https://doi.org/10.1111/caje.12661
%  [4] Cameron and Miller (2015) A Practitioner’s Guide to Cluster-Robust
%        Inference. J Hum Resour. 50(2):317-372
%  [5] Hall and Wilson (1991) Two Guidelines for Bootstrap Hypothesis Testing.
%        Biometrics, 47(2), 757-762
%  [6] Colquhoun (2019) The False Positive Risk: A Proposal Concerning What
%        to Do About p-Values, Am Stat. 73:sup1, 192-201
%  [7] Sellke, Bayarri and Berger (2001) Calibration of p-values for Testing
%        Precise Null Hypotheses. Am Stat. 55(1), 62-71
%
%  bootwild (version 2023.07.05)
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


function [stats, bootstat, bootsse, bootfit] = bootwild (y, X, ...
                                          dep, nboot, alpha, seed, L, ISOCTAVE)

  % Check the number of function arguments
  if (nargin < 1)
    error ('bootwild: y must be provided')
  end
  if (nargin > 8)
    error ('bootwild: Too many input arguments')
  end
  if (nargout > 4)
    error ('bootwild: Too many output arguments')
  end

  % Check if running in Octave (else assume Matlab)
  if (nargin < 8)
    info = ver; 
    ISOCTAVE = any (ismember ({info.Name}, 'Octave'));
  else
    if (~ islogical (ISOCTAVE))
      error ('bootwild: ISOCTAVE must be a logical scalar.')
    end
  end

  % Calculate the length of y
  if (nargin < 1)
    error ('bootwild: DATA must be provided')
  end
  sz = size (y);
  if ( (sz(1) < 2) || (sz (2) > 1) )
    error ('bootwild: y must be a column vector')
  end
  n = numel (y);

  % Evaluate the design matrix
  if ( (nargin < 2) || (isempty (X)) )
    X = ones (n, 1);
  elseif (size (X, 1) ~= n)
    error ('bootwild: X must have the same number of rows as y')
  end

  % Remove rows of the data whose outcome or value of any predictor is NaN or Inf
  excl = any ([isnan([y, X]), isinf([y, X])], 2);
  y(excl) = [];
  X(excl, :) = [];
  n = n - sum (excl);

  % Calculate the number of parameters
  k = size (X, 2);
  if ((k == 1) && (all (X == 1)) )
    p = 1;
    L = 1;
  else
    % Evaluate hypothesis matrix (L)
    if ( (nargin < 8) || isempty (L) )
      % If L is not provided, set L to 1
      L = 1;
      p = k;
    else
      % Calculate the number of parameters
      [m, p] = size (L);
      if (m ~= k)
        error (cat (2, 'bootwild: the number rows in L must be the same', ...
                       ' as the number of columns in X'))
      end
    end
  end

  % Evaluate cluster IDs or block size
  if ( (nargin > 2) && (~ isempty (dep)) )
    if (isscalar (dep))
      % Prepare for wild block bootstrap
      blocksz = dep;
      G = fix (n / blocksz);
      IC = (G + 1) * ones (n, 1);
      IC(1 : blocksz * G, :) = reshape (ones (blocksz, 1) * (1:G), [], 1);
      G = IC(end);
      method = 'block ';
    else
      % Prepare for wild cluster bootstrap
      dep(excl) = [];
      clustid = dep;
      if (bsxfun (@ne, size (clustid), sz))
        error ('bootwild: clustid must be the same size as y')
      end
      [C, IA, IC] = unique (clustid);
      G = numel (C); % Number of clusters
      method = 'cluster ';
    end
    UC = unique (IC);
    clusters = cell2mat (cellfun (@(i) IC == UC(i), ...
                                  num2cell (1:G), 'UniformOutput', false));
  else
    G = n;
    IC = [];
    clusters = [];
    method = '';
  end

  % Evaluate number of bootstrap resamples
  if ( (nargin < 4) || (isempty (nboot)) )
    nboot = 1999;
  else
    if (~ isa (nboot, 'numeric'))
      error ('bootwild: NBOOT must be numeric')
    end
    if (numel (nboot) > 1)
      error ('bootwild: NBOOT must be scalar')
    end
    if (nboot ~= abs (fix (nboot)))
      error ('bootwild: NBOOT must be a positive integer')
    end
    if (nboot < 999)
      error ('bootwild: NBOOT must be >= 999')
    end
  end
  % Compute resolution limit of the p-values as determined by resampling
  % with nboot resamples
  res_lim = 1 / (nboot + 1);

  % Evaluate alpha
  if ( (nargin < 5) || isempty (alpha) )
    alpha = 0.05;
    nalpha = 1;
  else
    nalpha = numel (alpha);
    if (~ isa (alpha, 'numeric') || (nalpha > 2))
      error (cat (2, 'bootwild: ALPHA must be a scalar (two-tailed', ...
                     ' probability) or a vector (pair of probabilities)'))
    end
    if (size (alpha, 1) > 1)
      alpha = alpha.';
    end
    if (any ((alpha < 0) | (alpha > 1)))
      error ('bootwild: Value(s) in ALPHA must be between 0 and 1')
    end
    if (nalpha > 1)
      % alpha is a pair of probabilities
      % Make sure probabilities are in the correct order
      if (alpha(1) > alpha(2) )
        error (cat (2, 'bootwild: The pair of probabilities must be in', ...
                       ' ascending numeric order'))
      end
    end
  end

  % Set random seed
  if ( (nargin > 5) && (~ isempty (seed)) )
    rand ('seed', seed);
  end

  % Compute unscaled covariance matrix
  [Q, R] = qr (X, 0);        % Economy-sized QR decomposition
  if ISOCTAVE
    ucov = chol2inv (R);     % Instead of pinv (X' * X). Faster
  else
    ucov = inv (R' * R);     % Instead of pinv (X' * X). Slower
  end

  % Create least squares anonymous function for bootstrap
  bootfun = @(y) lmfit (X, y, ucov, clusters, L, ISOCTAVE);

  % Calculate estimate(s)
  S = bootfun (y);
  original = S.b;
  std_err = S.se;
  sse = S.sse;
  t = original ./ std_err;

  % Wild bootstrap resampling (Webb's 6-point distribution)
  s = sign (rand (G, nboot) - 0.5) .* ...
      sqrt (0.5 * (fix (rand (G, nboot) * 3) + 1));
  if (~ isempty (IC))
    s = s(IC, :);  % Enforce clustering/blocking
  end
  yf = X * (X \ y);
  r = y - yf;
  rs = bsxfun (@times, r, s);
  Y = bsxfun (@plus, yf, rs);

  % Compute bootstap statistics
  bootout = cell2mat (cellfun (bootfun, num2cell (Y, 1), ...
                              'UniformOutput', false));
  bootstat = [bootout.b];
  bootse = [bootout.se];
  bootsse = [bootout.sse];
  bootfit = [bootout.fit];

  % Studentize the bootstrap statistics and compute two-tailed confidence
  % intervals and p-values following both guidelines described in Hall and
  % Wilson (1991) Biometrics, 47(2), 757-762
  T = bsxfun (@minus, bootstat, original) ./ bootse;
  unstable = any (or (lt (bootse, eps), isnan (T)), 2);
  ci = nan (p, 2);
  pval = nan (p, 1);
  if (any (~ isnan (alpha)))
    for j = 1:p
      [x, F, P] = bootcdf (abs (T(j,:)), true, 1);
      if (abs (t(j)) < x(1))
        pval(j) = interp1 (x, P, abs (t(j)), 'linear', 1);
      else
        pval(j) = interp1 (x, P, abs (t(j)), 'linear', res_lim);
      end
      if ( (~ isnan (std_err(j))) && (~ unstable(j)) )
        switch nalpha
          case 1
            ci(j, 1) = original(j) - std_err(j) * ...
                                    interp1 (F, x, 1 - alpha, 'linear', max (x));
            ci(j, 2) = original(j) + std_err(j) * ...
                                    interp1 (F, x, 1 - alpha, 'linear', max (x));
          case 2
            [x, F] = bootcdf (T(j,:), true, 1);
            ci(j, 1) = original(j) - std_err(j) * ...
                                    interp1 (F, x, alpha(2), 'linear', max (x));
            ci(j, 2) = original(j) - std_err(j) * ...
                                    interp1 (F, x, alpha(1), 'linear', min (x));
        end
      end
    end
  end

  % Compute minimum false positive risk
  fpr = pval2fpr (pval);

  % Prepare output arguments
  stats = struct;
  stats.original = original;
  stats.std_err = std_err;
  stats.CI_lower = ci(:,1);
  stats.CI_upper = ci(:,2);
  stats.tstat = t;
  stats.pval = pval;
  stats.fpr = fpr;
  stats.sse = sse;

  % Print output if no output arguments are requested
  if (nargout == 0) 
    print_output (stats, nboot, alpha, p, L, method);
  end

end

%--------------------------------------------------------------------------

% FUNCTION TO FIT THE LINEAR MODEL

function S = lmfit (X, y, ucov, clusters, L, ISOCTAVE)

  % Get model coefficients by solving the linear equation by matrix arithmetic

  % Solve linear equation to minimize least squares and compute the
  % regression coefficients (b) 
  b = X \ y;                 % Instead of inv (X' * X) * (X' * y);

  % Calculate heteroscedasticity-consistent (HC) or cluster robust (CR) standard 
  % errors for the regression coefficients. When the number of observations
  % equals the number of clusters, the calculations for CR reduce to HC.
  % References: 
  %   Long and Ervin (2000) Am. Stat, 54(3), 217-224
  %   Cameron, Gelbach and Miller (2008) Rev Econ Stat. 90(3), 414-427
  %   MacKinnon & Webb (2020) QED Working Paper Number 1421
  yf = X * b;
  u = y - yf;
  if ( (nargin < 3) || isempty (clusters) )
    % Heteroscedasticity-Consistent (HC0) standard errors
    meat = X' * diag (u.^2) * X;
  else
    % Cluster Robust (CR0) standard errors
    Sigma = cellfun (@(g) X(g,:)' * u(g) * u(g)' * X(g,:), ...
                  num2cell (clusters, 1), 'UniformOutput', false);
    meat = sum (cat (3, Sigma{:}), 3);
  end
  vcov = ucov * meat * ucov;
  S = struct; 
  if ( (nargin < 4) || isempty (L) )
    S.b = b;
    S.se = sqrt (max (diag (vcov), 0));
  else
    S.b = L' * b;
    S.se = sqrt (max (diag (L' * vcov * L), 0));
  end
  S.sse = sum (u.^2);
  S.fit = yf;

end

%--------------------------------------------------------------------------

% FUNCTION TO COMPUTE FALSE POSITIVE RISK (FPR)

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

% FUNCTION TO PRINT OUTPUT

function print_output (stats, nboot, alpha, p, L, method)

    fprintf (cat (2, '\nSummary of wild bootstrap null hypothesis', ...
                     ' significance tests for linear models\n', ...
                     '*******************************************', ...
                     '************************************\n\n'));
    fprintf ('Bootstrap settings: \n');
    if ( (numel(L) > 1) || (L ~= 1) )
      fprintf (' Function: L'' * X \\ y\n');
    else
      fprintf (' Function: X \\ y\n');
    end
    fprintf (' Resampling method: Wild %sbootstrap-t\n', method)
    fprintf (' Number of resamples: %u \n', nboot)
    fprintf (' Standard error calculations:');
    if (isempty (method))
      fprintf (' Heteroscedasticity-Consistent (HC0)\n');
    else
      fprintf (' Cluster Robust (CR0)\n');
    end
    nalpha = numel (alpha);
    if (nalpha > 1)
      % prob is a vector of probabilities
      fprintf (cat (2, ' Confidence interval (CI) type: Asymmetric', ...
                       ' bootstrap-t interval\n'));
      coverage = 100 * abs (alpha(2) - alpha(1));
      fprintf (cat (2, ' Nominal coverage (and the percentiles used):', ...
                       ' %.3g%% (%.1f%%, %.1f%%)\n', coverage, 100 * alpha(:)));
    else
      % prob is a two-tailed probability
      fprintf (cat (2, ' Confidence interval (CI) type: Symmetric', ...
                       ' bootstrap-t interval\n'));
      coverage = 100 * (1 - alpha);
      fprintf (' Nominal central coverage: %.3g%%\n', coverage);
    end
    fprintf (' Null value (H0) used for hypothesis testing (p-values): 0 \n')
    fprintf ('\nTest Statistics: \n');
    fprintf (cat (2, ' original     std_err      CI_lower     CI_upper', ...
                     '     t-stat      p-val     FPR\n'));
    for j = 1:p
      fprintf (cat (2, ' %#-+10.4g   %#-+10.4g   %#-+10.4g', ...
                       '   %#-+10.4g   %#-+9.3g'), ...
                       [stats.original(j), stats.std_err(j), ...
                        stats.CI_lower(j), stats.CI_upper(j), stats.tstat(j)]);
      if (stats.pval(j) <= 0.001)
        fprintf ('   <.001');
      elseif (stats.pval(j) < 0.9995)
        fprintf ('    .%03u', round (stats.pval(j) * 1e+03));
      elseif (isnan (stats.pval(j)))
        fprintf ('     NaN');
      else
        fprintf ('   1.000');
      end
      if (stats.fpr(j) <= 0.001)
        fprintf ('   <.001\n');
      elseif (stats.fpr(j) < 0.9995)
        fprintf ('    .%03u\n', round (stats.fpr(j) * 1e+03));
      elseif (isnan (stats.fpr(j)))
        fprintf ('     NaN\n');
      else
        fprintf ('   1.000\n');
      end
    end
    fprintf ('\n');

end

%--------------------------------------------------------------------------

%!demo
%!
%! % Input univariate dataset
%! heights = [183, 192, 182, 183, 177, 185, 188, 188, 182, 185].';
%!
%! % Test statistics and p-values (H0 = 0)
%! bootwild (heights);
%!
%! % Please be patient, the calculations will be completed soon...

%!demo
%!
%! % Input bivariate dataset
%! X = [ones(43,1),...
%!     [01,02,03,04,05,06,07,08,09,10,11,...
%!      12,13,14,15,16,17,18,19,20,21,22,...
%!      23,25,26,27,28,29,30,31,32,33,34,...
%!      35,36,37,38,39,40,41,42,43,44]'];
%! y = [188.0,170.0,189.0,163.0,183.0,171.0,185.0,168.0,173.0,183.0,173.0,...
%!     173.0,175.0,178.0,183.0,192.4,178.0,173.0,174.0,183.0,188.0,180.0,...
%!     168.0,170.0,178.0,182.0,180.0,183.0,178.0,182.0,188.0,175.0,179.0,...
%!     183.0,192.0,182.0,183.0,177.0,185.0,188.0,188.0,182.0,185.0]';
%!
%! % Compute test statistics and p-values
%! bootwild (y, X);
%!
%! % Please be patient, the calculations will be completed soon...

%!test
%! % Test if the mean is equal to a population value of 181.5 (one-tailed test)
%!
%! % Input univariate dataset
%! H0 = 181.5;
%! heights = [183, 192, 182, 183, 177, 185, 188, 188, 182, 185].';
%!
%! % Compute test statistics and p-values
%! [stats,bootstat] = bootwild(heights);
%! stats = bootwild(heights-H0);
%! stats = bootwild(heights-H0,ones(10,1));
%! stats = bootwild(heights-H0,[],2);
%! stats = bootwild(heights-H0,[],[1;1;2;2;3;3;4;4;5;5]);
%! stats = bootwild(heights-H0,[],[],1999);
%! stats = bootwild(heights-H0,[],[],[],0.05);
%! stats = bootwild(heights-H0,[],[],[],[0.025,0.975]);
%! stats = bootwild(heights-H0,[],[],[],[],1);
%! stats = bootwild(heights-H0,[],[],[],[],[]);
%! stats = bootwild(heights-H0,[],[],[],[],[],1);
%! stats = bootwild(heights-H0,[],[],[],[],[],[]);
%! stats = bootwild(heights-H0,[],[],[],0.05,1);
%! assert (stats.original, 3.0, 1e-06);
%! assert (stats.std_err, 1.242980289465605, 1e-06);
%! assert (stats.CI_lower, 0.1338910532454438, 1e-06);
%! assert (stats.CI_upper, 5.866108946754555, 1e-06);
%! assert (stats.tstat, 2.413553960127389, 1e-06);
%! assert (stats.pval, 0.04363142391272781, 1e-06);
%! assert (stats.fpr, 0.2708502563156392, 1e-06);
%! % ttest gives a p-value of 0.0478
%! stats = bootwild(heights-H0,[],[],[],[0.025,0.975],1);
%! assert (stats.original, 3.0, 1e-06);
%! assert (stats.std_err, 1.242980289465605, 1e-06);
%! assert (stats.CI_lower, 0.01340070207731392, 1e-06);
%! assert (stats.CI_upper, 5.801890495593613, 1e-06);
%! assert (stats.tstat, 2.413553960127389, 1e-06);
%! assert (stats.pval, 0.04363142391272781, 1e-06);
%! assert (stats.fpr, 0.2708502563156392, 1e-06);
%! stats = bootwild(heights-H0,[],2,[],0.05,1);
%! assert (stats.original, 3.0, 1e-06);
%! assert (stats.std_err, 1.240967364599086, 1e-06);
%! assert (stats.CI_lower, -2.816060435625108, 1e-06);
%! assert (stats.CI_upper, 8.816060435625108, 1e-06);
%! assert (stats.tstat, 2.41746889207614, 1e-06);
%! assert (stats.pval, 0.1297336562664251, 1e-06);
%! assert (stats.fpr, 0.4186764774953166, 1e-06);
%! stats = bootwild(heights-H0,[],[1;1;2;2;3;3;4;4;5;5],[],0.05,1);
%! assert (stats.original, 3.0, 1e-06);
%! assert (stats.std_err, 1.240967364599086, 1e-06);
%! assert (stats.CI_lower, -2.816060435625108, 1e-06);
%! assert (stats.CI_upper, 8.816060435625108, 1e-06);
%! assert (stats.tstat, 2.41746889207614, 1e-06);
%! assert (stats.pval, 0.1297336562664251, 1e-06);
%! assert (stats.fpr, 0.4186764774953166, 1e-06);

%!test
%! % Test if the regression coefficients equal 0
%!
%! % Input bivariate dataset
%! X = [ones(43,1),...
%!     [01,02,03,04,05,06,07,08,09,10,11,...
%!      12,13,14,15,16,17,18,19,20,21,22,...
%!      23,25,26,27,28,29,30,31,32,33,34,...
%!      35,36,37,38,39,40,41,42,43,44]'];
%! y = [188.0,170.0,189.0,163.0,183.0,171.0,185.0,168.0,173.0,183.0,173.0,...
%!     173.0,175.0,178.0,183.0,192.4,178.0,173.0,174.0,183.0,188.0,180.0,...
%!     168.0,170.0,178.0,182.0,180.0,183.0,178.0,182.0,188.0,175.0,179.0,...
%!     183.0,192.0,182.0,183.0,177.0,185.0,188.0,188.0,182.0,185.0]';
%!
%! % Compute test statistics and p-values
%! [stats,bootstat] = bootwild(y,X);
%! stats = bootwild(y,X);
%! stats = bootwild(y,X,3);
%! stats = bootwild(y,X,[],1999);
%! stats = bootwild(y,X,[],[],0.05);
%! stats = bootwild(y,X,[],[],[0.025,0.975]);
%! stats = bootwild(y,X,[],[],[],1);
%! stats = bootwild(y,X,[],[],[],[]);
%! stats = bootwild(y,X,[],[],[],[],1);
%! stats = bootwild(y,X,[],[],[],[],[]);
%! stats = bootwild(y,X,[],[],0.05,1);
%! assert (stats.original(2), 0.1904211996616223, 1e-06);
%! assert (stats.std_err(2), 0.08261019852213342, 1e-06);
%! assert (stats.CI_lower(2), -0.0008180267754335224, 1e-06);
%! assert (stats.CI_upper(2), 0.3816604260986781, 1e-06);
%! assert (stats.tstat(2), 2.305056797685863, 1e-06);
%! assert (stats.pval(2), 0.05133180571394391, 1e-06);
%! assert (stats.fpr(2), 0.2929561493533854, 1e-06);
%! % fitlm gives a CI of [0.0333, 0.34753] and a p-value of 0.018743
%! stats = bootwild(y,X,[],[],[0.025,0.975],1);
%! assert (stats.original(2), 0.1904211996616223, 1e-06);
%! assert (stats.std_err(2), 0.08261019852213342, 1e-06);
%! assert (stats.CI_lower(2), -0.0008180267754335224, 1e-06);
%! assert (stats.CI_upper(2), 0.3831445183919875, 1e-06);
%! assert (stats.tstat(2), 2.305056797685863, 1e-06);
%! assert (stats.pval(2), 0.05133180571394391, 1e-06);
%! assert (stats.fpr(2), 0.2929561493533854, 1e-06);
%! stats = bootwild(y,X,3,[],0.05,1);
%! assert (stats.original(2), 0.1904211996616223, 1e-06);
%! assert (stats.std_err(2), 0.07170701459665534, 1e-06);
%! assert (stats.CI_lower(2), -0.0242385109696315, 1e-06);
%! assert (stats.CI_upper(2), 0.4050809102928761, 1e-06);
%! assert (stats.tstat(2), 2.655544938423697, 1e-06);
%! assert (stats.pval(2), 0.07716525989964551, 1e-06);
%! assert (stats.fpr(2), 0.3495327983695262, 1e-06);
