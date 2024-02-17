% Performs resampling under the null hypothesis and computes p-values for
% (multiple) comparisons among independent samples in a one-way layout.
%
% -- Function File: boot1way (DATA, GROUP)5
% -- Function File: boot1way (..., NAME, VALUE)
% -- Function File: boot1way (..., 'bootfun', BOOTFUN)
% -- Function File: boot1way (..., 'nboot', NBOOT)
% -- Function File: boot1way (..., 'ref', REF)
% -- Function File: boot1way (..., 'alpha', ALPHA)
% -- Function File: boot1way (..., 'Options', PAROPT)
% -- Function File: PVAL = boot1way (DATA, GROUP, ...)
% -- Function File: [PVAL, C] = boot1way (DATA, GROUP, ...)
% -- Function File: [PVAL, C, STATS] = boot1way (DATA, GROUP, ...)
% -- Function File: [...] = boot1way (..., 'display', DISPLAYOPT)
%
%     'boot1way (DATA, GROUP)' performs a bootstrap version of a randomization
%     test [1] for comparing independent samples of data in a one-way layout.
%     Pairwise multiple comparison tests are computed by the single-step
%     maximum absolute t-statistic (maxT) procedure, which controls the family-
%     wise error rate (FWER) in a manner analagous to the Tukey-Kramer Honest
%     Significance Difference test. The results are displayed as a pretty table
%     and the differences between groups are plotted along with the symmetric
%     95% bootstrap-t confidence intervals (CI). The colours of the markers and
%     error bars depend on the value of the multiplicity-adjusted p-values:
%     red if p < .05, or blue if p > .05. All of the p-values reported represent
%     the outcome of two-tailed tests. DATA must be a numeric column vector or
%     matrix, where categorization of the DATA rows is achieved by labels in
%     GROUP. GROUP must be a vector or cell array with the same number of
%     rows as DATA.  
%
%     boot1way can take a number of optional parameters as NAME-VALUE pairs:
%
%     'boot1way (..., 'bootfun', BOOTFUN)' also specifies BOOTFUN: the function
%     calculated on the original sample and the bootstrap resamples. BOOTFUN
%     must be either a:
%        o function handle or anonymous function,
%        o string of function name, or
%        o a cell array where the first cell is one of the above function
%          definitions and the remaining cells are (additional) input arguments 
%          to that function (other than the data arguments).
%        In all cases, BOOTFUN must take DATA for the initial input argument(s).
%        BOOTFUN must calculate a statistic representative of the finite data
%        sample; it should NOT be an estimate of a population parameter (unless
%        they are one of the same). By default, BOOTFUN is @mean. If a robust
%        alternative to the mean is required, set BOOTFUN to 'robust' to
%        implement a smoothed version of the median (a.k.a. @smoothmedian). 
%
%     'boot1way (..., 'nboot', NBOOT)' is a scalar or a vector of upto two
%     positive integers indicating the number of resamples for the first
%     (bootstrap) and second (bootknife) levels of iterated resampling. If NBOOT
%     is a scalar value, or if NBOOT(2) is set to 0, then standard errors are
%     calculated either without resampling (if BOOTFUN @mean) or using Tukey's
%     jackknife. This implementation of jackknife requires the Statistics
%     package/toolbox. The default value of NBOOT is the vector: [999,99].
%
%     'boot1way (..., 'ref', REF)' sets the GROUP to use as the reference group
%     for the multiple comparison tests. If REF is a recognised member of GROUP,
%     then the maxT procedure for treatment versus reference controls the
%     family-wise error rate (FWER) in a manner analagous to Dunnett's multiple
%     comparison tests.
%
%     'boot1way (..., 'alpha', ALPHA)' specifies the two-tailed significance
%     level for CI coverage. The default value of ALPHA is 0.05 for 95%
%     confidence intervals.
%
%     'boot1way (..., 'Options', PAROPT)' specifies options that govern if
%     and how to perform bootstrap iterations using multiple processors (if the
%     Parallel Computing Toolbox or Octave Parallel package is available). This
%     argument is a structure with the following recognised fields:
%        o 'UseParallel':  If true, use parallel processes to accelerate
%                          bootstrap computations on multicore machines,
%                          specifically non-vectorized function evaluations,
%                          double bootstrap resampling and jackknife function
%                          evaluations. Default is false for serial computation.
%                          In MATLAB, the default is true if a parallel pool
%                          has already been started. 
%        o 'nproc':        nproc sets the number of parallel processes
%
%     'PVAL = boot1way (DATA, GROUP, ...)' returns the p-value(s) for the 
%     (multiple) two-tailed test(s). Note that the p-value(s) returned are
%     already adjusted to control the family-wise, type I error rate and 
%     truncated at the resolution limit determined by the number of bootstrap
%     replicates, specifically 1/NBOOT(1)  
%
%     '[PVAL, C] = boot1way (DATA, GROUP, ...)' also returns a 9 column matrix
%     that summarises multiple comparison test results. The columns of C are:
%       - column 1:  test GROUP number
%       - column 2:  reference GROUP number
%       - column 3:  value of BOOTFUN evaluated for the test GROUP
%       - column 4:  value of BOOTFUN evaluated for the reference GROUP
%       - column 5:  the difference between the groups (column 3 minus column 4)
%       - column 6:  LOWER bound of the 100*(1-ALPHA)% bootstrap-t CI
%       - column 7:  UPPER bound of the 100*(1-ALPHA)% bootstrap-t CI
%       - column 8:  t-ratio
%       - column 9:  multiplicity-adjusted p-value
%       - column 10: minimum false positive risk for the p-value
%
%     '[PVAL, C, STATS] = boot1way (DATA, GROUP, ...)' also returns a structure 
%     containing additional statistics. The stats structure contains the 
%     following fields:
%
%       gnames   - group names used in the GROUP input argument. The index of 
%                  gnames corresponds to the numbers used to identify GROUPs
%                  in columns 1 and 2 of the output argument C
%       ref      - index of the reference group
%       groups   - group index and BOOTFUN for each group with sample size,
%                  standard error and CI, which start to overlap at a
%                  multiplicity-adjusted p-value of approximately 0.05
%       Var      - weighted mean (pooled) sampling variance
%       nboot    - number of bootstrap resamples (1st and 2nd resampling layers)
%       alpha    - two-tailed significance level for the CI reported in C.
%
%     '[PVAL, C, STATS, BOOTSTAT] = boot1way (DATA, GROUP, ...)' also returns
%     the maximum test statistic computed for each bootstrap resample
%
%     '[...] = boot1way (..., 'display', DISPLAYOPT)' a logical value (true 
%      or false) used to specify whether to display the results and plot the
%      graph in addition to creating the output arguments. The default is true.
%
%  BIBLIOGRAPHY:
%  [1] Efron, and Tibshirani (1993) An Introduction to the Bootstrap. 
%        New York, NY: Chapman & Hall
%
%  boot1way (version 2023.10.04)
%  Bootstrap tests for comparing independent groups in a one-way layout
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


function [pval, c, stats, Q] = boot1way (data, group, varargin)

  % Evaluate the number of function arguments
  if (nargin < 2)
    error (cat (2, 'boot1way usage: ''boot1way (DATA, GROUP, varargin)'';', ...
                   ' atleast 2 input arguments required'))
  end

  % Store local functions in a stucture for parallel processes
  localfunc = struct ('maxstat',@maxstat, ...
                      'bootcdf',@bootcdf);

  % Check if running in Octave (else assume Matlab)
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Apply defaults
  bootfun = 'mean';
  nboot = [999,99];
  ref = [];
  alpha = 0.05;
  DisplayOpt = true;
  paropt = struct;
  paropt.UseParallel = false;
  if (ISOCTAVE)
    paropt.nproc = nproc;
  else
    paropt.nproc = feature ('numcores');
  end

  % Fetch extra input arguments
  argin3 = varargin;
  narg = numel (argin3);
  if (narg > 1)
    while ischar (argin3{end-1})
      if (strcmpi (argin3{end-1},'bootfun'))
        bootfun = argin3{end};
      elseif (strcmpi (argin3{end-1},'nboot'))
        nboot = argin3{end};
      elseif (strcmpi (argin3{end-1},'ref'))
        ref = argin3{end};
      elseif (any (strcmpi (argin3{end-1},{'Options','Option'})))
        paropt = argin3{end};
      elseif (strcmpi (argin3{end-1},'alpha'))
        alpha = argin3{end};
      elseif (any (strcmpi (argin3{end-1},{'DisplayOpt','Display'})))
        DisplayOpt = argin3{end};
      else
        error ('boot1way: Unrecognised input argument to boot1way')
      end
      argin3 = {argin3{1:end-2}};
      narg = numel (argin3);
      if (narg < 1)
        break
      end
    end
  end

  % Error checking
  % Check and process boot1way input arguments
  nvar = size (data,2);
  if (nargin < 2)
    error ('boot1way: Requires atleast two input arguments');
  end
  if (ischar (group))
    group = cellstr (group);
  end
  if ((size (group, 1)>1) && (size (data, 1) ~= size (group, 1)))
    error ('boot1way: DATA and GROUP must have the same number of rows')
  end
  if (iscell (group))
    if (~ iscellstr (group))
      group = cell2mat (group);
    end
  end
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
    if (strcmpi (bootfun_str, 'robust'))
      bootfun = @smoothmedian;
    else
      bootfun = str2func (bootfun);
    end
  elseif (isa (bootfun, 'function_handle'))
    bootfun_str = func2str (bootfun);
  else
    error ('boot1way: BOOTFUN must be a function name or function handle')
  end
  if (~ isa (nboot, 'numeric'))
    error ('boot1way: NBOOT must be numeric');
  end
  if (any (nboot ~= abs (fix (nboot))))
    error ('boot1way: NBOOT must contain positive integers')
  end
  if (numel (nboot) > 2)
    error ('boot1way: A vector for NBOOT cannot have length > 2')
  elseif (numel (nboot) < 2)
    nboot = cat (2, nboot, 0);
  end
  if (nboot(1) < 999)
    error ('boot1way: The minimum allowable value of NBOOT(1) is 999')
  end 
  if (nboot(2) == 0) && (nvar > 1)
    error (cat (2, 'boot1way: Jackknife currently only available for', ...
                   ' analysis of univariate data.'))
  end
  if ((nboot(2) == 0) && (~ strcmp (func2str (bootfun), 'mean')))
    if (~ exist ('jackknife','file'))
      if (ISOCTAVE)
        warning ('boot1way:jackfail', cat (2, '''jackknife'' function from', ...
                 ' statistics package not found. nboot(2) set to 100.'))
      else
        warning ('boot1way:jackfail',cat (2, '''jackknife'' function from', ...
                 ' statistics toolbox not found. nboot(2) set to 99.'))
      end
      nboot(2) = 99;
    end
  end
  
  % Error checking
  if (~ isempty (ref) && strcmpi (ref,'pairwise'))
    ref = [];
  end
  if (nargout > 4)
    error ('boot1way: Only supports up to 4 output arguments')
  end
  if (~ islogical (DisplayOpt) || (numel (DisplayOpt) > 1))
    error ('boot1way: The value DISPLAYOPT must be a logical scalar value')
  end

  % Data or group exclusion using NaN 
  if (isnumeric (group))
    if (any (isnan (group)))
      data(isnan (group),:) = [];
      group(isnan (group)) = [];
    end
  end
  if (any (any (isnan (data), 2)))
    group(any (isnan (data), 2)) = [];
    data(any (isnan (data), 2),:) = [];
  end

  % Assign non-zero numbers to group labels
  [gnames,jnk,g] = unique (group);
  clear jnk;
  gk = unique (g);
  k = numel (gk);
  if (k > 1)
    if (~ isempty (ref))
      if (isnumeric (ref))
        ref = gk(ismember (gnames, ref));
      else
        ref = gk(strcmp (gnames, ref));
      end
    end
  else
    error ('boot1way: GROUP must define atleast two groups')
  end
  N = numel (g);
  
  % If applicable, check we have parallel computing capabilities
  if (paropt.UseParallel)
    if (ISOCTAVE)
      software = pkg ('list');
      names = cellfun (@(S) S.name, software, 'UniformOutput', false);
      status = cellfun (@(S) S.loaded, software, 'UniformOutput', false);
      index = find (~ cellfun (@isempty, regexpi (names, '^parallel')));
      if (~ isempty (index))
        if (logical (status{index}))
          PARALLEL = true;
        else
          PARALLEL = false;
        end
      else
        PARALLEL = false;
      end
    else
      if (ismember ('Parallel Computing Toolbox', {info.Name}))
        PARALLEL = true;
      else
        PARALLEL = false;
      end
    end
  end
  
  % If applicable, setup a parallel pool (required for MATLAB)
  if (~ ISOCTAVE)
    % MATLAB
    if (paropt.UseParallel)
      % PARALLEL
      if (paropt.nproc > 0) 
        % MANUAL
        try 
          pool = gcp ('nocreate'); 
          if (isempty (pool))
            if (paropt.nproc > 1)
              % Start parallel pool with nproc workers
              parpool (paropt.nproc);
            else
              % Parallel pool is not running and nproc is 1 so run function
              % evaluations in serial
              paropt.UseParallel = false;
            end
          else
            if (pool.NumWorkers ~= paropt.nproc)
              % Check if number of workers matches nproc and correct it
              % accordingly if not
              delete (pool);
              if (paropt.nproc > 1)
                parpool (paropt.nproc);
              end
            end
          end
        catch
          % MATLAB Parallel Computing Toolbox is not installed
          warning ('boot1way:parallel', cat (2, 'Parallel Computing', ...
               ' Toolbox is not installed. Falling back to serial processing.'))
          paropt.UseParallel = false;
          paropt.nproc = 1;
        end
      else
        % AUTOMATIC
        try 
          pool = gcp ('nocreate'); 
          if (isempty (pool))
            % Parallel pool not running, start parallel pool using all
            % available workers
            parpool;
          else
            % Parallel pool is already running, set nproc to the
            % number of workers
            paropt.nproc = pool.NumWorkers;
          end
        catch
          % Parallel toolbox not installed, run function evaluations in serial
          paropt.UseParallel = false;
        end
      end
    end
  else
    if (paropt.UseParallel && (paropt.nproc > 1) && ~PARALLEL)
      if (ISOCTAVE)
        % OCTAVE Parallel Computing Package is not installed or loaded
        warning ('boot1way:parallel', cat (2, 'Parallel Computing Package', ...
         ' is not installed and/or loaded. Falling back to serial processing.'))
      else
        % MATLAB Parallel Computing Toolbox is not installed or loaded
        warning ('boot1way:parallel', cat (2, 'Parallel Computing Toolbox', ...
         ' is not installed and/or loaded. Falling back to serial processing.'))
      end
      paropt.UseParallel = false;
      paropt.nproc = 0;
    end
  end

  % Create maxstat anonymous function for bootstrap
  func = @(data) localfunc.maxstat (data, g, nboot(2), bootfun, ref, ISOCTAVE);

  % Perform resampling and calculate bootstrap statistics to estimate the 
  % sampling distribution under the null hypothesis.
  % LOO set to false for bootstrap (instead of bootknife) resampling.
  boot (1, 1, false, 1); % set random seed to make resampling deterministic
  if (paropt.UseParallel)
    [jnk, Q] = bootknife (data, nboot(1), func, NaN, [], paropt.nproc, ...
                          [], [], ISOCTAVE, true, false);
  else
    [jnk, Q] = bootknife (data, nboot(1), func, NaN, [], 0, ...
                          [], [], ISOCTAVE, true, false);
  end
  
  % Compute the estimate (theta) and it's pooled (weighted mean) sampling
  % variance 
  theta = zeros (k, 1);
  SE = zeros (k, 1);
  Var = zeros (k, 1);
  nk = zeros (size(gk));
  for j = 1:k
    if (nboot(2) == 0)
      nk(j) = sum (g == gk(j));
      if (strcmp (func2str (bootfun), 'mean'))
        theta(j) = mean (data(g == gk(j), :));
        % Quick analytical calculation for the standard error of the mean
        SE(j) = std (data(g == gk(j), :), 0) / sqrt (nk(j));
        if (j == 1); se_method = 'Calculated without resampling'; end
      else
        theta(j) = bootfun (data(g == gk(j), :));
        % If requested, compute unbiased estimates of the standard error using
        % jackknife resampling
        jackstat = jackknife (bootfun, data(g == gk(j), :));
        SE(j) = sqrt ((nk(j) - 1) / nk(j) * ...
                sum (((mean (jackstat) - jackstat)).^2));
        if (j == 1); se_method = 'Leave-one-out jackknife'; end
      end
    else
      % Compute unbiased estimate of the standard error by balanced bootknife
      % resampling. Bootknife resampling involves less computation than
      % Jackknife when sample sizes get larger
      theta(j) = bootfun (data(g == gk(j), :));
      nk(j) = sum (g == gk(j));
      bootout = bootknife (data(g == gk(j), :), [nboot(2), 0], bootfun, ...
                           NaN, [], 0, [], [], ISOCTAVE, false, true);
      SE(j) = bootout.std_error;
      if (j==1); se_method = 'Balanced, bootknife resampling'; end
    end
    Var(j) = ((nk(j) - 1) / (N - k)) * SE(j)^2;
  end
  if (any (SE == 0))
    error ('boot1way: Samples must have non-zero standard error')
  end
  if (any (isnan (SE)))
    error (cat (2, 'boot1way: Evaluating bootfun on the bootknife', ...
                   ' resamples created NaN values for the standard error'))
  end
  nk_bar = sum (nk.^2) ./ sum (nk); % weighted mean sample size
  Var = sum (Var .* nk / nk_bar);   % weighted pooled sampling variance
  %df = sum (nk) - k;                % degrees of freedom

  % Calculate weights to correct for unequal sample size  
  % when calculating standard error of the difference
  w = nk_bar ./ nk;

  % Prepare to make symmetrical bootstrap-t confidence intervals and
  % 2-tailed p-values, Create empirical distribution function
  [Q, F, P] = bootcdf (Q, true, 1);

  % Compute resolution limit of the p-values as determined by resampling
  % with nboot(1) resamples
  res_lim = 1 / (nboot(1) + 1);

  % Calculate p-values for comparisons adjusted to simultaneously
  % control the Family-Wise Error Rate (FWER)
  if (isempty (ref))
    % Single-step maxT procedure for pairwise comparisons is a resampling
    % version of Tukey-Kramer Honest Significant Difference (HSD) test
    A = tril (gk * ones (1, k), -1);
    B = ones (k, 1) * gk';
    M = [A(:) B(:)];
    ridx = (M(:,1) == 0);
    M(ridx, :) = [];
    n = size (M, 1);
    c = zeros (n, 10);
    c(:,1:2) = M;
    for i = 1:n
      c(i,3) = theta(c(i,1));
      c(i,4) = theta(c(i,2));
      c(i,5) = c(i,3) - c(i,4);
      SED = sqrt (Var * (w(c(i,1)) + w(c(i,2))));
      c(i,6) = c(i,5) - SED * interp1 (F, Q, 1 - alpha, 'linear', max (Q));
      c(i,7) = c(i,5) + SED * interp1 (F, Q, 1 - alpha, 'linear', max (Q));
      c(i,8) = c(i,5) / SED;
      if (c(i,8) < Q(1))
        c(i,9) = interp1 (Q, P, abs (c(i,8)), 'linear', 1);
      else
        c(i,9) = interp1 (Q, P, abs (c(i,8)), 'linear', res_lim);
      end
    end
    c(:,10) = pval2fpr (c(:,9));
  else
    % Single-step maxT procedure for treatment vs control comparisons is
    % a resampling version of Dunnett's test
    c = zeros (k, 10);
    c(:,2) = ref;
    c(:,4) = theta (ref);
    for j = 1:k
      c(j,1) = gk(j);
      c(j,3) = theta(c(j,1));
      c(j,5) = c(j,3) - c(j,4); 
      SED = sqrt (Var * (w(c(j,1)) + w(c(j,2))));
      c(j,6) = c(j,5) - SED * interp1 (F, Q, 1 - alpha, 'linear', max (Q));
      c(j,7) = c(j,5) + SED * interp1 (F, Q, 1 - alpha, 'linear', max (Q));
      c(j,8) = c(j,5) / SED;
      if (c(j,8) < Q(1))
        c(j,9) = interp1 (Q, P, abs (c(j,8)), 'linear', 1);
      else
        c(j,9) = interp1 (Q, P, abs (c(j,8)), 'linear', res_lim);
      end
    end
    c(:,10) = pval2fpr (c(:,9));
    c(ref,:) = [];
  end

  % Assign the calculated p-values to the return value
  pval = c(:,9);

  % Prepare stats output structure
  stats = struct;
  stats.gnames = gnames;
  stats.ref = ref;
  stats.groups = zeros (k,6);
  stats.groups = zeros (k,6);
  stats.groups(:,1) = gk;
  stats.groups(:,2) = theta;
  stats.groups(:,3) = nk;
  stats.groups(:,4) = SE;
  stats.groups(:,5) = theta - sqrt ((0.5 * (w + 1)) .* Var / 2) * ...
                              interp1 (F, Q, 1 - alpha, 'linear', max (Q));
  stats.groups(:,6) = theta + sqrt ((0.5 * (w + 1)) .* Var / 2) * ...
                              interp1 (F, Q, 1 - alpha, 'linear', max (Q));
  stats.Var = Var;
  stats.nboot = nboot;
  stats.alpha = alpha;

  % Print output and plot graph with confidence intervals if no output
  % arguments are requested
  cols = [1,2,5,8,9]; % columns in c that we want to print data for
  if ((nargout == 0) || (DisplayOpt == true))
    if (~ iscellstr (gnames))
      gnames = cellstr (num2str (gnames));
    end
    fprintf (cat (2, '\nSummary of bootstrap multiple comparison tests in', ...
                     ' a one-way layout\n', ...
                     '*******************************************', ...
                     '**********************************\n'));
    fprintf ('Bootstrap settings: \n');
    fprintf (' Function: %s\n', bootfun_str);
    fprintf (' Bootstrap resampling method: Balanced, bootstrap resampling\n')
    fprintf (' Number of bootstrap resamples: %u \n', nboot(1));
    fprintf (' Method for estimating standard errors: %s\n', se_method)
    if (nboot(2) > 0)
      fprintf (cat (2, ' Number of bootknife resamples used to estimate', ...
                       ' standard errors: %u \n'), nboot(2));
    end
    if (isempty (ref))
      fprintf (' Multiple comparison method:%s \n', ... 
               ' Single-step maxT procedure based on Tukey-Kramer');
    else
      fprintf (' Multiple comparison method:%s \n', ...
               ' Single-step maxT procedure based on Dunnett');
      fprintf (' Reference group used for comparisons: %u (%s) \n', ...
               gk(ref), gnames{ref});
    end
    if (size (c,1) >= 1)
      fprintf (cat (2, '\n----------------------------------------------', ...
                       '-------------------------------\n', ...
                       '| Comparison |     Test # |      Ref # |', ...
                       ' Difference |          t |        p |\n', ...
                       '|------------|------------|------------|', ...
                       '------------|------------|----------|\n'));
      if (isempty (ref))
        for i = 1:n
          tmp = num2cell (c(i, cols));
          tmp{end} = round (tmp{end} * 1000);
          if (c(i,9) <= 0.001)
            tmp(end) = [];
            fprintf (cat (2, '| %10u | %10u | %10u | %#+10.4g | %+10.2f |', ...
                             '    <.001 |***\n'), i, tmp{:});
          elseif (c(i,9) > 0.999)
            tmp(end) = [];
            fprintf (cat (2, '| %10u | %10u | %10u | %#+10.4g | %+10.2f |', ...
                             '    1.000 |\n'), i, tmp{:});
          else
            fprintf (cat (2, '| %10u | %10u | %10u | %#+10.4g | %+10.2f |', ...
                             '     .%03u |'), i, tmp{:});
            if (c(i,9) < 0.01)
              fprintf ('**\n')
            elseif (c(i,9) < 0.05)
              fprintf ('*\n')
            else
              fprintf ('\n')
            end
          end
        end
      else
        for j = 1:k-1
          tmp = num2cell (c(j, cols));
          tmp{end} = round (tmp{end} * 1000);
          if (c(j,9) <= 0.001)
            tmp(end) = [];
            fprintf (cat (2, '| %10u | %10u | %10u | %#+10.4g | %+10.2f |', ...
                             '    <.001 |***\n'), j, tmp{:});
          elseif (c(j,9) > 0.999)
            tmp(end) = [];
            fprintf (cat (2, '| %10u | %10u | %10u | %#+10.4g | %+10.2f |', ...
                             '    1.000 |\n'), j, tmp{:});
          else
            fprintf (cat (2, '| %10u | %10u | %10u | %#+10.4g | %+10.2f |', ...
                             '     .%03u |'), j, tmp{:});
            if (c(j,9) < 0.01)
              fprintf ('**\n')
            elseif (c(j,9) < 0.05)
              fprintf ('*\n')
            else
              fprintf ('\n')
            end
          end
        end
      end
      fprintf (cat (2, '\n-----------------------------------------------', ...
                      '------------------------------\n', ...
                      '|    GROUP # |                                   ', ...
                      '    GROUP label |        N |\n', ...
                      '|------------|-----------------------------------', ...
                      '----------------|----------|\n'));
      for j = 1:k
        fprintf ('| %10.6g | %49s | %8.3g |\n', gk(j), gnames{j}, nk(j));
      end
      fprintf ('\n')
    end

    % Plot graph of the difference in bootfun for each comparison with
    % 100*(1-alpha)% confidence intervals
    figure;
    nc = size(c,1);                    % Calculate number of comparisons to plot
    plot ([0; 0], [0; nc + 1]', 'k:'); % Plot vertical dashed line at 0 effect
    set (gca, 'Ydir', 'reverse')       % Flip y-axis direction
    ylim ([0.5, nc + 0.5]);            % Set y-axis limits
    hold on
    for i = 1 : nc
      if (c(i,9) < 0.05)
        % Plot marker for the difference estimate
        plot (c(i, 5), i, 'or', 'MarkerFaceColor', 'r');
        % Plot line for each confidence interval
        plot ([c(i, 6), c(i, 7)], i * ones (2, 1), 'r-');
      else
        % Plot marker for the difference estimate
        plot (c(i,5), i, 'ob', 'MarkerFaceColor', 'b');
        % Plot line for each confidence interval 
        plot ([c(i,6), c(i,7)], i * ones (2, 1), 'b-');   
      end
    end
    hold off
    xlabel (sprintf (cat (2, '%g%% bootstrap-t confidence interval for the', ...
                             ' difference'), 100 * (1 - alpha)));
    ylabel ('Comparison number (Test - Reference)');   

  end

end

%--------------------------------------------------------------------------

function maxT = maxstat (Y, g, nboot, bootfun, ref, ISOCTAVE)

  % Helper function file required for boot1way
  % Calculate maximum test statistic
  
  % maxstat cannot be a subfunction or nested function since 
  % Octave parallel threads won't be able to find it

  % Calculate the size of the data (N) and the number (k) of unique groups
  N = size (g, 1);
  gk = unique (g);
  k = numel (gk);

  % Compute the estimate (theta) and it's pooled (weighted mean) sampling
  % variance 
  theta = zeros (k, 1);
  SE = zeros (k, 1);
  Var = zeros (k, 1);
  nk = zeros (size (gk));
  for j = 1:k
    if (nboot == 0)
      nk(j) = sum (g == gk(j));
      if strcmp (func2str (bootfun), 'mean')
        theta(j) = mean (Y(g == gk(j), :));
        % Quick calculation for the standard error of the mean
        SE(j) = std (Y(g == gk(j), :), 0) / sqrt (nk(j));
      else
        theta(j) = bootfun (Y(g == gk(j), :));
        % If requested, compute unbiased estimates of the standard error
        % using jackknife resampling
        jackstat = jackknife (bootfun, Y(g == gk(j), :));
        SE(j) = sqrt ((nk(j) - 1) / nk(j) ...
                * sum (((mean (jackstat) - jackstat)).^2));
      end
    else
      % Compute unbiased estimate of the standard error by balanced bootknife
      % resampling. Bootknife resampling involves less computation than
      % Jackknife when sample sizes get larger
      theta(j) = bootfun (Y(g == gk(j), :));
      nk(j) = sum (g == gk(j));
      bootout = bootknife (Y(g == gk(j), :), [nboot, 0], bootfun, ...
                           NaN, [], 0, [], [], ISOCTAVE, false, true);
      SE(j) = bootout.std_error;
    end
    Var(j) = ((nk(j) - 1) / (N - k)) * SE(j)^2;
  end
  if (any (isnan (SE)))
    error (cat (2, 'boot1way:maxstat: Evaluating bootfun on the bootknife', ...
                   ' resamples created NaN values for the standard error'))
  end
  nk_bar = sum (nk.^2) ./ sum (nk);  % weighted mean sample size
  Var = sum (Var .* nk / nk_bar);    % weighted pooled sampling variance

  % Calculate weights to correct for unequal sample size  
  % when calculating standard error of the difference
  w = nk_bar ./ nk;

  % Calculate the maximum test statistic 
  if (isempty (ref))
    % Calculate Tukey-Kramer test statistic (without sqrt(2) factor)
    %
    % Bibliography:
    %  [1] en.wikipedia.org/wiki/Tukey%27s_range_test
    %  [2] cdn.graphpad.com/faq/1688/file/MulitpleComparisonAlgorithmsPrism8.pdf
    %  [3] www.graphpad.com/guides/prism/
    %          latest/statistics/stat_the_methods_of_tukey_and_dunne.htm
    idx = logical (triu (ones (k, k), 1));
    i = (1 : k)' * ones (1, k);
    j = ones (k, 1) * (1 : k);
    t = abs (theta(i(idx)) - theta(j(idx))) ./ ...
        sqrt (Var * (w(i(idx)) + w(j(idx))));
  else
    % Calculate Dunnett's test statistic 
    t = abs ((theta - theta(ref))) ./ sqrt (Var * (w + w(ref)));
  end
  maxT = max(t);
  
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
%! % COMPARISON OF TWO INDEPENDENT GROUPS WITH UNEQUAL SAMPLE SIZES 
%! % (analagous to Student's t-test)
%!
%! y =    [54       43
%!         23       34
%!         45       65
%!         54       77
%!         45       46
%!        NaN       65];
%! g = {'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'};
%!
%! boot1way (y(:), g(:), 'ref', 'male', 'nboot', 4999);
%!
%! % Please be patient, the calculations will be completed soon...

%!demo
%!
%! % COMPARISON OF TWO INDEPENDENT GROUPS WITH UNEQUAL SAMPLE SIZES 
%! % (a robust version of Student's t-test)
%!
%! y =    [54       43
%!         23       34
%!         45       65
%!         54       77
%!         45       46
%!        NaN       65];
%! g = {'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'};
%!
%! boot1way (y(:), g(:), 'ref', 'male', 'nboot', 4999, 'bootfun', 'robust');
%!
%! % Please be patient, the calculations will be completed soon...

%!demo
%!
%! % ONE-WAY ANOVA WITH EQUAL SAMPLE SIZES: Treatment vs. Control (1)
%!
%! y = [111.39 110.21  89.21  76.64  95.35  90.97  62.78;
%!      112.93  60.36  92.29  59.54  98.93  97.03  79.65;
%!       85.24 109.63  64.93  75.69  95.28  57.41  75.83;
%!      111.96 103.40  75.49  76.69  77.95  93.32  78.70];
%! g = [1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7];
%!
%! boot1way (y(:),g(:),'ref',1,'nboot',4999);
%!
%! % Please be patient, the calculations will be completed soon...

%!demo
%!
%! % ROBUST ONE-WAY ANOVA WITH EQUAL SAMPLE SIZES: Treatment vs. Control (1)
%!
%! y = [111.39 110.21  89.21  76.64  95.35  90.97  62.78;
%!      112.93  60.36  92.29  59.54  98.93  97.03  79.65;
%!       85.24 109.63  64.93  75.69  95.28  57.41  75.83;
%!      111.96 103.40  75.49  76.69  77.95  93.32  78.70];
%! g = [1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7];
%!
%! boot1way (y(:), g(:), 'ref', 1, 'nboot', 4999, 'bootfun', 'robust');
%!
%! % Please be patient, the calculations will be completed soon...

%!demo
%!
%! % ONE-WAY ANOVA WITH UNEQUAL SAMPLE SIZES: pairwise comparisons
%!
%! y = [54  87  45
%!      23  98  39
%!      45  64  51
%!      54  77  49
%!      45  89  50
%!      47 NaN  55];
%! g = [ 1   2   3
%!       1   2   3
%!       1   2   3 
%!       1   2   3
%!       1   2   3
%!       1   2   3];
%!
%! boot1way (y(:), g(:), 'nboot', 4999);
%!
%! % Please be patient, the calculations will be completed soon...

%!demo
%!
%! % COMPARE STANDARD DEVIATIONS BETWEEN 3 GROUPS: pairwise comparisons
%!
%! y = [54  87  45
%!      23  98  39
%!      45  64  51
%!      54  77  49
%!      45  89  50
%!      47 NaN  55];
%! g = [ 1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3];
%! p = boot1way (y(:),g(:),'bootfun',{@std,1});

%!demo
%!
%! % COMPARE CORRELATION COEFFICIENTS BETWEEN 2 DATA SETS
%! Y = randn (20, 2); g = [zeros(10, 1); ones(10, 1)];
%! func = @(M) cor (M(:,1), M(:,2));
%!
%! boot1way (Y, g, 'bootfun', func);
%!
%! % Please be patient, the calculations will be completed soon...

%!demo
%!
%! % COMPARE SLOPES FROM LINEAR REGRESSION ON 2 DATA SETS
%! y = randn (20, 1); x = randn (20, 1); X = [ones(20, 1), x];
%! g = [zeros(10, 1); ones(10, 1)];
%! func = @(M) subsref (M(:,2:end) \ M(:,1), ...
%!                      struct ('type', '()', 'subs', {{2}}));
%!
%! boot1way ([y, X], g, 'bootfun', func);
%!
%! % Please be patient, the calculations will be completed soon...

%!test
%! y = [111.39 110.21  89.21  76.64  95.35  90.97  62.78;
%!      112.93  60.36  92.29  59.54  98.93  97.03  79.65;
%!       85.24 109.63  64.93  75.69  95.28  57.41  75.83;
%!      111.96 103.40  75.49  76.69  77.95  93.32  78.70];
%! g = [1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7;
%!      1 2 3 4 5 6 7];
%! p = boot1way (y(:),g(:),'ref',1,'nboot',[999,0],'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (min(p), 0.01570455476054196, 1e-06);
%! end
%! p = boot1way (y(:),g(:),'nboot',[999,0],'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (min(p), 0.05374259279858003, 1e-06);
%! end
%! % Result from anova1 is 0.0387

%!test
%! y = [54       43
%!      23       34
%!      45       65
%!      54       77
%!      45       46
%!     NaN       65];
%! g = {'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'
%!      'male' 'female'};
%! p = boot1way (y(:),g(:),'ref','male','nboot',[999,0],'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (p, 0.281036161800208, 1e-06);
%! end
%! p = boot1way (y(:),g(:),'nboot',[999,0],'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (p, 0.281036161800208, 1e-06);
%! end
%! % Result from anova1 is 0.2613

%!test
%! y = [54  87  45
%!      23  98  39
%!      45  64  51
%!      54  77  49
%!      45  89  50
%!      47 NaN  55];
%! g = [ 1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3];
%! p = boot1way (y(:),g(:),'nboot',[999,0],'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (min(p), 0.001, 1e-06); % truncated at 0.001
%! end
%! % Result from anova1 is 4.162704768129188e-05

%!test
%! y = [54  87  45
%!      23  98  39
%!      45  64  51
%!      54  77  49
%!      45  89  50
%!      47 NaN  55];
%! g = [ 1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3
%!       1   2   3];
%! p = boot1way (y(:),g(:),'bootfun',@(y)std(y,1),'DisplayOpt',false);
%! p = boot1way (y(:),g(:),'bootfun',{@std,1},'DisplayOpt',false);
%! if (isempty (regexp (which ('boot'), 'mex$')))
%!   % test boot m-file result
%!   assert (min(p), 0.4523926257950379, 1e-06);
%! end

%!test
%! % Compare correlation coefficients
%! Y = randn (20, 2); g = [zeros(10, 1); ones(10, 1)];
%! func = @(M) cor (M(:,1), M(:,2));
%! p = boot1way (Y, g, 'bootfun', func, 'DisplayOpt', false);

%!test
%! % Compare slopes from linear regression
%! y = randn (20, 1); x = randn (20, 1); X = [ones(20, 1), x];
%! g = [zeros(10, 1); ones(10, 1)];
%! func = @(M) subsref (M(:,2:end) \ M(:,1), ...
%!                      struct ('type', '()', 'subs', {{2}}));
%! p = boot1way ([y, X], g, 'bootfun', func, 'DisplayOpt', false);
