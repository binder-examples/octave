% Balanced bootstrap resampling.
%
%
% -- Function File: BOOTSTAT = bootstrp (NBOOT, BOOTFUN, D)
% -- Function File: BOOTSTAT = bootstrp (NBOOT, BOOTFUN, D1, ..., DN)
% -- Function File: BOOTSTAT = bootstrp (..., 'seed', SEED)
% -- Function File: BOOTSTAT = bootstrp (..., 'Options', PAROPT)
% -- Function File: [BOOTSTAT, BOOTSAM] = bootstrp (...) 
%
%     BOOTSTAT = bootstrp (NBOOT, BOOTFUN, D) draws NBOOT bootstrap resamples
%     from the data D and returns the statistic computed by BOOTFUN in BOOTSTAT
%     [1]. bootstrp resamples from the rows of a data sample D (column vector
%     or a matrix). BOOTFUN is a function handle (e.g. specified with @), or a
%     string indicating the function name. The third input argument is data
%     (column vector or a matrix), that is used to create inputs for BOOTFUN.
%     The resampling method used throughout is balanced bootstrap resampling
%     [2-3].
%
%     BOOTSTAT = bootstrp (NBOOT, BOOTFUN, D1,...,DN) is as above except that 
%     the third and subsequent numeric input arguments are data vectors that
%     are used to create inputs for BOOTFUN.
%
%     BOOTSTAT = bootstrp (..., 'seed', SEED) initialises the Mersenne Twister
%     random number generator using an integer SEED value so that bootci results
%     are reproducible.
%
%     BOOTSTAT = bootstrp (..., 'Options', PAROPT) specifies options that govern
%     if and how to perform bootstrap iterations using multiple processors (if
%     the Parallel Computing Toolbox or Octave Parallel package is available).
%     This argument is a structure with the following recognised fields:
%        o 'UseParallel':  If true, use parallel processes to accelerate
%                          bootstrap computations on multicore machines. 
%                          Default is false for serial computation. In MATLAB,
%                          the default is true if a parallel pool
%                          has already been started. 
%        o 'nproc':        nproc sets the number of parallel processes
%
%     [BOOTSTAT, BOOTSAM] = bootstrp (...) also returns BOOTSAM, a matrix of
%     indices from the bootstrap. Each column in BOOTSAM corresponds to one
%     bootstrap sample and contains the row indices of the values drawn from
%     the nonscalar data argument to create that sample.
%
%  Bibliography:
%  [1] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [2] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [3] Booth, Hall and Wood (1993) Balanced Importance Resampling
%        for the Bootstrap. The Annals of Statistics. 21(1):286-298
%
%  bootstrp (version 2023.06.20)
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


function [bootstat, bootsam] = bootstrp (argin1, argin2, varargin)

  % Evaluate the number of function arguments
  if (nargin < 2)
    error (cat (2, 'bootstrp usage: ''bootstrp (nboot, {bootfun, data},', ...
                   ' varargin)''; atleast 2 input arguments required'))
  end

  % Check if using MATLAB or Octave
  info = ver; 
  ISOCTAVE = any (ismember ({info.Name}, 'Octave'));

  % Apply defaults
  paropt = struct;
  paropt.UseParallel = false;
  if (~ ISOCTAVE)
    ncpus = feature('numcores');
  else
    ncpus = nproc;
  end
  paropt.nproc = ncpus;

  % Assign input arguments to function variables
  nboot = argin1;
  bootfun = argin2;
  argin3 = varargin;
  narg = numel (argin3);
  if (narg > 1)
    while ischar(argin3{end-1})
      if (any(strcmpi({'Options','Option'},argin3{end-1})))
        paropt = argin3{end};
      elseif (any(strcmpi('seed',argin3{end-1})))
        seed = argin3{end};
        % Initialise the random number generator with the seed
        boot (1, 1, false, seed);
      else
        error ('bootstrp: Unrecognised input argument to bootstrp')
      end
      argin3 = {argin3{1:end-2}};
      narg = numel (argin3);
      if (narg < 2)
        break
      end
    end
  end
  if (numel (argin3) > 1)
    x = argin3;
  else
    x = argin3{1};
  end
  if (paropt.UseParallel)
    ncpus = paropt.nproc;
  else
    ncpus = 0;
  end

  % Error checking
  % nboot input argument
  if ((nargin < 2) || isempty (nboot))
    nboot = 2000;
  else
    if (~ isa (nboot, 'numeric'))
      error ('bootstrp: NBOOT must be numeric');
    end
    if (numel (nboot) > 1)
      error ('bootstrp: NBOOT cannot contain more than 1 value');
    end
    if (nboot ~= abs (fix (nboot)))
      error ('bootstrp: NBOOT must contain positive integers');
    end    
  end
  if (~ all (size (nboot) == [1, 1]))
    error ('bootstrp: NBOOT must be a scalar value')
  end

  % bootfun input argument
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
      error ('bootstrp: BOOTFUN must be a function name or function handle')
    end
  end

  % Determine properties of the DATA (x)
  szx = size (x);
  n = szx(1);
  nvar = szx(2);
  if (n < 2)
    error ('bootstrp: DATA must be numeric and contain > 1 row')
  end

  % If applicable, check we have parallel computing capabilities
  if (ncpus > 1)
    if (ISOCTAVE)
      pat = '^parallel';
      software = pkg ('list');
      names = cellfun (@(S) S.name, software, 'UniformOutput', false);
      status = cellfun (@(S) S.loaded, software, 'UniformOutput', false);
      index = find (~ cellfun (@isempty, regexpi (names,pat)));
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
      info = ver; 
      if (ismember ('Parallel Computing Toolbox', {info.Name}))
        PARALLEL = true;
      else
        PARALLEL = false;
      end
    end
  end

  % Evaluate bootfun on the DATA
  T0 = bootfun (x);
  if (any (isnan (T0)))
    error ('bootstrp: BOOTFUN returned NaN with the DATA provided')
  end

  % Check whether bootfun is vectorized
  if (nvar > 1)
    M = cell2mat (cellfun (@(i) repmat (x(:, i), 1, 2), ...
                  num2cell (1:nvar), 'UniformOutput', false));
  else
    M = repmat (x, 1, 2);
  end
  if (any (szx > 1))
    vectorized = false;
  else
    try
      chk = bootfun (M);
      if (all (size (chk) == [size(T0, 1), 2]) && all (chk == bootfun (x)))
        vectorized = true;
      else
        vectorized = false;
      end
    catch
      vectorized = false;
    end
  end

  % If applicable, setup a parallel pool (required for MATLAB)
  if (~ ISOCTAVE)
    % MATLAB
    % bootfun is not vectorized
    if (ncpus > 0) 
      % MANUAL
      try 
        pool = gcp ('nocreate'); 
        if isempty (pool)
          if (ncpus > 1)
            % Start parallel pool with ncpus workers
            parpool (ncpus);
          else
            % Parallel pool is not running and ncpus is 1 so run function
            % evaluations in serial
            ncpus = 1;
          end
        else
          if (pool.NumWorkers ~= ncpus)
            % Check if number of workers matches ncpus and correct it
            % accordingly if not
            delete (pool);
            if (ncpus > 1)
              parpool (ncpus);
            end
          end
        end
      catch
        % MATLAB Parallel Computing Toolbox is not installed
        warning ('bootstrp:parallel', ...
                 cat (2, 'Parallel Computing Toolbox not installed or', ...
                         ' operational. Falling back to serial processing.'))
        ncpus = 1;
      end
    end
  else
    if ((ncpus > 1) && ~ PARALLEL)
      if (ISOCTAVE)
        % OCTAVE Parallel Computing Package is not installed or loaded
        warning ('bootstrp:parallel', ...
                 cat (2, 'Parallel Computing Package not installed and/or', ...
                         ' loaded. Falling back to serial processing.'))
      else
        % MATLAB Parallel Computing Toolbox is not installed or loaded
        warning ('bootstrp:parallel', ...
                 cat (2, 'Parallel Computing Toolbox not installed and/or', ...
                         ' loaded. Falling back to serial processing.'))
      end
      ncpus = 0;
    end
  end

  % Calculate the number of elements in the return value of bootfun 
  m = numel (T0);
  if (m > 1)
    % Vectorized along the dimension of the return values of bootfun so
    % reshape the output to be a column vector before proceeding with bootstrap
    if (size (T0, 2) > 1)
      bootfun = @(x) reshape (bootfun (x), [], 1);
      T0 = reshape (T0, [], 1);
      vectorized = false;
    end
  end

  % Perform balanced bootstrap resampling
  unbiased = false;  % Set to false for bootstrap resampling
  if (nvar > 1) || (nargout > 1)
    % We can save some memory by making bootsam an int32 datatype
    bootsam = zeros (n, nboot, 'int32');
    bootsam(:, :) = boot (n, nboot, unbiased);
  else
    % For more efficiency, if we don't need bootsam, we can directly resample
    % values of x
    bootsam = [];
    X = boot (x, nboot, unbiased);
  end

  % Evaluate bootfun each bootstrap resample
  if (isempty (bootsam))
    if (vectorized)
      % Vectorized evaluation of bootfun on the DATA resamples
      bootstat = bootfun (X);
    else
      if (ncpus > 1)
        % Evaluate bootfun on each bootstrap resample in PARALLEL
        if (ISOCTAVE)
          % OCTAVE
          bootstat = parcellfun (ncpus, bootfun, num2cell (X, 1), ...
                                 'UniformOutput', false);
        else
          % MATLAB
          bootstat = cell (1, nboot);
          parfor b = 1:nboot; bootstat{b} = bootfun (X(:, b)); end
        end
      else
        bootstat = cellfun (bootfun, num2cell (X, 1), 'UniformOutput', false);
      end
    end
  else
    if (vectorized)
      % DATA resampling (using bootsam) and vectorized evaluation of bootfun on 
      % the DATA resamples 
      if (nvar > 1)
        % Multivariate
        % Perform DATA sampling
        X = cell2mat (cellfun (@(i) reshape (x(bootsam, i), n, nboot), ...
                      num2cell (1:nvar, 1), 'UniformOutput', false));
      else
        % Univariate
        % Perform DATA sampling
        X = x(bootsam);
      end
      % Function evaluation on bootstrap samples
      bootstat = bootfun (X);
    else 
      cellfunc = @(bootsam) bootfun (x(bootsam, :));
      if (ncpus > 1)
        % Evaluate bootfun on each bootstrap resample in PARALLEL
        if (ISOCTAVE)
          % OCTAVE
          bootstat = parcellfun (ncpus, cellfunc, num2cell (bootsam, 1), ...
                                 'UniformOutput', false);
        else
          % MATLAB
          bootstat = cell (1, nboot);
          parfor b = 1:nboot; bootstat{b} = cellfunc (bootsam(:, b)); end
        end
      else
        % Evaluate bootfun on each bootstrap resample in SERIAL
        bootstat = cellfun (cellfunc, num2cell (bootsam, 1), ...
                            'UniformOutput', false);
      end
    end
  end
  if (iscell (bootstat))
    bootstat = cell2mat (bootstat);
  end

  % Format output to be consistent with MATLAB's bootstrp
  bootstat = bootstat.';
  bootsam = double (bootsam);

end

%!demo
%!
%! % Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! % Compute 50 bootstrap statistics for the mean and calculate the bootstrap
%! % standard arror
%! bootstat = bootstrp (50, @mean, data)
%! std (bootstat)

%!test
%!
%! % Input univariate dataset
%! data = [48 36 20 29 42 42 20 42 22 41 45 14 6 ...
%!         0 33 28 34 4 32 24 47 41 24 26 30 41]';
%!
%! % Compute 50 bootstrap statistics for the mean and calculate the bootstrap
%! % standard arror
%! bootstat = bootstrp (50, @mean, data);
