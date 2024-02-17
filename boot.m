% This function returns resampled data or indices created by balanced bootstrap 
% or bootknife resampling.
%
% -- Function File: BOOTSAM = boot (N, NBOOT)
% -- Function File: BOOTSAM = boot (X, NBOOT)
% -- Function File: BOOTSAM = boot (..., NBOOT, LOO)
% -- Function File: BOOTSAM = boot (..., NBOOT, LOO, SEED)
% -- Function File: BOOTSAM = boot (..., NBOOT, LOO, SEED, WEIGHTS)
%
%     'BOOTSAM = boot (N, NBOOT)' generates NBOOT bootstrap samples of length N.
%     The samples generated are composed of indices within the range 1:N, which
%     are chosen by random resampling with replacement [1]. N and NBOOT must be
%     positive integers. The returned value, BOOTSAM, is a matrix of indices,
%     with N rows and NBOOT columns. The efficiency of the bootstrap simulation
%     is ensured by sampling each of the indices exactly NBOOT times, for first-
%     order balance [2-3]. Balanced resampling only applies when NBOOT > 1.
%
%     'BOOTSAM = boot (X, NBOOT)' generates NBOOT bootstrap samples, each the
%     same length as X (N). X must be a numeric vector, and NBOOT must be
%     positive integer. BOOTSAM is a matrix of values from X, with N rows
%     and NBOOT columns. The samples generated contains values of X, which
%     are chosen by balanced bootstrap resampling as described above [1-3].
%     Balanced resampling only applies when NBOOT > 1.
%
%     Note that the values of N and NBOOT map onto int32 data types in the 
%     boot MEX file. Therefore, these values must never exceed (2^31)-1.
%
%     'BOOTSAM = boot (..., NBOOT, LOO)' sets the resampling method. If LOO
%     is false, the resampling method used is balanced bootstrap resampling.
%     If LOO is true, the resampling method used is balanced bootknife
%     resampling [4]. The latter involves creating leave-one-out jackknife
%     samples of size N - 1, and then drawing resamples of size N with
%     replacement from the jackknife samples, thereby incorporating Bessel's
%     correction into the resampling procedure. LOO must be a scalar logical
%     value. The default value of LOO is false.
%
%     'BOOTSAM = boot (..., NBOOT, LOO, SEED)' sets a seed to initialize
%     the pseudo-random number generator to make resampling reproducible between
%     calls to the boot function. Note that the mex function compiled from the
%     source code boot.cpp is not thread-safe. Below is an example of a line of
%     code one can run in Octave/Matlab before attempting parallel operation of
%     boot.mex in order to ensure that the initial random seeds of each thread
%     are unique:
%       • In Octave:
%            pararrayfun (nproc, @boot, 1, 1, false, 1:nproc)
%       • In Matlab:
%            ncpus = feature('numcores'); 
%            parfor i = 1:ncpus; boot (1, 1, false, i); end;
%
%     'BOOTSAM = boot (..., NBOOT, LOO, SEED, WEIGHTS)' sets a weight
%     vector of length N. If WEIGHTS is empty or not provided, the default 
%     is a vector of length N, with each element equal to NBOOT (i.e. uniform
%     weighting). Each element of WEIGHTS is the number of times that the
%     corresponding index (or element in X) is represented in BOOTSAM.
%     Therefore, the sum of WEIGHTS must equal N * NBOOT. 
%
%  Bibliography:
%  [1] Efron, and Tibshirani (1993) An Introduction to the
%        Bootstrap. New York, NY: Chapman & Hall
%  [2] Davison et al. (1986) Efficient Bootstrap Simulation.
%        Biometrika, 73: 555-66
%  [3] Booth, Hall and Wood (1993) Balanced Importance Resampling
%        for the Bootstrap. The Annals of Statistics. 21(1):286-298
%  [4] Hesterberg T.C. (2004) Unbiasing the Bootstrap—Bootknife Sampling 
%        vs. Smoothing; Proceedings of the Section on Statistics & the 
%        Environment. Alexandria, VA: American Statistical Association.
%
%  boot (version 2023.01.28)
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

%!demo
%!
%! % N as input; balanced bootstrap resampling with replacement
%! boot(3, 20, false)

%!demo
%!
%! % N as input; (unbiased) balanced bootknife resampling with replacement
%! boot(3, 20, true)

%!demo
%! 
%! % N as input; balanced resampling with replacement; setting the random seed
%! boot(3, 20, false, 1) % Set random seed
%! boot(3, 20, false, 1) % Reset random seed, BOOTSAM is the same
%! boot(3, 20, false)    % Without setting random seed, BOOTSAM is different

%!demo
%!
%! % Vector (X) as input; balanced resampling with replacement; setting weights
%! x = [23; 44; 36];
%! boot(x, 10, false, 1)            % equal weighting
%! boot(x, 10, false, 1, [20;0;10]) % unequal weighting, no x(2) in BOOTSAM 

%!demo
%!
%! % N as input; resampling without replacement; 3 trials
%! boot(6, 1, false, 1) % Sample 1; Set random seed for first sample only
%! boot(6, 1, false)    % Sample 2
%! boot(6, 1, false)    % Sample 3

%!test
%! % Test that random samples vary between calls to boot.
%! I1 = boot (3, 20);
%! I2 = boot (3, 20);
%! assert (all (I1(:) == I2(:)), false);

%!test
%! % Test that random seed gives identical resamples when LOO is false.
%! I1 = boot (3, 20, false, 1);
%! I2 = boot (3, 20, false, 1);
%! assert (all (I1(:) == I2(:)), true);

%!test
%! % Test that random seed gives identical resamples when LOO is true.
%! I1 = boot (3, 20, true, 1);
%! I2 = boot (3, 20, true, 1);
%! assert (all (I1(:) == I2(:)), true);

%!test
%! % Test that default setting for LOO is false.
%! I1 = boot (3, 20, [], 1);
%! I2 = boot (3, 20, false, 1);
%! assert (all (I1(:) == I2(:)), true);

%!test
%! % Test that resampling is balanced when LOO is false.
%! I = boot (3, 20, false, 1);
%! assert (sum (I(:) == 1), 20, 1e-03);
%! assert (sum (I(:) == 2), 20, 1e-03);
%! assert (sum (I(:) == 3), 20, 1e-03);

%!test
%! % Test that resampling is balanced when LOO is true.
%! I = boot (3, 20, true, 1);
%! assert (sum (I(:) == 1), 20, 1e-03);
%! assert (sum (I(:) == 2), 20, 1e-03);
%! assert (sum (I(:) == 3), 20, 1e-03);

%! % Test for unbiased sampling (except for last sample).
%! % The exception is a requirement to maintain balance.
%! I = boot (3, 20, true, 1);
%! assert (all (diff (sort (I(1:end-1)))), false);

%!test
%! % Test feature for changing resampling weights when LOO is false
%! I = boot (3, 20, false, 1, [30,30,0]);
%! assert (any (I(:) == 3), false);

%!test
%! % Test feature for changing resampling weights when LOO is true
%! I = boot (3, 20, true, 1, [30,30,0]);
%! assert (any (I(:) == 3), false);
