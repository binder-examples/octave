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
