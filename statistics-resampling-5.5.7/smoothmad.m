% Calculates a smoothed version of the median absolute deviation (MAD).
%
% -- Function File: MAD = smoothmad (X)
% -- Function File: MAD = smoothmad (X, GROUP)
% -- Function File: MAD = smoothmad (X, GROUP, CONSTANT)
%
%     'MAD = smoothmad (X)' calculates a smoothed version of the median
%     absolute deviation (MAD) for each column of the data in x. The
%     statistics are scaled by a constant of 1.41 to make the estimator
%     consistent with the standard deviation for normally distributed data.
% 
%     'MAD = smoothmad (X, GROUP)' defines group membership of the rows
%     of x and returns the pooled smooth MAD. GROUP must be a numeric
%     vector with the same number of rows as x.
%
%     'MAD = smoothmad (X, GROUP, CONSTANT)' sets the CONSTANT to scale
%     the value of the MAD. (Default is 1.41).
%
%  smoothmad (version 2023.05.02)
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


function PMAD = smoothmad (x, group, constant)

  % Evaluate input arguments
  if ((nargin < 1) || (nargin > 3))
    error ('smoothmad: Invalid number of input arguments')
  end

  if (nargout > 1)
    error ('smoothmad: Invalid number of output arguments')
  end

  if (numel (size (x)) > 2)
    error ('smoothmad: Arrays of more than 2 dimensions are not supported')
  end

  [m, n] = size (x);
  if (nargin < 2 || isempty (group))
    group = ones(m,1);
  else
    if size (group, 1) ~= size (x, 1)
      error (cat (2, 'smoothmad: GROUP must be a column vector with the', ...
                     ' same number of rows as the data'))
    end
  end

  if (nargin < 3)
    constant = 1.41;
  end

  % Initialize variables
  g = unique (group);
  l = numel (g);
  M = zeros (l, n);
  MAD = zeros (l, n);
  nk = zeros (l, 1);
  PMAD = zeros (1, n);

  % Perform calculations on each data group
  for k = 1:l

    % Collect the data for group k
    I = (group == g(k));

    % Calculate the size of the data group
    nk(k) = sum(I);

    % Calculate the smoothed median of the data group 
    M(k,:) = smoothmedian (x(I,:));

    % Calculate the smoothed median absolute deviation of the data group
    MAD(k,:) = smoothmedian (abs (x(I, :) - M(k, :))) * constant;

    % Begin pooling the smoothed median absolute deviations
    PMAD = PMAD + (nk(k) - 1) * MAD(k, :).^2;

  end

  % Calculate the pooled, smoothed median absolute deviation
  PMAD = sqrt (PMAD / (sum (nk) - l));
