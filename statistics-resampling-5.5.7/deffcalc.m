% Computes the design effect (DEFF), which can subsequently be used to correct
% sample size calculations using the 'sampszcalc' function.
%
% -- Function File: DEFF = deffcalc (BOOTSTAT, BOOTSTAT_SRS)
%
%     'DEFF = deff_calc (BOOTSTAT, BOOTSTAT_SRS)' computes the design effect
%     (DEFF) by taking the ratio of the variance of the bootstrap statistics
%     from a complex design over the variance of bootstrap statistics from
%     simple random sampling with replacement:
%
%            DEFF = var (BOOTSTAT, 0, 2) ./ var (BOOTSTAT_SRS, 0, 2);
%
%     BOOTSTAT and BOOTSTAT_SRS must be row vectors, or matrices with dimenions
%     of size P * NBOOT, where P is the number of parameters being estimated
%     and NBOOT is the number of bootstrap statistics. The number of parameters
%     being estimated (but not the number of bootstrap resamples) must be the
%     same to compute DEFF using this function.
%
%  deffcalc (version 2023.09.17)
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

function DEFF = deffcalc (BOOTSTAT, BOOTSTAT_SRS)

  % Check input and output arguments
  if (nargin < 2)
    error ('deffcalc: Too few input arguments.')
  end
  if (nargin > 2)
    error ('deffcalc: Too many input arguments.')
  end
  if (nargout > 1)
    error ('deffcalc: Too many output arguments.')
  end

  % Check that the input arguments have consistent dimenions
  szA = size (BOOTSTAT);
  szB = size (BOOTSTAT_SRS);
  if (szA(1) ~= szB(1))
    error ('deffcalc: input arguments must have the same number of rows')
  end
  DEFF = var (BOOTSTAT, 0, 2) ./ var (BOOTSTAT_SRS, 0, 2);
