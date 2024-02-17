% Performs sample size calculations, with optional correction for the design
% effect deviating from unity.
%
% -- Function File: N = sampszcalc (TESTTYPE, EFFSZ)
% -- Function File: N = sampszcalc (TESTTYPE, EFFSZ, POW)
% -- Function File: N = sampszcalc (TESTTYPE, EFFSZ, POW, ALPHA)
% -- Function File: N = sampszcalc (TESTTYPE, EFFSZ, POW, ALPHA, TAILS)
% -- Function File: N = sampszcalc (TESTTYPE, EFFSZ, POW, ALPHA, TAILS, DEFF)
%
%      'N = sampszcalc (TESTTYPE, EFFSZ)' returns the required sample size to
%      reach the significance level (alpha) of 0.05 in a two-tailed version of
%      the test specified in TESTTYPE for the specified effect size, EFFSZ,
%      with a power of 0.8 (i.e. a type II error rate of 1 - 0.8 = 0.2)
%
%        TESTTYPE can be:
%
%          't2' (default) : two-sample unpaired t-test
%
%          't' : paired t-test or one-sample t-test
%
%          'z2' (default) : two-sample unpaired z-test (Normal approximation)
%
%          'z' : paired z-test or one-sample z-test (Normal approximation)
%
%          'r' : significance test for correlation
%
%        EFFSZ can be numeric value corresponding to the standardized effect
%        size: Cohen's d or h (when TESTTYPE is 't2', 't', 'z' or 'z'), or 
%        Pearson's correlation coefficient (when TESTTYPE is 'r'). For
%        convenience, EFFSZ can also be one of the following strings:
%
%          'small' : which is 0.2 for Cohen's d (or h), or 0.1 for Pearson's r.
%
%          'medium' : which is 0.5 for Cohen's d (or h), or 0.3 for Pearson's r.
%
%          'large' : which is 0.8 for Cohen's d (or h), or 0.5 for Pearson's r.
%
%       'N = sampszcalc (TESTTYPE, EFFSZ, POW)' also sets the desired power of
%       the test. The power corresponds to 1 - beta, where beta is the type II
%       error rate (i.e. the probability of not rejecting the null hypothesis
%       when it is actually false). (Default is 0.8)
%
%       'N = sampszcalc (TESTTYPE, EFFSZ, POW, ALPHA)' also sets the desired
%       significance level, ALPHA, of the test. ALPHA corresponds to the type I
%       error rate (i.e. the probability of rejecting the null hypothesis when
%       it is actually true). (Default is 0.05)
%
%       HINT: If the test is expected to be among a family of tests, divide
%       ALPHA by the number of tests so that the sample size calculations will
%       maintain the desired power after correction for multiple comparisons.
%
%       'N = sampszcalc (TESTTYPE, EFFSZ, POW, ALPHA, TAILS)' also sets whether
%       the test is one-sided or two-sided. (Default is 2)
%
%       'N = sampszcalc (TESTTYPE, EFFSZ, POW, ALPHA, TAILS, DEFF)' also sets
%       the design effect to correct the sample size calculation. (Default is 1)
%       DEFF can be estimated by dividing the sampling variance of the parameter
%       of interest from a complex experimental design by the equivalent
%       statistic computed using simple random sampling with replacement.
%
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
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see http://www.gnu.org/licenses/

function n = sampszcalc (testtype, effsz, power, alpha, tails, deff)

    % Set default values
    if (nargin < 2)
      error ('sampszcalc: atleast 2 input arguments required')
    end
    if ( (nargin < 3) || isempty (power) )
      power = 0.8; 
    end
    if ( (nargin < 4) || isempty (alpha) )
      alpha = 0.05;
    end
    if ( (nargin < 5) || isempty (tails) )
      tails = 2;
    end
    if ( (nargin < 6) || isempty (deff) )
      deff = 1;
    end

    % Error checking
    if (ischar (testtype))
      if (~ismember (testtype, {'t2', 't', 'z2', 'z', 'r'}))
        error ('sampszcalc: TESTTYPE not supported')
      end
    else
      error ('sampszcalc: TESTTYPE must be a character string')
    end
    if (~ any (isa (power, {'single','double'})))
      error ('sampszcalc: POW must be single or double precision')
    end
    if (numel (power) > 1)
      error ('sampszcalc: POW must be scalar')
    end
    if ( (power <= 0) || (power >= 1) )
      error ('sampszcalc: the value of POW must be: 0 < POW < 1')
    end
    if (~ any (isa (alpha, {'single','double'})))
      error ('sampszcalc: ALPHA must be single or double precision')
    end
    if (numel (alpha) > 1)
      error ('sampszcalc: ALPHA must be scalar')
    end
    if ( (alpha <= 0) || (alpha >= 1) )
      error ('sampszcalc: the value of ALPHA must be: 0 < ALPHA < 1')
    end
    if (~ ismember (tails, 1:2))
      error ('sampszcalc: tails must be either 1 or 2')
    end
    if (~ any (isa (deff, {'single','double'})))
      error ('sampszcalc: DEFF must be single or double precision')
    end
    if (numel (deff) > 1)
      error ('sampszcalc: DEFF must be scalar')
    end

    % Perform sample size calculation
    if strcmpi (testtype, 'r')
      if (ischar (effsz))
        switch (lower (effsz))
          case 'small'
            STAT = 0.1;
          case 'medium'
            STAT = 0.3;
          case 'large'
            STAT = 0.5;
          otherwise
            error ('sampszcalc: string description for EFFSIZE not recognised')
        end
      else 
        STAT = abs (effsz);
      end
      STAT = atanh (STAT);
      testtype = 'z';
      c = 3;
    else 
      if (ischar (effsz))
        switch (lower (effsz))
          case 'small'
            STAT = 0.2;
          case 'medium'
            STAT = 0.5;
          case 'large'
            STAT = 0.8;
          otherwise
            error ('sampszcalc: string description for EFFSIZE not recognised')
        end
      else 
        STAT = abs (effsz);
      end
      c = 0;
    end

    % Sample size calculations for the difference between means
    % Assume effect size is Cohen's d
    k = numel (testtype);
    % Calculate group sample size based on Normal approximation
    stdnorminv = @(p) sqrt (2) * erfinv (2 * p-1);
    n0 = k * (((stdnorminv (power) + ...
                stdnorminv (1 - alpha / tails)) / STAT)^2 + c);
    switch ( lower (testtype) )
      case {'z','z2'}
        n = ceil (n0 * deff);
      case {'t','t2'}
        % Create function to optimize sample size based on Student-t
        % distribution and n * k - k degrees of freedom
        if (exist ('betaincinv', 'file'))
          studinv = @(P, DF) sign (P - 0.5) * ...
                sqrt ( DF ./ betaincinv (2 * min (P, 1 - P), DF / 2, 0.5) - DF);
        else
          studinv = @(P, DF) sign (P - 0.5) * ...
                   sqrt ( DF ./ betainv (2 * min (P, 1 - P), DF / 2, 0.5) - DF);
        end
        func = @(n) n - k * ...
                 (((studinv (power, n * k - k) + ...
                    studinv (1 - alpha / tails, n * k - k)) / STAT)^2 + c);
        n = ceil (fzero (func, n0) * deff); % Find the root using fzero
    end

end

%!demo
%!
%! % The difference between a sample mean from a zero constant (one sample test)
%! % or the difference between two dependent means (matched pair). Sample size
%! % determined for Cohen's d = 0.8.
%! % d effect size
%!
%! n = sampszcalc ('t', 0.8)

%!demo
%!
%! % The difference between two independent means (two groups). Sample size
%! % determined for Cohen's d = 0.8.
%!
%! n = sampszcalc ('t2', 0.8)

%!demo
%!
%! % The difference between two independent means (two groups). Sample size
%! % determined for Cohen's d = 0.8 and a design effect of 1.5
%!
%! n = sampszcalc ('t2', 0.8, [], [], [], 1.5)

%!demo
%!
%! % The difference between two independent proportions (two sample test). 
%! % Sample size determined for Cohen's h = 0.8 using Normal approximation.
%!
%! n = sampszcalc ('z2', 0.8)

%!demo
%!
%! % The test for Pearson's correlation coefficient (r) equal to 0 (constant),
%! % Sample size determined for r effect size = 0.5.
%!
%! n = sampszcalc ('r', 0.5)

%!demo
%!
%! % Sample size calculation for nested two-sample test using the design effect
%! % from a pilot experiment. N below corresponds to the number of independent
%! % sampling units (i.e. clusters).
%! % See also the help documentation for functions bootlm and deffcalc.
%! score = [21, 26, 33, 22, 18, 25, 26, 24, 21, 25, 35, 28, 32, 36, 38, ...
%!          26, 34, 27, 38, 44, 34, 45, 38, 31, 41, 34, 35, 38, 46]';
%! method = {'A','A','A','A','A','A','A','A','A','A','A','A','A','A','A', ...
%!           'B','B','B','B','B','B','B','B','B','B','B','B','B','B'}';
%! room = [1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, ...
%!         1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3]';
%!
%! [STATS_STD] = bootlm (score, {method}, 'clustid', room, ...
%!                             'seed', 1, 'display', 'off', 'dim', 1, ...
%!                             'posthoc', 'trt_vs_ctrl', ...
%!                             'method', 'bayesian', 'prior', 'auto', ...
%!                             'standardize', true);
%!
%! [STATS, BOOTSTAT] = bootlm (score, {method}, 'clustid', room, ...
%!                             'seed', 1, 'display', 'off', 'dim', 1, ...
%!                             'posthoc', 'trt_vs_ctrl', ...
%!                             'method', 'bayesian', 'prior', 'auto');
%!
%! [STATS_SRS, BOOTSTAT_SRS] = bootlm (score, {method}, 'clustid', [], ...
%!                             'seed', 1, 'display', 'off', 'dim', 1, ...
%!                             'posthoc', 'trt_vs_ctrl', ...
%!                             'method', 'bayesian', 'prior', 'auto');
%!
%! fprintf('Cohen''s d = %.2f\n', STATS_STD.estimate)
%!
%! N = sampszcalc ('t2', STATS_STD.estimate, 0.80, 0.05, 2)
%!
%! DEFF = deffcalc (BOOTSTAT, BOOTSTAT_SRS)
%!
%! N_corrected = sampszcalc ('t2', STATS_STD.estimate, 0.80, 0.05, 2, DEFF)

%!test
%! % The difference between a sample mean from a zero constant (one sample
%! % test) or the difference between two dependent means (matched pair).
%! % Required sample sizes for small, medium and large effects with power, 
%! % alpha and the number of tails at 0.8, 0.05 and 2 respectively (defaults)
%! % The standardized effect size corresponds to Cohen's d
%! % Results compared to G*Power 3.1
%! ns = sampszcalc ('t', 0.20, 0.80, 0.05, 2);
%! assert (ns, 199, 1);
%! nm = sampszcalc ('t', 0.50, 0.80, 0.05, 2);
%! assert (nm, 34, 1);
%! nl = sampszcalc ('t', 0.80, 0.80, 0.05, 2);
%! assert (nl, 15, 1);
%! ns = sampszcalc ('t', 'small');
%! assert (ns, 199, 1);
%! nm = sampszcalc ('t', 'medium');
%! assert (nm, 34, 1);
%! nl = sampszcalc ('t', 'large');
%! assert (nl, 15, 1);

%!test
%! % The difference between two independent means (two groups).
%! % Required sample sizes for small, medium and large effects with power, 
%! % alpha and the number of tails at 0.8, 0.05 and 2 respectively (defaults)
%! % The standardized effect size corresponds to Cohen's d
%! % Results compared to G*Power 3.1
%! ns = sampszcalc ('t2', 0.20, 0.80, 0.05, 2);
%! assert (ns, 394, 1);
%! nm = sampszcalc ('t2', 0.50, 0.80, 0.05, 2);
%! assert (nm, 64, 1);
%! nl = sampszcalc ('t2', 0.80, 0.80, 0.05, 2);
%! assert (nl, 26, 1);
%! ns = sampszcalc ('t2', 'small');
%! assert (ns, 394, 1);
%! nm = sampszcalc ('t2', 'medium');
%! assert (nm, 64, 1);
%! nl = sampszcalc ('t2', 'large');
%! assert (nl, 26, 1);

%!test
%! % The difference between a proportion and a constant (one sample test)
%! % or the difference between two dependent proportions (matched pair).
%! % Required sample sizes for small, medium and large effects with power, 
%! % alpha and the number of tails at 0.8, 0.05 and 2 respectively (defaults)
%! % The standardized effect size corresponds to Cohen's h
%! % Results compared to NCSS PASS
%! ns = sampszcalc ('z', 0.20, 0.80, 0.05, 2);
%! assert (ns, 197, 1);
%! nm = sampszcalc ('z', 0.50, 0.80, 0.05, 2);
%! assert (nm, 32, 1);
%! nl = sampszcalc ('z', 0.80, 0.80, 0.05, 2);
%! assert (nl, 13, 1);
%! ns = sampszcalc ('z', 'small');
%! assert (ns, 197, 1);
%! nm = sampszcalc ('z', 'medium');
%! assert (nm, 32, 1);
%! nl = sampszcalc ('z', 'large');
%! assert (nl, 13, 1);

%!test
%! % The difference between two independent proportions (two-sample test).
%! % Required sample sizes for small, medium and large effects with power, 
%! % alpha and the number of tails at 0.8, 0.05 and 2 respectively (defaults)
%! % The standardized effect size corresponds to Cohen's h
%! % Results compared to NCSS PASS
%! ns = sampszcalc ('z2', 0.20, 0.80, 0.05, 2);
%! assert (ns, 393, 1);
%! nm = sampszcalc ('z2', 0.50, 0.80, 0.05, 2);
%! assert (nm, 63, 1);
%! nl = sampszcalc ('z2', 0.80, 0.80, 0.05, 2);
%! assert (nl, 25, 1);
%! ns = sampszcalc ('z2', 'small');
%! assert (ns, 393, 1);
%! nm = sampszcalc ('z2', 'medium');
%! assert (nm, 63, 1);
%! nl = sampszcalc ('z2', 'large');
%! assert (nl, 25, 1);

%!test
%! % The test for Pearson's correlation coefficient equal to 0.
%! % Required sample sizes for small, medium and large effects with power, 
%! % alpha and the number of tails at 0.8, 0.05 and 2 respectively (defaults)
%! % Results compared to G*Power 3.1
%! ns = sampszcalc ('r', 0.10, 0.80, 0.05, 2);
%! assert (ns, 783, 1);
%! nm = sampszcalc ('r', 0.30, 0.80, 0.05, 2);
%! assert (nm, 85, 1);
%! nl = sampszcalc ('r', 0.50, 0.80, 0.05, 2);
%! assert (nl, 30, 1);
%! ns = sampszcalc ('r', 'small');
%! assert (ns, 783, 1);
%! nm = sampszcalc ('r', 'medium');
%! assert (nm, 85, 1);
%! nl = sampszcalc ('r', 'large');
%! assert (nl, 30, 1);

