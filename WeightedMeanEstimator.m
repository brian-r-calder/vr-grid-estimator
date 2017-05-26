% Compute the precision weighted mean of observations, the precision of
% the input data, and the standard deviation of the observations about
% the mean.  The code accumulates data as it's presented, and therefore
% can provide mean, precision, and standard deviation estimates at any
% point, and continue to accumulate data afterwards.

% Copyright (c) 2017, University of New Hampshire, Center for Coastal and
% Ocean Mapping & NOAA-UNH Joint Hydrographic Center.
%
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, write to the Free Software Foundation, Inc.,
% 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
%
% (You might also be able to get a copy of the license electronically if
% required, from <http://www.gnu.org/licenses/>)

classdef WeightedMeanEstimator < Estimator
    properties (GetAccess = public, SetAccess = private)
        n_obs = 0           % Number of observations added
    end
    properties (GetAccess = private)
        sum_depths = 0.0    % Sum of observations (m)
        sum_sq = 0.0        % Sum of the squares of the obs (m^2) (for std. dev.)
        prec = 0.0          % Running precision of all added obs. (m^-2)
    end
    methods
        function obj = Update(obj, obs)
            % Add another observation to the estimator, accumulating
            % statistics.  The _obs_ variable must be an Observation class,
            % or something that provides the same properties.  Note that
            % this method cannot be used with an array of objects.
            obj.sum_depths = obj.sum_depths + (1.0/obs.uncrt)*obs.depth;
            obj.prec = obj.prec + 1.0/obs.uncrt;
            obj.sum_sq = obj.sum_sq + obs.depth*obs.depth;
            obj.n_obs = obj.n_obs + 1;
        end
        function r = GetCurrentDepth(obj)
            % Return the current depth estimate for the estimator,
            % returning NaN if there are not enough observations.  This
            % method can be used for an array of objects.
            r = [obj.sum_depths] ./ [obj.prec];
        end
        function prec = GetCurrentPrecision(obj)
            % Return the current estimate precision for the estimator,
            % returning zero if there are no observations.  This method can
            % be used for an array of objects.
            prec = [obj.prec];
        end
        function stddev = GetCurrentStdDev(obj)
            % Return the current estimate of sample standard deviation of
            % the observations about the current mean, returning NaN if
            % there are not at least two observations.  This method can be
            % used for an array of objects.
            mean_depth = GetCurrentDepth(obj);
            stddev = sqrt([obj.sum_sq]./([obj.n_obs] - 1.0) - mean_depth.*mean_depth.*[obj.n_obs]./([obj.n_obs] - 1.0));
        end
        function nobs = GetObservationCount(obj)
            % Inspector for the current count of observations (which can
            % also be read directly if required).  This can be applied to
            % an array of objects
            nobs = [obj.n_obs];
        end
    end
end
