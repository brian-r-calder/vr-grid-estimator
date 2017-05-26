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

classdef Estimator < handle
    % Abstract superclass for any estimator that the code can use
    properties
    end
    methods (Abstract = true)
        % Add new observation (class Observation) from obs to estimator
        obj = Update(obj, obs)
        % Return current estimate of depth
        r = GetCurrentDepth(obj)
        % Return current estimate of depth precision
        r = GetCurrentPrecision(obj)
        % Return current estimate of sample standard deviation
        r = GetCurrentStdDev(obj)
        % Return the number of observations used
        r = GetObservationCount(obj)
    end
    methods(Static)
        function obj = Create(id)
            % Create an instance (in this case a handle) for a given class
            % of estimator object.
            
            switch id
                case 1
                    obj = MeanEstimator;
                case 2
                    obj = WeightedMeanEstimator;
                otherwise
                    error('Estimator ID not recognised.');
            end
        end
    end
end
