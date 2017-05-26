% Encapsulation of a single observation of the depth, and its
% uncertainty value and observation beamwidth.

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

classdef Observation
    properties
        depth       % Field value (m)
        uncrt       % Uncertainty as variance (m^2)
    end
    properties (GetAccess = private)
        beamwidth   % Beamwidth of observation system (rad)
    end
    methods
        function obj = Observation(depth, uncrt, beamwidth)
            % Default constructor for observation at given beamwidth
            if nargin == 0
                obj.depth = NaN;
                obj.uncrt = NaN;
                obj.beamwidth = NaN;
                return
            end
            
            obj.depth = depth;
            if uncrt <= 0
                error('uncertainty must be positive.');
            else
                obj.uncrt = uncrt;
            end
            if beamwidth <= 0 || beamwidth > pi/2.0
                error('beamwidth must be positive.');
            else
                obj.beamwidth = beamwidth;
            end
        end
        function r = FootprintWidth(obj)
            % Compute the beam footprint (m) for the observation
            r = 2.0*abs(obj.depth)*tan(obj.beamwidth/2.0);
        end
    end
end
