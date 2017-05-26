% This is a base class for tiles to estimate data density and
% estimation resolution, or to do the high-resolution estimation of
% depth.

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

classdef Tile < handle & matlab.mixin.Heterogeneous
    properties (Access = private)
        Nc      % Number of SuperCells in the tile
        S       % Width of SuperCells (m)
    end
    methods (Abstract = true)
        % Add an observation to the tile at offset _x_ w.r.t. the tile's
        % left edge (m).
        AddObservation(obj, x, obs)
        % Report the number of estimates that are going to be generated
        % from the tile
        n = EstimateCount(obj)
        % Extract the depth estimates (and associated statistics) for the
        % whole tile (3xN array)
        d = GetTileDepthEstimates(obj)
        % Extract the positions of the depth estimates for the whole tile
        x = GetTileEstimateLocations(obj)
        % Extract the number of observations for the depth estimates for
        % the whole tile
        n = GetTileObsCounts(obj)
    end
    methods (Access = protected)
        function r = SuperCellCount(obj)
            % Inspector for the count of SuperCells in the tile
            r = obj.Nc;
        end
        function r = SuperCellWidth(obj)
            % Inspector for the width (in metres) of the SuperCells
            r = obj.S;
        end
    end
    methods
        function obj = Tile(ncells, cellwidth)
            % Core constructor for the class
            obj.Nc = ncells;
            obj.S = cellwidth;
        end
        function w = TileWidth(obj)
            % Compute the width (m) of the tile span.
            w = [obj.S] .* [obj.Nc];
        end
    end
    methods (Static, Sealed, Access = protected)
        function defaultObject = getDefaultScalarElement
            defaultObject = DummyTile;
        end
    end
end
