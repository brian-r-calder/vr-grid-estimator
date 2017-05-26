% This object encapsulates the information required to do estimation of
% data density, and refinement, across a single tile.

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

classdef EstimatorTile < Tile
    properties (GetAccess = public, SetAccess = private)
        nobservations   % Array for number of observations per cell
        occupancy       % 2D array for occupancy of each cell by observations
        lowres          % Object array for the low-res estimates in each cell
        
        lengthused      % Array for estimate of the super-cell length (m) used by observations
        obsdensity      % Array for estimate of observation density (m^-1)
        nodespacing     % Array for estimate of refined node spacing (m)
        nrefinednodes   % Array for count of refined nodes in each super-cell
    end
    properties (Access = private)
        B               % Number of bins in occupancy array for each super-cell
        maxGridScalar   % Maximum proportion of cell that can be used for a grid
        minResolution   % Minimum resolution at which to refine (m)
    end
    methods
        function obj = EstimatorTile(ncells,cellwidth,oversampling, min_resolution, estimator)
            % Default constructor for a fixed number of super-cells in the
            % Tile, with given cell width.  The oversampling parameter is
            % the number of bins used to stencil in the footprints of data
            % within each super-cell in order to determine how much of the
            % super-cell is considered "observed".
            
            obj@Tile(ncells, cellwidth);
            obj.B = oversampling;
            obj.minResolution = min_resolution;
            e(obj.SuperCellCount) = Estimator.Create(estimator);
            obj.lowres = e;
            obj.nobservations = zeros(1, obj.SuperCellCount);
            obj.occupancy = zeros(obj.B, obj.SuperCellCount);
            obj.maxGridScalar = ComputeGridScalar(obj);
        end
        
        function AddObservation(obj, x, obs)
            % Add a new observation of the field, with uncertainty, at the
            % specified location (which must be relative to the left hand
            % side of the tile - no georeferencing is stored in the tile).
            
            % Since we're adding new observations, the estimate for length
            % in cells used, the data density, and refinement node spacing
            % are about to be invalidated.
            if ~isempty(obj.lengthused)
                obj.lengthused = [];
                obj.obsdensity = [];
                obj.nodespacing = [];
                obj.nrefinednodes = [];
            end
            
            % We need to count the observation, and update the occupancy
            % first, then update the low-resolution estimate
            c = floor(x/obj.SuperCellWidth) + 1; % +1 for index into MATLAB arrays
            obj.nobservations(c) = obj.nobservations(c) + 1;
            incell_x = x - (c - 1)*obj.SuperCellWidth; % Position within the cell
            obj.UpdateOccupancy(c, incell_x, obs);
            obj.lowres(c).Update(obs);
            
            footprint_width = obs.FootprintWidth;
            if incell_x - footprint_width/2.0 < 0
                % Propagate left for occupancy
                prop_cell = c - 1;
                if prop_cell >= 1
                    prop_incell_x = x - (prop_cell - 1)*obj.SuperCellWidth;
                    obj.UpdateOccupancy(prop_cell, prop_incell_x, obs);
                end
            end
            if incell_x + footprint_width/2.0 > obj.SuperCellWidth
                % Propagate right for occupancy
                prop_cell = c + 1;
                if prop_cell <= obj.SuperCellCount
                    prop_incell_x = x - (prop_cell - 1)*obj.SuperCellWidth;
                    obj.UpdateOccupancy(prop_cell, prop_incell_x, obs);
                end
            end
        end
        
        function n = EstimateCount(obj)
            % Estimate the number of results that will be returned for an
            % extraction (i.e., the number of cells)
            
            n = obj.SuperCellCount;
        end
        
        function ComputeRefinements(obj, mean_obs, noise_factor)
            % Compute the length of cell observed in each cell, the data
            % density, the refinement node spacing, and the number of
            % nodes in each cell.
            
            obj.lengthused = sum(obj.occupancy)*(obj.SuperCellWidth/obj.B);
            obj.obsdensity = obj.nobservations ./ obj.lengthused;
            obj.nodespacing = mean_obs*(1 + noise_factor)./ obj.obsdensity;
            obj.nodespacing(~isfinite(obj.nodespacing)) = NaN;
            obj.nodespacing(isfinite(obj.nodespacing)) = min(obj.SuperCellWidth, obj.nodespacing(isfinite(obj.nodespacing)));
            obj.nrefinednodes = floor(obj.SuperCellWidth ./ obj.nodespacing) + 1;
            
            % We now need to determine whether the grid can fit within the
            % cell without touching the edges (which can lead to
            % overlapping nodes, causing difficulties later in the
            % estimation process).
            refined_spaces = obj.nrefinednodes - 1;
            grid_width = refined_spaces .* obj.nodespacing;
            for c = 1:obj.SuperCellCount
                if grid_width(c) > obj.SuperCellWidth * obj.maxGridScalar
                    % The grid in this cell would go into the boundary
                    % region around the edge of the cell, so we adjust the
                    % refined resolution to avoid this
                    obj.nodespacing(c) = obj.maxGridScalar * obj.SuperCellWidth / refined_spaces(c);
                    if obj.nodespacing(c) < obj.MinimumResolution
                        % This means we dropped through the minimum
                        % resolution allowed, so we need to adjust the node
                        % count if we can: the only way to fit in the
                        % number of nodes that we'd prefer is to drop
                        % through the min. res. (which we can't do) so the
                        % only solution is to reduce the number of nodes.
                        if refined_spaces(c) > 1
                            % We have nodes to remove
                            refined_spaces(c) = floor(obj.maxGridScalar*obj.SuperCellWidth/obj.MinimumResolution);
                            obj.nodespacing(c) = obj.maxGridScalar * obj.SuperCellWidth / refined_spaces(c);
                        else
                            % This means that we have no further moves to
                            % make, and should be rare.  We report as an
                            % error
                            error('cannot refine grid (is min_res = cell_width?');
                        end
                    end
                    % Getting to here means we changed the spacing count,
                    % so we need to recompute the number of refined nodes
                    % to match
                    obj.nrefinednodes(c) = floor(obj.SuperCellWidth / obj.nodespacing(c)) + 1;
                end
            end
        end
        
        function b = IsRefined(obj)
            % Determine whether the tile has had the refinement
            % computations done or not; returns true if so, otherwise,
            % false.
            
            if isempty(obj.nodespacing)
                b = false;
            else
                b = true;
            end
        end
        
        function r = MinimumResolution(obj)
            % Inspector for the minimum refinement resolution allowed
            
            r = obj.minResolution;
        end
        
        function g = MaxGridScalar(obj)
            % Inspector for the grid placement scalar
            
            g = obj.maxGridScalar;
        end
        
        function d = GetTileDepthEstimates(obj)
            % Extract the current depth, precision, and standard deviation
            % for all of the low-resolution estimators (one per cell) in
            % the tile.  The rows in the output matrix (3xN) are:
            %  1. Depth estimate (m)
            %  2. Depth precision (m^{-2})
            %  3. Depth sample std. dev. (m)
            
            d = zeros(3, obj.EstimateCount);
            d(1, :) = [obj.lowres.GetCurrentDepth];
            d(2, :) = [obj.lowres.GetCurrentPrecision];
            d(3, :) = [obj.lowres.GetCurrentStdDev];
        end
        
        function x = GetTileEstimateLocations(obj)
            % Extract the locations of all of the depth estimates being
            % generated (i.e., the centres of the cells) in the tile.
            
            x = obj.SuperCellWidth*(0:obj.SuperCellCount-1) + 0.5*obj.SuperCellWidth;
        end
        
        function n = GetTileObsCounts(obj)
            % Extract the count of observations per low-resolution depth
            % estimate (i,.e., one per cell) in the tile.
            
            n = obj.nobservations;
        end
    end
    
	methods(Access = private)
        function UpdateOccupancy(obj, cell, x, obs)
            % Compute the footprint of the observation on the super-cell,
            % and accumulate occupancy within the cell's over-sampled
            % accumluator bins.  Anti-aliasing (partial bin coverage) is
            % implemented with saturation accumulation.
            
            footprint_width = obs.FootprintWidth;
            
            b_left = (x - footprint_width/2.0)/(obj.SuperCellWidth/obj.B) + 1;
            b_right = (x + footprint_width/2.0)/(obj.SuperCellWidth/obj.B) + 1;
            
            left_bin = floor(b_left);
            right_bin = floor(b_right);
            
            partial_length_left = 1.0 - (b_left - left_bin);
            partial_length_right = b_right - right_bin;
            
            if left_bin == right_bin
                % All in the same bin, so update by the difference; since
                % the observation has to be in the cell somewhere, if this
                % happens, we're guaranteed that the bin's valid
                delta = b_right - b_left;
                obj.occupancy(left_bin, cell) = min(1.0, obj.occupancy(left_bin, cell) + delta);
            else
                if left_bin >= 1
                    obj.occupancy(left_bin, cell) = min(1.0, obj.occupancy(left_bin, cell) + partial_length_left);
                end
                if right_bin <= obj.B
                    obj.occupancy(right_bin, cell) = min(1.0, obj.occupancy(right_bin, cell) + partial_length_right);
                end
                if right_bin - left_bin - 2 >= 0
                    left_limit = max(1, left_bin+1);
                    right_limit = min(obj.B, right_bin-1);
                    obj.occupancy(left_limit:right_limit, cell) = 1.0;
                end
            end
        end
        
        function s = ComputeGridScalar(obj)
            % Compute a proportion of the cell that can be safely used for
            % refinement grids without getting too close to the edges.
            scale_constant = obj.ComputeScaleConstant;
            guard_fraction = scale_constant * obj.MinimumResolution / (2.0 * obj.SuperCellWidth);
            s = 1.0 - 2.0*guard_fraction;
        end
        
        function c = ComputeScaleConstant(obj)
            % Compute a safe scale constant to drive re-assessment of
            % resolution in the case where the ideal resolution would cause
            % the refinement grid to be too close to the edge of the cell
            % it's refining.
            n_probes = 1000;
            probe_quantum = 1.0 / n_probes;
            dx = ((obj.SuperCellWidth - obj.MinimumResolution) - obj.MinimumResolution)/(n_probes - 1);
            
            for probe = n_probes-1:-1:1
                k = probe * probe_quantum;
                test_const = obj.SuperCellWidth / (obj.SuperCellWidth - k * obj.MinimumResolution);
                failed = false;
                for x = 0:n_probes-1
                    p = test_const*floor(obj.SuperCellWidth/(x*dx + obj.MinimumResolution));
                    if p - floor(p) < 1.0e-6
                        failed = true;
                    end
                end
                if ~failed
                    c = k;
                    return;
                end
            end
            % If we get to here, it means that we went through all of the
            % probes and available constants and didn't find anything.
            % This should be vanishingly rare, so simply report an error
            error('failed to find scale constant for grid placement.');
        end
    end
end
