% This implements the operations of a tile of refined piecewise
% constant grids, starting from the EstimatorTile containing the
% estimates of node spacing and node count for the various cells in the
% same tile space.
    
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

classdef RefinedTile < Tile
    properties (Access = private)
        samplespacing       % Estimate node spacing (m) per cell
        nodecount           % Count of nodes in each cell
        nodes               % Array for all estimator nodes for the tile
        startnode           % Index array for where to find the start of each SuperCell
        
        default_mean_obs = 5.0      % Default mean observation count
        default_noise_factor = 0.0  % Default noise factor
    end
    methods
        function obj = RefinedTile(esttile, estimator)
            % Default constructor for a refined tile, based on an
            % EstimatorTile input structure.  The basic shape of the tile
            % is copied, but the refined grids are created instead of the
            % occupancy estimates.
            obj@Tile(esttile.SuperCellCount, esttile.SuperCellWidth);
            if ~esttile.IsRefined
                esttile.ComputeRefinements(obj.default_mean_obs, obj.default_noise_factor);
            end
            obj.samplespacing = esttile.nodespacing;
            obj.nodecount = esttile.nrefinednodes;
            total_nodes = sum(obj.nodecount(~isnan(obj.nodecount)));
            if total_nodes < 1
                error('no refined nodes specified in source EstimatorTile.');
            end
            estimator_array(total_nodes) = Estimator.Create(estimator);
            obj.nodes = estimator_array;
            counts = zeros(size(obj.nodecount));
            counts(~isnan(obj.nodecount)) = obj.nodecount(~isnan(obj.nodecount));
            obj.startnode = [1 cumsum(counts(1:end-1))+1];
            obj.startnode(isnan(obj.nodecount)) = NaN;
        end
        
        function AddObservation(obj, x, obs)
            % Add an observation to the refined grids in the tile,
            % propagating between SuperCells if required.  Location of the
            % observation (class Observation) must be relative to the left
            % edge of the tile.
            target_cell = floor(x / obj.SuperCellWidth);
                % Cell in which the observation falls.  Note, 0-based.
                
            if isnan(obj.nodecount(target_cell+1))
                % If there's data there, the cell should have been refined
                % in the EstimateTile stage; this is therefore a user
                % error.
                error('attempt to add an observation to an un-refined cell.');
            end
            
            [target_index, in_cell_offset] = ComputeTargetIndex(obj, target_cell, x);
            if target_index < 0
                return
            end
            
            obj.nodes(target_index).Update(obs);
                % Carries out the update for the primary node.
            
            % We now need to check whether the observation was sufficiently
            % close to the edge of the cell, and if so, offer it to the
            % neighbouring cells too.  If we don't do this, we'll get edge
            % effects on the cell edges.  Note that we only make the offer
            % if (a) there is a neighbouring cell, and (b) the cell has
            % been refined.  Note that propagation between tiles is done at
            % the TileManager level, so we don't propagate left from the
            % first cell in the tile, or right from the last cell.
            if target_cell > 0 && ...
                    in_cell_offset < 0.5*obj.samplespacing(target_cell) && ...
                    ~isnan(obj.nodecount(target_cell))
                [target_index, ~] = ComputeTargetIndex(obj, target_cell-1, x);
                if target_index > 0
                    obj.nodes(target_index).Update(obs);
                end
            end
            if target_cell < obj.SuperCellCount - 1 && ...
                    in_cell_offset > obj.SuperCellWidth - 0.5*obj.samplespacing(target_cell+2) && ...
                    ~isnan(obj.nodecount(target_cell+2))
                [target_index, ~] = ComputeTargetIndex(obj, target_cell+1, x);
                if target_index > 0
                    obj.nodes(target_index).Update(obs);
                end
            end
        end
        
        function n = EstimateCount(obj)
            % Count the number of estimates that will be provided by the
            % tile on extraction
            
            n = length(obj.nodes);
        end
        
        function [z, p, s] = GetCellEstimates(obj, target_cell)
            % Extract the depth, precision, and sample std. dev. estimates
            % for all of the nodes within a given cell in the tile.
            
            if isnan(obj.nodecount(target_cell+1))
                error('that cell contains no refined nodes.');
            end
            
            [start_node, end_node] = GetCellNodeBounds(obj, target_cell);
            z = [obj.nodes(start_node:end_node).GetCurrentDepth];
            p = [obj.nodes(start_node:end_node).GetCurrentPrecision];
            s = [obj.nodes(start_node:end_node).GetCurrentStdDev];
        end
        
        function d = GetTileDepthEstimates(obj)
            % Extract the depth, precision, and sample std. dev. estimates
            % for all of the nodes in all cells within the tile.
            
            d = zeros(3, obj.EstimateCount);
            d(1,:) = [obj.nodes.GetCurrentDepth];
            d(2,:) = [obj.nodes.GetCurrentPrecision];
            d(3,:) = [obj.nodes.GetCurrentStdDev];
        end
        
        function n = GetCellObsCounts(obj, target_cell)
            % Extract the observation counts for each node within a given
            % cell in the tile.
            
            if isnan(obj.nodecount(target_cell+1))
                error('that cell contains no refined nodes.');
            end
            
            [start_node, end_node] = GetCellNodeBounds(obj, target_cell);
            n = [obj.nodes(start_node:end_node).GetObservationCount];
        end
        
        function n = GetTileObsCounts(obj)
            % Extract the observation counts for all nodes within all cells
            % in the tile.
            
            n = [obj.nodes.GetObservationCount];
        end
        
        function x = GetCellEstimateLocations(obj, target_cell)
            % Compute the locations, relative to the left edge of the tile,
            % for all of the nodes in a given cell within the tile.
            
            if isnan(obj.nodecount(target_cell+1))
                error('that cell contains no refined nodes.');
            end
            
            cell_centre = (target_cell + 0.5) * obj.SuperCellWidth;
            grid_width = (obj.nodecount(target_cell + 1) - 1)*obj.samplespacing(target_cell + 1);
            grid_left = cell_centre - grid_width*0.5;
            x = grid_left + (0:obj.nodecount(target_cell + 1)-1)*obj.samplespacing(target_cell + 1);
        end
        
        function x = GetTileEstimateLocations(obj)
            % Compute the locations, relative to the left edge of the tile,
            % for all nodes in all cells within the tile.
            
            x = zeros(1, length(obj.nodes));
            for target_cell = 0:obj.SuperCellCount-1
                if isnan(obj.nodecount(target_cell+1))
                    continue;
                end
                [start_node, end_node] = GetCellNodeBounds(obj, target_cell);
                x(start_node:end_node) = GetCellEstimateLocations(obj, target_cell);
            end
        end
        
        function w = CellWidth(obj)
            w = obj.SuperCellWidth;
        end
    end
    methods (Access = private)
        function [target_index, in_cell_offset] = ComputeTargetIndex(obj, target_cell, x)
            % Compute the index within the composite node array for a
            % particular offset within a target cell.  Note that the offset
            % given here should be relative to the tile's left edge.
            
            cell_left = obj.SuperCellWidth * target_cell;
                % Offset to the left edge of the target cell (m).
            in_cell_offset = x - cell_left;
                % Offset with respect to the left edge of the cell holding
                % the observation
            target_node = floor((in_cell_offset - 0.5*obj.SuperCellWidth)/obj.samplespacing(target_cell+1) + 0.5*obj.nodecount(target_cell+1));
                % Target node within the refined grid in the SuperCell
            if target_node < 0 || target_node >= obj.nodecount(target_cell+1)
                % Nothing to do here: it's outside the refined grid.  This
                % can sometimes happen when data are propagated from other
                % cells to here.
                target_index = -1;
                return
            end
            target_index = target_node + obj.startnode(target_cell+1);
                % Index into the storage array (note that this is already
                % set to be 1-based, and doesn't need incremented for index
                % into the computation array.
        end
        
        function [s,e] = GetCellNodeBounds(obj, target_cell)
            % Compute the start and end index for a given taget cell's
            % nodes in the composite array, so they can be manipulated as a
            % group.
            
            s = obj.startnode(target_cell+1);
            e = s + obj.nodecount(target_cell+1) - 1;
        end
    end
end
