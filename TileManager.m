% Manage an array of Tile objects to allow for low-resource estimation
% across large areas.

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

classdef TileManager < handle
    properties (GetAccess = public, SetAccess = private)
        cellsPerTile    % Number of cells in a tile
        cellWidth       % Width of each cell (m)
        lengthBinCount  % Number of bins used to compute length observed in cells
        minResolution   % Minimum refinement resolution (m)
        depthEstimator  % Estimator ID for depth estimation
        
        reqd_obs_mean = 5;      % Mean number of observations required at each node
        reqd_obs_noise = 0.0;  	% Noise factor for translation to sample spacing
    end
    properties (Access = private)
        startTile = []  % First tile that's active
        endTile = []    % Last tile that's active
        tiles = []      % Array of Tile objects, expanding as required
    end
    methods
        function obj = TileManager(ncells, cellwidth, oversampling, min_resolution, estimator)
            % Default constructor, which just stores the parameters for
            % creating tiles until they're required.
            
            obj.cellsPerTile = ncells;
            obj.cellWidth = cellwidth;
            obj.lengthBinCount = oversampling;
            obj.minResolution = min_resolution;
            obj.depthEstimator = estimator;
        end
        
        function SetRefinementParameters(obj, mean, noise)
            % Set parameters for conversion from density to sample spacing
            
            obj.reqd_obs_mean = mean;
            obj.reqd_obs_noise = noise;
        end
        
        function AddObservation(obj, x, obs)
            % Add an observation to a tile, making new tiles if required in
            % order to have the space to do it.  Newly created tiles are
            % always made as EstimateTile objects.
            
            if x < 0
                error('Along-track observation positions must be positive.');
            end
            
            target_tile = floor(x / obj.TileWidth) + 1;
                % The first tile has its left edge at the zero of the
                % coordinate system; note tile numbers are 1-based to avoid
                % having to offset for MATLAB's indexing.
            
            if isempty(obj.startTile) || target_tile < obj.startTile || target_tile > obj.endTile
                % Outside of the current active area, so we need to
                % re-arrange the tiles and extend
                if isempty(obj.startTile)
                    new_left = target_tile;
                    new_right = target_tile;
                else
                    new_left = min(target_tile, obj.startTile);
                    new_right = max(target_tile, obj.endTile);
                end
                new_size = new_right - new_left + 1;
                
                % Allocate an array of the appropriate size using a dummy
                % EstimateTile with the right parameters.
                dummy_tile = DummyTile;
                t(1, new_size) = dummy_tile;
                
                % Copy over the old tiles in the appropriate place (if
                % there are old tiles to copy)
                if ~isempty(obj.startTile)
                    first_tile = obj.startTile - new_left + 1;
                    last_tile = first_tile + length(obj.tiles) - 1;
                    t(first_tile:last_tile) = obj.tiles;
                
                    % Make sure any new tiles to left or right are initialised
                    % correctly to an EstimatorTile of appropriate
                    % dimensions.  Because all Tile sub-classes are handle
                    % classes, we need to make a new Estimator tile for
                    % each entry, rather than just copying a single entity
                    % into each entry.
                    if new_left < obj.startTile
                        t0 = 1;
                        t1 = first_tile-1;
                    else
                        t0 = last_tile+1;
                        t1 = length(t);
                    end
                    for i = t0:t1
                        t(i) = EstimatorTile(obj.cellsPerTile, obj.cellWidth, obj.lengthBinCount, obj.minResolution, obj.depthEstimator);
                    end
                else
                    % This is the first tile being added, so we initialise
                    % it to an EstimateTile.
                    t(1) = EstimatorTile(obj.cellsPerTile, obj.cellWidth, obj.lengthBinCount, obj.minResolution, obj.depthEstimator);
                end
                
                % Copy the array back, and update markers
                obj.tiles = t;
                obj.startTile = new_left;
                obj.endTile = new_right;
            end
            
            % Now that we have the array arranged, all we need to do is to
            % pass the observation to the appropriate tile, offsetting for
            % the geo-refencing of the left side of the target tile, since
            % the tiles don't store their own georeferencing information,
            % and estimate w.r.t. the left hand edge.
            tile_left = obj.TileGeoref(target_tile);
            idx = target_tile - obj.startTile + 1;
            obj.tiles(idx).AddObservation(x - tile_left, obs);
        end
        
        function [s,e] = GetTileRange(obj)
            % Inspector for the range of the currently defined tiles.
            
            s = obj.startTile;
            e = obj.endTile;
        end
        
        function t = GetTile(obj, target_tile)
            % Extract a handle for a given tile in the array.  Since this
            % is a handle, destructive operations should be avoided in
            % order to preserve the integrity of the TileManager data
            % structures.  Extracting estimates works fine, however.
            
            obj.ValidateTile(target_tile);
            t = obj.tiles(target_tile - obj.startTile + 1);
        end
        
        function b = IsRefined(obj, target_tile)
            % Test whether the target tile has been refined, or not.
            
            obj.ValidateTile(target_tile);
            b = isa(obj.tiles(target_tile - obj.startTile + 1), 'RefinedTile');
        end
        
        function RefineTileRange(obj, start_tile, end_tile)
            % Refine any EstimatorTile object in the given range
            
            obj.ValidateTileRange(start_tile, end_tile);
            
            first_tile = start_tile - obj.startTile + 1;
            last_tile = end_tile - obj.startTile + 1;
            
            for t = first_tile:last_tile
                if isa(obj.tiles(t), 'EstimatorTile')
                    if ~obj.tiles(t).IsRefined
                        obj.tiles(t).ComputeRefinements(obj.reqd_obs_mean, obj.reqd_obs_noise);
                    end
                    r = RefinedTile(obj.tiles(t), obj.depthEstimator);
                    obj.tiles(t) = r;
                end
            end
        end
        
        function RefineAllTiles(obj)
            % Convenience method to refine all of the tiles in the cache.
            
            obj.RefineTileRange(obj.startTile, obj.endTile);
        end
        
        function n = EstimateCount(obj)
            % Count the number of estimates from all tiles in the cache.
            
            n = 0;
            for t = 1:length(obj.tiles)
                n = n + obj.tiles(t).EstimateCount;
            end
        end
        
        function r = GetNodeSpacingEstimates(obj, start_tile, end_tile)
            % Get all of the EstimatorTile sample-spacing estimates in the
            % given tile range.  Columns in the result array are:
            %  1. Position of the estimate.
            %  2. Number of observations used in the estimate.
            %  3. Length of cell covered by observations (m).
            %  4. Observation density (m^{-1}).
            %  5. Refinement node spacing (m).
            %  6. Number of refined nodes.
            %  7. Depth estimate for all observations in the cell (m).
            %  8. Depth precision estimate for all observations (m^{-2}).
            %  9. Depth sample std. dev. for all observations (m).
            
            if nargin < 2
                start_tile = obj.startTile;
                end_tile = obj.endTile;
            end
            
            obj.ValidateTileRange(start_tile, end_tile);
            first_tile = start_tile - obj.startTile + 1;
            last_tile = end_tile - obj.startTile + 1;
            tile_range = first_tile:last_tile;
            
            % Run through all of the tiles and find the total number of
            % counts in all the EstimatorTiles
            n_total = 0;
            for t = tile_range
                if isa(obj.tiles(t), 'EstimatorTile')
                    n_total = n_total + obj.tiles(t).EstimateCount;
                    if ~obj.tiles(t).IsRefined
                        obj.tiles(t).ComputeRefinements(obj.reqd_obs_mean, obj.reqd_obs_noise);
                    end
                end
            end
            
            r = zeros(n_total, 9);
            
            out_row = 1;
            for t = tile_range
                if ~isa(obj.tiles(t), 'EstimatorTile')
                    continue;
                end
                count = obj.tiles(t).EstimateCount;
                out_range = out_row:out_row+count-1;
                r(out_range, 1) = obj.tiles(t).GetTileEstimateLocations' + obj.TileGeoref(t + obj.startTile - 1);
                r(out_range, 2) = obj.tiles(t).nobservations';
                r(out_range, 3) = obj.tiles(t).lengthused';
                r(out_range, 4) = obj.tiles(t).obsdensity';
                r(out_range, 5) = obj.tiles(t).nodespacing';
                r(out_range, 6) = obj.tiles(t).nrefinednodes';
                r(out_range, 7:9) = obj.tiles(t).GetTileDepthEstimates';
                out_row = out_row + count;
            end
        end
        
        function r = GetDepthEstimates(obj, start_tile, end_tile)
            % Extract depth estimates from all tiles in the current cache.
            % For EstimatorTile objects, this is the low-resolution (fixed
            % spacing) one-per-cell estimates; for RefinedTile objects,
            % this is the high-resolution (variable spacing) estimate set.
            % Columns in the output matrix are:
            %  1. Estimate absolute location (m).
            %  2. Count of observations used in the estimate
            %  3. Depth estimate (m).
            %  4. Depth precision (m^{-2}).
            %  5. Depth sample std. dev. (m).
            
            if nargin < 2
                start_tile = obj.startTile;
                end_tile = obj.endTile;
            end
            
            obj.ValidateTileRange(start_tile, end_tile);
            first_tile = start_tile - obj.startTile + 1;
            last_tile = end_tile - obj.startTile + 1;
            tile_range = first_tile:last_tile;

            n_est = obj.EstimateCount;
            r = zeros(n_est, 5);
            
            out_row = 1;
            for t = tile_range
                count = obj.tiles(t).EstimateCount;
                out_range = out_row:out_row+count-1;
                r(out_range, 1) = obj.tiles(t).GetTileEstimateLocations' + obj.TileGeoref(t + obj.startTile - 1);
                r(out_range, 2) = obj.tiles(t).GetTileObsCounts';
                r(out_range, 3:5) = obj.tiles(t).GetTileDepthEstimates';
                out_row = out_row + count;
            end
        end
        
        function Dimensions(obj)
            % Report some statistics about the current configuration of the
            % TileManager's internal tile cache.
            
            if isempty(obj.startTile)
                disp('No tiles have been initialised in cache.');
                return;
            end
            
            disp(['Currently with tile range [', ...
                  num2str(obj.startTile, '%d'), ', ', num2str(obj.endTile, '%d'), ...
                  '] (', num2str(length(obj.tiles), '%d'), ' tiles).']);
            n_estimators = 0;
            n_refined = 0;
            n_dummy = 0;
            n_unknown = 0;
            for i = 1:length(obj.tiles)
                if isa(obj.tiles(i), 'EstimatorTile')
                    n_estimators = n_estimators + 1;
                else
                    if isa(obj.tiles(i), 'RefinedTile')
                        n_refined = n_refined + 1;
                    else
                        if isa(obj.tiles(i),'DummyTile')
                            n_dummy = n_dummy + 1;
                        else
                            n_unknown = n_unknown + 1;
                        end
                    end
                end
            end
            disp(['With ', num2str(n_estimators, '%d'), ' EstimatorTiles, ', ...
                 num2str(n_refined, '%d'), ' RefinedTiles, ', ...
                 num2str(n_dummy, '%d'), ' DummyTiles, and ', ...
                 num2str(n_unknown, '%d'), ' Unrecognised Tile sub-classes.']);
        end
        
        function x = GetTileGeoref(obj, target_tile)
            % Parameter-checked method to extract the georeferencing point
            % for a given tile (i.e., the west edge of the western-most
            % cell in the tile).
            
            obj.ValidateTile(target_tile);
            x = obj.TileGeoref(target_tile);
        end
    end
    
    methods (Access = private)
        function w = TileWidth(obj)
            % Compute the width (m) of any tiles that are going to be made
            w = obj.cellsPerTile * obj.cellWidth;
        end
        
        function x = TileGeoref(obj, target_tile)
            % Compute the location (m) along-track for the given target
            % tile number (not that target tile numbers are 1-based, but
            % the geo-referencing requires the first tile at zero offset).
            x = obj.TileWidth * (target_tile - 1);
        end
        
        function ValidateTile(obj, target_tile)
            % Test whether a target tile is present in the cache
            
            if target_tile < obj.startTile || target_tile > obj.endTile
                error('that tile is not currently defined.');
            end
        end
        
        function ValidateTileRange(obj, start_tile, end_tile)
            % Test whether a target tile range is all present in cache
            
            if start_tile > end_tile
                error('tile range is invalid.');
            end
            if start_tile < obj.startTile || end_tile > obj.endTile
                error('that tile range is not currently wholy defined.');
            end
        end
    end
end
