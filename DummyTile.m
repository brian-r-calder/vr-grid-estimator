classdef DummyTile < Tile
    methods
        function AddObservation(obj, x, obs)
            error('you cannot add observations to a DummyTile.');
        end
        function n = EstimateCount(obj)
            n = obj.SuperCellCount;
        end
        function d = GetTileDepthEstimates(obj)
            d = [];
        end
        function x = GetTileEstimateLocations(obj)
            x = [];
        end
        function n = GetTileObsCounts(obj)
            n = [];
        end
        function obj = DummyTile
            obj@Tile(1, 1);
        end
    end
end
