% FUNCTION summarise_reftile:
%
%   summarise_refttile(tile)
%
% Generate a figure that summarises the content of a RefinedTile,
% comprising the estimated refined depths, precisions, sample std. dev.,
% and observation count.

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

function summarise_reftile(tile)
    figure('WindowStyle', 'docked');
    d = tile.GetTileDepthEstimates;
    n = tile.GetTileObsCounts;
    x = tile.GetTileEstimateLocations;
    W = tile.TileWidth;
    S = tile.CellWidth;
    subplot(411);
    plot(x, d(1,:));
    DecorateSubplot('Depth (m)', S, 0, W);
    subplot(412);
    plot(x, d(2,:));
    DecorateSubplot('Precision (m^{-2})', S, 0, W);
    subplot(413);
    plot(x, d(3,:));
    DecorateSubplot('Sample Std. Dev. (m)', S, 0, W);
    subplot(414);
    plot(x, n);
    DecorateSubplot('Observation Count', S, 0, W);
end

function DecorateSubplot(yaxis_label, cell_width, min_x, max_x)
    axis([min_x max_x, -inf inf]);
    grid;
    set(gca, 'FontSize', 14);
    xlabel('Distance Within Tile (m)');
    n_divisions = (max_x - min_x) / cell_width;
    if n_divisions > 10
        multiplier = floor((max_x - min_x)/(10*cell_width));
        tick_width = multiplier*cell_width;
    else
        tick_width = cell_width;
    end
    set(gca, 'XTick', 0:tick_width:max_x);
    ylabel(yaxis_label);
end