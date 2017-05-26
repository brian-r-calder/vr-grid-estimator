% FUNCTION summarise_esttile:
%
%   summarise_esttile(tile)
%
% Generate a figure that summarises the content of a (refined)
% EstimateTile, comprising the length used per tile, observation density,
% node count and number of observations.

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

function summarise_esttile(tile)
    figure('WindowStyle', 'docked');
    subplot(411);
    plot([tile.nobservations]);
    DecorateSubplot('Number of Observations');
    subplot(412);
    plot([tile.lengthused]);
    DecorateSubplot('Length Used (m)');
    subplot(413);
    plot([tile.obsdensity]);
    DecorateSubplot('Observation Density (m^{-1})');
    subplot(414);
    plot([tile.nodespacing]);
    DecorateSubplot('Refined Node Spacing (m)');
   
    figure('WindowStyle', 'docked');
    d = tile.GetTileDepthEstimates;
    subplot(411);
    plot([tile.nrefinednodes]);
    DecorateSubplot('Number of Refined Nodes');
    subplot(412);
    plot(d(1,:));
    DecorateSubplot('Low-Res Depth (m)');
    subplot(413);
    plot(d(3,:));
    DecorateSubplot('Low-Res Std. Dev. (m)');
    subplot(414);
    plot(d(2,:));
    DecorateSubplot('Low-Res Precision (m^{-2})');
end

function DecorateSubplot(yaxis_label)
    axis([0 100 -inf inf]);
    grid;
    set(gca, 'FontSize', 14);
    xlabel('Cell Index Number');
    ylabel(yaxis_label);
end