% FUNCTION display_refined_estimates:
%
%   display_refined_estimates(d)
%
% Generate a figure to summarise the high-resolution (variable spacing)
% depth estimates generated by a call to TileManager.GetDepthEstimates.

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

function display_refined_estimates(d)
    figure('WindowStyle', 'docked');
    subplot(411);
    plot(d(:,1), d(:,2));
    DecorateSubplot('Observation Count', d(1,1), d(end,1));
    subplot(412);
    plot(d(:,1), d(:,3));
    DecorateSubplot('Depth (m)', d(1,1), d(end,1));
    subplot(413);
    plot(d(:,1), d(:,4));
    DecorateSubplot('Precision (m^{-2})', d(1,1), d(end,1));
    subplot(414);
    plot(d(:,1), d(:,5));
    DecorateSubplot('Sample Std. Dev. (m)', d(1,1), d(end,1));
end

function DecorateSubplot(ylabel_text, x_min, x_max)
    grid;
    axis([x_min x_max -inf inf]);
    set(gca, 'FontSize', 14);
    xlabel('Along-track Distance (m)');
    ylabel(ylabel_text);
end