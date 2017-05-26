% FUNCTION display_density_estimates:
%
%   display_density_estimates(d)
%
% Generate figures to illustrate the summary of density estimation results
% (from a call to TileManager.GetNodeSpacingEstimates).

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

function display_density_estimates(d)

figure('WindowStyle', 'docked');
subplot(311);
plot(d(:,1), d(:,2));
DecorateSubplot('Observations/Cell', d(1,1), d(end,1));
subplot(312);
plot(d(:,1), d(:,3));
DecorateSubplot('Length Used (m)', d(1,1), d(end,1));
subplot(313);
plot(d(:,1), d(:,4));
DecorateSubplot('Observation Density (m^{-1})', d(1,1), d(end,1));

figure('WindowStyle', 'docked');
subplot(211);
plot(d(:,1), d(:,5));
DecorateSubplot('Refined Node Spacing (m)', d(1,1), d(end,1));
subplot(212);
plot(d(:,1), d(:,6));
DecorateSubplot('Refined Node Count (m)', d(1,1), d(end,1));

figure('WindowStyle', 'docked');
subplot(311);
plot(d(:,1), d(:,7));
DecorateSubplot('Low-Resolution Depth (m)', d(1,1), d(end,1));
subplot(312);
plot(d(:,1), d(:,8));
DecorateSubplot('Low-Resolution Precision (m^{-2})', d(1,1), d(end,1));
subplot(313);
plot(d(:,1), d(:,9));
DecorateSubplot({'Low-Resolution Sample', 'Standard Deviation (m)'}, d(1,1), d(end,1));

end

function DecorateSubplot(ylabel_text, x_min, x_max)
    grid;
    axis([x_min x_max -inf inf]);
    set(gca, 'FontSize', 14);
    xlabel('Cell Absolute Location (m)');
    ylabel(ylabel_text);
end
