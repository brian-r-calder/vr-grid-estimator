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

% This file provides an example script to use the code provided in this
% archive, illustrating the intended usage as a jumping off point for
% further development.  Many other configurations of code are naturally
% possible.

%% Simulator Parameters
% This section provides parameters for setting up the scenario used to test
% the estimators, and is the most likely place for modification.

trajectory_length = 40000.0;        % Length of bathymetry to simulate (m)
min_depth = -300.0;                 % Minimum depth in the region (m)
max_depth = -5.0;                   % Maximum depth in the region (m)
bathy_roughness = 1.8;              % Roughness in range (1,2)
bathy_min_wavelength = 1.0;         % Minimum wavelength for the bathymetry (m)
bathy_sampling = 0.5;               % Bathymetric model resolution (m)
sound_speed = 1500.0;               % Geometric mean sound speed (m/s)
ship_speed = 5.0;                   % Ship speed through water (m/s)
sounder_beamwidth = deg2rad(10);    % Sounder beamwidth (rad)

%% Estimator Parameters
% This section provides parameters for setting up the estimator
cell_width = 32.0;      % Width of each low-resolution estimation cell (m)
min_spacing = 0.25;     % Minimum sample spacing to allow in refinement (m)
length_bins = 10;       % Number of bins used to estimate length of cell covered
cells_per_tile = 100;   % Number of cells in a single tile of estimators
base_estimator = 1;     % Type of estimator (1 = Mean, 2 = Weighted Mean)
mean_obs_reqd = 5.0;    % Number of observations, on average, to require
mean_obs_noise = 0.05;  % Noise factor to allow for observation counts

%% Debugging for sequences
if ~exist('rerun', 'var')
    % Setting 'rerun' true in the global space after a run of this script
    % allows the script to be re-run with exactly the same sequence of data
    % points, which can allow for easier debugging (or examination of
    % behaviour).
    rerun = false;
end

%% Debugging for observations
if ~exist('record', 'var')
    % Normally, the code discards the observations as they are used.
    % Setting 'record' true in the global space will cause the script to
    % record them instead.  Depending on the bathymetry generated, and the
    % simulation parameters, this could be a lot of data.
    record = false;
end

%% Debugging for EstimatorTiles
if ~exist('preserve', 'var')
    % Normally, the EstimatorTile objects used by the code to estimate the
    % observation density and thus the node spacing are replaced with the
    % RefinedTile objects when the refinement occurs.  Setting 'preserve'
    % to true in the global space copies them before they're replaced.
    preserve = false;
end

%% Simulator
% This section sets up the bathymetric model, and its observation
% simulator, which are used to generate plausible 'observations' for a
% fixed-beamwidth echosounder moving over the model.

if ~rerun
    rng('shuffle');     % Generate a different sequence on each run
    bathygen = GenerateBathymetry(trajectory_length, bathy_sampling);
    [bathy_x, bathy_z] = bathygen.Sample(bathy_roughness, min_depth, max_depth, bathy_min_wavelength);
    base_position = rand(1)*trajectory_length;
        % This generates a random starting point for the sequence,
        % demonstrating the estimator's ability to generate tiles on demand,
        % rather than having to pre-define the estimation domain.

    rng_state = rng;    % Store the RNG state before we start generating observations
    obssim = SimulateObservations(bathy_sampling, bathy_z, sound_speed, ship_speed, sounder_beamwidth);
end

%% Estimator Configuration
% This section sets up the basic estimator configuration, essentially
% making a TileManager and then configuring the refinement parameters
% required to translate the density of observations into a sample spacing.

estimator = TileManager(cells_per_tile, cell_width, length_bins, min_spacing, base_estimator);
estimator.SetRefinementParameters(mean_obs_reqd, mean_obs_noise);

%% Simulation
% This generates the data corresponding to the trajectory configured above,
% recording the outputs at each stage.  Note that since the data is
% generated randomly, but we need to replay this data after refinement, the
% code records the state of the PRNG before simulating, and resets to this
% after refinement to ensure that exactly the same data is generated.

if rerun
    rng(rng_state);
    obssim.Reset;
end
disp('Executing first pass over observations, estimating sample spacing ...');
n_observations = 0;
while ~obssim.IsComplete
    estimator.AddObservation(obssim.GetCurrentPosition + base_position, ...
                       obssim.GetCurrentObservation);
    n_observations = n_observations + 1;
    obssim.Step;
end
disp(['Added ', num2str(n_observations), ' observations to estimator.']);

% If requested, copy out the EstimatorTile objects before they're replaced
% with RefinedTile objects.
if preserve
    [start_tile, end_tile] = estimator.GetTileRange;
    n_tiles = end_tile - start_tile + 1;
    estimator_tiles(1,n_tiles) = estimator.GetTile(end_tile);
    for t = start_tile:end_tile-1
        estimator_tiles(t-start_tile+1) = estimator.GetTile(t);
    end
end

% If requested, set up arrays to record all of the observations that are
% going into the data sequence.
if record
    obs_positions = zeros(1, n_observations);
    obs_values(1,n_observations) = Observation(0, 1, 1);
end

% Get the results of the first stage of the process: the density estimates
% (and plot them)
disp('Getting sample spacing estimates ...');
density_ests = estimator.GetNodeSpacingEstimates;
display_density_estimates(density_ests);

%%
% Refine the tiles to allow for high-resolution estimation; this removes
% all of the old EstimatorTiles and replaces them with RefinedTiles.
disp('Refining all active tiles to give variable-spacing estimate slots ...');
estimator.RefineAllTiles;

% Reset the PRNG state so that we repeat the same data sequence, and
% re-apply the observations
rng(rng_state);
obssim.Reset;
disp('Executing second pass over observations, estimating depth ...');
if record
    obs = 1;
end
while ~obssim.IsComplete
    if record
        obs_positions(obs) = obssim.GetCurrentPosition + base_position;
        obs_values(obs) = obssim.GetCurrentObservation;
        obs = obs + 1;
    end
    estimator.AddObservation(obssim.GetCurrentPosition + base_position, ...
                       obssim.GetCurrentObservation);
    obssim.Step;
end

% Extract the high-resolution variable-spacing results, and display
disp('Getting high-resolution variable-spacing estimates ...');
refined_ests = estimator.GetDepthEstimates;
display_refined_estimates(refined_ests);

%% Display
% This plots a pair of auxiliary information windows to examine how well
% the estimator has matched the model bathymetry used to generate the
% simulated observations.

interp_bathy = interp1(bathy_x + base_position, bathy_z, ...
                       refined_ests(:,1), 'spline');
interp_bathy(:,2) = sqrt(max(0, interp_bathy(:,2)));
    % Convert to std. dev. and avoid any interpolation overshoot
max_uncrt = 1.05*1.96*max(interp_bathy(:,2));

figure('WindowStyle', 'docked');
subplot(211);
plot(refined_ests(:,1), refined_ests(:,3), ...
     refined_ests(:,1), interp_bathy(:,1), ...
     'LineWidth', 1.0);
grid; axis([refined_ests(1,1) refined_ests(end,1) min_depth max_depth]);
set(gca, 'FontSize', 14);
xlabel('Along-track Distance (m)');
ylabel('Depth (m)');
legend('Estimated', 'Model');
subplot(212);
plot(refined_ests(:,1), refined_ests(:,5), ...
     refined_ests(:,1), interp_bathy(:,2), ...
     'LineWidth', 1.0);
grid; axis([refined_ests(1,1) refined_ests(end,1) 0 max_uncrt]);
set(gca, 'FontSize', 14);
xlabel('Along-track Distance (m)');
ylabel('Standard Deviation (m)');
legend('Estimated', 'Model');

figure('WindowStyle', 'docked');
subplot(211);
plot(refined_ests(:,1), refined_ests(:,3), ...
     refined_ests(:,1), interp_bathy(:,1), ...
     'LineWidth', 1.0);
grid; axis([refined_ests(1,1) refined_ests(end,1) min_depth max_depth]);
set(gca, 'FontSize', 14);
xlabel('Along-track Distance (m)');
ylabel('Depth (m)');
legend('Estimated', 'Model');
subplot(212);
plot(refined_ests(:,1), refined_ests(:,3) - interp_bathy(:,1), ...
     refined_ests(:,1), 1.96*interp_bathy(:,2), 'r-', ...
     refined_ests(:,1), -1.96*interp_bathy(:,2), 'r-', ...
     'LineWidth', 1.0);
grid; axis([refined_ests(1,1) refined_ests(end,1) -max_uncrt max_uncrt]);
set(gca, 'FontSize', 14);
xlabel('Along-track Distance (m)');
ylabel('Depth Difference (m)');
legend('Difference', 'Model 95% CI');
