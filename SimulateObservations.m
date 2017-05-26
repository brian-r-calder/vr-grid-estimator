% Simulate a sequence of observations that might be taken by a
% single-beam echosounder traversing at fixed speed-through-water over
% a given bathymetry, assuming a constant (geometric mean) sound speed
% in the water.

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

classdef SimulateObservations < handle
    properties (GetAccess = private, SetAccess = private)
        sampleSpacing   % Sample spacing in the input bathymetry (m)
        z               % Depth/Variance pairs along-track (m, m^2)
        soundSpeed      % Geometric sound speed (m/s) for the simluation
        shipSpeed       % Ship speed through water (m/s) for the simulation
        beamwidth       % Beamwidth of the echosounder (rad)
        
        currentPosition = 0.0   % Position along track at current time (m)
        currentObservation = [] % Current observation (when valid)
    end
    methods
        function obj = SimulateObservations(sample_spacing, bathy_z, sound_speed, ship_speed, beamwidth)
            % Default constructor with basic sanity checks for parameters
            
            if length(sample_spacing) > 1 || sample_spacing < 0.01
                error('Inter-sample spacing in bathymetry must be above 0.01m.');
            end
            if size(bathy_z, 2) ~= 2
                error('Z must have two columns (depth, variance).');
            end
            if sound_speed < 1300 || sound_speed > 2000
                error('sound speed must be in range [1300, 2000] m/s.');
            end
            if ship_speed < 0.1
                error('ship speed must be above 0.1 m/s');
            end
            if beamwidth < 0.001
                error('beamwidth must be above 0.001 rad.');
            end
            
            obj.sampleSpacing = sample_spacing;
            obj.z = bathy_z;
            obj.soundSpeed = sound_speed;
            obj.shipSpeed = ship_speed;
            obj.beamwidth = beamwidth;
            obj.Reset;
        end
        function Reset(obj)
            % Reset the state of the simulator so that the sequence can be
            % restarted.  Note that the sequence uses the random number
            % generator, and therefore a different sequence will be
            % simulated unless the RNG seed is reset before starting to
            % generate samples again.
            
            obj.currentPosition = 0.0;
            depth = obj.SimulateDepth(obj.z(1,1), obj.z(1,2));
            obj.currentObservation = Observation(depth, obj.z(1,2), obj.beamwidth);
        end
        function Step(obj)
            % Step the simulator object forward one sample, generating a
            % plausible observation of depth, and its uncertainty estimate,
            % given the current sample.
            
            if isempty(obj.currentObservation)
                % Already generated all the samples, so the state stays the
                % same, and we just return.
                return
            end
            two_way_time = 2.0 * abs(obj.currentObservation.depth) / obj.soundSpeed;
            if two_way_time < 1/40.0
                two_way_time = 1/40.0; % Limit above at 40Hz ping rate
            end
            obj.currentPosition = obj.currentPosition + obj.shipSpeed * two_way_time;
            ind = floor(obj.currentPosition / obj.sampleSpacing) + 1;
            if ind > size(obj.z, 1)
                obj.currentPosition = -1.0;
                obj.currentObservation = [];
            else
                depth = obj.SimulateDepth(obj.z(ind,1), obj.z(ind, 2));
                obj.currentObservation = Observation(depth, obj.z(ind,2), obj.beamwidth);
            end
        end
        function b = IsComplete(obj)
            % Test whether all samples have been generated from the
            % simulation object (returns true if so).
            
            if obj.currentPosition < 0
                b = true;
            else
                b = false;
            end
        end
        function obs = GetCurrentObservation(obj)
            % Inspector for the current observation value.  An empty object
            % indicates that all samples for the simulation have been
            % generated (and you should call the Reset() method).
            
            obs = obj.currentObservation;
        end
        function x = GetCurrentPosition(obj)
            % Inspector for the current position along the bathymetry
            % signal used for the simulator object.  A value less than zero
            % indicates that all samples for the simulation have been
            % generated (and you should call the Reset() method).
            
            x = obj.currentPosition;
        end
    end
    methods (Access = private)
        function z = SimulateDepth(obj, mean, variance)
            % Generate one sample from a Gaussian distribution with given
            % mean and variance.
            
            z = randn(1)*sqrt(variance) + mean;
        end
    end
end