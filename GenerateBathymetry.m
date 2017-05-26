% Generator class for fractal-based synthetic bathymetry, used as test
% data.  The class generates samples of bathymetry with the same basic
% properties on each call of the Sample method.  Uncertainty is
% computed using the IHO S.44 model.

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

classdef GenerateBathymetry
    
    properties (GetAccess = public, SetAccess = private)
        nSamples            % Number of samples to generate in a trajectory
        sampleWidth         % Width of each sample (m)
        fixedUncert = 0.5   % Fixed portion of IHO-style uncertainty (m)
        varUncert = 0.013   % Depth-scaling portion of IHO-style uncertainty
   end
    properties (GetAccess = private, SetAccess = private)
        oversampleRate = 20	% Over-sampling rate for signal generation
    end
    methods
        function obj = GenerateBathymetry(trajectory_len, sample_rate)
            % Default constructor, computing the number of samples to
            % generate, and storing the sample spacing required.
            
            obj.nSamples = floor(trajectory_len / sample_rate) + 1;
            obj.sampleWidth = sample_rate;
        end
        function obj = SetUncertaintyParameters(obj, fixed, variable)
            % Set the IHO S.44 style parameters used to compute the
            % bathymetric uncertainty given depth.
            
            obj.fixedUncert = fixed;
            obj.varUncert = variable;
        end
        function len = TrajectoryLength(obj)
            % Inspector for the length of the trajectory to be generated (m).
            
            len = (obj.nSamples - 1)*obj.sampleWidth;
        end
        function [x, z] = Sample(obj, roughness, z_min, z_max, min_wavelength)
            % Generate a sample of the bathymetry of the appropriate
            % length, given a target roughness, dynamic range, and minimum
            % wavelength in the output signal.
            
            if roughness < 1 || roughness >= 2.0
                error('roughness must be in range [1,2).');
            end
            if z_min >= z_max
                error('dynamic range must have z_min < z_max.');
            end
            
            % We can get away with generating at twice the frequency of the
            % highest required frequency (which is determined from the
            % minimum wavelength to allow)
            max_reqd_frequency = 1.0 / min_wavelength;
            reqd_samp_frequency = 2.0*max_reqd_frequency;
            reqd_samples = floor(obj.nSamples*reqd_samp_frequency);
            
            % Check if we're going to be aliased (i.e., the output rate can
            % support the frequency demanded)
            if obj.sampleWidth > 1.0/reqd_samp_frequency
                error('bathymetry would be aliased with that minimum wavelength.');
            end
            
            % If we generate the power-spectrum fractal at exactly the
            % right number of samples, we'll end up with the first harmonic
            % dominating.  So we generate over a wider area (at the lower
            % resolution) to spread this out, and then sub-sample and
            % interpolate back to the output sample spacing
            gen_samples = floor(reqd_samples * obj.oversampleRate);
            signal = obj.compute_fractal(gen_samples, roughness);
            
            % Now we need to sub-sample out a suitable length of the signal
            % (avoiding the first harmonic problem) and interpolate up to
            % the right number of output samples
            sub_loc = floor(rand(1)*(gen_samples - reqd_samples)) + 1;
            z = zeros(obj.nSamples, 2);
            z(:,1) = interp1((1.0/reqd_samples)*(0:reqd_samples-1), signal(sub_loc:sub_loc+reqd_samples-1), ...
                             (1.0/obj.nSamples)*(0:obj.nSamples-1), 'spline');
            
            % Scale output in vertical to match the required dynamic range.
            z_drange = max(z(:,1)) - min(z(:,1));
            op_drange = z_max - z_min;
            z(:,1) = z_max - (op_drange/z_drange)*(z(:,1) - min(z(:,1)));
            
            % Bathymetric uncertainty is reported as variance, since this
            % is what's required for the Observation object and subsequent
            % estimation.  This avoids repeated square and square root
            % operations.  The 1.96^2 factor is due to IHO S.44 reporting
            % the 95% CI (assuming Gaussian statistics).
            z(:,2) = (obj.fixedUncert.^2 + (obj.varUncert*z(:,1)).^2)/1.96.^2;
            
            % Provide a location vector for the samples, to assist with
            % plotting and lookup.
            x = obj.sampleWidth*(0:obj.nSamples-1);
        end
    end
    methods (Access = private)
        function f = compute_fractal(obj, n, d)
            % Generate a power-spectrum fractal signal to the appropriate
            % sample count (n) and fractal dimension (d).
            
            beta = (5.0 - d)/2.0;
            
            spectrum = zeros(n,1);
            spectral_half_length = obj.compute_spectral_length(n);
            
            mag(1:spectral_half_length) = (1:spectral_half_length).^(-beta);
            ph(1:spectral_half_length) = -pi + 2.0*pi*rand(1,spectral_half_length);
            
            spectrum(2:2+spectral_half_length-1) = mag.*(cos(ph) + 1i*sin(ph));
            spectrum(end:-1:end-spectral_half_length) = conj(spectrum(2:2+spectral_half_length));
            
            f = ifft(spectrum);
        end
        function len = compute_spectral_length(obj, N)
            % Determine the half-way point for an FFT spectrum of length N.
            
            if mod(N,2) == 0
                len = (N-2)/2;
            else
                len = (N-1)/2;
            end
        end
    end
end
