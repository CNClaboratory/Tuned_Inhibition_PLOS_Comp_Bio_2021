function [perf, traces] = TI_sim_trial(param)
% [perf, traces] = TI_sim_trial(param)
%
% Simulate a trial using the Differential Tuned Inhibition model of 
% Maniscalco, Odegaard, Grimaldi, Cho, Basso, Lau, & Peters 2021 PLOS Comp Bio
%
% INPUT
% -----
%
% The input struct "param" must have the following fields:
%
% * param.T     - decision threshold
% * param.sigma - standard deviation of instantaneous noise
% * param.S_i   - a (2 x 1) vector containing instantaneous stimulus evidence 
%                 for stimulus alternatives 1 and 2
% * param.tau   - number of time steps for post-decision evidence accumulation
% 
% "param" may also have the following optional fields:
%
% * param.tmax [default value = 1e4]
% maximum # of time steps the simulation will run before returning
% 
% * param.exit_on_response [default value = 1]
% if 1, simulation returns once perceptual decision and confidence are determined
% if 0, simulation runs for tmax timesteps
%
% * param.save_traces [default value = 0]
% if 1, the "traces" output contains data
% if 0, the "traces" output is empty
%
%
% OUTPUT
% ------
%
% The "perf" output struct contains the following fields:
%
% * perf.resp
% the perceptual decision
% 0 --> chose stimulus alternative 1
% 1 --> chose stimulus alternative 2
% 
% * perf.rt
% the time step at which the perceptual decision was made
% (by convention, the trial begins at t = 1)
% 
% * perf.conf_Cx
% confidence in the perceptual decision as computed from the Cx model.
% this is computed as the value of the absolute evidence accumulator x
% corresponding to the perceptual decision at time t = rt + tau
% 
% * perf.conf_Cx_i
% a (2 x 1) vector holding the values of both absolute evidence accumulators 
% x_i at time t = rt + tau, such that conf_Cx = conf_Cx_i(resp+1).
% this output allows for inspection of the state of x_i for the non-selected
% stimulus alternative at the time of confidence rating.
% 
% * perf.conf_Cd
% confidence in the perceptual decision as computed from the Cd model.
% this is computed as the value of the relative evidence accumulator delta
% corresponding to the perceptual decision at time t = rt + tau
% 
% * perf.conf_Cd_i
% a (2 x 1) vector holding the values of both relative evidence accumulators 
% delta_i at time t = rt + tau, such that conf_Cd = conf_Cd_i(resp+1).
% this output allows for inspection of the state of delta_i for the non-selected
% stimulus alternative at the time of confidence rating.
%
% NOTE: all output fields listed above default to values of -1 if the
% corresponding perceptual decision or confidence value was not registered
% within tmax time steps.
% 
% * perf.rt2
% the time step at which confidence was determined. by definition, rt2 = rt + tau
%
% * perf.responded
% 0 --> no perceptual decision was registered within tmax time steps
% 1 --> a perceptual decision was successfully registered
%
% * perf.responded2
% 0 --> no confidence value was registered within tmax time steps
% 1 --> a confidence value was successfully registered
%
% The "traces" output struct contains the following fields:
%
% * traces.x_it
% a (2 x tmax) matrix containing the state of x_i at each time step t.
%
% * traces.delta_it
% a (2 x tmax) matrix containing the state of delta_i at each time step t.
%
% if param.exit_on_response = 1, then time steps after rt+tau are not
% simulated and are thus registered in these outputs as NaN.
% if param.save_traces = 0, then traces.x_it and traces.delta_it are empty.
%
% 2/23/2021 Brian Maniscalco & Megan Peters

%% handle inputs

% check existence of mandatory inputs
if ~all(isfield(param,{'T','sigma','S_i','tau'}))
    error('param is missing a mandatory field; see help TI_sim_trial')
end

% check size of param.S_i
if ~all(size(param.S_i) == [2,1])
    error('param.S_i must be a 2 x 1 vector')
end

% set default values for optional inputs
if ~isfield(param, 'tmax')
    param.tmax = 1e4;
end

if ~isfield(param, 'exit_on_response')
    param.exit_on_response = 1;
end

if ~isfield(param, 'save_traces')
    param.save_traces = 0;
end


%% initialize simulation

% initialize evidence accumulation units
x_i     = zeros(2, 1); % absolute evidence accumulators
delta_i = zeros(2, 1); % instantaneous evidence difference

% initialize performance output
perf.resp       = -1;
perf.rt         = -1;
perf.conf_Cx    = -1;
perf.conf_Cx_i  = [-1, -1];
perf.conf_Cd    = -1;
perf.conf_Cd_i  = [-1, -1];
perf.rt2        = -1;
perf.responded  = 0;
perf.responded2 = 0;

% initialize traces
if param.save_traces
    traces.x_it     = [x_i,     nan(2, param.tmax-1)];
    traces.delta_it = [delta_i, nan(2, param.tmax-1)];
else
    traces.x_it     = [];
    traces.delta_it = [];
end


%% run simulation

t  = 2;
while t <= param.tmax && ~(perf.responded && perf.responded2 && param.exit_on_response)
    
    % update accumulators
    x_noise     = normrnd(0, param.sigma, [2, 1]);
    delta_noise = normrnd(0, param.sigma, [2, 1]);
 
    dx_i = param.S_i + x_noise;
    x_i  = max( x_i + dx_i, 0 );
    
    % compute instantaneous evidence difference
    delta_i(1) = max(x_i(1) - x_i(2) + delta_noise(1), 0 );
    delta_i(2) = max(x_i(2) - x_i(1) + delta_noise(2), 0 );

    if param.save_traces
        traces.x_it(:,t)     = x_i;
        traces.delta_it(:,t) = delta_i;
    end
    
    % check for type 1 response
    if ~perf.responded
        
        if delta_i(1) > param.T && delta_i(1) > delta_i(2)
            perf.resp      = 0;
            perf.rt        = t;
            perf.responded = 1;
            
        elseif delta_i(2) > param.T && delta_i(2) > delta_i(1)
            perf.resp      = 1;
            perf.rt        = t;
            perf.responded = 1;
            
        end
    end
    
    % check for type 2 response
    if perf.responded
        
        t2 = t - perf.rt;
        
        if t2 == param.tau
            
            % conf from x units
            perf.conf_Cx   = x_i(perf.resp+1);
            perf.conf_Cx_i = x_i;

            % conf from delta units
            perf.conf_Cd   = delta_i(perf.resp+1);
            perf.conf_Cd_i = delta_i;

            perf.rt2        = t;
            perf.responded2 = 1;
        end
        
    end
    
    t = t + 1;
end

end