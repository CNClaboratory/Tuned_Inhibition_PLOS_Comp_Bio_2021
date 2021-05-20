clear

addpath('../../');
addpath(genpath('../../toolbox/'));

%% set up simulation

% options for saving simulation results
saveplot = 1;
savedata = 1;
savedir  = 'results\';
mkdir(savedir);

% model parameters
param.T     = 1;   % decision threshold
param.sigma = 0.1; % standard deviation of accumulation noise
param.tau   = 0;   % number of timesteps for post-decision evidence accumulation
param.tmax  = 1e4; % max # of timesteps before simulation exits

% simulation setup
sim.S_list  = linspace(param.T/1000, param.T/40, 10); % values of S_i to probe in simulation 
sim.ntrials = 1e5;                                    % number of trials per simulation
sim.nreps   = 10;                                     % number of simulation repetitions

% % speed up simulations with the parallel processing toolbox
% parpool;


%% perform the simulation

tic

for i_rep = 1:sim.nreps
    perf = TI_MPL2016_d_search(param, sim);
    
    d(:, i_rep)     = perf.d';
    presp(:, i_rep) = perf.presp';
end

sim_runtime_in_minutes = toc / 60

% delete(gcp) 

% average across simulation repetitions
d     = mean(d, 2);
presp = mean(presp, 2);


%% regress d' onto PE

d_target = 1.5282;

B_PE2d = polyfit(sim.S_list', d, 2);

a = B_PE2d(1);
b = B_PE2d(2);
c = (B_PE2d(3) - d_target);

PE_fit(1) = (-b + sqrt(b^2 - (4*a*c)) ) / (2*a);
PE_fit(2) = (-b - sqrt(b^2 - (4*a*c)) ) / (2*a);


%% save data

filename = 'TI_MPL2016_d_search';

if savedata
    results  = v2struct(d, presp, B_PE2d, d_target, PE_fit);
    savefile = [savedir filename  '.mat'];
    save(savefile, 'param', 'sim', 'results', 'sim_runtime_in_minutes')
end


%% plot results

fs = 11.5;

figure; hold on;
plot(sim.S_list, d, 'bo-')

plot([PE_fit(1), PE_fit(1)], [0, d_target], 'k-');
plot([0, PE_fit(1)], [d_target, d_target], 'k-');

plot(sim.S_list, a*sim.S_list.^2 + b*sim.S_list + c + d_target, 'k-');

xlabel('PE strength')
ylabel('d''')
set(gca, 'FontSize', fs);

title(['d'' target = ' num2str(d_target) ', PE fit = ' num2str(PE_fit)]);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end
