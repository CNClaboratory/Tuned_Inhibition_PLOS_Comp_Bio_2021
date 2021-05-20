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
sim.S1      = 0.010923;                  % derived from d_search fit
sim.S2_list = linspace(0,sim.S1*2.5,10);
sim.ntrials = 1e5;                       % number of trials per simulation
sim.nreps   = 10;                        % number of simulation repetitions

% % speed up simulations with the parallel processing toolbox
% parpool;


%% perform the simulation

tic

for i_rep = 1:sim.nreps
    perf = TI_MPL2016_d_search2(param, sim);

    d(:, i_rep)     = perf.d';
    presp(:, i_rep) = perf.presp';
end

sim_runtime_in_minutes = toc / 60
    
% delete(gcp) 

% average across simulation repetitions
d     = mean(d, 2);
presp = mean(presp, 2);


%% regress d' onto S2 stim strength

d_target = [0.8763, 1.1007, 1.9535, 2.3411];

% % B_S22d = regress(d, param.S2_list');
% % S2_fit = d_target / B_S22d;
% 
B_S22d = polyfit(sim.S2_list',d,2);
% B_S22d1 = polyfit(param.S2_list',d,1);

a = B_S22d(1);
b = B_S22d(2);

for i_d = 1:4
    c = B_S22d(3) - d_target(i_d);
    S2_fit(i_d) = (-b + sqrt(b^2 - (4*a*c)) ) / (2*a);
end

%% save data

filename = 'TI_MPL2016_d_search2';

if savedata
    results = v2struct(d, presp, B_S22d, d_target, S2_fit);
    savefile = [savedir filename  '.mat'];
    save(savefile, 'param', 'sim', 'results', 'sim_runtime_in_minutes')    
end


%% plot results

fs = 11.5;
figure; hold on;
plot(sim.S2_list, d, 'bo-')

for i_d = 1:4
    plot(S2_fit(i_d)*[1,1], d_target(i_d)*[0,1], 'k-');
    plot(S2_fit(i_d)*[0,1], d_target(i_d)*[1,1], 'k-');
end

plot(sim.S2_list, B_S22d(1)*sim.S2_list.^2 + B_S22d(2)*sim.S2_list + B_S22d(3), 'k-');

xlabel('S2 strength')
ylabel('d''')
set(gca, 'FontSize', fs);

title({['d'' target = ' num2str(d_target)], ['S2 fit = ' num2str(S2_fit)]});
% title(['NE multipler = ' num2str(param.NE_mult)])

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end
