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
param.tmax  = 1e4; % max # of timesteps before simulation exits

% simulation setup
sim.S        = 0.010923;   % derived from d_search fit
sim.tau_list = 5:5:50;
sim.ntrials  = 1e5;        % number of trials per simulation
sim.nreps    = 10;         % number of simulation repetitions

% % speed up simulations with the parallel processing toolbox
% parpool;


%% perform the simulation

tic

for i_rep = 1:sim.nreps
    perf = TI_MPL2016_md_Cd_search(param, sim);

    d(:, i_rep)     = perf.d';
    rt(:, i_rep)    = perf.rt_median';
    md_Cd(:, i_rep) = perf.md_Cd';
    presp(:, i_rep) = perf.presp';
end

sim_runtime_in_minutes = toc / 60

% delete(gcp) 

% average across simulation repetitions
d     = mean(d, 2);
md_Cd = mean(md_Cd, 2);
rt    = mean(rt, 2);
presp = mean(presp, 2);


%% regress meta-d' onto tau

d_target  = 1.5282;
md_target = 1.1421;

B_tau2md = polyfit(sim.tau_list', md_Cd, 2);

a = B_tau2md(1);
b = B_tau2md(2);
c = (B_tau2md(3) - md_target);

tau_fit(1) = (-b + sqrt(b^2 - (4*a*c)) ) / (2*a);
tau_fit(2) = (-b - sqrt(b^2 - (4*a*c)) ) / (2*a);


%% save data

filename = 'TI_MPL2016_md_Cd_search';

if savedata
    results  = v2struct(d, md_Cd, presp, B_tau2md, md_target, tau_fit);
    savefile = [savedir filename  '.mat'];
    save(savefile, 'param', 'sim', 'results', 'sim_runtime_in_minutes')
end


%% plot results

fs = 11.5;
figure; hold on;
plot(sim.tau_list, md_Cd, 'bo-')

plot(tau_fit(1)*[1,1], md_target*[0,1], 'k-');
plot(tau_fit(1)*[0,1], md_target*[1,1], 'k-');

plot(sim.tau_list, B_tau2md(1)*sim.tau_list.^2 + B_tau2md(2)*sim.tau_list + B_tau2md(3), 'k-');
% plot(sim.tau_list, d_target*ones(size(sim.tau_list)), 'k--');
% plot(mean(rt)*[1,1], ylim, 'k--');

xlabel('tau')
ylabel('meta-d''')
set(gca, 'FontSize', fs);

title(['meta-d'' target = ' num2str(md_target) ', tau fit = ' num2str(tau_fit(1))]);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end
