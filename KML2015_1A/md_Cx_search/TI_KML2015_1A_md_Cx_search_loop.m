function TI_KML2015_1A_md_Cx_search_loop(alpha)

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
sim.ntrials  = 1e5; % number of trials per simulation
sim.nreps    = 10;  % number of simulation repetitions
sim.tau_list = 10:10:100;
sim.alpha    = alpha;

% for each value of alpha, set the value of S that fits mean d' in the low 
% PE condition; derived from d_search fit
switch sim.alpha
    case 0.1
        sim.S = 0.013779;
    case 0.3
        sim.S = 0.018091;
    case 0.5
        sim.S = 0.026095;
end


%% perform the simulation

tic

for i_rep = 1:sim.nreps
    perf = TI_KML2015_1A_md_Cx_search(param, sim);

    d(:, i_rep)     = perf.d';
    rt(:, i_rep)    = perf.rt_median';
    md_Cx(:, i_rep) = perf.md_Cx';
    presp(:, i_rep) = perf.presp';
end

sim_runtime_in_minutes = toc / 60
    
% delete(gcp) 

% average across simulation repetitions
d     = mean(d, 2);
md_Cx = mean(md_Cx, 2);
rt    = mean(rt, 2);

presp = mean(presp, 2);


%% regress meta-d' onto tau

% mean d' and meta-d' in the low PE condition
d_target  = 1.7090;
md_target = 1.1031;

B_tau2md = polyfit(sim.tau_list', md_Cx, 2);

a = B_tau2md(1);
b = B_tau2md(2);
c = (B_tau2md(3) - md_target);

tau_fit(1) = (-b + sqrt(b^2 - (4*a*c)) ) / (2*a);
tau_fit(2) = (-b - sqrt(b^2 - (4*a*c)) ) / (2*a);


%% save data

filename = ['TI_KML2015_1A_md_Cx_search_alpha=' num2str(sim.alpha)];

if savedata
    results  = v2struct(d, rt, md_Cx, presp, B_tau2md, md_target, tau_fit);
    savefile = [savedir filename  '.mat'];
    save(savefile, 'param', 'sim', 'results', 'sim_runtime_in_minutes')
end


%% plot results

fs = 11.5;
figure; hold on;
plot(sim.tau_list, md_Cx, 'bo-')

plot(tau_fit(1)*[1,1], md_target*[0,1], 'k-');
plot(tau_fit(1)*[0,1], md_target*[1,1], 'k-');

plot(sim.tau_list, B_tau2md(1)*sim.tau_list.^2 + B_tau2md(2)*sim.tau_list + B_tau2md(3), 'k-');
% plot(sim.tau_list, d_target*ones(size(sim.tau_list)), 'k--');
% plot(mean(rt)*[1,1], ylim, 'k--');

xlabel(['tau (\alpha = ' num2str(sim.alpha) ')'])
ylabel('meta-d''')
set(gca, 'FontSize', fs);

title(['meta-d'' target = ' num2str(md_target) ', tau fit = ' num2str(tau_fit(1))]);

set(gcf, 'position', [488  193  900  570])



if saveplot
    savefile = [savedir filename '.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end

end