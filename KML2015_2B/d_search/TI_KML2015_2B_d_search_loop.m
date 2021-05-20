function TI_KML2015_2B_d_search_loop(sigma)

addpath('../../');
addpath(genpath('../../toolbox/'));

%% set up simulation

% options for saving simulation results
saveplot = 1;
savedata = 1;
savedir  = 'results\';
mkdir(savedir);

% model parameters
param.T     = 1;     % decision threshold
param.sigma = sigma; % standard deviation of accumulation noise
param.tmax  = 1e4;   % max # of timesteps before simulation exits
param.tau   = 0;     % number of timesteps for post-decision evidence accumulation

% simulation setup
sim.ntrials = 1e5; % number of trials per simulation
sim.nreps   = 10;  % number of simulation repetitions

sim.alpha   = 0;

% PE search space is tailored to sigma on the basis of prior preliminary
% simulations
switch param.sigma
    case 0.1
        sim.PE_list = linspace(param.T/1000, 0.015, 10);
    case 0.11
        sim.PE_list = linspace(param.T/1000, 0.015, 10);
    case 0.12
        sim.PE_list = linspace(param.T/1000, 0.02, 10);        
    case 0.13
        sim.PE_list = linspace(param.T/1000, 0.02, 10);        
    case 0.14
        sim.PE_list = linspace(param.T/1000, 0.025, 10);        
    case 0.15
        sim.PE_list = linspace(param.T/1000, 0.03,10);
    case 0.16
        sim.PE_list = linspace(param.T/1000, 0.03,10);
    case 0.17
        sim.PE_list = linspace(param.T/1000, 0.04,10);
    case 0.18
        sim.PE_list = linspace(param.T/1000, 0.04,10);
    case 0.19
        sim.PE_list = linspace(param.T/1000, 0.05,10);
    case 0.2
        sim.PE_list = linspace(param.T/1000, 0.05,10);
    case 0.25
        sim.PE_list = linspace(param.T/1000, 0.08,10);        
    case 0.3
        sim.PE_list = linspace(param.T/1000, 0.12,10);
end


%% perform the simulation

tic

for i_rep = 1:sim.nreps
    perf = TI_KML2015_2B_d_search(param, sim);

    d(:, i_rep)     = perf.d';
    presp(:, i_rep) = perf.presp';
end

sim_runtime_in_minutes = toc / 60
    
% average across simulation repetitions
d     = mean(d, 2);
presp = mean(presp, 2);


%% regress d' onto PE

% low PE  -> 0.6896, 0.9952, 1.3008
% high PE -> 0.8407, 1.1756, 1.5106
% overall mean -> 1.0854
d_target = [0.6896, 0.8407, 0.9952, 1.0854, 1.1756, 1.3008, 1.5106];

B_PE2d = polyfit(sim.PE_list',d,2);

for i_d = 1:length(d_target)
    a = B_PE2d(1);
    b = B_PE2d(2);
    c = (B_PE2d(3) - d_target(i_d));
    cc = B_PE2d(3);
    
    PE_fit(1,i_d) = (-b + sqrt(b^2 - (4*a*c)) ) / (2*a);
    PE_fit(2,i_d) = (-b - sqrt(b^2 - (4*a*c)) ) / (2*a);
end


%% save data

filename = ['TI_KML2015_2B_d_search_sigma=' num2str(sigma)];

if savedata
    results  = v2struct(d, presp, B_PE2d, d_target, PE_fit);
    savefile = [savedir filename  '.mat'];
    save(savefile, 'param', 'sim', 'results', 'sim_runtime_in_minutes')
end


%% plot results

fs = 11.5;

figure; hold on;
plot(sim.PE_list, d, 'bo-')

for i_d = 1:length(d_target)
    plot(PE_fit(1,i_d)*[1,1], d_target(i_d)*[0,1], 'k-');
    plot(PE_fit(1,i_d)*[0,1], d_target(i_d)*[1,1], 'k-');
end

plot(sim.PE_list, a*sim.PE_list.^2 + b*sim.PE_list + cc, 'k-');

xlabel(['PE strength (\sigma = ' num2str(sigma) ')'])
ylabel('d''')
set(gca, 'FontSize', fs);

title({['d'' target = ' num2str(d_target)], ['PE fit = ' num2str(PE_fit(1,:))]});

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end

end