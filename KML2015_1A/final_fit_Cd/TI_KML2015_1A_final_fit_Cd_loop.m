function TI_KML2015_1A_final_fit_Cd_loop(alpha_lowPE, alpha_highPE)

addpath('../../');
addpath('../d_search/');
addpath(genpath('../../toolbox/'));

%% set up simulation

saveplot = 1;
savedata = 1;
savedir  = 'results\';
mkdir(savedir);

% model parameters
param.T     = 1;   % decision threshold
param.sigma = 0.1; % standard deviation of accumulation noise
param.tmax  = 1e4; % max # of timesteps before simulation exits

% values for tau derived from md_Cd_search fit
switch alpha_lowPE
    case 0.1
        param.tau = 18;
    case 0.3
        param.tau = 16;
    case 0.5
        param.tau = 15;
end

% simulation setup
sim.ntrials = 3*1e5; % number of trials per simulation
sim.nreps   = 10;    % number of simulation repetitions

sim.alpha_lowPE  = alpha_lowPE;
sim.alpha_highPE = alpha_highPE;

% use d_search fits to determine PE levels for each condition and alpha
% level, using the helper function PE_fit_vs_alpha
sim.S_lowPE(1) = PE_fit_vs_alpha(alpha_lowPE, 'lowPE hard');
sim.S_lowPE(2) = PE_fit_vs_alpha(alpha_lowPE, 'lowPE mean');
sim.S_lowPE(3) = PE_fit_vs_alpha(alpha_lowPE, 'lowPE easy');

sim.S_highPE(1) = PE_fit_vs_alpha(alpha_highPE, 'highPE hard');
sim.S_highPE(2) = PE_fit_vs_alpha(alpha_highPE, 'lowPE mean'); % this is set to lowPE mean so d' in lowPE and highPE mean conditions are identical
sim.S_highPE(3) = PE_fit_vs_alpha(alpha_highPE, 'highPE easy');

sim.S_lowNE  = alpha_lowPE  * sim.S_lowPE;
sim.S_highNE = alpha_highPE * sim.S_highPE;

% % speed up simulations with the parallel processing toolbox
% parpool;


%% perform the simulation

tic

for i_rep = 1:sim.nreps
    perf = TI_KML2015_1A_final_fit_Cd(param, sim);

    d(:, :, i_rep)         = perf.d;              % d'
    md_Cd(:, :, i_rep)     = perf.md_Cd;          % meta-d'
    rating_Cd(:, :, i_rep) = perf.rating_Cd_mean; % mean conf rating
    
    % RT
    rt(:, :, i_rep)        = perf.rt_median;
    rt_corr(:, :, i_rep)   = perf.rt_median_corr;
    rt_incorr(:, :, i_rep) = perf.rt_median_incorr;
    rtmin(:, :, i_rep)     = perf.rt_min;
    rtmax(:, :, i_rep)     = perf.rt_max;
   
    presp(:, :, i_rep)     = perf.presp;
    
end

sim_runtime_in_minutes = toc / 60
    
% delete(gcp) 

% average across simulation repetitions
d         = squeeze( mean(d, 3) );
md_Cd     = squeeze( mean(md_Cd, 3) );
rating_Cd = squeeze( mean(rating_Cd, 3) );
rt        = squeeze( mean(rt, 3) );
rt_corr   = squeeze( mean(rt_corr, 3) );
rt_incorr = squeeze( mean(rt_incorr, 3) );
rtmin     = squeeze( mean(rtmin, 3) );
rtmax     = squeeze( mean(rtmax, 3) );
presp     = squeeze( mean(presp, 3) );


%% save data

filename = ['TI_KML2015_1A_final_fit_Cd_alpha_lowPE=' num2str(alpha_lowPE) '_alpha_highPE=' num2str(alpha_highPE)];

if savedata
    results = v2struct(d, md_Cd, rating_Cd, ...
                       rt, rt_corr, rt_incorr, rtmin, rtmax, ...
                       presp);    
    savefile = [savedir filename  '.mat'];
    save(savefile, 'param', 'sim', 'results', 'perf', 'sim_runtime_in_minutes')
end


%% plot confidence rating vs d'

% interpolate high PE conf at d' = 1.709 (mean d' for low PE condition)
d_for_highPE    = [1.4760, 1.9491];
conf_for_highPE = [2.3111, 2.5250];
m = (conf_for_highPE(2) - conf_for_highPE(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = conf_for_highPE(1) - m*d_for_highPE(1);
conf_interp = m*1.709+b;

% define empirical d' and conf for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
d_target      = [1.3225, 1.7090, 2.0955; 1.4760, 1.7090, 1.9491];
rating_target = [2.0345, 2.1728, 2.3111; 2.3111, conf_interp, 2.5250];

fs = 11.5;
lw = 2;
figure; hold on;

style = {'bv-', 'r^-'};
for i_PE = 1:2
    plot(d(i_PE,:), rating_Cd(i_PE,:), style{i_PE})
end

style = {'b*--', 'r*--'};
for i_PE = 1:2
    plot(d_target(i_PE,:), rating_target(i_PE,:), style{i_PE})
end


xlabel('d''')
ylabel('confidence')
legend('low PE', 'high PE', 'location', 'northwest')
set(gca, 'FontSize', fs);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '_conf.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end



%% plot meta-d' vs d'

% interpolate high PE meta-d' at d' = 1.709 (mean d' for low PE condition)
d_for_highPE  = [1.4760, 1.9491];
md_for_highpe = [0.9641, 1.1873];
m = (md_for_highpe(2) - md_for_highpe(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = md_for_highpe(1) - m*d_for_highPE(1);
md_interp = m*1.709+b;

% define empirical d' and meta-d' for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
d_target = [1.3225, 1.7090, 2.0955; 1.4760, 1.7090, 1.9491];
md_target = [0.9165, 1.1031, 1.2896; 0.9641, md_interp, 1.1873];

fs = 11.5;
lw = 2;
figure; hold on;

style = {'bv-', 'r^-'};
for i_PE = 1:2
    plot(d(i_PE,:), md_Cd(i_PE,:), style{i_PE})
end

style = {'b*--', 'r*--'};
for i_PE = 1:2
    plot(d_target(i_PE,:), md_target(i_PE,:), style{i_PE})
end

plot(xlim, xlim, 'k--')

xlabel('d''')
ylabel('meta-d''')
legend('low PE', 'high PE', 'location', 'northwest')

set(gca, 'FontSize', fs);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '_md_Cx.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end



%% plot RT vs d'

% interpolate high PE conf at d' = 1.709 (mean d' for low PE condition)
d_for_highPE  = [1.4760, 1.9491];
rt_for_highpe = [0.6112, 0.5946];
m = (rt_for_highpe(2) - rt_for_highpe(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = rt_for_highpe(1) - m*d_for_highPE(1);
rt_interp = m*1.709+b;

% define empirical d' and RT for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
d_target  = [1.3225, 1.7090, 2.0955; 1.4760, 1.7090, 1.9491];
rt_target = [0.6196, 0.60435, 0.5891; 0.6112, rt_interp, 0.5946];

fs = 11.5;
lw = 2;
figure; 
subplot(1,2,1); hold on;

style = {'bv-', 'r^-'};
for i_PE = 1:2
    plot(d(i_PE,:), rt(i_PE,:), style{i_PE})
end
xlabel('d''')
ylabel('RT')
legend('low PE', 'high PE')
title('simulation')

subplot(1,2,2); hold on;
style = {'b*--', 'r*--'};
for i_PE = 1:2
    plot(d_target(i_PE,:), rt_target(i_PE,:), style{i_PE})
end
xlabel('d''')
ylabel('RT')
title('data')


set(gca, 'FontSize', fs);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '_rt.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end