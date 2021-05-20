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
param.tau   = 59;  % derived from md_Cx_search fit

% simulation setup
sim.S1      = 0.010923;                                             % derived from d_search fit
sim.S2_list = [0.0015652, 0.0046458, 0.010923, 0.017813, 0.024899]; % derived from d_search2 fit
sim.ntrials = 3*1e5;                                                % number of trials per simulation
sim.nreps   = 10;                                                   % number of simulation repetitions

% % speed up simulations with the parallel processing toolbox
% parpool;


%% perform the simulation

tic

for i_rep = 1:sim.nreps
    perf = TI_MPL2016_final_fit_Cx(param, sim);

    % d'
    d(:, i_rep)        = perf.d';
    
    % meta-d'
    md_Cx(:, i_rep)     = perf.md_Cx';
    md_Cx_rS1(:, i_rep) = perf.md_Cx_rS1';
    md_Cx_rS2(:, i_rep) = perf.md_Cx_rS2';

    % mean conf rating
    rating_Cx(:, i_rep)            = perf.rating_Cx_mean';
    rating_Cx_rS1_corr(:, i_rep)   = perf.rating_Cx_rS1_corr';
    rating_Cx_rS1_incorr(:, i_rep) = perf.rating_Cx_rS1_incorr';
    rating_Cx_rS2_corr(:, i_rep)   = perf.rating_Cx_rS2_corr';
    rating_Cx_rS2_incorr(:, i_rep) = perf.rating_Cx_rS2_incorr';
    
    % RT
    rt(:, i_rep)        = perf.rt_median';
    rt_corr(:, i_rep)   = perf.rt_median_corr';
    rt_incorr(:, i_rep) = perf.rt_median_incorr';
    rtmin(:, i_rep)     = perf.rt_min';
    rtmax(:, i_rep)     = perf.rt_max';

    rt_rS1_corr(:, i_rep)   = perf.rt_rS1_corr';
    rt_rS1_incorr(:, i_rep) = perf.rt_rS1_incorr';
    rt_rS2_corr(:, i_rep)   = perf.rt_rS2_corr';
    rt_rS2_incorr(:, i_rep) = perf.rt_rS2_incorr';
    
    presp(:, i_rep)     = perf.presp';

end

sim_runtime_in_minutes = toc / 60
    
% delete(gcp) 

% average across simulation repetitions
d         = mean(d, 2);

md_Cx     = mean(md_Cx, 2);
md_Cx_rS1 = mean(md_Cx_rS1, 2);
md_Cx_rS2 = mean(md_Cx_rS2, 2);

rt       = mean(rt, 2);
rtmin    = mean(rtmin, 2);
rtmax    = mean(rtmax, 2);

rating_Cx            = mean(rating_Cx, 2);
rating_Cx_rS1_corr   = mean(rating_Cx_rS1_corr, 2);
rating_Cx_rS1_incorr = mean(rating_Cx_rS1_incorr, 2);
rating_Cx_rS2_corr   = mean(rating_Cx_rS2_corr, 2);
rating_Cx_rS2_incorr = mean(rating_Cx_rS2_incorr, 2);

rt_rS1_corr   = mean(rt_rS1_corr, 2);
rt_rS1_incorr = mean(rt_rS1_incorr, 2);
rt_rS2_corr   = mean(rt_rS2_corr, 2);
rt_rS2_incorr = mean(rt_rS2_incorr, 2);

presp         = mean(presp, 2);


%% save data

filename = 'TI_MPL2016_final_fit_Cx';

if savedata
    results = v2struct(d, md_Cx, md_Cx_rS1, md_Cx_rS2, rt, rtmin, rtmax, ...
                       rt_rS1_corr, rt_rS1_incorr, rt_rS2_corr, rt_rS2_incorr, ...
                       rating_Cx_rS1_corr, rating_Cx_rS1_incorr, rating_Cx_rS2_corr, rating_Cx_rS2_incorr, ...
                       rating_Cx, presp);    
    savefile = [savedir filename  '.mat'];
    save(savefile, 'param', 'sim', 'results', 'perf', 'sim_runtime_in_minutes')
end


%% plot response-specific meta-d'

d_target      = [0.8763    1.1007    1.5282    1.9535    2.3411];
md_target     = [0.6689    0.7588    1.1421    1.4588    1.8415];
md_rS1_target = [1.0339    1.0929    1.1319    0.7670    0.4621];
md_rS2_target = [0.3226    0.5407    1.1804    1.5802    1.9274];

fs = 11.5;
lw = 2;

figure; hold on;
plot(d, md_Cx, 'ko-', 'LineWidth',lw)
plot(d, md_Cx_rS1, 'r^-', 'LineWidth',lw)
plot(d, md_Cx_rS2, 'bv-', 'LineWidth',lw)

plot(d_target, md_target, 'ko--')
plot(d_target, md_rS1_target, 'r^--')
plot(d_target, md_rS2_target, 'bv--')

plot(d, d, 'k--', 'LineWidth',lw)

xlabel('d''')
ylabel('meta-d''')
set(gca, 'FontSize', fs);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end


%% response-specific RT

rt_rS1_corr_target   = [0.5475    0.5475    0.5475    0.5475    0.5475];
rt_rS1_incorr_target = [0.5794    0.5352    0.5274    0.5349    0.5335];
rt_rS2_corr_target   = [0.5334    0.5288    0.5088    0.4829    0.4765];
rt_rS2_incorr_target = [0.5193    0.5193    0.5193    0.5193    0.5193];

fs = 11.5;

figure; 
subplot(2,1,1); hold on;
plot(1:5, rt_rS1_corr, 'ro-')
plot(1:5, rt_rS1_incorr, 'rx-')
xlabel('S2 stimulus strength')
ylabel('RT')
title('response = "S1"')

subplot(2,1,2); hold on;
plot(1:5, rt_rS2_corr, 'bo-')
plot(1:5, rt_rS2_incorr, 'bx-')
xlabel('S2 stimulus strength')
ylabel('RT')
title('response = "S2"')

set(gca, 'FontSize', fs);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '_rt.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end

%% response-specific rating

rating_rS1_corr_target   = [2.0171    2.0171    2.0171    2.0171    2.0171];
rating_rS1_incorr_target = [1.6041    1.5853    1.5791    1.7086    1.7804];
rating_rS2_corr_target   = [1.7237    1.8420    2.1363    2.3546    2.5412];
rating_rS2_incorr_target = [1.5788    1.5788    1.5788    1.5788    1.5788];


fs = 11.5;
lw = 2;

figure; 
subplot(2,1,1); hold on;
plot(1:5, rating_Cx_rS1_corr, 'ro-', 'LineWidth', lw)
plot(1:5, rating_Cx_rS1_incorr, 'rx-', 'LineWidth', lw)
plot(1:5, rating_rS1_corr_target, 'ro--', 'LineWidth', 1)
plot(1:5, rating_rS1_incorr_target, 'rx--', 'LineWidth', 1)
xlabel('S2 stimulus strength')
ylabel('conf')
title('response = "S1"')

legend('correct (sim)', 'incorrect (sim)', 'correct (data)', 'incorrect (data)')

subplot(2,1,2); hold on;
plot(1:5, rating_Cx_rS2_corr, 'bo-', 'LineWidth', lw)
plot(1:5, rating_Cx_rS2_incorr, 'bx-', 'LineWidth', lw)
plot(1:5, rating_rS2_corr_target, 'bo--', 'LineWidth', 1)
plot(1:5, rating_rS2_incorr_target, 'bx--', 'LineWidth', 1)
xlabel('S2 stimulus strength')
ylabel('conf')
title('response = "S2"')

set(gca, 'FontSize', fs);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '_rating.png'];
    saveas(gcf, savefile, 'png')
    % delete(gcf)
end
