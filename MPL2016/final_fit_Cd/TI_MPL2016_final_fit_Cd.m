function perf = TI_MPL2016_final_fit_Cd(param, sim)

%% initialize simulation

stimID  = [];
S_PE    = [];
S_NE    = [];
cond_S2 = [];
for i_cond_S2 = 1:length(sim.S2_list)
    cond_S2 = [cond_S2, i_cond_S2*ones(1,sim.ntrials/5)];

    stimID = [stimID, zeros(1,sim.ntrials/10),         ones(1,sim.ntrials/10)];
    S_PE   = [S_PE,   sim.S1*ones(1,sim.ntrials/10), sim.S2_list(i_cond_S2)*ones(1,sim.ntrials/10)];
    S_NE   = [S_NE,   zeros(1,sim.ntrials/5)];
end

resp      = -1*ones(1,sim.ntrials);
conf_Cd   = -1*ones(1,sim.ntrials);
rt        = -1*ones(1,sim.ntrials);
responded = -1*ones(1,sim.ntrials);
conf_Cd_i = -1*ones(sim.ntrials,2);


%% run simulation

parfor i_trial = 1:sim.ntrials
% for i_trial = 1:sim.ntrials
    
    S_i = zeros(2,1);
    if stimID(i_trial)==0
        S_i(1) = S_PE(i_trial);
        S_i(2) = S_NE(i_trial);
    else
        S_i(1) = S_NE(i_trial);
        S_i(2) = S_PE(i_trial);
    end
    
    trial_param     = param;
    trial_param.S_i = S_i;
    trial_perf      = TI_sim_trial(trial_param);

    resp(i_trial)       = trial_perf.resp;
    rt(i_trial)         = trial_perf.rt;
    responded(i_trial)  = trial_perf.responded;

    conf_Cx(i_trial)    = trial_perf.conf_Cx;
    conf_Cd(i_trial)    = trial_perf.conf_Cd;
    rt2(i_trial)        = trial_perf.rt2;
    responded2(i_trial) = trial_perf.responded2;

    conf_Cd_i(i_trial, :) = trial_perf.conf_Cd_i;

end


%% analysis

% prating(i) corresponds to mean p(conf=i) across subjects and conditions
% in MPL 2016 data set
prating = [0.3497, 0.3713, 0.2267, 0.0524];

% compute confidence thresholds so as to match p(conf=i) in simulated and real data
U2 = quantile(conf_Cd(resp >= 0), prating(1));
U3 = quantile(conf_Cd(resp >= 0), prating(1)+prating(2));
U4 = quantile(conf_Cd(resp >= 0), prating(1)+prating(2)+prating(3));

for i_cond_S2 = 1:length(sim.S2_list)
    
    % select trials
    f_cond = (cond_S2==i_cond_S2) | stimID==0;
    presp(i_cond_S2) = mean(responded(f_cond));

    f = responded==1 & f_cond;

    stimIDf   = stimID(f);
    respf     = resp(f);
    conff     = conf_Cd(f);
    rtf       = rt(f);
    rtcorrf   = rt(f & stimID==resp);
    rtincorrf = rt(f & stimID~=resp);

    % compute d'
    hr  = (sum(stimIDf==1 & respf==1) + 1/2) / (sum(stimIDf==1) + 1);
    far = (sum(stimIDf==0 & respf==1) + 1/2) / (sum(stimIDf==0) + 1);
    d(i_cond_S2) = norminv(hr) - norminv(far);
    
    
    % compute accuracy and RT
    accf                        = respf == stimIDf;
    pcorr(i_cond_S2)            = mean(accf);
    rt_mean(i_cond_S2)          = mean(rtf);
    rt_median(i_cond_S2)        = median(rtf);
    rt_min(i_cond_S2)           = min(rtf);
    rt_max(i_cond_S2)           = max(rtf);
    rt_median_corr(i_cond_S2)   = median(rtcorrf);
    rt_median_incorr(i_cond_S2) = median(rtincorrf);

    rt_rS1_corr(i_cond_S2)   = median(rtf(respf==0 & accf==1));
    rt_rS1_incorr(i_cond_S2) = median(rtf(respf==0 & accf==0));
    rt_rS2_corr(i_cond_S2)   = median(rtf(respf==1 & accf==1));
    rt_rS2_incorr(i_cond_S2) = median(rtf(respf==1 & accf==0));
    
    
    % compute conf rating
    rating_Cd = ones(size(conf_Cd));
    rating_Cd(conf_Cd >= U2) = 2;
    rating_Cd(conf_Cd >= U3) = 3;
    rating_Cd(conf_Cd >= U4) = 4;
    rating_Cd(conf_Cd == -1) = -1;
    
    ratingf = rating_Cd(f);

    conf_Cd_mean(i_cond_S2)   = mean(conff);
    rating_Cd_mean(i_cond_S2) = mean(ratingf);

    rating_Cd_rS1_corr(i_cond_S2)   = mean(ratingf(accf==1 & respf==0));
    rating_Cd_rS1_incorr(i_cond_S2) = mean(ratingf(accf==0 & respf==0));
    rating_Cd_rS2_corr(i_cond_S2)   = mean(ratingf(accf==1 & respf==1));
    rating_Cd_rS2_incorr(i_cond_S2) = mean(ratingf(accf==0 & respf==1));

    % compute meta-d'
    [nR_S1, nR_S2] = trials2counts(stimIDf, respf, ratingf, 4, 1);
    fit = fit_meta_d_MLE(nR_S1, nR_S2);

    md_Cd(i_cond_S2) = fit.meta_da;

    fit = fit_rs_meta_d_MLE(nR_S1, nR_S2);

    md_Cd_rS1(i_cond_S2) = fit.meta_da_rS1;
    md_Cd_rS2(i_cond_S2) = fit.meta_da_rS2;

end


%% package output

perf   = v2struct(stimID, S_PE, S_NE, cond_S2, responded, responded2, resp, rt, ...
                  conf_Cd, conf_Cd, rating_Cd, presp, pcorr, ...
                  rt_mean, rt_median, rt_min, rt_max, rt_median_corr, rt_median_incorr, ...
                  rt_rS1_corr, rt_rS1_incorr, rt_rS2_corr, rt_rS2_incorr, ...
                  conf_Cd_mean, rating_Cd_mean, ...
                  rating_Cd_rS1_corr, rating_Cd_rS1_incorr, rating_Cd_rS2_corr, rating_Cd_rS2_incorr, ... 
                  conf_Cd_i, ...
                  d, md_Cd, md_Cd_rS1, md_Cd_rS2);

end