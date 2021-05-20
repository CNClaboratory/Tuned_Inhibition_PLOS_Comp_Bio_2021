function perf = TI_KML2015_2B_final_fit_Cx(param, sim)

%% initialize variables

stimID    = [];
S_PE      = [];
S_NE      = [];
cond_PE   = [];
cond_diff = [];

for i_cond_diff = 1:3
    cond_diff = [cond_diff, i_cond_diff*ones(1,sim.ntrials/3)];
    cond_PE   = [cond_PE, ones(1,sim.ntrials/6), 2*ones(1,sim.ntrials/6)];

    stimID = [stimID, repmat( [zeros(1,sim.ntrials/12), ones(1,sim.ntrials/12)], 1, 2)];
    S_PE   = [S_PE,   sim.S_lowPE(i_cond_diff)*ones(1,sim.ntrials/6), sim.S_highPE(i_cond_diff)*ones(1,sim.ntrials/6)];
    S_NE   = [S_NE,   sim.S_lowNE(i_cond_diff)*ones(1,sim.ntrials/6), sim.S_highNE(i_cond_diff)*ones(1,sim.ntrials/6)];
end

resp      = -1*ones(1,sim.ntrials);
conf_Cx   = -1*ones(1,sim.ntrials);
rt        = -1*ones(1,sim.ntrials);
responded = -1*ones(1,sim.ntrials);
conf_Cx_i = -1*ones(sim.ntrials,2);


%% run simulation

parfor i_trial = 1:sim.ntrials
% for i_trial = 1:ntrials
    
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
    
    % for the low PE condition, sigma is fixed to 0.1
    if cond_PE(i_trial) == 1
        trial_param.sigma = 0.1;
    end
    
    trial_perf      = TI_sim_trial(trial_param);

    resp(i_trial)       = trial_perf.resp;
    rt(i_trial)         = trial_perf.rt;
    responded(i_trial)  = trial_perf.responded;

    conf_Cx(i_trial)    = trial_perf.conf_Cx;
    conf_Cd(i_trial)    = trial_perf.conf_Cd;
    rt2(i_trial)        = trial_perf.rt2;
    responded2(i_trial) = trial_perf.responded2;

    conf_Cx_i(i_trial, :) = trial_perf.conf_Cx_i;

end


%% analysis

% prating(i) corresponds to mean p(conf=i) across subjects and conditions
% in KML 2015 expt 2B data set
prating = [0.3140, 0.2745, 0.2752, 0.1363];

% compute confidence thresholds so as to match p(conf=i) in simulated and real data
U2 = quantile(conf_Cx(resp >= 0), prating(1));
U3 = quantile(conf_Cx(resp >= 0), prating(1)+prating(2));
U4 = quantile(conf_Cx(resp >= 0), prating(1)+prating(2)+prating(3));

for i_cond_PE = 1:2
    for i_cond_diff = 1:3
        
        % select trials
        f_cond = (cond_PE==i_cond_PE) & (cond_diff==i_cond_diff);
        presp(i_cond_PE, i_cond_diff) = mean(responded(f_cond));

        f = responded==1 & responded2==1 & f_cond;

        stimIDf   = stimID(f);
        respf     = resp(f);
        conff     = conf_Cx(f);
        rtf       = rt(f);
        rtcorrf   = rt(f & stimID==resp);
        rtincorrf = rt(f & stimID~=resp);

        % compute d'
        hr  = (sum(stimIDf==1 & respf==1) + 1/2) / (sum(stimIDf==1) + 1);
        far = (sum(stimIDf==0 & respf==1) + 1/2) / (sum(stimIDf==0) + 1);
        d(i_cond_PE, i_cond_diff) = norminv(hr) - norminv(far);

       % compute accuracy and RT
        accf                                     = respf == stimIDf;
        pcorr(i_cond_PE, i_cond_diff)            = mean(accf);
        rt_mean(i_cond_PE, i_cond_diff)          = mean(rtf);
        rt_median(i_cond_PE, i_cond_diff)        = median(rtf);
        rt_min(i_cond_PE, i_cond_diff)           = min(rtf);
        rt_max(i_cond_PE, i_cond_diff)           = max(rtf);
        rt_median_corr(i_cond_PE, i_cond_diff)   = median(rtcorrf);
        rt_median_incorr(i_cond_PE, i_cond_diff) = median(rtincorrf);

        % compute conf rating
        rating_Cx = ones(size(conf_Cx));
        rating_Cx(conf_Cx >= U2) = 2;
        rating_Cx(conf_Cx >= U3) = 3;
        rating_Cx(conf_Cx >= U4) = 4;
        rating_Cx(conf_Cx == -1) = -1;

        ratingf = rating_Cx(f);

        conf_Cx_mean(i_cond_PE, i_cond_diff)   = mean(conff);
        rating_Cx_mean(i_cond_PE, i_cond_diff) = mean(ratingf);

        % compute meta-d'
        [nR_S1, nR_S2] = trials2counts(stimIDf, respf, ratingf, 4, 1);
        fit = fit_meta_d_MLE(nR_S1, nR_S2);

        md_Cx(i_cond_PE, i_cond_diff) = fit.meta_da;

    end
end


%% package output

perf   = v2struct(stimID, S_PE, S_NE, cond_PE, cond_diff, responded, responded2, resp, rt, ...
                  conf_Cx, conf_Cd, rating_Cx, presp, pcorr, ...
                  rt_mean, rt_median, rt_min, rt_max, rt_median_corr, rt_median_incorr, ...
                  conf_Cx_mean, rating_Cx_mean, ...
                  conf_Cx_i, ...
                  d, md_Cx);

end