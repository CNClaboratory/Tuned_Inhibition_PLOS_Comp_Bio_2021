function perf = TI_MPL2016_md_Cx_search(param, sim)

%% initialize simulation

stimID    = [];
S_PE      = [];
S_NE      = [];
cond_tau  = [];
tau_trial = [];
for i_cond_tau = 1:length(sim.tau_list)
    cond_tau = [cond_tau, i_cond_tau*ones(1,sim.ntrials/10)];

    stimID = [stimID, zeros(1,sim.ntrials/20), ones(1,sim.ntrials/20)];
    S_PE   = [S_PE, sim.S*ones(1,sim.ntrials/10)];
    S_NE   = [S_NE, zeros(1,sim.ntrials/10)];
    
    tau_trial = [tau_trial, sim.tau_list(i_cond_tau)*ones(1,sim.ntrials/10)];
end

resp      = -1*ones(1,sim.ntrials);
conf_Cx   = -1*ones(1,sim.ntrials);
responded = -1*ones(1,sim.ntrials);


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
    trial_param.tau = tau_trial(i_trial);
    trial_perf      = TI_sim_trial(trial_param);

    resp(i_trial)       = trial_perf.resp;
    rt(i_trial)         = trial_perf.rt;
    responded(i_trial)  = trial_perf.responded;

    conf_Cx(i_trial)    = trial_perf.conf_Cx;
    responded2(i_trial) = trial_perf.responded;

end


%% analysis

% prating(i) corresponds to mean p(conf=i) across subjects and conditions
% in MPL 2016 data set
prating = [0.3497, 0.3713, 0.2267, 0.0524];

% compute confidence thresholds so as to match p(conf=i) in simulated and real data
U2 = quantile(conf_Cx(resp >= 0), prating(1));
U3 = quantile(conf_Cx(resp >= 0), prating(1)+prating(2));
U4 = quantile(conf_Cx(resp >= 0), prating(1)+prating(2)+prating(3));

for i_cond_tau = 1:length(sim.tau_list)
    
    % select trials
    f_cond = cond_tau == i_cond_tau;
    presp(i_cond_tau) = mean(responded(f_cond));

    f = responded==1 & responded2==1 & f_cond;

    stimIDf   = stimID(f);
    respf     = resp(f);
    rtf       = rt(f);    
    conff     = conf_Cx(f);

    % compute d'
    hr  = (sum(stimIDf==1 & respf==1) + 1/2) / (sum(stimIDf==1) + 1);
    far = (sum(stimIDf==0 & respf==1) + 1/2) / (sum(stimIDf==0) + 1);
    d(i_cond_tau) = norminv(hr) - norminv(far);
    
    rt_median(i_cond_tau)   = median(rtf);    
    
    % compute conf rating
    rating = ones(size(conf_Cx));
    rating(conf_Cx >= U2) = 2;
    rating(conf_Cx >= U3) = 3;
    rating(conf_Cx >= U4) = 4;
    rating(conf_Cx == -1) = -1;
    
    ratingf = rating(f);
    
    % compute meta-d'
    [nR_S1, nR_S2] = trials2counts(stimIDf, respf, ratingf, 4, 1);
    fit = fit_meta_d_MLE(nR_S1, nR_S2);

    md_Cx(i_cond_tau) = fit.meta_da;

end


%% package output

perf = v2struct(stimID, S_PE, S_NE, cond_tau, responded, responded2, d, rt_median, md_Cx, presp);

end