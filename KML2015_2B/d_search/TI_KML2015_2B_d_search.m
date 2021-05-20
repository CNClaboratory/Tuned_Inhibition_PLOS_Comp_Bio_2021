function perf = TI_KML2015_2B_d_search(param, sim)

%% initialize simulation

stimID    = [];
S_PE      = [];
S_NE      = [];
cond_stim = [];
for i_cond_stim = 1:length(sim.PE_list)
    cond_stim = [cond_stim, i_cond_stim*ones(1,sim.ntrials/10)];

    stimID = [stimID, zeros(1,sim.ntrials/20), ones(1,sim.ntrials/20)];
    S_PE   = [S_PE, sim.PE_list(i_cond_stim)*ones(1,sim.ntrials/10)];
    S_NE   = [S_NE, sim.alpha*sim.PE_list(i_cond_stim)*ones(1,sim.ntrials/10)];
end

resp      = -1*ones(1,sim.ntrials);
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
    trial_perf      = TI_sim_trial(trial_param);

    resp(i_trial)       = trial_perf.resp;
    responded(i_trial)  = trial_perf.responded;
    
end


%% analysis

for i_cond_stim = 1:length(sim.PE_list)
    
    f_cond = cond_stim==i_cond_stim;
    presp(i_cond_stim) = mean(responded(f_cond));

    f = responded==1 & f_cond;

    stimIDf   = stimID(f);
    respf     = resp(f);

    % compute d'
    hr  = (sum(stimIDf==1 & respf==1) + 1/2) / (sum(stimIDf==1) + 1);
    far = (sum(stimIDf==0 & respf==1) + 1/2) / (sum(stimIDf==0) + 1);
    d(i_cond_stim) = norminv(hr) - norminv(far);

end


%% package output

perf = v2struct(stimID, S_PE, S_NE, cond_stim, responded, presp, d);

end