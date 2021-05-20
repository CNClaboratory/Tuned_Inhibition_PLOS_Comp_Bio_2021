

%% get simulation results

sigma = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2];

for i_sigma = 1:length(sigma)

    if i_sigma == 1
        diff_conf_Cx(i_sigma) = 0;
        diff_md_Cx(i_sigma)   = 0;
        diff_rt_Cx(i_sigma)   = 0;
    else
        load(['..\final_fit_Cx\results\TI_KML2015_2B_final_fit_Cx_sigma=' num2str(sigma(i_sigma)) '.mat']);

        diff_conf_Cx(i_sigma) = results.rating_Cx(2,2) - results.rating_Cx(1,2);
        diff_md_Cx(i_sigma)   = results.md_Cx(2,2) - results.md_Cx(1,2);
        diff_rt_Cx(i_sigma)   = (results.rt(2,2) - results.rt(1,2)) / (results.rt(1,3) - results.rt(1,1));
    end

    if i_sigma == 1
        diff_conf_Cd(i_sigma) = 0;
        diff_md_Cd(i_sigma)   = 0;
        diff_rt_Cd(i_sigma)   = 0;
    else
        load(['..\final_fit_Cd\results\TI_KML2015_2B_final_fit_Cd_sigma=' num2str(sigma(i_sigma)) '.mat']);

        diff_conf_Cd(i_sigma) = results.rating_Cd(2,2) - results.rating_Cd(1,2);
        diff_md_Cd(i_sigma)   = results.md_Cd(2,2) - results.md_Cd(1,2);
        diff_rt_Cd(i_sigma)   = (results.rt(2,2) - results.rt(1,2)) / (results.rt(1,3) - results.rt(1,1));
    end

end


%% get targets

% interpolate high PE conf at d' = 0.9952 (mean d' for low PE condition)
d_for_highPE = [0.8407, 1.5106];
c_for_highPE = [2.1769, 2.4311];
m = (c_for_highPE(2) - c_for_highPE(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = c_for_highPE(1) - m*d_for_highPE(1);
c_interp = m*0.9952+b;

% define empirical d' and conf for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
d_target      = [0.6896, 0.9952, 1.3008; 0.8407, 0.9952, 1.5106];
rating_target = [2.0558, 2.1638, 2.2718; 2.1769, c_interp, 2.4311];


% interpolate high PE meta-d' at d' = 0.9952 (mean d' for low PE condition)
md_for_highPE = [0.4990, 1.0779];
m = (md_for_highPE(2) - md_for_highPE(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = md_for_highPE(1) - m*d_for_highPE(1);
md_interp = m*0.9952+b;

% define empirical meta-d' for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
md_target = [0.5035, 0.7619, 1.0203; 0.4990, md_interp, 1.0779];


% interpolate high PE RT at d' = 0.9952 (mean d' for low PE condition)
rt_for_highPE = [0.3834, 0.3761];
m = (rt_for_highPE(2) - rt_for_highPE(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = rt_for_highPE(1) - m*d_for_highPE(1);
rt_interp = m*0.9952+b;

% define empirical RT for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
rt_target = [0.3943, 0.3898, 0.3853; 0.3834, rt_interp, 0.3761];


% compute target values for the (high PE - low PE) difference at the value of d'
% corresponding the mean d' of the low PE condition
diff_conf_target = rating_target(2,2) - rating_target(1,2);
diff_md_target   = md_target(2,2) - md_target(1,2);
diff_rt_target   = (rt_target(2,2) - rt_target(1,2)) / (rt_target(1,3) - rt_target(1,1));
