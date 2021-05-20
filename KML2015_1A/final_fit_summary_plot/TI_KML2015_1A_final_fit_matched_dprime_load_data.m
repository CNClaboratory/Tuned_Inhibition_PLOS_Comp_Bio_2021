

%% get simulation results

alpha_lowPE{1}  = 0.1;
alpha_highPE{1} = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];

alpha_lowPE{2}  = 0.3;
alpha_highPE{2} = [0.3 0.4 0.5 0.6 0.7 0.8 0.9];

alpha_lowPE{3}  = 0.5;
alpha_highPE{3} = [0.5 0.6 0.7 0.8 0.9];

for i_alpha_lowPE = 1:3
    for i_alpha_highPE = 1:length(alpha_highPE{i_alpha_lowPE})

        if i_alpha_highPE == 1
            diff_conf_Cx{i_alpha_lowPE}(i_alpha_highPE) = 0;
            diff_md_Cx{i_alpha_lowPE}(i_alpha_highPE)   = 0;
            diff_rt_Cx{i_alpha_lowPE}(i_alpha_highPE)   = 0;
        else
            load(['..\final_fit_Cx\results\TI_KML2015_1A_final_fit_Cx_alpha_lowPE=' num2str(alpha_lowPE{i_alpha_lowPE}) '_alpha_highPE=' num2str(alpha_highPE{i_alpha_lowPE}(i_alpha_highPE)) '.mat']);

            % compute (high PE - low PE) difference at the value of d'
            % corresponding the mean d' of the low PE condition
            diff_conf_Cx{i_alpha_lowPE}(i_alpha_highPE) = results.rating_Cx(2,2) - results.rating_Cx(1,2);
            diff_md_Cx{i_alpha_lowPE}(i_alpha_highPE)   = results.md_Cx(2,2) - results.md_Cx(1,2);
            diff_rt_Cx{i_alpha_lowPE}(i_alpha_highPE)   = (results.rt(2,2) - results.rt(1,2)) / (results.rt(1,3) - results.rt(1,1));
        end
        
        if i_alpha_highPE == 1
            diff_conf_Cd{i_alpha_lowPE}(i_alpha_highPE) = 0;
            diff_md_Cd{i_alpha_lowPE}(i_alpha_highPE)   = 0;
            diff_rt_Cd{i_alpha_lowPE}(i_alpha_highPE)   = 0;
        else
            load(['..\final_fit_Cd\results\TI_KML2015_1A_final_fit_Cd_alpha_lowPE=' num2str(alpha_lowPE{i_alpha_lowPE}) '_alpha_highPE=' num2str(alpha_highPE{i_alpha_lowPE}(i_alpha_highPE)) '.mat']);

            % compute (high PE - low PE) difference at the value of d'
            % corresponding the mean d' of the low PE condition
            diff_conf_Cd{i_alpha_lowPE}(i_alpha_highPE) = results.rating_Cd(2,2) - results.rating_Cd(1,2);
            diff_md_Cd{i_alpha_lowPE}(i_alpha_highPE)   = results.md_Cd(2,2) - results.md_Cd(1,2);
            diff_rt_Cd{i_alpha_lowPE}(i_alpha_highPE)   = (results.rt(2,2) - results.rt(1,2)) / (results.rt(1,3) - results.rt(1,1));
        end
    
    end
end


%% get targets

% interpolate high PE conf at d' = 1.709 (mean d' for low PE condition)
d_for_highPE = [1.4760, 1.9491];
c_for_highPE = [2.3111, 2.5250];
m = (c_for_highPE(2) - c_for_highPE(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = c_for_highPE(1) - m*d_for_highPE(1);
c_interp = m*1.709+b;

% define empirical d' and conf for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
d_target      = [1.3225, 1.7090, 2.0955; 1.4760, 1.7090, 1.9491];
rating_target = [2.0345, 2.1728, 2.3111; 2.3111, c_interp, 2.5250];


% interpolate high PE meta-d' at d' = 1.709 (mean d' for low PE condition)
md_for_highPE = [0.9641, 1.1873];
m = (md_for_highPE(2) - md_for_highPE(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = md_for_highPE(1) - m*d_for_highPE(1);
md_interp = m*1.709+b;

% define empirical meta-d' for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
md_target = [0.9165, 1.1031, 1.2896; 0.9641, md_interp, 1.1873];


% interpolate high PE RT at d' = 1.709 (mean d' for low PE condition)
rt_for_highPE = [0.6112, 0.5946];
m = (rt_for_highPE(2) - rt_for_highPE(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = rt_for_highPE(1) - m*d_for_highPE(1);
rt_interp = m*1.709+b;

% define empirical RT for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
rt_target = [0.6196, 0.60435, 0.5891; 0.6112, rt_interp, 0.5946];


% compute target values for the (high PE - low PE) difference at the value of d'
% corresponding the mean d' of the low PE condition
diff_conf_target = rating_target(2,2) - rating_target(1,2);
diff_md_target   = md_target(2,2) - md_target(1,2);
diff_rt_target   = (rt_target(2,2) - rt_target(1,2)) / (rt_target(1,3) - rt_target(1,1));
