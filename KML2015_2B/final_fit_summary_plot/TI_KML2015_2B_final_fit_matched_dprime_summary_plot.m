clear

%% load target and simulated matched-d' data

TI_KML2015_2B_final_fit_matched_dprime_load_data


%% interpolate conf, meta-d', and RT at fitted sigma values

% given the level of sigma that yields a fit to the d'-matched (high PE - low PE) 
% confidence difference in the data, find the corresponding d'-matched meta-d' 
% and RT values

% sigma values yielding fits to d'-matched conf differences
% manually enter values obtained from 
% TI_KML2015_2B_final_fit_matched_dprime_interpolate_sigma.m
sigma_Cx = 0.10867;
sigma_Cd = 0.16795;

% consecutive indeces for sigma containing the nearest sigma values to the 
% fitted sigma. used for interpolating meta-d' and RT values in the loop below
ind_Cx = [1, 2];
ind_Cd = [7, 8];

% values of sigma at ind_C?(1). 
% used for interpolating meta-d' and RT values in the loop below
b_Cx = .1;
b_Cd = .16;

   
% interpolate Cx conf
m = (diff_conf_Cx(ind_Cx(2)) - diff_conf_Cx(ind_Cx(1))) / .01;
b = diff_conf_Cx(ind_Cx(1)) - m*b_Cx;
conf_Cx_interp = m*sigma_Cx + b;

% interpolate Cx meta-d'
m = (diff_md_Cx(ind_Cx(2)) - diff_md_Cx(ind_Cx(1))) / .01;
b = diff_md_Cx(ind_Cx(1)) - m*b_Cx;
md_Cx_interp = m*sigma_Cx + b;

% interpolate Cx RT
m = (diff_rt_Cx(ind_Cx(2)) - diff_rt_Cx(ind_Cx(1))) / .01;
b = diff_rt_Cx(ind_Cx(1)) - m*b_Cx;
rt_Cx_interp = m*sigma_Cx + b;    

    
% interpolate Cd conf
m = (diff_conf_Cd(ind_Cd(2)) - diff_conf_Cd(ind_Cd(1))) / .01;
b = diff_conf_Cd(ind_Cd(1)) - m*b_Cd;
conf_Cd_interp = m*sigma_Cd + b;

% interpolate Cd meta-d'
m = (diff_md_Cd(ind_Cd(2)) - diff_md_Cd(ind_Cd(1))) / .01;
b = diff_md_Cd(ind_Cd(1)) - m*b_Cd;
md_Cd_interp = m*sigma_Cd + b;

% interpolate Cd RT
m = (diff_rt_Cd(ind_Cd(2)) - diff_rt_Cd(ind_Cd(1))) / .01;
b = diff_rt_Cd(ind_Cd(1)) - m*b_Cd;
rt_Cd_interp = m*sigma_Cd + b;

    
%% compute error for interpolated and actual conf, meta-d', and RT values

% at the value of sigma for which the model fits d'-matched confidence, 
% compute the error for the corresponding conf, meta-d' and RT values. 
% (note, we expect the conf error to be near-zero since that is the objecive 
% of the fitting procedure)

% Cx errror
conf_Cx_error = conf_Cx_interp - diff_conf_target;
md_Cx_error   = md_Cx_interp - diff_md_target;
rt_Cx_error   = rt_Cx_interp - diff_rt_target;

% Cd errror
conf_Cd_error = conf_Cd_interp - diff_conf_target;
md_Cd_error   = md_Cd_interp - diff_md_target;
rt_Cd_error   = rt_Cd_interp - diff_rt_target;



%% plot the results

% plot options
ms = 20;
fs = 15;
lw = 1.5;

xcol = [255 200 0] / 255;
dcol = [100 200 100] / 255;


figure;
   
%%% plot d'-matched confidence
subplot(2,2,1); hold on;

% simulated and target values
plot(sigma, diff_conf_Cx, 'ro-', 'color', xcol, 'linewidth', lw);
plot(sigma, diff_conf_Cd, 'bo-', 'color', dcol, 'linewidth', lw);
plot(xlim, diff_conf_target*[1,1], 'k--', 'linewidth', lw);

% interpolated Cx values
plot(sigma_Cx*[1,1], [diff_conf_target, conf_Cx_interp], 'r--', 'color', xcol, 'linewidth', lw);
plot(sigma_Cx, conf_Cx_interp, 'r.', 'MarkerSize', ms, 'color', xcol)

% interpolated Cd values
plot(sigma_Cd*[1,1], [diff_conf_target, conf_Cd_interp], 'b--', 'color', dcol, 'linewidth', lw);
plot(sigma_Cd, conf_Cd_interp, 'b.', 'MarkerSize', ms, 'color', dcol)

xlim([min(sigma), max(sigma)])

title('conf_{high PE} - conf_{low PE}')
legend('C_x','C_\delta','data','location','northwest')
legend('boxoff')
xlabel('\sigma_{high PE}')

% ylim([0 0.8])
% set(gca,'YTick', 0:.2:.8)
% set(gca,'XTick', .1:.5:.2)
set(gca,'fontsize',fs)


%%% plot d'-matched meta-d'
subplot(2,2,2); hold on;

% simulated and target values
plot(sigma, diff_md_Cx, 'ro-', 'color', xcol, 'linewidth', lw);
plot(sigma, diff_md_Cd, 'bo-', 'color', dcol, 'linewidth', lw);
plot(xlim, diff_md_target*[1,1], 'k--', 'linewidth', lw);

% interpolated Cx values
plot(sigma_Cx*[1,1], [diff_md_target, md_Cx_interp], 'r--', 'color', xcol, 'linewidth', lw);
plot(sigma_Cx, md_Cx_interp, 'r.', 'MarkerSize', ms, 'color', xcol);

% interpolated Cd values
plot(sigma_Cd*[1,1], [diff_md_target, md_Cd_interp], 'b--', 'color', dcol, 'linewidth', lw);
plot(sigma_Cd, md_Cd_interp, 'b.', 'MarkerSize', ms, 'color', dcol);

xlim([min(sigma), max(sigma)])

title('meta-d''_{high PE} - meta-d''_{low PE}')   
xlabel('\sigma_{high PE}')

%     ylim([-0.8 0.3])
%     set(gca,'YTick', -0.8:.2:.2)
% set(gca,'XTick', .1:.5:.2)
set(gca,'fontsize',fs);


%%% plot d'-matched RT'
subplot(2,2,3); hold on;

% simulated and target values
plot(sigma, diff_rt_Cx, 'ro-', 'color', xcol, 'linewidth', lw);
plot(sigma, diff_rt_Cd, 'bo-', 'color', dcol, 'linewidth', lw);
plot(xlim, diff_rt_target*[1,1], 'k--', 'linewidth', lw);

% interpolated Cx values
plot(sigma_Cx*[1,1], [diff_rt_target, rt_Cx_interp], 'r--', 'color', xcol, 'linewidth', lw);
plot(sigma_Cx, rt_Cx_interp, 'r.', 'MarkerSize', ms, 'color', xcol);

% interpolated Cd values
plot(sigma_Cd*[1,1], [diff_rt_target, rt_Cd_interp], 'b--', 'color', dcol, 'linewidth', lw);
plot(sigma_Cd, rt_Cd_interp, 'b.', 'MarkerSize', ms, 'color', dcol);

xlim([min(sigma), max(sigma)])

title('\eta_{RT}')
xlabel('\sigma_{high PE}')

%     ylim([0 2.2])
%     set(gca,'YTick', 0:.5:2)
% set(gca,'XTick', .1:.5:.2)
set(gca,'fontsize',fs);


% plot fitting errors
subplot(2,2,4); hold on;
y = [conf_Cx_error, md_Cx_error, rt_Cx_error; ... 
     conf_Cd_error, md_Cd_error, rt_Cd_error];
b = bar(y, 'FaceColor', 'flat');
colormap gray
for k = 1:size(y,2)
    b(k).CData = k;
end

set(gca, 'XTick', [1 2])
set(gca, 'XTickLabel', {'C_x', 'C_\delta'});
title('fitting error')
legend('confidence', 'meta-d''','\eta_{RT}','location','northwest')
legend('boxoff')

%     ylim([-.5, 2])
%     set(gca,'YTick', -.5:.5:2)
set(gca,'fontsize',fs);
