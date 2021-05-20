clear

%% load target and simulated matched-d' data

TI_KML2015_1A_final_fit_matched_dprime_load_data


%% interpolate conf, meta-d', and RT at fitted high PE alpha values

% for each level of low PE alpha, given the level of high PE alpha that
% yields a fit to the d'-matched (high PE - low PE) confidence difference 
% in the data, find the corresponding d'-matched meta-d' and RT values

% highPE alpha values yielding fits to d'-matched conf differences
% manually enter values obtained from 
% TI_KML2015_1A_final_fit_matched_dprime_interpolate_alpha.m
alpha_highPE_Cx = [0.25483, 0.39663, 0.55252];
alpha_highPE_Cd = [0.8928, .9, .9];

% consecutive indeces for alpha_highPE{i_alpha_lowPE} containing the
% nearest high PE alpha values to the fitted high PE alpha. used for
% interpolating meta-d' and RT values in the loop below
ind_Cx = [2, 3; 1, 2; 1, 2];
ind_Cd = [8, 9; 6, 7; 4, 5];

% values of alpha_highPE{i_alpha_lowPE} at ind_C?(i_alpha_lowPE,1). 
% used for interpolating meta-d' and RT values in the loop below
b_Cx = [.2, .3, .5];
b_Cd = [.8, .8, .8];


for i_alpha_lowPE = 1:3
    
    % interpolate Cx conf
    m = (diff_conf_Cx{i_alpha_lowPE}(ind_Cx(i_alpha_lowPE,2)) - diff_conf_Cx{i_alpha_lowPE}(ind_Cx(i_alpha_lowPE,1))) / .1;
    b = diff_conf_Cx{i_alpha_lowPE}(ind_Cx(i_alpha_lowPE,1)) - m*b_Cx(i_alpha_lowPE);
    conf_Cx_interp(i_alpha_lowPE) = m*alpha_highPE_Cx(i_alpha_lowPE) + b;
    
    % interpolate Cx meta-d'
    m = (diff_md_Cx{i_alpha_lowPE}(ind_Cx(i_alpha_lowPE,2)) - diff_md_Cx{i_alpha_lowPE}(ind_Cx(i_alpha_lowPE,1))) / .1;
    b = diff_md_Cx{i_alpha_lowPE}(ind_Cx(i_alpha_lowPE,1)) - m*b_Cx(i_alpha_lowPE);
    md_Cx_interp(i_alpha_lowPE) = m*alpha_highPE_Cx(i_alpha_lowPE) + b;
    
    % interpolate Cx RT
    m = (diff_rt_Cx{i_alpha_lowPE}(ind_Cx(i_alpha_lowPE,2)) - diff_rt_Cx{i_alpha_lowPE}(ind_Cx(i_alpha_lowPE,1))) / .1;
    b = diff_rt_Cx{i_alpha_lowPE}(ind_Cx(i_alpha_lowPE,1)) - m*b_Cx(i_alpha_lowPE);
    rt_Cx_interp(i_alpha_lowPE) = m*alpha_highPE_Cx(i_alpha_lowPE) + b;    

    
    % interpolate Cd conf
    m = (diff_conf_Cd{i_alpha_lowPE}(ind_Cd(i_alpha_lowPE,2)) - diff_conf_Cd{i_alpha_lowPE}(ind_Cd(i_alpha_lowPE,1))) / .1;
    b = diff_conf_Cd{i_alpha_lowPE}(ind_Cd(i_alpha_lowPE,1)) - m*b_Cd(i_alpha_lowPE);
    conf_Cd_interp(i_alpha_lowPE) = m*alpha_highPE_Cd(i_alpha_lowPE) + b;    
    
    % interpolate Cd meta-d'
    m = (diff_md_Cd{i_alpha_lowPE}(ind_Cd(i_alpha_lowPE,2)) - diff_md_Cd{i_alpha_lowPE}(ind_Cd(i_alpha_lowPE,1))) / .1;
    b = diff_md_Cd{i_alpha_lowPE}(ind_Cd(i_alpha_lowPE,1)) - m*b_Cd(i_alpha_lowPE);
    md_Cd_interp(i_alpha_lowPE) = m*alpha_highPE_Cd(i_alpha_lowPE) + b;
    
    % interpolate Cd RT
    m = (diff_rt_Cd{i_alpha_lowPE}(ind_Cd(i_alpha_lowPE,2)) - diff_rt_Cd{i_alpha_lowPE}(ind_Cd(i_alpha_lowPE,1))) / .1;
    b = diff_rt_Cd{i_alpha_lowPE}(ind_Cd(i_alpha_lowPE,1)) - m*b_Cd(i_alpha_lowPE);
    rt_Cd_interp(i_alpha_lowPE) = m*alpha_highPE_Cd(i_alpha_lowPE) + b;        
end

%% compute error for interpolated and actual conf, meta-d', and RT values

% for each level of low PE alpha, at the value of high PE alpha for which
% the model fits d'-matched confidence, compute the error for the corresponding 
% conf, meta-d' and RT values. (note, we expect the conf error to be
% near-zero since that is the objecive of the fitting procedure)

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
i = 0;
for i_alpha_lowPE = 1:3
   
    %%% plot d'-matched confidence
    i = i + 1;
    subplot(3,4,i); hold on;
    
    % simulated and target values
    plot(alpha_highPE{i_alpha_lowPE}, diff_conf_Cx{i_alpha_lowPE}, 'ro-', 'color', xcol, 'linewidth', lw);
    plot(alpha_highPE{i_alpha_lowPE}, diff_conf_Cd{i_alpha_lowPE}, 'bo-', 'color', dcol, 'linewidth', lw);
    plot(xlim, diff_conf_target*[1,1], 'k--', 'linewidth', lw);
        
    % interpolated Cx values
    plot(alpha_highPE_Cx(i_alpha_lowPE)*[1,1], [diff_conf_target, conf_Cx_interp(i_alpha_lowPE)], 'r--', 'color', xcol, 'linewidth', lw);
    plot(alpha_highPE_Cx(i_alpha_lowPE), conf_Cx_interp(i_alpha_lowPE), 'r.', 'MarkerSize', ms, 'color', xcol)
    
    % interpolated Cd values
    plot(alpha_highPE_Cd(i_alpha_lowPE)*[1,1], [diff_conf_target, conf_Cd_interp(i_alpha_lowPE)], 'b--', 'color', dcol, 'linewidth', lw);
    plot(alpha_highPE_Cd(i_alpha_lowPE), conf_Cd_interp(i_alpha_lowPE), 'b.', 'MarkerSize', ms, 'color', dcol)
    
    xlim([min(alpha_highPE{i_alpha_lowPE}), max(alpha_highPE{i_alpha_lowPE})])
    ylabel(['\alpha_{low PE} = ' num2str(alpha_lowPE{i_alpha_lowPE})])

    if i_alpha_lowPE == 1        
        title('conf_{high PE} - conf_{low PE}')
        legend('C_x','C_\delta','data','location','northwest')
        legend('boxoff')
    elseif i_alpha_lowPE == 3
        xlabel('\alpha_{high PE}')
    end
    
    ylim([0 2])
    set(gca,'YTick', 0:.5:2)
    if i_alpha_lowPE == 1
        set(gca,'XTick', .1:.2:.9)
    elseif i_alpha_lowPE == 2
        set(gca,'XTick', .3:.2:.9)
    else
        set(gca,'XTick', .5:.1:.9)
    end
    set(gca,'fontsize',fs)
    
    
    %%% plot d'-matched meta-d'
    i = i + 1;
    subplot(3,4,i); hold on;

    % simulated and target values
    plot(alpha_highPE{i_alpha_lowPE}, diff_md_Cx{i_alpha_lowPE}, 'ro-', 'color', xcol, 'linewidth', lw);
    plot(alpha_highPE{i_alpha_lowPE}, diff_md_Cd{i_alpha_lowPE}, 'bo-', 'color', dcol, 'linewidth', lw);
    plot(xlim, diff_md_target*[1,1], 'k--', 'linewidth', lw);

    % interpolated Cx values
    plot(alpha_highPE_Cx(i_alpha_lowPE)*[1,1], [diff_md_target, md_Cx_interp(i_alpha_lowPE)], 'r--', 'color', xcol, 'linewidth', lw);
    plot(alpha_highPE_Cx(i_alpha_lowPE), md_Cx_interp(i_alpha_lowPE), 'r.', 'MarkerSize', ms, 'color', xcol);
    
    % interpolated Cd values
    plot(alpha_highPE_Cd(i_alpha_lowPE)*[1,1], [diff_md_target, md_Cd_interp(i_alpha_lowPE)], 'b--', 'color', dcol, 'linewidth', lw);
    plot(alpha_highPE_Cd(i_alpha_lowPE), md_Cd_interp(i_alpha_lowPE), 'b.', 'MarkerSize', ms, 'color', dcol);
    
    xlim([min(alpha_highPE{i_alpha_lowPE}), max(alpha_highPE{i_alpha_lowPE})])
    
    if i_alpha_lowPE == 1
        title('meta-d''_{high PE} - meta-d''_{low PE}')   
    elseif i_alpha_lowPE == 3
        xlabel('\alpha_{high PE}')
    end
    
    ylim([-0.8 0.3])
    set(gca,'YTick', -0.8:.2:.2)
    if i_alpha_lowPE == 1
        set(gca,'XTick', .1:.2:.9)
    elseif i_alpha_lowPE == 2
        set(gca,'XTick', .3:.2:.9)
    else
        set(gca,'XTick', .5:.1:.9)
    end
    set(gca,'fontsize',fs);
    
    
    %%% plot d'-matched RT'
    i = i + 1;
    subplot(3,4,i); hold on;

    % simulated and target values
    plot(alpha_highPE{i_alpha_lowPE}, diff_rt_Cx{i_alpha_lowPE}, 'ro-', 'color', xcol, 'linewidth', lw);
    plot(alpha_highPE{i_alpha_lowPE}, diff_rt_Cd{i_alpha_lowPE}, 'bo-', 'color', dcol, 'linewidth', lw);
    plot(xlim, diff_rt_target*[1,1], 'k--', 'linewidth', lw);
    
    % interpolated Cx values
    plot(alpha_highPE_Cx(i_alpha_lowPE)*[1,1], [diff_rt_target, rt_Cx_interp(i_alpha_lowPE)], 'r--', 'color', xcol, 'linewidth', lw);
    plot(alpha_highPE_Cx(i_alpha_lowPE), rt_Cx_interp(i_alpha_lowPE), 'r.', 'MarkerSize', ms, 'color', xcol);
    
    % interpolated Cd values
    plot(alpha_highPE_Cd(i_alpha_lowPE)*[1,1], [diff_rt_target, rt_Cd_interp(i_alpha_lowPE)], 'b--', 'color', dcol, 'linewidth', lw);
    plot(alpha_highPE_Cd(i_alpha_lowPE), rt_Cd_interp(i_alpha_lowPE), 'b.', 'MarkerSize', ms, 'color', dcol);
    
    xlim([min(alpha_highPE{i_alpha_lowPE}), max(alpha_highPE{i_alpha_lowPE})])
    
    if i_alpha_lowPE == 1
        title('relative PE : difficulty RT effect')
        title('\eta_{RT}')
    elseif i_alpha_lowPE == 3
        xlabel('\alpha_{high PE}')
    end
    
    ylim([0 2.2])
    set(gca,'YTick', 0:.5:2)
    if i_alpha_lowPE == 1
        set(gca,'XTick', .1:.2:.9)
    elseif i_alpha_lowPE == 2
        set(gca,'XTick', .3:.2:.9)
    else
        set(gca,'XTick', .5:.1:.9)
    end
    set(gca,'fontsize',fs);
    
    
    % plot fitting errors
    i = i + 1;
    subplot(3,4,i); hold on;
    y = [conf_Cx_error(i_alpha_lowPE), md_Cx_error(i_alpha_lowPE), rt_Cx_error(i_alpha_lowPE); ... 
         conf_Cd_error(i_alpha_lowPE), md_Cd_error(i_alpha_lowPE), rt_Cd_error(i_alpha_lowPE)];
    b = bar(y, 'FaceColor', 'flat');
    colormap gray
    for k = 1:size(y,2)
        b(k).CData = k;
    end
    
    set(gca, 'XTick', [1 2])
    set(gca, 'XTickLabel', {'C_x', 'C_\delta'});
    if i_alpha_lowPE == 1
        title('fitting error')
        legend('confidence', 'meta-d''','\eta_{RT}','location','northwest')
        legend('boxoff')
    end
    ylim([-.5, 2])
    set(gca,'YTick', -.5:.5:2)
    set(gca,'fontsize',fs);
    
end
