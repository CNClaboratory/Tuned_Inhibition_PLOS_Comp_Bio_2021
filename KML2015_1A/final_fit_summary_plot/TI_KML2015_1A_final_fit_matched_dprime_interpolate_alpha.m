clear

% fits are stored in the first column of the variables
% alpha_highPE_fit_Cx
% alpha_highPE_fit_Cd

%% load target and simulated matched-d' data

TI_KML2015_1A_final_fit_matched_dprime_load_data


%% interpolate high PE alpha

% for each level of low PE alpha, find the value of high PE alpha that
% yields simulated mean confidence difference equal to the target mean 
% confidence difference

for i_alpha_lowPE = 1:3
    
    % for the Cx model, apply a quadratic fit to the first 3 data points of
    % the conf vs high PE alpha curve since the function is approximately
    % quadratic over this range, and 
    B_alpha2dconf_Cx{i_alpha_lowPE} = polyfit(alpha_highPE{i_alpha_lowPE}(1:3),diff_conf_Cx{i_alpha_lowPE}(1:3),2);

    a = B_alpha2dconf_Cx{i_alpha_lowPE}(1);
    b = B_alpha2dconf_Cx{i_alpha_lowPE}(2);
    c = (B_alpha2dconf_Cx{i_alpha_lowPE}(3) - diff_conf_target);

    alpha_highPE_fit_Cx(i_alpha_lowPE, 1) = (-b + sqrt(b^2 - (4*a*c)) ) / (2*a);
    alpha_highPE_fit_Cx(i_alpha_lowPE, 2) = (-b - sqrt(b^2 - (4*a*c)) ) / (2*a);


    % for the Cd model, apply the quadratic fit over all data points, since
    % the entire curve is approximately quadratic
    B_alpha2dconf_Cd{i_alpha_lowPE} = polyfit(alpha_highPE{i_alpha_lowPE},diff_conf_Cd{i_alpha_lowPE},2);

    a = B_alpha2dconf_Cd{i_alpha_lowPE}(1);
    b = B_alpha2dconf_Cd{i_alpha_lowPE}(2);
    c = (B_alpha2dconf_Cd{i_alpha_lowPE}(3) - diff_conf_target);

    alpha_highPE_fit_Cd(i_alpha_lowPE, 1) = (-b + sqrt(b^2 - (4*a*c)) ) / (2*a);
    alpha_highPE_fit_Cd(i_alpha_lowPE, 2) = (-b - sqrt(b^2 - (4*a*c)) ) / (2*a);
    
    % imaginary results occur for the Cd model when the simulated diff_conf 
    % curve never reaches the target value. in this case, manually replace 
    % the imaginary result with the closest fit available, alpha_highPE = 0.9
    if ~isreal( alpha_highPE_fit_Cd(i_alpha_lowPE, 1) )
        alpha_highPE_fit_Cd(i_alpha_lowPE, 1) = 0.9;
    end
end

%% plot interpolation results

markerSize = 20;
fontSize = 12;

figure;

for i_alpha_lowPE = 1:3

    subplot(3,1,i_alpha_lowPE); hold on;

    % plot simulated data and target value
    plot(alpha_highPE{i_alpha_lowPE}, diff_conf_Cx{i_alpha_lowPE}, 'bo');
    plot(alpha_highPE{i_alpha_lowPE}, diff_conf_Cd{i_alpha_lowPE}, 'ro');
    plot(xlim, diff_conf_target*[1,1], 'k--')
    
    % plot fits
    plot(alpha_highPE_fit_Cx(i_alpha_lowPE,1), diff_conf_target, 'b.', 'MarkerSize', markerSize)
    plot(alpha_highPE_fit_Cd(i_alpha_lowPE,1), diff_conf_target, 'r.', 'MarkerSize', markerSize)

    plot(alpha_highPE_fit_Cx(i_alpha_lowPE,1)*[1 1], [0 diff_conf_target], 'b-')
    plot(alpha_highPE_fit_Cd(i_alpha_lowPE,1)*[1 1], [0 diff_conf_target], 'r-')
    
    
    % plot Cx interpolation curve
    plot(alpha_highPE{i_alpha_lowPE}, B_alpha2dconf_Cx{i_alpha_lowPE}(1) * alpha_highPE{i_alpha_lowPE}.^2 + ...
                                      B_alpha2dconf_Cx{i_alpha_lowPE}(2) * alpha_highPE{i_alpha_lowPE} + ...
                                      B_alpha2dconf_Cx{i_alpha_lowPE}(3), 'b-');


    % plot Cd interpolation curve
    plot(alpha_highPE{i_alpha_lowPE}, B_alpha2dconf_Cd{i_alpha_lowPE}(1) * alpha_highPE{i_alpha_lowPE}.^2 + ...
                                      B_alpha2dconf_Cd{i_alpha_lowPE}(2) * alpha_highPE{i_alpha_lowPE} + ...
                                      B_alpha2dconf_Cd{i_alpha_lowPE}(3), 'r-');

    xlim([min(alpha_highPE{i_alpha_lowPE}), max(alpha_highPE{i_alpha_lowPE})])
    ylabel(['\alpha_{low PE} = ' num2str(alpha_lowPE{i_alpha_lowPE})])
    
    if i_alpha_lowPE == 1
        title('conf_{high PE} - conf_{low PE}')
        legend('C_x simulation','C_\delta simulation','target','C_x fit','C_\delta fit','location','northwest')
    end
    
    if i_alpha_lowPE == 3
        xlabel('\alpha_{high PE}')
    end
    
    set(gca, 'FontSize', fontSize);
end