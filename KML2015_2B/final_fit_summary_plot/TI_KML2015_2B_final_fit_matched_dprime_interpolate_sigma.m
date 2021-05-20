clear

% fits are stored in the first column of the variables
% sigma_fit_Cx
% sigma_fit_Cd

%% load target and simulated matched-d' data

TI_KML2015_2B_final_fit_matched_dprime_load_data


%% interpolate sigma
% find the value of sigma that yields simulated mean confidence difference
% equal to the target mean confidence difference

% for the Cx model, apply a quadratic fit to the conf diff vs sigma curve 
% since the function is approximately quadratic
B_sigma2dconf_Cx = polyfit(sigma, diff_conf_Cx, 2);

a = B_sigma2dconf_Cx(1);
b = B_sigma2dconf_Cx(2);
c = (B_sigma2dconf_Cx(3) - diff_conf_target);

sigma_fit_Cx(1) = (-b + sqrt(b^2 - (4*a*c)) ) / (2*a);
sigma_fit_Cx(2) = (-b - sqrt(b^2 - (4*a*c)) ) / (2*a); 


% for the Cd model, apply a quadratic fit to the conf diff vs sigma curve 
% since the function is approximately quadratic
B_sigma2dconf_Cd = polyfit(sigma, diff_conf_Cd, 2);

a = B_sigma2dconf_Cd(1);
b = B_sigma2dconf_Cd(2);
c = (B_sigma2dconf_Cd(3) - diff_conf_target);

sigma_fit_Cd(1) = (-b + sqrt(b^2 - (4*a*c)) ) / (2*a);
sigma_fit_Cd(2) = (-b - sqrt(b^2 - (4*a*c)) ) / (2*a); 
    
    
%% plot interpolation results

markerSize = 20;
fontSize = 12;

figure; hold on;

% plot simulated data and target value
plot(sigma, diff_conf_Cx, 'bo');
plot(sigma, diff_conf_Cd, 'ro');
plot(xlim, diff_conf_target*[1,1], 'k--')

% plot fits
plot(sigma_fit_Cx(1), diff_conf_target, 'b.', 'MarkerSize', markerSize)
plot(sigma_fit_Cd(1), diff_conf_target, 'r.', 'MarkerSize', markerSize)

plot(sigma_fit_Cx(1)*[1 1], [0 diff_conf_target], 'b-')
plot(sigma_fit_Cd(1)*[1 1], [0 diff_conf_target], 'r-')


% plot Cx interpolation curve
plot(sigma, B_sigma2dconf_Cx(1) * sigma.^2 + ...
            B_sigma2dconf_Cx(2) * sigma + ...
            B_sigma2dconf_Cx(3), 'b-');


% plot Cd interpolation curve
plot(sigma, B_sigma2dconf_Cd(1) * sigma.^2 + ...
            B_sigma2dconf_Cd(2) * sigma + ...
            B_sigma2dconf_Cd(3), 'r-');

xlim([min(sigma), max(sigma)])
title('conf_{high PE} - conf_{low PE}')
legend('C_x simulation','C_\delta simulation','target','C_x fit','C_\delta fit','location','northwest')
xlabel('\sigma')

set(gca, 'FontSize', fontSize);
