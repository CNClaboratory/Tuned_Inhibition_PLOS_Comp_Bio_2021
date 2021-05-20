clear

savedir  = 'results\';
saveplot = 1;

alpha_lowPE  = 0.1;
alpha_highPE_list = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];

alpha_lowPE  = 0.3;
alpha_highPE_list = [0.4 0.5 0.6 0.7 0.8 0.9];

alpha_lowPE  = 0.5;
alpha_highPE_list = [0.6 0.7 0.8 0.9];

for i = 1:length(alpha_highPE_list)
    
    alpha_highPE = alpha_highPE_list(i);
    
    load(['results/TI_KML2015_1A_final_fit_Cd_alpha_lowPE=' num2str(alpha_lowPE) '_alpha_highPE=' num2str(alpha_highPE) '.mat']);
    v2struct(results);
    
    filename = ['TI_KML2015_1A_final_fit_Cd_alpha_lowPE=' num2str(alpha_lowPE) '_alpha_highPE=' num2str(alpha_highPE)];

%% plot confidence rating vs d'

% interpolate high PE conf at d' = 1.709 (mean d' for low PE condition)
d_for_highPE    = [1.4760, 1.9491];
conf_for_highPE = [2.3111, 2.5250];
m = (conf_for_highPE(2) - conf_for_highPE(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = conf_for_highPE(1) - m*d_for_highPE(1);
conf_interp = m*1.709+b;

% define empirical d' and conf for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
d_target      = [1.3225, 1.7090, 2.0955; 1.4760, 1.7090, 1.9491];
rating_target = [2.0345, 2.1728, 2.3111; 2.3111, conf_interp, 2.5250];

fs = 11.5;
lw = 2;
figure; hold on;

style = {'bv-', 'r^-'};
for i_PE = 1:2
    plot(d(i_PE,:), rating_Cd(i_PE,:), style{i_PE})
end

style = {'b*--', 'r*--'};
for i_PE = 1:2
    plot(d_target(i_PE,:), rating_target(i_PE,:), style{i_PE})
end


xlabel('d''')
ylabel('confidence')
legend('low PE', 'high PE', 'location', 'northwest')
set(gca, 'FontSize', fs);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '_conf.png'];
    saveas(gcf, savefile, 'png')
    delete(gcf)
end



%% plot meta-d' vs d'

% interpolate high PE meta-d' at d' = 1.709 (mean d' for low PE condition)
d_for_highPE  = [1.4760, 1.9491];
md_for_highpe = [0.9641, 1.1873];
m = (md_for_highpe(2) - md_for_highpe(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = md_for_highpe(1) - m*d_for_highPE(1);
md_interp = m*1.709+b;

% define empirical d' and meta-d' for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
d_target = [1.3225, 1.7090, 2.0955; 1.4760, 1.7090, 1.9491];
md_target = [0.9165, 1.1031, 1.2896; 0.9641, md_interp, 1.1873];

fs = 11.5;
lw = 2;
figure; hold on;

style = {'bv-', 'r^-'};
for i_PE = 1:2
    plot(d(i_PE,:), md_Cd(i_PE,:), style{i_PE})
end

style = {'b*--', 'r*--'};
for i_PE = 1:2
    plot(d_target(i_PE,:), md_target(i_PE,:), style{i_PE})
end

plot(xlim, xlim, 'k--')

xlabel('d''')
ylabel('meta-d''')
legend('low PE', 'high PE', 'location', 'northwest')

set(gca, 'FontSize', fs);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '_md_Cx.png'];
    saveas(gcf, savefile, 'png')
    delete(gcf)
end



%% plot RT vs d'

% interpolate high PE conf at d' = 1.709 (mean d' for low PE condition)
d_for_highPE  = [1.4760, 1.9491];
rt_for_highpe = [0.6112, 0.5946];
m = (rt_for_highpe(2) - rt_for_highpe(1)) / (d_for_highPE(2) - d_for_highPE(1));
b = rt_for_highpe(1) - m*d_for_highPE(1);
rt_interp = m*1.709+b;

% define empirical d' and RT for PE (row 1 = low PE, row 2 = high PE) and
% difficulty (col 1 = hard, col 2 = mean of low PE easy/difficult, col 3 = easy)
d_target  = [1.3225, 1.7090, 2.0955; 1.4760, 1.7090, 1.9491];
rt_target = [0.6196, 0.60435, 0.5891; 0.6112, rt_interp, 0.5946];

fs = 11.5;
lw = 2;
figure; 
subplot(1,2,1); hold on;

style = {'bv-', 'r^-'};
for i_PE = 1:2
    plot(d(i_PE,:), rt(i_PE,:), style{i_PE})
end
xlabel('d''')
ylabel('RT')
legend('low PE', 'high PE')
title('simulation')

subplot(1,2,2); hold on;
style = {'b*--', 'r*--'};
for i_PE = 1:2
    plot(d_target(i_PE,:), rt_target(i_PE,:), style{i_PE})
end
xlabel('d''')
ylabel('RT')
title('data')


set(gca, 'FontSize', fs);

set(gcf, 'position', [488  193  900  570])


if saveplot
    savefile = [savedir filename '_rt.png'];
    saveas(gcf, savefile, 'png')
    delete(gcf)
end

end
    