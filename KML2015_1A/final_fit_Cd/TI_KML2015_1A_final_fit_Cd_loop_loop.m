clear

sim_time = [];

alpha_lowPE  = 0.1;
alpha_highPE = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];

for i = 1:length(alpha_highPE)
    tic
    TI_KML2015_1A_final_fit_Cd_loop(alpha_lowPE, alpha_highPE(i))
    sim_time(end+1) = toc;
end


alpha_lowPE  = 0.3;
alpha_highPE = [0.4 0.5 0.6 0.7 0.8 0.9];

for i = 1:length(alpha_highPE)
    tic
    TI_KML2015_1A_final_fit_Cd_loop(alpha_lowPE, alpha_highPE(i))
    sim_time(end+1) = toc;
end


alpha_lowPE  = 0.5;
alpha_highPE = [0.6 0.7 0.8 0.9];

for i = 1:length(alpha_highPE)
    tic
    TI_KML2015_1A_final_fit_Cd_loop(alpha_lowPE, alpha_highPE(i))
    sim_time(end+1) = toc;
end