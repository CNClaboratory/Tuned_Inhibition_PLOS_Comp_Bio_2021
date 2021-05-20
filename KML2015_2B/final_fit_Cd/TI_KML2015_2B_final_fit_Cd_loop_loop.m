clear

sim_time = [];

sigma = [0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2];

for i = 1:length(sigma)
    tic
    TI_KML2015_2B_final_fit_Cd_loop(sigma(i))
    sim_time(end+1) = toc;
end