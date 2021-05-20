clear

% sigma = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2];

sigma = 0.1;

for i = 1:length(sigma)
    TI_KML2015_2B_md_Cx_search_loop(sigma(i))
end