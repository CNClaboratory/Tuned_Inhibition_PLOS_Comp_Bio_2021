function PE_output = PE_fit_vs_sigma(sigma_input, cond, makePlot)

if ~exist('makePlot','var')
    makePlot = 0;
end

% low PE  -> 0.6896, 0.9952, 1.3008
% high PE -> 0.8407, 1.1756, 1.5106
% overall mean -> 1.0854
d_target = [0.6896, 0.8407, 0.9952, 1.0854, 1.1756, 1.3008, 1.5106];

switch cond
    case 'lowPE hard',   ind = 1;
    case 'highPE hard',  ind = 2;
    case 'lowPE mean',   ind = 3;
    case 'overall mean', ind = 4;
    case 'highPE mean',  ind = 5;
    case 'lowPE easy',   ind = 6;
    case 'highPE easy',  ind = 7;
end

sigma = [0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2];

for i = 1:length(sigma)
    load(['..\d_search\results\TI_KML2015_2B_d_search_sigma=' num2str(sigma(i)) '.mat']);
    PE_fit(i) = results.PE_fit(1,ind);
end


B = polyfit(sigma,PE_fit,4);

if any(sigma == sigma_input)
    PE_output = PE_fit(sigma == sigma_input);
else
    x = sigma_input;
    PE_output = B(1)*x.^4 + B(2)*x.^3 + B(3)*x.^2 + B(4)*x.^1 + B(5)*x.^0;
end

if makePlot
    x = 0.1:.01:.3;
    % s = 6;
    % m = (PE_fit(end)-PE_fit(1)) / exp(s*NE_mult(end));
    % y_exp = (m*(exp(s*x)-1)) + PE_fit(1);
    y_poly = B(1)*x.^4 + B(2)*x.^3 + B(3)*x.^2 + B(4)*x.^1 + B(5)*x.^0;

    figure; hold on;
    plot(sigma, PE_fit, 'bo-')
    plot(x, y_poly, 'r-')
    % plot(x, y_exp, 'r-')
    xlabel('sigma')
    plot(sigma_input*[1,1], PE_output*[0,1], 'k-')
    plot(sigma_input*[0,1], PE_output*[1,1], 'k-')
    ylabel(['PE yielding d'' = ' num2str(d_target(ind), 3)])
    title(cond)
    xlim([0.1,.3])
end

end