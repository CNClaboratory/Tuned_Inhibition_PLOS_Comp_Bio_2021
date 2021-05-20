function PE_output = PE_fit_vs_alpha(alpha_input, cond, makePlot)

if ~exist('makePlot','var')
    makePlot = 0;
end

% low PE d'       -> 1.3225, 1.7090, 2.0955
% high PE d'      -> 1.4760, 1.7126, 1.9491
% overall mean d' -> 1.7108
d_target = [1.3225, 1.4760, 1.7090, 1.7108, 1.7126, 1.9491, 2.0955]; 
switch cond
    case 'lowPE hard',   ind = 1;
    case 'highPE hard',  ind = 2;
    case 'lowPE mean',   ind = 3;
    case 'overall mean', ind = 4;
    case 'highPE mean',  ind = 5;
    case 'highPE easy',  ind = 6;
    case 'lowPE easy',   ind = 7;
end

alpha = [0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9];

for i = 1:length(alpha)
    load(['..\d_search\results\TI_KML2015_1A_d_search_alpha=' num2str(alpha(i)) '.mat']);
    PE_fit(i) = results.PE_fit(1,ind);
end


B = polyfit(alpha,PE_fit,6);

if any(alpha == alpha_input)
    PE_output = PE_fit(alpha == alpha_input);
else
    x = alpha_input;
    PE_output = B(1)*x.^6 + B(2)*x.^5 + B(3)*x.^4 + B(4)*x.^3 + B(5)*x.^2 + B(6)*x + B(7);
end

if makePlot
    x = 0:.01:.9;
    % s = 6;
    % m = (PE_fit(end)-PE_fit(1)) / exp(s*NE_mult(end));
    % y_exp = (m*(exp(s*x)-1)) + PE_fit(1);
    y_poly = B(1)*x.^6 + B(2)*x.^5 + B(3)*x.^4 + B(4)*x.^3 + B(5)*x.^2 + B(6)*x + B(7);

    figure; hold on;
    plot(alpha, PE_fit, 'bo-')
    plot(x, y_poly, 'r-')
    % plot(x, y_exp, 'r-')
    xlabel('NE mult')
    plot(alpha_input*[1,1], PE_output*[0,1], 'k-')
    plot(alpha_input*[0,1], PE_output*[1,1], 'k-')
    ylabel(['PE yielding d'' = ' num2str(d_target(ind), 3)])
    title(cond)
end

end