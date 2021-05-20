clear

param.T     = 1;
param.sigma = 0.1;
param.S_i   = [.01; 0];
param.tau   = 20;

param.save_traces = 1;

[perf, traces] = TI_sim_trial(param);

figure; 
subplot(1,2,1); hold on;
plot(traces.x_it(1,:), 'b-')
plot(traces.x_it(2,:), 'r-')
xl = xlim; 
yl = ylim;
plot(xl, param.T*[1,1], 'k--')
plot(perf.rt*[1,1], yl, 'k--')
plot(perf.rt2*[1,1], yl, 'k--')
xlabel('time')
ylabel('evidence')
title('absolute evidence units x_i')
legend('x_1', 'x_2','location','northwest')

subplot(1,2,2); hold on;
plot(traces.delta_it(1,:), 'b-')
plot(traces.delta_it(2,:), 'r-')
plot(xl, param.T*[1,1], 'k--')
plot(perf.rt*[1,1], yl, 'k--')
plot(perf.rt2*[1,1], yl, 'k--')
xlabel('time')
ylabel('evidence')
title('relative evidence units \delta_i')
legend('\delta_1', '\delta_2','location','northwest')