% This script takes ice sheet model output at the GISP2 and JAK-CR1 bedrock
% core sites and determines the erosion rate that yields the best fit to
% the cosmogenic-nuclide data. Then, the erosion law e = Cu^b is fit using
% the average ice sheet velocity during the relevant glacial eriods over the 
% model run and the erosion rate from the cosmogenic-nuclide data. 

function [b C gisp jak] = erosion_law_calibration_wrapper(jak_hist, gisp_hist, plot_flag)

jak = erosion_jak_wrapper(jak_hist, plot_flag);
gisp = erosion_gisp2_wrapper(gisp_hist, plot_flag);


if isnan(gisp.u_mean.recent)
    gisp.u_mean.recent = gisp.u_mean.tot;
end

%% Solve for C and b

% calculate all in m

b = log((jak.longterm_e.*10)/(gisp.recent_e.*10))/log(jak.u_mean.long_term/gisp.u_mean.recent);

C = (jak.longterm_e.*10)/(jak.u_mean.long_term.^b);

%% plot

if plot_flag == 1

plot_x = linspace(0.001, 3500, 1500);
plot_y = (C.*(plot_x.^b));
plot_b1 = (1e-4.*(plot_x));
plot_b2 = (1e-3.*(plot_x));

figure
hold on

plot(plot_x, plot_y, 'k')
plot(plot_x, plot_b1, 'k:')
plot(plot_x, plot_b2, 'k:')
plot(gisp.u_mean.recent, gisp.recent_e.*10, 'k*')
plot(jak.u_mean.long_term, jak.longterm_e.*10, 'ks')

set(gca, 'yscale', 'log')
grid on
xlabel('Sliding velocity (m yr^{-1})')
ylabel('Erosion rate (mm yr^{-1})')
else 
end

end

%% calibration with fit

% x(1) = 1e-5; % C, or bedrock constant
% x(2) = 2; % b, or exponent
% 
% x0 = [x(1) x(2)];
% 
% if isnan(gisp.u_mean.recent)
%     gisp.u_mean.recent = gisp.u_mean.tot;
% end
% 
% erosion_results = [(jak.longterm_e.*0.01) gisp.recent_e.*0.01]; % modeled erosion results in m/yr
% velocity_results = [jak.u_mean.long_term gisp.u_mean.recent]; % relevent mean velocities
% 
% [optx fval] = fmincon(@(x) erosion_const_calib(x, velocity_results, erosion_results), x0, [], []);
% 
% %%
% 
% function out = erosion_const_calib(x, velocity_results, erosion_results)
%    
%     e_pred = x(1) .* (velocity_results.^x(2));
% 
%     miss = sum((e_pred - erosion_results).^2);
% 
%     out = miss;
% end
% 
