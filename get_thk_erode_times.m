function [model_times, model_mask, erode] = get_thk_erode_times(data, fix_time); 

% erode vector will come out in cm. 

if nargin < 2
    fix_time = 0;
end

% if there is a cutoff time, truncate datset

data = data(data(:, 1).*-1 >= fix_time, :);

thickness = data(:, 4);

time = data(:, 1).*-1;

etot = data(:, 2) .* 100; % total erosion per timestep in cm

etot(find(~thickness > 0)) = 0;

%% glaciation mask on full dataset

glaciation_mask = thickness > 0; % this is when the site is glaciated

    
%% find model times and ice-cover mask

switch_diff = diff(glaciation_mask); % if difference is -1, then going from ice cover to no ice cover; if difference is 1, going from no ice cover to cover. 
switch_indx = find(or(switch_diff == 1, switch_diff == -1));
switch_indx_full = [1; switch_indx; length(time)];

switch_val = switch_diff(switch_indx);

switch_times = time(switch_indx_full);

model_times = abs(diff(switch_times));

model_mask = zeros(length(model_times), 1);
model_mask(1) = glaciation_mask(1);

for a = 1:length(model_mask(2:end))
    if switch_val(a) == 1
        model_mask(a+1) = 1;
    elseif switch_val(a) == -1
        model_mask(a+1) = 0;
    end
end

%% erosion vector

erode = zeros(length(model_times), 1);

for a = 1:length(erode)
if a == 1
erode(a) = sum(etot(switch_indx_full(a):switch_indx_full(a+1)));
else 
erode(a) = sum(etot((switch_indx_full(a)+1):switch_indx_full(a+1)));
end
end

%%
%to check if times translated correctly. 
% figure
% hold on
% plot(time(glaciation_mask == 1), 3.*ones(length(time(glaciation_mask == 1)), 1), 'bo')
% plot(time(glaciation_mask == 0), 3.*ones(length(time(glaciation_mask == 0)), 1), 'ro')
% 
% for a = 1:length(model_times)
%     if model_mask(a) == 1
%         plot([switch_times(a) (switch_times(a)-model_times(a))], [5 5], 'b', 'LineWidth', 4)
%     elseif model_mask(a) == 0
%         plot([switch_times(a) (switch_times(a)-model_times(a))], [5 5], 'r', 'LineWidth', 4)
%     end
% end
% 
% ylim([0 8])
end
