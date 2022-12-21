function [model_times, model_mask, switch_times, time, u_mean] = get_thk_times(data, fix_time, plot_flag); 

if nargin < 2
    fix_time = 0;
end

data = load(data);

% if there is a cutoff time, truncate datset

data = data(data(:, 1).*-1 >= fix_time, :);

thickness = data(:, 4);

time = data(:, 1).*-1;

u_tot = data(:, 3); % velocity in m/yr 

u_tot(thickness<=1) = 0;


%% glaciation mask on full dataset

glaciation_mask = thickness > 1; % this is when the site is glaciated
    
%% find model times and ice-cover mask

switch_diff = diff(glaciation_mask); % if difference is -1, then going from ice cover to no ice cover; if difference is 1, going from no ice cover to cover. 
switch_indx = find(or(switch_diff == 1, switch_diff == -1))+1; % shift by one so switch occurs at end of timestep

switch_val = switch_diff(switch_indx-1);

switch_times = [time(1); time(switch_indx); time(end)];

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

%% get velocity averages

if u_tot(end) == 0
    print('ERROR: need to update script to accomodate histories ending in exposure')
end

u_tot_recent = u_tot((find(glaciation_mask == 0, 1, 'last')+1):end); % find values for u_tot during last period of ice cover
u_tot_longterm = u_tot(1:(find(glaciation_mask == 0, 1, 'last'))); % find values for u_tot prior to last period of ice cover

u_mean.tot = mean(u_tot(u_tot > 0));
u_mean.recent = mean(u_tot_recent(u_tot_recent > 0));
u_mean.long_term = mean(u_tot_longterm(u_tot_longterm > 0));

%%
%to check if times translated correctly. 

if plot_flag == 1
figure
hold on
plot(time(glaciation_mask == 1), 3.*ones(length(time(glaciation_mask == 1)), 1), 'bo')
plot(time(glaciation_mask == 0), 3.*ones(length(time(glaciation_mask == 0)), 1), 'ro')
plot(time, glaciation_mask)

for a = 1:length(model_times)
    if model_mask(a) == 1
        plot([switch_times(a) (switch_times(a)-model_times(a))], [5 5], 'b', 'LineWidth', 4)
    elseif model_mask(a) == 0
        plot([switch_times(a) (switch_times(a)-model_times(a))], [5 5], 'r', 'LineWidth', 4)
    end
end

ylim([0 8])
else 
end
end
