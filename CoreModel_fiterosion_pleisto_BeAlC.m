function [out] = CoreModel_fiterosion_pleisto_BeAlC(X, p, hist, data, consts, flag)

%% Unwrap variables 

rho = 2.65;

dz = 0.1; % [cm] needs to be same precision as your sample thickness measurements.
profile_cm = [0:dz:500]; % 5 meters - actual core is 4 m. 
profile_gcm2 = profile_cm .* rho;

P10z = p.P10z;
P26z = p.P26z;
P14z = p.P14z;

erate_recent = X(1) .* rho; % g cm-2/yr; average subglacial erosion rate during most recent period of ice cover
erate_long = X(2) .* rho; % g cm-2/yr; average subglacial erosion rate prior to most recent period of ice cover
xexp = X(3); % extra years of exposure prior to ice model start


%% Finish setting up model times

timesteps = [xexp hist.model_times' hist.times_append'];

ice_mask = [0 hist.model_mask' hist.mask_append'];


%% Set up erosion 
erosion_rate = ones(size(ice_mask)) .* ice_mask .* erate_long; % vector with erosion rate in g cm-2 yr-1

% find index for last period of ice cover

recent_cover = find(ice_mask == 1, 1, 'last');

erosion_rate(recent_cover) = erate_recent;

erode = erosion_rate .* timesteps; 
%% Model Set Up

starttime = sum(timesteps); % define model start time

N_final_10 = zeros(length(profile_gcm2), 1); % create matrix for final nuclide concentrations
N_final_26 = zeros(length(profile_gcm2), 1);
N_final_14 = zeros(length(profile_gcm2), 1);


% define erosion

erode = round(erode, 2); % round to nearest 0.01 g cm.^2
startdepth = sum(erode); % find starting depth of modern rock surface

start_index = ceil(((startdepth./rho)./dz))+1; % find index of the starting depth [turn to cm so can use dz]

%% Run Model

% intitialize

topdepth_new = startdepth; 
topindex_new = start_index;
bottomindex_new = ceil(start_index + (length(profile_gcm2)-1));

topdepth_old = topdepth_new;
topindex_old = topindex_new;

N_old_10 = zeros(length(profile_gcm2), 1);
N_old_26 = zeros(length(profile_gcm2), 1);
N_old_14 = zeros(length(profile_gcm2), 1);

% loop through model time

for a = 1:length(timesteps)

    % expose for even a (interglacials)

    if ice_mask(a) == 0
        % calculate cosmogenic nuclide accumulation during exposure
        N_new_10 = N_old_10.*exp(-consts.l10.*timesteps(a)) + P10z(topindex_new:bottomindex_new)'./consts.l10 .* (1-exp(-consts.l10.*timesteps(a)));
        N_new_26 = N_old_26.*exp(-consts.l26.*timesteps(a)) + P26z(topindex_new:bottomindex_new)'./consts.l26 .* (1-exp(-consts.l26.*timesteps(a)));
        N_new_14 = N_old_14.*exp(-consts.l14.*timesteps(a)) + P14z(topindex_new:bottomindex_new)'./consts.l14 .* (1-exp(-consts.l14.*timesteps(a)));
        % no erosion takes place
        topindex_new = topindex_old;
    elseif ice_mask(a) == 1
    % bury and erode for odd a (glacials)
        % cosmogenic nuclide decay
        N_new_10 = N_old_10.*exp(-consts.l10.*timesteps(a));
        N_new_26 = N_old_26.*exp(-consts.l26.*timesteps(a));
        N_new_14 = N_old_14.*exp(-consts.l14.*timesteps(a));
        
        % shift depth profile up so that next exposure period begins at
        % new production rate
        topdepth_new = round(topdepth_old - erode(a), 2); % depth at the end of the glacial; round to nearest 0.01 g cm-2
        topindex_new = ceil(((topdepth_new./rho)./dz)+1); % find index for top of depth profile. will be used to find production rate during next exposure period. 
        bottomindex_new = ceil(topindex_new + (length(profile_gcm2)-1)); % find index for bottom of depth profile.
    end

    topdepth_old = topdepth_new;
    topindex_old = topindex_new;
    N_old_10 = N_new_10;
    N_old_26 = N_new_26;
    N_old_14 = N_new_14;

end

% save results

N_final_10 = N_new_10; % save final depth profiles at 0.1 cm spacing for each combo of erosion rates.
N_final_26 = N_new_26; 
N_final_14 = N_new_14; 
depths_final = profile_gcm2(topindex_new:bottomindex_new)'; % all final depths should be zero. this passes through test in model run script.

if depths_final(1) ~= 0
    disp('WARNING: depths not equal to zero at end of model run')
end


%% Calculate Concentration at Sample Depths

for a = 1:length(data)
    for i = 1:length(data{a}.td)
        [~, min_td] = min(abs(profile_gcm2 - data{a}.td(i))); % maps sample depths to depth profile in case precision is different
        [~, min_bd] = min(abs(profile_gcm2 - data{a}.bd(i))); % sample depths in g cm^-2
    
        N_temp_10 = N_final_10(min_td:min_bd);
        N_temp_26 = N_final_26(min_td:min_bd);
        N_temp_14 = N_final_14(min_td:min_bd);

        depth_temp = profile_gcm2(min_td:min_bd);

        this_N_10(i) = trapz(depth_temp, N_temp_10)./(data{a}.bd(i)-data{a}.td(i)); % sum modeled nuclide concentrations over sample depth
        this_N_26(i) = trapz(depth_temp, N_temp_26)./(data{a}.bd(i)-data{a}.td(i));
        this_N_14(i) = trapz(depth_temp, N_temp_14)./(data{a}.bd(i)-data{a}.td(i));

        clear N_temp_10 % probably not necessary to do this clearing here.
        clear N_temp_26
        clear N_temp_14
    end
   
    if data{a}.nuclide == 10
        N_model(a) = sum(this_N_10.*data{a}.mq)./sum(data{a}.mq);
    elseif data{a}.nuclide == 26
        N_model(a) = sum(this_N_26.*data{a}.mq)./sum(data{a}.mq);
    elseif data{a}.nuclide == 14
        N_model(a) = sum(this_N_14.*data{a}.mq)./sum(data{a}.mq);
    end

    sample_conc(a) = data{a}.N;
    sample_error(a) = data{a}.dN;

    clear this_N_10
    clear this_N_26
    clear this_N_14

end

result.N_model = N_model;

result.N_final_10 = N_new_10;
result.N_final_26 = N_new_26; 
result.N_final_14 = N_new_14;
result.z_cm = profile_cm;
result.z_gcm2 = profile_gcm2;

%% Calculate misfit

miss = sum(((sample_conc - N_model)./sample_error).^2);

result.chi2 = miss;

if flag == 0
    out = miss;
else out = result;
end
end

