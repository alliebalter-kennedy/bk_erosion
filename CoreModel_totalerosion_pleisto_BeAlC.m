function [out] = CoreModel_totalerosion_pleisto_BeAlC(lat, lon, elv, time, etot, vel, thk)

data = [time, etot, vel, thk];

[history.model_times, history.model_mask, erosion.erode] = get_thk_erode_times(data);

% this version is for sites that are currently covered by ice. 

%% Load Constants

consts = bedrock_constants();

% should update this so it is a stand-alone script. 

%% Inputs to calculate production 

% Define a rock core depth -- make it longer than your actual core. 

erosion.dz = 0.1; % [cm] needs to be same precision as your sample thickness measurements.
erosion.profile_vec = [0:erosion.dz:500]; % 5 meters - actual core is 4 m. 


dz = erosion.dz; 

rho = 2.65;            % [g cm^-3]; rock density
z_cm = [0:dz:4000000]; % [cm]; length of starting production profile.
                       % needs to be longer than total erosion over model.
z_gcm2 = z_cm.*rho;    % [g cm^-2]; depth


profile_vec = erosion.profile_vec'.*rho; % [g cm^-2];

% Atmospheric pressure at site
site_p = ERA40atm(lat,lon,elv); % site air pressure


% Build and load muon profile
% build_muon_profile.m builds a production rate profile defined on a grid
% for efficient integration later after Balco, 2017. 
m = build_muon_profile_w14c(site_p,consts,0);

% Define production rate info
SFsp = stone2000(lat,site_p,1); % scaling factor

% Build a data structure with production rate information
p.P10sp = consts.P10q_St .* SFsp; % Be-10 spallation production rate at surface
p.P26sp = p.P10sp.*consts.R2610q; % Al-26 spallation production rate at surface
p.P14sp = consts.P14q_St.*SFsp; % C-14 spallation production rate at surface

% Attenuation
p.Lsp = 140; % g/cm2.

% Define total production

P10z = PofZ(z_gcm2, m, p, 10); % sum of production by spallation and muons; Balco (2017)
P26z = PofZ(z_gcm2, m, p, 26);
P14z = PofZ(z_gcm2, m, p, 14);
 
%% Unwrap variables 

model_times = history.model_times'; % pre-deglaciation exposure history as determined by the d18O threshold

%% Model Set Up

% Create gl_int vectors with d18O data

% add Holocene history information
timesteps = [model_times]; % concatenate all exposure history information; in this version, just from d18O threshold

ice_mask = history.model_mask;

starttime = max(fliplr(cumsum(timesteps))); % define model start time

N_final_10 = zeros(length(profile_vec)); % create matrix for final nuclide concentrations
N_final_26 = zeros(length(profile_vec));
N_final_14 = zeros(length(profile_vec));


% define erosion

erode = erosion.erode .* rho; % [g cm^-2]
erode = round(erode, 2); % round to nearest 0.1 g cm^-2.

startdepth = max(fliplr(cumsum(erode))); % find starting depth of modern rock surface

start_index = ceil(((startdepth./rho)./dz)+1); % find index of the starting depth

%% Run Model

% intitialize

topdepth_new = startdepth; 
topindex_new = start_index;
bottomindex_new = ceil(start_index + (length(profile_vec)-1));

topdepth_old = topdepth_new;
topindex_old = topindex_new;

N_old_10 = zeros(length(profile_vec), 1);
N_old_26 = zeros(length(profile_vec), 1);
N_old_14 = zeros(length(profile_vec), 1);

% loop through model time

for a = 1:length(timesteps)

    % expose for interglacials

    if ice_mask(a) == 0
        % calculate cosmogenic nuclide accumulation during exposure
        N_new_10 = N_old_10.*exp(-consts.l10.*timesteps(a)) + P10z(topindex_new:bottomindex_new)'./consts.l10 .* (1-exp(-consts.l10.*timesteps(a)));
        N_new_26 = N_old_26.*exp(-consts.l26.*timesteps(a)) + P26z(topindex_new:bottomindex_new)'./consts.l26 .* (1-exp(-consts.l26.*timesteps(a)));
        N_new_14 = N_old_14.*exp(-consts.l14.*timesteps(a)) + P14z(topindex_new:bottomindex_new)'./consts.l14 .* (1-exp(-consts.l14.*timesteps(a)));
        % no erosion takes place
        topindex_new = topindex_old;
    elseif ice_mask(a) == 1
    % bury and erode for glacials
        % cosmogenic nuclide decay
        N_new_10 = N_old_10.*exp(-consts.l10.*timesteps(a));
        N_new_26 = N_old_26.*exp(-consts.l26.*timesteps(a));
        N_new_14 = N_old_14.*exp(-consts.l14.*timesteps(a));
        
        % shift depth profile up so that next exposure period begins at
        % new production rate
        topdepth_new = round(topdepth_old - erode(a), 2); % depth at the end of the glacial; round to nearest 0.1 cm
        topindex_new = ceil(((topdepth_new./rho)./dz)+1); % find index for top of depth profile. will be used to find production rate during next exposure period. 
        bottomindex_new = ceil(topindex_new + (length(profile_vec)-1)); % find index for bottom of depth profile.
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
depths_final = z_gcm2(topindex_new:bottomindex_new)'; % all final depths should be zero. this passes through test in model run script.

out = [N_final_10(1) N_final_26(1)];
end

