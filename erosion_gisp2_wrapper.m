% This script takes ice sheet model output at the GISP2 and gisp-CR1 bedrock
% core sites and determines the erosion rate that yields the best fit to
% the cosmogenic-nuclide data. Then, the erosion law e = Cu^b is fit using
% the average ice sheet velocity over the model run and the erosion rate
% from the cosmogenic-nuclide data. 

function out = erosion_gisp2_wrapper(history_filename, plot_flag);
%% Set up depth profile

%% Inputs to calculate production 

% Define a rock core depth -- make it longer than your actual core. 

dz = 0.1; % [cm] needs to be same precision as your sample thickness measurements.

rho = 2.65;            % [g cm^-3]; rock density
z_cm = 0:dz:4000000;   % [cm]; length of starting production profile.
                       % needs to be longer than total erosion over model.
z_gcm2 = z_cm.*rho;    % [g cm^-2]; depth


consts = bedrock_constants();
%% Site and sample info for gisp-CR1

datafile = 'data/gisp2_conc.xlsx';

gisp_data = readtable(datafile);

% get data in structure

for a = 1:height(gisp_data)
    gisp.data{a}.id = table2array(gisp_data(a, 'id'));
    gisp.data{a}.nuclide = table2array(gisp_data(a, 'nuclide'));
    gisp.data{a}.td = cell2mat(cellfun(@times, cellfun(@str2num, table2array(gisp_data(a, 'td')), 'UniformOutput', false), {rho}, 'UniformOutput', false));
    gisp.data{a}.bd = cell2mat(cellfun(@times, cellfun(@str2num, table2array(gisp_data(a, 'bd')), 'UniformOutput', false), {rho}, 'UniformOutput', false));
    gisp.data{a}.mq = cell2mat(cellfun(@str2num,table2array(gisp_data(a, 'mq')), 'UniformOutput', false));
    gisp.data{a}.N = table2array(gisp_data(a, 'N'));
    gisp.data{a}.dN = table2array(gisp_data(a, 'dN'));
end

gisp.loc.lat = 72.5796;
gisp.loc.lon = -38.4592;
gisp.loc.elv = 0;

%% Make model time 

% pull in data and run get_thk_times

gisp.filename = history_filename;

gisp.fix_time = 0; % years, where model stops and "reality" takes over 

[gisp.hist.model_times, gisp.hist.model_mask, gisp.plot.switch_times, gisp.time, gisp.u_mean] = get_thk_times(gisp.filename, gisp.fix_time, plot_flag); 

%% build time input for gisp

gisp.hist.times_append = [];

%% make mask_append 

gisp.hist.mask_append = []; % append for ice mask. 1 is ice cover, 0 is no ice

% make -1 if no LIA
gisp.hist.historical_erosion = -1; % cm/yr if relevant

%% build depth profile inputs

% Atmospheric pressure at site
gisp.site_p = ERA40atm(gisp.loc.lat,gisp.loc.lon,gisp.loc.elv); % site air pressure


% Build and load muon profile
% build_muon_profile.m builds a production rate profile defined on a grid
% for efficient integration later after Balco, 2017. 
gisp.m = build_muon_profile_w14c(gisp.site_p,consts,0);

% Define production rate info
gisp.SFsp = stone2000(gisp.loc.lat,gisp.site_p,1); % scaling factor

% Build a data structure with production rate information
gisp.p.P10sp = consts.P10q_St .* gisp.SFsp; % Be-10 spallation production rate at surface
gisp.p.P26sp = gisp.p.P10sp.*consts.R2610q; % Al-26 spallation production rate at surface
gisp.p.P14sp = consts.P14q_St.*gisp.SFsp; % C-14 spallation production rate at surface

% Attenuation
gisp.p.Lsp = 160; % g/cm2.

% Define total production

gisp.p.P10z = PofZ(z_gcm2, gisp.m, gisp.p, 10); % sum of production by spallation and muons; Balco (2017)
gisp.p.P26z = PofZ(z_gcm2, gisp.m, gisp.p, 26);
gisp.p.P14z = PofZ(z_gcm2, gisp.m, gisp.p, 14);
 
%% Find best e and prior exposure

% ee_recent = [0:1e-5:2e-4]; % erosion in cm/yr
% ee_longterm = [0:1e-5:1e-4]; % erosion in cm/yr
% xexp = [2e5:1e4:3e5];
% xbur = [0]; % within 1.1 Myr burial found in Schaefer 2016
% 
% for a = 1:length(ee_recent)
%     for b = 1:length(ee_longterm)
%         for c = 1:length(xexp)
%             for d = 1:length(xbur)
%                 chi2(a,b,c,d) = CoreModel_fiterosion_pleisto_BeAlC([ee_recent(a) ee_recent(a) xexp(c) xbur(d)], gisp.p, gisp.hist, gisp.data, consts, 0);
%         
%             end
%         end
%     end
% end
% 
% [v,loc] = min(chi2(:));
% [ii,jj,kk,ll] = ind2sub(size(chi2),loc);
% 
%% Best fit fminsearch

x(1) = 2e-4; % recent erosion (cm)
x(2) = 0; % long-term erosion
x(3) = 1e2; % exposure prior to modeled history

x0 = [x(1) x(2) x(3)]; % initial guess

lb = [0 0 0];
ub = [Inf 0 Inf];

[optx, fval] = fminsearchbnd(@(x) CoreModel_fiterosion_pleisto_BeAlC(x, gisp.p, gisp.hist, gisp.data, consts, 0), x0, lb,ub);

% [optx, fval] = fminsearch(@(x) CoreModel_fiterosion_pleisto_BeAlC(x, gisp.p, gisp.hist, gisp.data, consts, 0), x0);
%% best 

% opt_ee_recent = ee_recent(ii);
% opt_ee_longterm = ee_longterm(jj);
% opt_xexp = xexp(kk);
% opt_xbur = xbur(ll);

%% plot

% result = CoreModel_fiterosion_pleisto_BeAlC([ee_recent(ii) ee_recent(ii) xexp(kk) xbur(ll)], gisp.p, gisp.hist, gisp.data, consts, 1);

if plot_flag == 1
result = CoreModel_fiterosion_pleisto_BeAlC(optx, gisp.p, gisp.hist, gisp.data, consts, 1);

N_model = result.N_model;

figure
for a = 1:length(N_model)
    if gisp.data{a}.nuclide == 10
        plot([N_model(a) N_model(a)], [gisp.data{a}.td(1) gisp.data{a}.bd(end)], 'r')
        hold on
        plot([gisp.data{a}.N gisp.data{a}.N], [gisp.data{a}.td(1) gisp.data{a}.bd(end)], 'k')    
    elseif gisp.data{a}.nuclide == 26
        plot([N_model(a) N_model(a)], [gisp.data{a}.td(1) gisp.data{a}.bd(end)], 'b')
        hold on
        plot([gisp.data{a}.N gisp.data{a}.N], [gisp.data{a}.td(1) gisp.data{a}.bd(end)], 'g')
    end
end
hold on
plot(result.N_final_10, result.z_gcm2,'r:')
plot(result.N_final_26, result.z_gcm2, 'b:')

    set(gca, 'ydir', 'reverse', 'ylim', [0 500])
    grid on
else
end
%% outputs for erosion fitting

gisp.out.u_mean = gisp.u_mean;
gisp.out.recent_e = optx(1);
gisp.out.longterm_e = optx(2);
gisp.out.pre_exp = optx(3);
gisp.out.fval = fval;

out = gisp.out;
end