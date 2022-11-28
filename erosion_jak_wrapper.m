% This script takes ice sheet model output at the GISP2 and gisp-CR1 bedrock
% core sites and determines the erosion rate that yields the best fit to
% the cosmogenic-nuclide data. Then, the erosion law e = Cu^b is fit using
% the average ice sheet velocity over the model run and the erosion rate
% from the cosmogenic-nuclide data. 

function out = erosion_jak_wrapper(history_filename, plot_flag);

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

datafile = 'data/JAK-CR1_conc.xlsx';

gisp_data = readtable(datafile);

% get data in structure

for a = 1:height(gisp_data)
    jak.data{a}.id = table2array(gisp_data(a, 'id'));
    jak.data{a}.nuclide = table2array(gisp_data(a, 'nuclide'));
    jak.data{a}.td = table2array(gisp_data(a, 'td')).*rho;
    jak.data{a}.bd = table2array(gisp_data(a, 'bd')).*rho;
    jak.data{a}.mq = table2array(gisp_data(a, 'mq'));
    jak.data{a}.N = table2array(gisp_data(a, 'N'));
    jak.data{a}.dN = table2array(gisp_data(a, 'dN'));
end

jak.loc.lat = 69.2308;
jak.loc.lon = -49.8089;
jak.loc.elv = 93;

%% Make model time 

% pull in data and run get_thk_times

jak.filename = history_filename;

jak.fix_time = 30000; % years, where model stops and "reality" takes over 

[jak.hist.model_times, jak.hist.model_mask, jak.plot.switch_times, jak.time, jak.u_mean] = get_thk_times(jak.filename, jak.fix_time, plot_flag); 

%% build time input for jak

jak.gl_hist.deglac_t = 7508; % local deglaciation age
jak.gl_hist.bridge = jak.fix_time - jak.gl_hist.deglac_t; % time between model stop and deglac_t
jak.gl_hist.historical_cover = 213; % duration of historical cover
jak.gl_hist.historical_deglac = 10; % recent deglaciation (years before 2018)
jak.gl_hist.earlyholo = jak.gl_hist.deglac_t - jak.gl_hist.historical_cover - jak.gl_hist.historical_deglac;


jak.hist.times_append = [jak.gl_hist.bridge; jak.gl_hist.earlyholo; jak.gl_hist.historical_cover; jak.gl_hist.historical_deglac];

%% make mask_append 

jak.hist.mask_append = [1; 0; 1; 0]; % append for ice mask. 1 is ice cover, 0 is no ice

% make -1 if no LIA
jak.hist.historical_erosion = 0.02; % cm/yr, abrasion rate from Balter-Kennedy et al., 2021 with Lsp = 160

%% build depth profile inputs

% Atmospheric pressure at site
jak.site_p = ERA40atm(jak.loc.lat,jak.loc.lon,jak.loc.elv); % site air pressure

% Build and load muon profile
% build_muon_profile.m builds a production rate profile defined on a grid
% for efficient integration later after Balco, 2017. 
jak.m = build_muon_profile_w14c(jak.site_p,consts,0);

% Define production rate info
jak.SFsp = stone2000(jak.loc.lat,jak.site_p,1); % scaling factor

% Build a data structure with production rate information
jak.p.P10sp = consts.P10q_St .* jak.SFsp; % Be-10 spallation production rate at surface
jak.p.P26sp = jak.p.P10sp.*consts.R2610q; % Al-26 spallation production rate at surface
jak.p.P14sp = consts.P14q_St.*jak.SFsp; % C-14 spallation production rate at surface

% Attenuation
jak.p.Lsp = 160; % g/cm2.

% Define total production

jak.p.P10z = PofZ(z_gcm2, jak.m, jak.p, 10); % sum of production by spallation and muons; Balco (2017)
jak.p.P26z = PofZ(z_gcm2, jak.m, jak.p, 26);
jak.p.P14z = PofZ(z_gcm2, jak.m, jak.p, 14);
 
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
%                 chi2(a,b,c,d) = CoreModel_fiterosion_pleisto_BeAlC([ee_recent(a) ee_recent(a) xexp(c) xbur(d)], jak.p, jak.hist, jak.data, consts, 0);
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

x(1) = 0.02; % recent erosion (cm)
x(2) = 0.02; % long-term erosion (cm)
x(3) = 0; % exposure prior to modeled history

x0 = [x(1) x(2) x(3)]; % initial guess

lb = [0.02 0 0];
ub = [0.02 Inf 0];

[optx, fval] = fminsearchbnd(@(x) CoreModel_fiterosion_pleisto_BeAlC(x, jak.p, jak.hist, jak.data, consts, 0), x0, lb,ub);

% [optx, fval] = fminsearch(@(x) CoreModel_fiterosion_pleisto_BeAlC(x, jak.p, jak.hist, jak.data, consts, 0), x0);
%% best 

% opt_ee_recent = ee_recent(ii);
% opt_ee_longterm = ee_longterm(jj);
% opt_xexp = xexp(kk);
% opt_xbur = xbur(ll);

%% plot
if plot_flag == 1
% result = CoreModel_fiterosion_pleisto_BeAlC([ee_recent(ii) ee_recent(ii) xexp(kk) xbur(ll)], jak.p, jak.hist, jak.data, consts, 1);

result = CoreModel_fiterosion_pleisto_BeAlC(optx, jak.p, jak.hist, jak.data, consts, 1);

N_model = result.N_model;

figure
for a = 1:length(N_model)
    if jak.data{a}.nuclide == 10
        plot([N_model(a) N_model(a)], [jak.data{a}.td(1) jak.data{a}.bd(end)], 'r')
        hold on
        plot([jak.data{a}.N jak.data{a}.N], [jak.data{a}.td(1) jak.data{a}.bd(end)], 'k')    
    elseif jak.data{a}.nuclide == 26
        plot([N_model(a) N_model(a)], [jak.data{a}.td(1) jak.data{a}.bd(end)], 'b')
        hold on
        plot([jak.data{a}.N jak.data{a}.N], [jak.data{a}.td(1) jak.data{a}.bd(end)], 'g')
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

jak.out.u_mean = jak.u_mean;
jak.out.recent_e = optx(1);
jak.out.longterm_e = optx(2);
jak.out.pre_exp = optx(3);
jak.out.fval = fval;

out = jak.out;

end