function out = PofZ(zin,m,p,nuclide)

% Returns production rate at any depth z (g/cm2) at SLHL.
% m is structure with precalculated data; p is structure with production
% rate data P10sp and P26sp; Lsp is Lsp
% output is production rate
%
% Exponential approximation for muons below 20,000 g/cm2
% 
% The idea is to allow for reasonably fast time integration of muon
% production during erosion.
%
% 
% Greg Balco -- Berkeley Geochronology Center -- May 2016
% Updated by Allie Balter-Kennedy in Dec 2022 to include C-14


 
% Unwrap matrix input
z = reshape(zin,1,numel(zin));

inbounds = (z <= max(max(m.zz)));
outofbounds = ~inbounds;

if nuclide == 10
    P10mu(inbounds) = interp1(m.zz,m.P10mu,z(inbounds));
    P10mu(outofbounds) = m.P10mu(end).*exp(-(z(outofbounds)-max(m.zz))./m.L10bot);
    out = P10mu + p.P10sp.*exp(-z./p.Lsp);
elseif nuclide ==26
    P26mu(inbounds) = interp1(m.zz,m.P26mu,z(inbounds));
    P26mu(outofbounds) = m.P26mu(end).*exp(-(z(outofbounds)-max(m.zz))./m.L26bot);
    out = P26mu + p.P26sp.*exp(-z./p.Lsp);
elseif nuclide == 14
    P14mu(inbounds) = interp1(m.zz,m.P14mu,z(inbounds));
    P14mu(outofbounds) = m.P14mu(end).*exp(-(z(outofbounds)-max(m.zz))./m.L14bot);
    out = P14mu + p.P14sp.*exp(-z./p.Lsp);
else 
    print('Nuclide not supported')
end;

% Rewrap matrix input
out = reshape(out,size(zin,1),size(zin,2));
