% Headless driver for the CO2-brine PVT example (paper Fig. 2).
% Based on mrst-adblackoil-gcs/co2-brine_props/example_props.m.
% Changes for a headless Octave run:
%  - pvtBrineWithCO2BlackOil needs the MATLAB 'table' class, which Octave
%    lacks. The repository ships that function's output as text tables
%    (PVDG_*.txt, PVTO_*.txt), so we read those instead of regenerating.
%  - subplot replaces tiledlayout/nexttile (Octave has no tiledlayout).
set(0, 'defaultfigurevisible', 'off');
outdir = '/tmp/claude-1000/-home-adriano-codes-MRST/ad017599-1c83-453b-a05e-6b33ad5bc00f/scratchpad';
run('/home/adriano/codes/MRST/startup.m');
mrstModule add ad-blackoil deckformat ad-core ad-props
% Octave replacement of fastInterpTable: ad-core's version needs
% griddedInterpolant, which Octave lacks. Put the Octave version first.
addpath('/home/adriano/codes/MRST/core/utils/octave_only/mrst');

ddir = ['/home/adriano/codes/MRST/mrst-adblackoil-gcs/co2-brine_props/' ...
        'input_DATA_files/'];

%% Load .DATA input file and generate MRST fluid object
fn   = fullfile(ddir, 'example_co2brine_3reg_unsatVals.DATA');
deck = convertDeckUnits(readEclipseDeck(fn));
fluid = initDeckADIFluid(deck);

%% Values computed from the fluid object (dashed lines in paper Fig. 2)
np        = 100;
p_val     = linspace(1,600,np)'*barsa;
rho_co2   = fluid.rhoGS*fluid.bG(p_val);     % dry gas
mu_co2    = fluid.muG(p_val);
rss_val   = fluid.rsSat(p_val);              % live oil (brine with CO2)
rho_b_sat = fluid.bO(p_val,rss_val,true(np,1)) .* ...
            (rss_val.*fluid.rhoGS + fluid.rhoOS);
rho_b     = fluid.rhoOS*fluid.bO(p_val,zeros(np,1),false(np,1));
mu_b_sat  = fluid.muO(p_val,rss_val,true(np,1));
mu_b      = fluid.muO(p_val,zeros(np,1),false(np,1));

%% Tabulated output of pvtBrineWithCO2BlackOil (solid lines in Fig. 2)
% PVDG table: P [bar], B_g [m3/Sm3], mu_g [cP]
tg = dlmread(fullfile(ddir, 'PVDG_T353.15_P1P600_Sbrine100000ppm.txt'), '\t', 1, 0);
% PVTO table: Rs [Sm3/Sm3], P [bar], B_aq [m3/Sm3], mu_aq [cP]
ta = dlmread(fullfile(ddir, 'PVTO_T353.15_P1P600_Sbrine100000ppm.txt'), '\t', 1, 0);

rho_g_txt = fluid.rhoGS ./ tg(:,2);          % rho_g = rhoGS / B_g

% Saturated rows: each new Rs value marks a bubble point
isat = [true; diff(ta(:,1)) ~= 0];
Rs_s  = ta(isat,1);  P_s = ta(isat,2);  B_s = ta(isat,3);  mu_s = ta(isat,4);
rho_aq_sat_txt = (fluid.rhoOS + Rs_s*fluid.rhoGS) ./ B_s;

% Minimal-CO2 branch (first Rs block, undersaturated rows): approximates
% the "no CO2" function curve of Fig. 2
i1 = find(ta(:,1) == ta(1,1));
P_u = ta(i1,2);  B_u = ta(i1,3);  mu_u = ta(i1,4);
rho_b_txt = (fluid.rhoOS + ta(1,1)*fluid.rhoGS) ./ B_u;

%% Quantitative comparison at matching pressures
rho_g_fluid_at  = fluid.rhoGS*fluid.bG(tg(:,2:end)*0 + tg(:,1)*barsa);  % same P
rho_g_fluid_at  = rho_g_fluid_at(:,1);
mu_g_fluid_at   = fluid.muG(tg(:,1)*barsa);
ed_rho_g = max(abs(rho_g_fluid_at - rho_g_txt) ./ rho_g_txt);
ed_mu_g  = max(abs(mu_g_fluid_at*1e3 - tg(:,3)) ./ tg(:,3));

rss_at    = fluid.rsSat(P_s*barsa);
rho_aq_at = fluid.bO(P_s*barsa, rss_at, true(numel(P_s),1)) .* ...
            (rss_at*fluid.rhoGS + fluid.rhoOS);
mu_aq_at  = fluid.muO(P_s*barsa, rss_at, true(numel(P_s),1));
ed_rho_aq = max(abs(rho_aq_at - rho_aq_sat_txt) ./ rho_aq_sat_txt);
ed_mu_aq  = max(abs(mu_aq_at*1e3 - mu_s) ./ mu_s);

fprintf('\nSurface densities from deck: brine %.2f kg/m3, CO2 %.4f kg/m3\n', ...
        fluid.rhoOS, fluid.rhoGS);
fprintf('Max rel. difference fluid object vs table:\n');
fprintf('  CO2 density:            %.3e\n', ed_rho_g);
fprintf('  CO2 viscosity:          %.3e\n', ed_mu_g);
fprintf('  Sat. brine density:     %.3e\n', ed_rho_aq);
fprintf('  Sat. brine viscosity:   %.3e\n', ed_mu_aq);
fprintf('CO2 density at 100/300/600 bar [kg/m3]: %.1f / %.1f / %.1f\n', ...
        interp1(p_val/barsa, rho_co2, [100 300 600]));
fprintf('Sat. brine density at 100/300/600 bar [kg/m3]: %.1f / %.1f / %.1f\n', ...
        interp1(p_val/barsa, rho_b_sat, [100 300 600]));

%% Figure (paper Fig. 2 layout, subplot version)
figure('Position', [100 100 1300 320]);
subplot(1,4,1)
plot(p_val/barsa, rho_co2, '-', 'color', [80 0 0]/255, 'linewidth', 3.5);
hold on
plot(tg(:,1), rho_g_txt, '-', 'color', [220 0 0]/255, 'linewidth', 1.25);
hold off, grid on
xlabel('p [bar]'), ylabel('rho_g [kg/m^3]')
xlim([0 600]), ylim([0 1200])
title('CO2, T=80 C')
legend('fluid object', 'pvtBrine fcn', 'location', 'northwest')

subplot(1,4,2)
semilogy(p_val/barsa, mu_co2*1e3, '-', 'color', [80 0 0]/255, 'linewidth', 3.5);
hold on
semilogy(tg(:,1), tg(:,3), '-', 'color', [220 0 0]/255, 'linewidth', 1.25);
hold off, grid on
xlabel('p [bar]'), ylabel('mu_g [cP]')
xlim([0 600]), ylim([0.01 1])

subplot(1,4,3)
plot(p_val/barsa, rho_b, '-', 'color', [0 0 80]/255, 'linewidth', 3.5); hold on
plot(p_val/barsa, rho_b_sat, '-', 'color', [0 80 80]/255, 'linewidth', 3.5)
plot(P_u, rho_b_txt, '-', 'color', [150 150 255]/255, 'linewidth', 1.25)
plot(P_s, rho_aq_sat_txt, '-', 'color', [0 255 255]/255, 'linewidth', 1.25)
hold off, grid on
xlabel('p [bar]'), ylabel('rho_b [kg/m^3]')
xlim([0 600]), ylim([1000 1100])
title('Brine, T=80 C, S=1e5 ppm')
legend('no CO2', 'sat. CO2', 'no CO2, fcn', 'sat. CO2, fcn', ...
       'location', 'southeast')

subplot(1,4,4)
plot(p_val/barsa, mu_b*1e3, '-', 'color', [0 0 80]/255, 'linewidth', 3.5); hold on
plot(p_val/barsa, mu_b_sat*1e3, '-', 'color', [0 80 80]/255, 'linewidth', 3.5)
plot(P_u, mu_u, '-', 'color', [150 150 255]/255, 'linewidth', 1.25)
plot(P_s, mu_s, '-', 'color', [0 255 255]/255, 'linewidth', 1.25)
hold off, grid on
xlabel('p [bar]'), ylabel('mu_b [cP]')
xlim([0 600]), ylim([0.3 0.6])

print(gcf, fullfile(outdir, 'fig_pvt_props.png'), '-dpng', '-r130');
fprintf('\nFigure saved: fig_pvt_props.png\n');
