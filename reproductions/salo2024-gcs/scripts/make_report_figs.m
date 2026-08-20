% Extra figures for the LaTeX report:
%  1) Hysteretic gas relative permeability curves (drainage, bounding
%     imbibition, one scanning curve) from the PUNQ-S3 deck.
%  2) PUNQ-S3 plume metrics vs time for both cases.
set(0, 'defaultfigurevisible', 'off');
outdir = '/tmp/claude-1000/-home-adriano-codes-MRST/ad017599-1c83-453b-a05e-6b33ad5bc00f/scratchpad';

run('/home/adriano/codes/MRST/startup.m');
mrstModule add ad-props deckformat ad-core ad-blackoil
addpath('/home/adriano/codes/MRST/core/utils/octave_only/mrst');

%% 1. kr hysteresis curves
fn = ['/home/adriano/codes/MRST/mrst-adblackoil-gcs/kr-hysteresis/' ...
      'input_files_kr/punq-s3/CASE2.DATA'];
deck  = convertDeckUnits(readEclipseDeck(fn));
fluid = initDeckADIFluid(deck);
fluid.krHyst = 2;
% Scanning curves via the standalone repository's addScanKr
addpath('/home/adriano/codes/MRST/mrst-adblackoil-gcs/kr-hysteresis/code');
G    = computeGeometry(initEclipseGrid(deck));
rock = compressRock(initEclipseRock(deck), G.cells.indexMap);
fluid = addScanKr(fluid, rock.regions.imbibition, 0.02);

sgd = linspace(fluid.krPts.g(1,2), fluid.krPts.g(1,3), 50)';
krd = fluid.krG{1}(sgd);                       % primary drainage
sgi = linspace(fluid.krPts.g(2,2), fluid.krPts.g(2,3), 50)';
kri = fluid.krG{2}(sgi);                       % bounding imbibition
sgs = linspace(0.28, 0.4, 50)';
krs = fluid.krGi{1}(sgs, repelem(0.4, 50, 1)); % scanning curve, sg_max=0.4

figure('Position', [100 100 480 400]);
plot(sgd, krd, '-r', 'linewidth', 2); hold on
plot(sgi, kri, '--b', 'linewidth', 2);
plot(sgs, krs, '-.k', 'linewidth', 1.5); hold off
grid on
xlabel('S_g [-]'), ylabel('k_{rg} [-]')
legend('drainage k_{rg}^d', 'bounding imbibition k_{rg}^{ib}', ...
       'scanning curve (S_{gi}=0.4)', 'location', 'northwest')
xlim([0 1]), ylim([0 1])
print(gcf, fullfile(outdir, 'fig_kr_hysteresis.png'), '-dpng', '-r150');

%% 2. PUNQ-S3 plume metrics vs time
c1 = load(fullfile(outdir, 'punq_case1_nohyst.mat'));
c2 = load(fullfile(outdir, 'punq_case2_hyst.mat'));
ty = c1.tvec/365.25/86400;                     % years
n1 = sum(c1.sg_all > 0.05, 1);                 % plume extent (cells)
n2 = sum(c2.sg_all > 0.05, 1);
m1 = max(c1.sg_all, [], 1);                    % max saturation
m2 = max(c2.sg_all, [], 1);

figure('Position', [100 100 900 360]);
subplot(1,2,1)
semilogx(ty, n1, '-r', 'linewidth', 2); hold on
semilogx(ty, n2, '-b', 'linewidth', 2); hold off
grid on
xlabel('t [years]'), ylabel('cells with S_g > 0.05')
legend('no hysteresis', 'hysteresis', 'location', 'northwest')
subplot(1,2,2)
semilogx(ty, m1, '-r', 'linewidth', 2); hold on
semilogx(ty, m2, '-b', 'linewidth', 2); hold off
grid on
xlabel('t [years]'), ylabel('max S_g [-]')
print(gcf, fullfile(outdir, 'fig_punq_evolution.png'), '-dpng', '-r150');
fprintf('Report figures saved.\n');
