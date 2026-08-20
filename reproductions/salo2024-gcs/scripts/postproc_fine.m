% Fine-mesh diffusion comparison (paper Fig. 6 style): dissolved CO2 (rs)
% at t = 30 days, D = 0 vs D = 2e-11 m2/s, 105,279 cells.
set(0, 'defaultfigurevisible', 'off');
outdir = '/tmp/claude-1000/-home-adriano-codes-MRST/ad017599-1c83-453b-a05e-6b33ad5bc00f/scratchpad';

d0 = load(fullfile(outdir, 'diff_fine_D0.mat'));
d2 = load(fullfile(outdir, 'diff_fine_D2e11.mat'));

yq = linspace(min(d0.centroids(:,2)), max(d0.centroids(:,2)), 700);
zq = linspace(min(d0.centroids(:,3)), max(d0.centroids(:,3)), 460);
[Yq, Zq] = meshgrid(yq, zq);

cases = {d0, d2};
names = {'D = 0', 'D = 2e-11 m^2/s'};
figure('Position', [50 50 1200 700]);
for c = 1:2
    d = cases{c};
    rs_end = d.rs_all(:, end);
    F = griddata(d.centroids(:,2), d.centroids(:,3), rs_end, Yq, Zq, 'nearest');
    subplot(2,1,c);
    imagesc(yq, zq, F);
    axis tight
    set(gca, 'YDir', 'reverse');
    colormap(flipud(bone)); colorbar;
    xlabel('y [m]'), ylabel('depth [m]')
    title(sprintf('%s   (max rs = %.4f, mean rs = %.5f)', names{c}, ...
          max(rs_end), mean(rs_end)));
end
print(gcf, fullfile(outdir, 'fig_diffusion_fine30d.png'), '-dpng', '-r130');

fprintf('\n%-32s %-12s %-12s\n', 'Metric (t = 30 d, fine)', 'D = 0', 'D = 2e-11');
fprintf('%-32s %-12.4f %-12.4f\n', 'Max rs', ...
        max(d0.rs_all(:,end)), max(d2.rs_all(:,end)));
fprintf('%-32s %-12d %-12d\n', 'Cells with rs > 1e-3', ...
        sum(d0.rs_all(:,end) > 1e-3), sum(d2.rs_all(:,end) > 1e-3));
fprintf('%-32s %-12.5f %-12.5f\n', 'Mean rs', ...
        mean(d0.rs_all(:,end)), mean(d2.rs_all(:,end)));
fprintf('%-32s %-12.1f %-12.1f\n', 'Runtime [s]', d0.t_elapsed, d2.t_elapsed);
fprintf('\nFigure saved: fig_diffusion_fine30d.png\n');
