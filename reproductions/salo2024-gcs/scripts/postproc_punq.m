% Post-process the two PUNQ-S3 runs: CO2 saturation maps per layer at
% t = 500 years, one row per case (compare with paper Fig. 4, bottom row).
set(0, 'defaultfigurevisible', 'off');
outdir = '/tmp/claude-1000/-home-adriano-codes-MRST/ad017599-1c83-453b-a05e-6b33ad5bc00f/scratchpad';

c1 = load(fullfile(outdir, 'punq_case1_nohyst.mat'));
c2 = load(fullfile(outdir, 'punq_case2_hyst.mat'));
nx = c1.cartDims(1); ny = c1.cartDims(2); nz = c1.cartDims(3);

% ECLIPSE-like colormap from the example script
colors = [0 0 1; 0 1 1; 0 1 0; 1 1 0; 1 0 0];
cmap = [interp1([0 1], colors(1:2,:), linspace(0,1,30));
        interp1([0 1], colors(2:3,:), linspace(0,1,25));
        interp1([0 1], colors(3:4,:), linspace(0,1,25));
        interp1([0 1], colors(4:5,:), linspace(0,1,20))];

cases = {c1, c2};
names = {'No hysteresis', 'Hysteresis'};
figure('Position', [50 50 1500 550]);
for c = 1:2
    d = cases{c};
    % Inactive cells get a sentinel value, shown gray (gnuplot has no alpha)
    sgc = -0.011*ones(nx*ny*nz, 1);
    sgc(d.indexMap) = d.sg_all(:, end);      % final state, t = 500 y
    sg3 = reshape(sgc, nx, ny, nz);
    for k = 1:nz
        subplot(2, nz, (c-1)*nz + k);
        imagesc(squeeze(sg3(:,:,k))');
        axis equal tight; caxis([-0.011 1]); colormap([0.85 0.85 0.85; cmap]);
        set(gca, 'XTick', [], 'YTick', []);
        if c == 1, title(sprintf('Layer %d', k)); end
        if k == 1, ylabel(names{c}); end
    end
end
subplot(2, nz, 2*nz); colorbar;
print(gcf, fullfile(outdir, 'fig_punq_sg500y.png'), '-dpng', '-r130');

% Metrics
fprintf('\n%-28s %-14s %-14s\n', 'Metric (t = 500 y)', 'No hysteresis', 'Hysteresis');
fprintf('%-28s %-14.3f %-14.3f\n', 'Max Sg', ...
        max(c1.sg_all(:,end)), max(c2.sg_all(:,end)));
fprintf('%-28s %-14d %-14d\n', 'Cells with Sg > 0.05', ...
        sum(c1.sg_all(:,end) > 0.05), sum(c2.sg_all(:,end) > 0.05));
fprintf('%-28s %-14d %-14d\n', 'Cells with Sg > 0.3', ...
        sum(c1.sg_all(:,end) > 0.3), sum(c2.sg_all(:,end) > 0.3));
fprintf('%-28s %-14.3f %-14.3f\n', 'Mean Sg where Sg > 0.05', ...
        mean(c1.sg_all(c1.sg_all(:,end)>0.05, end)), ...
        mean(c2.sg_all(c2.sg_all(:,end)>0.05, end)));

% End of injection (step 10, t = 10 y)
fprintf('\n%-28s %-14s %-14s\n', 'Metric (t = 10 y)', 'No hysteresis', 'Hysteresis');
fprintf('%-28s %-14.3f %-14.3f\n', 'Max Sg', ...
        max(c1.sg_all(:,10)), max(c2.sg_all(:,10)));
fprintf('%-28s %-14d %-14d\n', 'Cells with Sg > 0.05', ...
        sum(c1.sg_all(:,10) > 0.05), sum(c2.sg_all(:,10) > 0.05));
fprintf('\nFigure saved: fig_punq_sg500y.png\n');
