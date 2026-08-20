% Headless driver for co2lab-spillpoint/examples/firstTrappingExample.m
% Changes relative to the original example:
%  - trapAnalysis uses the edge-based method (false). The cell-based method
%    needs matlab_bgl MEX binaries, which Octave cannot load.
%  - 3-D patch plots are replaced by 2-D maps. The gnuplot toolkit cannot
%    draw 3-D filled quad patches.
%  - The interactive viewer section is skipped (headless run).
set(0, 'defaultfigurevisible', 'off');
outdir = '/tmp/claude-1000/-home-adriano-codes-MRST/ad017599-1c83-453b-a05e-6b33ad5bc00f/scratchpad';

run('/home/adriano/codes/MRST/startup.m');
mrstModule add co2lab-spillpoint co2lab-common co2lab-legacy co2lab-demo coarsegrid

%% Make grid and rock data
[Lx,Ly,H] = deal(10000, 5000, 50);
[nx,ny]   = deal(100, 100);
G  = cartGrid([nx ny 1],[Lx Ly H]);
x  = G.nodes.coords(1:G.nodes.num/2,1)/Lx;
y  = G.nodes.coords(1:G.nodes.num/2,2)/Ly;
z  = G.nodes.coords(1:G.nodes.num/2,3)/H;
zt = z + x - 0.2*sin(5*pi*x).*sin(5*pi*y.^1.5) - 0.075*sin(1.25*pi*y) + 0.15*sin(x+y);
zb = 1 + x;
G.nodes.coords(:,3) = [zt; zb]*H+1000;
G = computeGeometry(G);

%% Construct top surface grid
Gt = topSurfaceGrid(G);

%% Geometric analysis of caprock (spill-point analysis)
res = trapAnalysis(Gt, false);
num_traps = max(res.traps);
fprintf('\nNumber of traps found: %d\n', num_traps);

% Map of top-surface depth with trap cells marked
xc = Gt.cells.centroids(:,1);
yc = Gt.cells.centroids(:,2);
zmap = reshape(Gt.cells.z, nx, ny)';
figure;
imagesc([min(xc) max(xc)], [min(yc) max(yc)], zmap);
axis xy equal tight; colormap(jet); colorbar;
title('Top surface depth (m); white = trap cells');
hold on
tmask = reshape(double(res.traps>0), nx, ny)';
[yy,xx] = find(tmask);
plot((xx-0.5)*Lx/nx, (yy-0.5)*Ly/ny, 'w.', 'MarkerSize', 4);
hold off
print(gcf, fullfile(outdir, 'fig1_depth_and_traps.png'), '-dpng', '-r120');

%% Show connection between traps and spill paths
trap_field = zeros(size(res.traps));
trap_field(res.traps>0) = 2;
for r = [res.cell_lines{:}]'
    for c = 1:numel(r)
        trap_field(r{c}) = 1;
    end
end

figure;
imagesc([min(xc) max(xc)], [min(yc) max(yc)], reshape(trap_field, nx, ny)');
axis xy equal tight; colormap(jet);
title('Traps (red), spill paths (green), other cells (blue)');
for i=1:num_traps
   ind = res.top(i);
   text(Gt.cells.centroids(ind,1), Gt.cells.centroids(ind,2), ...
      num2str(res.traps(ind)), 'Color', .99*[1 1 1], ...
      'FontSize', 12, 'HorizontalAlignment', 'center');
end
print(gcf, fullfile(outdir, 'fig2_traps_and_spillpaths.png'), '-dpng', '-r120');

%% Compute the total trapping capacity
rock2D.poro = 0.25 * ones(G.cells.num, 1);
trap_volumes = volumesOfTraps(Gt, res, 1:num_traps, 'poro', rock2D.poro);

total_trapping_capacity = sum(trap_volumes);
pv = sum(poreVolume(Gt,rock2D));
fprintf('Total trapping capacity is: %6.3e\n', total_trapping_capacity);
fprintf('This amounts to %.2f %% of a total pore volume of %6.2e\n\n', ...
   total_trapping_capacity/pv*100, pv);
[sorted_vols, sorted_ix] = sort(trap_volumes, 'descend');

fprintf('trap ix   | trap vol(m3)  | cells in trap\n');
fprintf('----------+---------------+--------------\n');
tcells = zeros(num_traps,1);
for i = 1:num_traps
   tcells(i) = sum(res.traps == sorted_ix(i));
   fprintf('%7d   |  %10.3e   | %5d\n', ...
       sorted_ix(i), sorted_vols(i), tcells(i));
end
fprintf(['\nTogether, the five largest traps cover %6.2e m3, which represents %3.1f%% of' ...
   '\nthe total trapping capacity of this grid.\n'], ...
   sum(sorted_vols(1:5)), sum(sorted_vols(1:5)) / total_trapping_capacity * 100)

figure;
subplot(2,1,1);
bar(1:num_traps, sorted_vols);
set(gca,'XTick',1:num_traps,'XTickLabel',sorted_ix);
title('Trap volume (m^3), sorted');
subplot(2,1,2);
bar(1:num_traps, tcells);
set(gca,'XTick',1:num_traps,'XTickLabel',sorted_ix);
title('Number of cells in trap');
print(gcf, fullfile(outdir, 'fig3_trap_volumes.png'), '-dpng', '-r120');

fprintf('\nDriver finished. Figures saved to %s\n', outdir);
