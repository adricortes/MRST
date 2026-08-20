% Rank the coarse-mesh fault-strip cells by bounding-box aspect ratio
% (Ry = max(dx/dy, dy/dx)) and pick an elongated one close to the paper's
% injector location, for the Peaceman-vs-Palagi-Aziz well-index comparison
% in notebooks/04_diffusion.ipynb, Appendix A.4.
run('/home/adriano/codes/MRST/startup.m');
mrstModule add upr coarsegrid co2lab-mit ad-props ad-blackoil deckformat ad-core
addpath('/home/adriano/codes/MRST/reproductions/salo2024-gcs/scripts');
addpath('/home/adriano/codes/MRST/reproductions/salo2024-gcs/scripts/shim');
addpath('/home/adriano/codes/MRST/core/utils/octave_only/mrst');
mrstVerbose off

G_dat = simpleExtrudedFluidFlowerMesh('coarse');
sealID = [1 3 5 8];
faultID = 6;
G = G_dat.G;
if isfield(G.faces, 'tag'), G.faces = rmfield(G.faces, 'tag'); end
[G, cellmap] = removeCells(G, ismember(G_dat.compartID, sealID));

cid = 1:G_dat.G.cells.num;
fid_orig = cid(G_dat.compartID == faultID);
fid = ismember(cellmap, fid_orig);
rock.poro = 0.3*ones(G.cells.num, 1);
rock.poro(fid) = 0.1;
rock.perm = 1000*ones(G.cells.num, 1);
rock.perm(fid) = 10;
rock.perm = rock.perm*(milli*darcy);

faultCells = find(fid);
fprintf('Number of fault cells: %d\n', numel(faultCells));

% Aspect ratio (bounding box) of every fault cell
[dx, dy, dz] = cellDims(G, faultCells);
Ry = max([dx./dy, dy./dx], [], 2);
[sortedRy, ix] = sort(Ry, 'descend');

fprintf('\nTop 10 most elongated fault cells (bounding-box aspect ratio):\n');
fprintf('%-8s %-8s %-10s %-10s %-10s %-14s\n', 'cell', 'Ry', 'dx[m]', 'dy[m]', 'dz[m]', 'centroid y,z');
for i = 1:10
    c = faultCells(ix(i));
    fprintf('%-8d %-8.2f %-10.4e %-10.4e %-10.4e (%.4f, %.4f)\n', ...
            c, sortedRy(i), dx(ix(i)), dy(ix(i)), dz(ix(i)), ...
            G.cells.centroids(c,2), G.cells.centroids(c,3));
end

% Among the most elongated cells (top quartile), pick the one closest to
% the paper's injector location, so the comparison isolates the effect of
% cell shape rather than also moving to a very different part of the domain.
nTop = max(10, round(numel(faultCells)*0.25));
candidates = faultCells(ix(1:nTop));
dist = sqrt(sum((G.cells.centroids(candidates,:) - [0.005, 0.903, 1000.58]).^2, 2));
[~, j] = min(dist);
fprintf('\nChosen fault cell (elongated, closest to original well y,z): %d\n', candidates(j));
fprintf('Ry = %.2f, centroid = (%.4f, %.4f, %.4f)\n', ...
        Ry(faultCells==candidates(j)), G.cells.centroids(candidates(j),:));
