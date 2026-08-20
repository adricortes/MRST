% Headless driver for the CO2 diffusion example (paper Fig. 5b/6).
% Based on co2lab-mit/examples/flowWithDiffusionExample.m with mesh='coarse'.
% Changes for a headless Octave run:
%  - Fluid deck comes from the local mrst-adblackoil-gcs repository (equal
%    to the published CO2LabMIT dataset file).
%  - ode23/deval replaced by explicit integration; pdist2 by direct norm.
%  - Default AD backend and BackslashSolverAD (no MEX in Octave).
%  - No plots; states saved to .mat.
% Env vars: DIFFD (pseudo-diffusivity, m2/s), OUTTAG (output file tag).

Denv   = str2double(getenv('DIFFD'));
outtag = getenv('OUTTAG');
meshres = getenv('MESHRES');
if isempty(meshres), meshres = 'coarse'; end
outdir = '/tmp/claude-1000/-home-adriano-codes-MRST/ad017599-1c83-453b-a05e-6b33ad5bc00f/scratchpad';

run('/home/adriano/codes/MRST/startup.m');
mrstModule add upr coarsegrid co2lab-mit ad-props ad-blackoil deckformat ad-core
addpath([outdir '/shim']);                                   % delaunayTriangulation
addpath('/home/adriano/codes/MRST/core/utils/octave_only/mrst'); % fastInterpTable
mrstVerbose off

%% 1. Mesh (cached; generated earlier)
G_dat = simpleExtrudedFluidFlowerMesh(meshres);
% Guard: an invalid mesh must stop the run before hours of simulation
assert(all(G_dat.G.cells.volumes > 0), 'Mesh has non-positive cell volumes');
assert(all(G_dat.G.faces.areas > 0), 'Mesh has non-positive face areas');
fprintf('Mesh %s: %d cells, min volume %.3e\n', meshres, ...
        G_dat.G.cells.num, min(G_dat.G.cells.volumes));
sealID = [1 3 5 8];
faultID = 6;
G = G_dat.G;
if isfield(G.faces, 'tag'), G.faces = rmfield(G.faces, 'tag'); end
[G, cellmap] = removeCells(G, ismember(G_dat.compartID, sealID));

%% 2. Rock
cid = 1:G_dat.G.cells.num;
fid = cid(G_dat.compartID == faultID);
fid = ismember(cellmap, fid);
rock.poro = 0.3*ones(G.cells.num, 1);
rock.poro(fid) = 0.1;
rock.perm = 1000*ones(G.cells.num, 1);
rock.perm(fid) = 10;
rock.perm = rock.perm*(milli*darcy);
rock.regions.saturation = ones(G.cells.num, 1);
rock.regions.saturation(fid) = 3;
rock.regions.rocknum = ones(G.cells.num, 1);

%% 3. Fluid
fn = ['/home/adriano/codes/MRST/mrst-adblackoil-gcs/diffusion/' ...
      'input_files_diff/fluid_props/example_co2brine_1kmDepth_3regions.DATA'];
deck = convertDeckUnits(readEclipseDeck(fn));
deck.REGIONS.ROCKNUM = rock.regions.rocknum;
fluid = initDeckADIFluid(deck);

%% 4. Initialize
gravity reset on
g = norm(gravity);
water_column = 1000;
p_r = 1*barsa + g*fluid.rhoOS*water_column;
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
nz = 2000;
zv = linspace(z_0, z_max, nz)';
dz = zv(2) - zv(1);
ph = zeros(nz,1);  ph(1) = p_r;
for k = 2:nz
    dpdz = g * fluid.bO(ph(k-1), 0, false) * fluid.rhoOS;
    pmid = ph(k-1) + 0.5*dz*dpdz;
    ph(k) = ph(k-1) + dz * g * fluid.bO(pmid, 0, false) * fluid.rhoOS;
end
p0 = interp1(zv, ph, G.cells.centroids(:,3));
s0  = repmat([1, 0], [G.cells.num, 1]);
rs0 = zeros(G.cells.num, 1);
rv0 = 0;
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

%% 5. Wells and times
t = [60*minute 1*day 30*day];
reportTimes = [(12:12:t(1)/minute)*minute, ...
               (2:1:24)*hour, ...
               (1440+5:5:1465)*minute, ...
               ([25 26 28 32 36 40 48 60 72 96 120])*hour, ...
                (6:30)*day];
wellno = 1;
rate   = 8;                                      % mL/min (surface conditions)
injrate = rate*(milli*litre)/(minute*wellno);    % Sm3/s
dist = sqrt(sum((G.cells.centroids - [0.005, 0.903, 1000.58]).^2, 2));
[~, wellInx] = min(dist);
W = addWell([ ], G, rock, wellInx, 'Name', 'I1', 'Dir', 'z', ...
            'Type', 'rate', 'Val', injrate, 'compi', [0, 1], ...
            'refDepth', G.cells.centroids(wellInx, G.griddim), ...
            'Radius', 1e-3);
timesteps = [reportTimes(1) diff(reportTimes)];
assert(sum(timesteps)==t(end), 'sum of timesteps must equal simTime')

%% 6. Model
model = GenericBlackOilModel(G, rock, fluid, 'disgas', true, 'water', false);
model.minimumPressure = min(state0.pressure);
model = model.validateModel();

if Denv > 0
    diffFlux = CO2TotalFluxWithDiffusion(model);
    diffFlux.componentDiffusion = [0 Denv];
    diffFlux.faceAverage = true;
    model.FlowDiscretization.ComponentTotalFlux = diffFlux;
end

nls = getNonLinearSolver(model, 'TimestepStrategy', 'iteration');
nls.LinearSolver = BackslashSolverAD();
nls.useLinesearch = true;
nls.maxIterations = 10;
nls.maxTimestepCuts = 12;
nls.acceptanceFactor = 2;

%% BCs: open aquifer via pore volume multipliers on lateral boundary cells
L = max(G.faces.centroids(:,2));
f = any([G.faces.centroids(:,2) == 0, ...
         G.faces.centroids(:,2) > L-1e-3], 2);
cellsext = unique(reshape(G.faces.neighbors(f, :), [], 1));
cellsext(cellsext==0) = [];
model.operators.pv(cellsext) = model.operators.pv(cellsext)*10^5;
bc = [];

%% Schedule with rate ramp up and down
schedule_inj = simpleSchedule(timesteps, 'W', W, 'bc', bc);
n_ramp = 5;
v = injrate;
injrates = [0.01*v 0.1*v 0.2*v 0.5*v v ...
            0.8*v 0.5*v 0.1*v 0.01*v 0];
tmp = cell(numel(injrates), 1);
schedule = struct('step', schedule_inj.step);
schedule.control = struct('W', tmp, 'bc', tmp, 'src', tmp);
for n=1:numel(injrates)
    schedule.control(n).W = W;
    schedule.control(n).W.val = injrates(n);
    schedule.control(n).bc = bc;
end
idStep = find(cumsum(schedule.step.val) < t(1));
schedule.step.control(idStep) = 1:max(idStep);
schedule.step.control(idStep(end)+1:end) = max(idStep)+1;
idStep2 = find(cumsum(schedule.step.val) > t(2), 1);
schedule.step.control(idStep2:idStep2+(n_ramp -1)) = (1:n_ramp) + n_ramp;
schedule.step.control(idStep2+n_ramp:end) = 2*n_ramp;

%% Simulation
t_start = tic;
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
                                                'NonLinearSolver', nls);
t_elapsed = toc(t_start);
fprintf('\nSimulation done in %.1f s (%d steps)\n', t_elapsed, numel(states));

%% Save results
ns = numel(states);
sg_all = zeros(G.cells.num, ns);
rs_all = zeros(G.cells.num, ns);
for n = 1:ns
    sg_all(:,n) = states{n}.s(:,2);
    rs_all(:,n) = states{n}.rs;
end
tvec = cumsum(schedule.step.val);
centroids = G.cells.centroids;
save('-v7', fullfile(outdir, ['diff_' outtag '.mat']), ...
     'sg_all', 'rs_all', 'tvec', 'centroids', 'fid', 'wellInx', 't_elapsed');
fprintf('Saved diff_%s.mat | max rs final = %.4f | max sg final = %.3f\n', ...
        outtag, max(rs_all(:,end)), max(sg_all(:,end)));
