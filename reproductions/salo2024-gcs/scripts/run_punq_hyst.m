% Headless driver for the PUNQ-S3 hysteresis validation (paper Fig. 4).
% Based on co2lab-mit/examples/relativePermeabilityHysteresisExample.m.
% Changes for a headless Octave run:
%  - Deck path points to the local mrst-adblackoil-gcs repository.
%  - The ode23/deval hydrostatic initialization is replaced by explicit
%    integration (Octave has no deval).
%  - The linear solver is forced to BackslashSolverAD (no MEX in Octave).
%  - No plots. States are saved to a .mat file for post-processing.
% Env vars: KRHYST (0/1), NSTEPS (0 = all steps), OUTTAG (file tag).

krhyst = str2double(getenv('KRHYST'));
nsteps = str2double(getenv('NSTEPS'));
outtag = getenv('OUTTAG');
outdir = '/tmp/claude-1000/-home-adriano-codes-MRST/ad017599-1c83-453b-a05e-6b33ad5bc00f/scratchpad';

run('/home/adriano/codes/MRST/startup.m');
mrstModule add ad-props deckformat ad-core ad-blackoil
gravity reset on;
mrstVerbose off

fn = ['/home/adriano/codes/MRST/mrst-adblackoil-gcs/kr-hysteresis/' ...
      'input_files_kr/punq-s3/CASE2.DATA'];

%% Input deck and generate grid, rock and fluid
deck = convertDeckUnits(readEclipseDeck(fn));
G = initEclipseGrid(deck);
shiftZup = 840 - min(G.nodes.coords(:, 3));
G.nodes.coords(:, 3) = G.nodes.coords(:, 3) + shiftZup;
G = computeGeometry(G);

rock  = initEclipseRock(deck);
rock  = compressRock(rock, G.cells.indexMap);

fluid = initDeckADIFluid(deck);
rock.regions.saturation = ones(G.cells.num, 1); % 1 saturation region
if krhyst == 1
    fluid.krHyst = 2;   % imbibition region where hysteresis is active
else
    fluid = rmfield(fluid, 'krHyst');           % no hysteresis
end

%% Model
model = selectModelFromDeck(G, rock, fluid, deck);
model = model.validateModel();  % Set up the state function groups
nls = getNonLinearSolver(model, 'TimestepStrategy', 'iteration');
nls.LinearSolver = BackslashSolverAD();
nls.useLinesearch = true;
nls.maxIterations = 15;
nls.maxTimestepCuts = 12;
nls.acceptanceFactor = 2;

% Pore volume multiplier of 1000 on boundary cells (Juanes et al., 2006)
idAct = find(deck.GRID.ACTNUM);
idG = (1:prod(deck.GRID.cartDims))';
idG_mult = idG(deck.GRID.PORV==1000);
cellsb = ismember(idAct, idG_mult);
model.operators.pv(cellsb) = model.operators.pv(cellsb)*1e3;
model.FlowPropertyFunctions = ...
model.FlowPropertyFunctions.setStateFunction('RelativePermeability', ...
                                             HystereticRelativePermeability(model));

%% Initialization (hydrostatic, explicit integration instead of ode23/deval)
g = norm(gravity);
[z_0, z_max] = deal(min(G.cells.centroids(:,3)), max(G.cells.centroids(:,3)));
p_r = g*fluid.rhoWS*z_0;
nz = 2000;
zv = linspace(z_0, z_max, nz)';
dz = zv(2) - zv(1);
pv_hydro = zeros(nz,1);  pv_hydro(1) = p_r;
for k = 2:nz
    dpdz = g * fluid.bW(pv_hydro(k-1)) * fluid.rhoWS;
    pmid = pv_hydro(k-1) + 0.5*dz*dpdz;          % midpoint rule
    pv_hydro(k) = pv_hydro(k-1) + dz * g * fluid.bW(pmid) * fluid.rhoWS;
end
p0 = interp1(zv, pv_hydro, G.cells.centroids(:,3));
s0  = repmat([1, 0], [G.cells.num, 1]);  % fully saturated in water
rs0 = 0;
rv0 = 0;
state0 = struct('s', s0, 'rs', rs0, 'rv', rv0, 'pressure', p0);

%% Deck schedule into a MRST schedule
schedule = convertDeckScheduleToMRST(model, deck);

% Well rate at surface conditions ('resv' wells are not supported here)
resv = 18/day;                                  % rm3/day per well
pavg = mean(p0([schedule.control(1).W.cells]));
rate = resv*fluid.bG(pavg);                     % rm3 to sm3
[schedule.control(1).W.val] = deal(rate);
[schedule.control(2).W.val] = deal(0);
[schedule.control(1).W.type] = deal('rate');
[schedule.control(2).W.type] = deal('rate');
nwell = numel(schedule.control(1).W);
for n=1:nwell
    schedule.control(1).W(n).lims.bhp = 160*barsa;
end

if nsteps > 0   % truncated schedule for the smoke test
    schedule.step.val = schedule.step.val(1:nsteps);
    schedule.step.control = schedule.step.control(1:nsteps);
end

%% Simulation
t_start = tic;
[wellSols, states, report] = simulateScheduleAD(state0, model, schedule, ...
                                                'NonLinearSolver', nls);
t_elapsed = toc(t_start);
fprintf('\nSimulation done in %.1f s (%d steps)\n', t_elapsed, numel(states));

%% Save results
ns = numel(states);
sg_all = zeros(G.cells.num, ns);
p_all  = zeros(G.cells.num, ns);
for n = 1:ns
    sg_all(:,n) = states{n}.s(:,2);
    p_all(:,n)  = states{n}.pressure;
end
tvec = cumsum(schedule.step.val);
indexMap = G.cells.indexMap;
cartDims = deck.GRID.cartDims;
save('-v7', fullfile(outdir, ['punq_' outtag '.mat']), ...
     'sg_all', 'p_all', 'tvec', 'indexMap', 'cartDims', 'p0', 't_elapsed');
fprintf('Results saved: punq_%s.mat | max Sg final = %.3f | sum Sg>0.05: %d cells\n', ...
        outtag, max(sg_all(:,end)), sum(sg_all(:,end) > 0.05));
