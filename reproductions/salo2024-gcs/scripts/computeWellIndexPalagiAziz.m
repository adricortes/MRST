function [WI, info] = computeWellIndexPalagiAziz(G, rock, wellCell, radius, varargin)
% Simplified Peaceman-type well index for an unstructured (PEBI/Voronoi)
% grid, following Palagi & Aziz, "Modeling Vertical and Horizontal Wells
% With Voronoi Grid," SPE 24072 (SPE Res. Eng., Feb. 1994), Eq. 17.
%
% Unlike MRST's own computeWellIndex (core/utils/computeWellIndex.m),
% which approximates the well cell as an axis-aligned Cartesian block
% using its node bounding box (dx, dy) regardless of its true shape, this
% derives the equivalent radius from the well cell's ACTUAL connections:
%
%   ln(req) = (sum_j T_ij*ln(L_ij) - theta*Kh) / sum_j T_ij
%   WI      = theta*Kh / (ln(req/rw) + skin)
%
% where T_ij is the two-point face transmissibility to neighbor j, L_ij
% the distance between cell centroids, and theta = 2*pi for an interior
% well (matching the implicit theta = 2*pi in MRST's own WI = 2*pi*Kh/...).
%
% Kh (permeability-thickness) is computed exactly as MRST's own
% computeWellIndex does, so the two well indices differ ONLY in how the
% equivalent radius is derived, not in Kh, skin, or rw.
%
% PARAMETERS:
%   G, rock   - As for computeWellIndex.
%   wellCell  - Single cell index (this model handles one well cell only).
%   radius    - Well radius (rw).
%
% OPTIONAL PARAMETERS:
%   'Dir'  - Well direction, 'x'/'y'/'z' (default 'z'), same convention as
%            computeWellIndex: picks which cell dimension is the "length"
%            (ell) contributing to Kh, vs. the transverse plane.
%   'Skin' - Skin factor (default 0).
%
% RETURNS:
%   WI   - Palagi-Aziz well index.
%   info - Struct with Kh, req, nConn (number of interior-face
%          connections used), Tij, Lij, for diagnostics.

    opt = struct('Dir', 'z', 'Skin', 0);
    opt = merge_options(opt, varargin{:});

    assert(numel(wellCell) == 1, 'This model handles one well cell only.');

    [dx, dy, dz] = cellDims(G, wellCell);
    k = rock.perm(wellCell, :);
    if size(k, 2) == 1
        k = k(:, [1 1 1]);
    elseif size(k, 2) == 2
        k = k(:, [1 2 2]);
    end

    switch lower(opt.Dir)
        case 'x', k1 = k(1); k2 = k(2); ell = dx;   % noqa: transverse ky,kz not needed here
        case 'y', k1 = k(1); k2 = k(3); ell = dy;
        otherwise, k1 = k(1); k2 = k(2); ell = dz;   % 'z'
    end
    Kh = ell * sqrt(k1 * k2);

    hT = computeTrans(G, rock);
    ixc   = G.cells.facePos;
    frng  = ixc(wellCell):(ixc(wellCell + 1) - 1);
    faces = G.cells.faces(frng, 1);
    hTwc  = hT(frng);

    nb = G.faces.neighbors(faces, :);
    isInterior = all(nb ~= 0, 2);
    faces = faces(isInterior);
    hTwc  = hTwc(isInterior);
    nb    = nb(isInterior, :);
    neigh = sum(nb, 2) - wellCell;

    nConn = numel(faces);
    Tij = zeros(nConn, 1);
    Lij = zeros(nConn, 1);
    for i = 1:nConn
        cN   = neigh(i);
        frngN = ixc(cN):(ixc(cN + 1) - 1);
        fN   = G.cells.faces(frngN, 1);
        hTN  = hT(frngN);
        hTo  = hTN(fN == faces(i));
        Tij(i) = 1 / (1/hTwc(i) + 1/hTo(1));
        Lij(i) = norm(G.cells.centroids(cN, :) - G.cells.centroids(wellCell, :));
    end

    theta = 2 * pi;
    ln_req = (sum(Tij .* log(Lij)) - theta * Kh) / sum(Tij);
    req = exp(ln_req);
    WI = theta * Kh / (log(req / radius) + opt.Skin);

    info = struct('Kh', Kh, 'req', req, 'nConn', nConn, 'Tij', Tij, 'Lij', Lij, ...
                   'dx', dx, 'dy', dy, 'dz', dz);
end
