# Full 500-year, two-phase PUNQ-S3 CO2 injection: TPFA vs MPFA (tensor
# assembly, with the new gravity operator), both MRST. Inputs:
# results/punq_fine_nohyst.mat (TPFA, existing fair2 reference) and
# results/punq_mpfa_pvdgfine.mat (MPFA, run_punq_pvdgfine_mpfa.m).
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

RESDIR = '/home/adriano/codes/MRST/reproductions/punqs3-geos-vs-mrst/results/'
YR = 365.0 * 24 * 3600

tp = loadmat(RESDIR + 'punq_fine_nohyst.mat')
mp = loadmat(RESDIR + 'punq_mpfa_pvdgfine.mat')

def unpack(d):
    sg = d['sg_all']            # (ncells, nsteps)
    p = d['p_all']
    t = d['tvec'].ravel() / YR  # years
    idx = d['indexMap'].ravel().astype(int)
    cd = tuple(int(v) for v in d['cartDims'].ravel())
    return sg, p, t, idx, cd

sg_t, p_t, t_t, idx_t, cd_t = unpack(tp)
sg_m, p_m, t_m, idx_m, cd_m = unpack(mp)

print(f'TPFA: {sg_t.shape[0]} cells, {sg_t.shape[1]} steps, cartDims={cd_t}')
print(f'MPFA: {sg_m.shape[0]} cells, {sg_m.shape[1]} steps, cartDims={cd_m}')
assert cd_t == cd_m, 'grids differ between runs'
assert np.array_equal(idx_t, idx_m), 'active-cell index maps differ between runs'
assert sg_t.shape[1] == sg_m.shape[1], 'different number of saved steps'
print(f'max |t_mpfa - t_tpfa| = {np.max(np.abs(t_t - t_m)):.6g} years (should be ~0)')

nx, ny, nz = cd_t
full = lambda: np.full(nx * ny * nz, np.nan)

# ---- summary metrics over time ----
print('\n=== max Sg and plume size (Sg > 0.05) over time ===')
print(f'{"t [y]":>7s} {"maxSg TPFA":>11s} {"maxSg MPFA":>11s} '
      f'{"plume TPFA":>11s} {"plume MPFA":>11s}')
for i in range(0, len(t_t), max(1, len(t_t)//10)):
    print(f'{t_t[i]:7.1f} {sg_t[:,i].max():11.4f} {sg_m[:,i].max():11.4f} '
          f'{int((sg_t[:,i]>0.05).sum()):11d} {int((sg_m[:,i]>0.05).sum()):11d}')

dsg_final = sg_m[:, -1] - sg_t[:, -1]
print(f'\n=== Final-time (t={t_t[-1]:.0f}y) cellwise Sg difference (MPFA - TPFA) ===')
print(f'mean |dSg| = {np.mean(np.abs(dsg_final)):.4f}, max |dSg| = {np.max(np.abs(dsg_final)):.4f}')
print(f'mean Sg  TPFA = {sg_t[:,-1].mean():.4f}, MPFA = {sg_m[:,-1].mean():.4f}')

# ---- plume-shape figure, final time, layer by layer ----
def to_grid(vals, idx):
    v = full()
    v[idx - 1] = vals
    return v.reshape(nz, ny, nx)

g_t = to_grid(sg_t[:, -1], idx_t)
g_m = to_grid(sg_m[:, -1], idx_m)
g_d = g_m - g_t

fig, axs = plt.subplots(3, nz, figsize=(2.6*nz, 7.2))
cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('0.85')
cmapd = plt.get_cmap('RdBu_r').copy(); cmapd.set_bad('0.85')
vmax_d = np.nanmax(np.abs(g_d))
im_d = None
for k in range(nz):
    axs[0, k].imshow(g_t[k], origin='lower', cmap=cmap, vmin=0, vmax=0.7)
    axs[0, k].set_title(f'layer {k+1}'); axs[0, k].set_xticks([]); axs[0, k].set_yticks([])
    axs[1, k].imshow(g_m[k], origin='lower', cmap=cmap, vmin=0, vmax=0.7)
    axs[1, k].set_xticks([]); axs[1, k].set_yticks([])
    im_d = axs[2, k].imshow(g_d[k], origin='lower', cmap=cmapd, vmin=-vmax_d, vmax=vmax_d)
    axs[2, k].set_xticks([]); axs[2, k].set_yticks([])
axs[0, 0].set_ylabel('TPFA $S_g$'); axs[1, 0].set_ylabel('MPFA $S_g$'); axs[2, 0].set_ylabel('MPFA $-$ TPFA')
fig.suptitle(f'PUNQ-S3 CO2 plume at t = {t_t[-1]:.0f} years: TPFA vs MPFA (tensor assembly)')
fig.colorbar(im_d, ax=axs[2, :], shrink=0.85, label='$\\Delta S_g$')
fig.tight_layout()
fig.savefig(RESDIR + '../figures/cmp_tpfa_mpfa_twophase_final.png', dpi=150)
print('\nsaved figures/cmp_tpfa_mpfa_twophase_final.png')

# ---- time series of max Sg and plume volume ----
fig, axs = plt.subplots(1, 2, figsize=(10, 4))
axs[0].plot(t_t, sg_t.max(axis=0), label='TPFA', color='#4269D0', lw=2)
axs[0].plot(t_m, sg_m.max(axis=0), label='MPFA', color='#C4520A', lw=2, ls='--')
axs[0].set_xlabel('time [years]'); axs[0].set_ylabel('max $S_g$'); axs[0].legend()
axs[1].plot(t_t, (sg_t > 0.05).sum(axis=0), label='TPFA', color='#4269D0', lw=2)
axs[1].plot(t_m, (sg_m > 0.05).sum(axis=0), label='MPFA', color='#C4520A', lw=2, ls='--')
axs[1].set_xlabel('time [years]'); axs[1].set_ylabel('cells with $S_g$ > 0.05'); axs[1].legend()
fig.tight_layout()
fig.savefig(RESDIR + '../figures/cmp_tpfa_mpfa_twophase_timeseries.png', dpi=150)
print('saved figures/cmp_tpfa_mpfa_twophase_timeseries.png')
