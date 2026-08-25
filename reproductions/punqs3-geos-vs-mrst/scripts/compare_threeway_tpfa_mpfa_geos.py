# Three-way comparison of the PUNQ-S3 no-hysteresis, full 500-year,
# two-phase CO2 case: TPFA (MRST), MPFA tensor-assembly with gravity
# (MRST), and GEOS. The core question: does the more consistent
# discretization (MPFA) bring MRST closer to GEOS than TPFA does, on the
# vertical/buoyancy-direction connections flagged as non-K-orthogonal in
# check_tpfa_consistency.m?
#
# Inputs (this directory): geos_fair2_data.npz, mrst_fine_data.npz
# (TPFA), mrst_mpfa_data.npz (MPFA). All three produced by scripts
# already in this directory (extract_geos_fair2.py, mrst_mat_to_npz.py).
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

YR = 365.25 * 86400.0
C_TPFA, C_MPFA, C_GEOS = '#4269D0', '#3CA951', '#C4520A'

geos = np.load('geos_fair2_data.npz')
tpfa = np.load('mrst_fine_data.npz')
mpfa = np.load('mrst_mpfa_data.npz')

key = lambda ijk: [tuple(r) for r in ijk]
gk, tk, mk = key(geos['ijk']), key(tpfa['ijk']), key(mpfa['ijk'])
assert set(gk) == set(tk) == set(mk), 'active cell sets differ between runs'
assert tk == mk, 'TPFA and MPFA cell ordering differs -- unexpected, same grid/deck'
g_of = {c: n for n, c in enumerate(gk)}
perm_g = np.array([g_of[c] for c in tk])   # geos row for each tpfa/mpfa row (same order)

t_g = geos['times'] / YR
t_m = tpfa['times'] / YR                    # same schedule for tpfa and mpfa
assert np.allclose(t_m, mpfa['times'] / YR), 'tpfa/mpfa step times differ'

sg_g, p_g, mg_g = geos['sg'][perm_g], geos['p'][perm_g], geos['mg'][perm_g]
sg_t, p_t = tpfa['sg'], tpfa['p']
sg_m, p_m = mpfa['sg'], mpfa['p']
vol, poro_ref = geos['vol'][perm_g], geos['poro_ref'][perm_g]
interior = poro_ref < 1.0
pv_ref = poro_ref * vol
ijk = tpfa['ijk']

# ---- interpolate MRST (both) onto GEOS snapshot times for cellwise diffs ----
def interp_rows(t_src, X, t_dst):
    out = np.empty((X.shape[0], len(t_dst)))
    for r in range(X.shape[0]):
        out[r] = np.interp(t_dst, np.r_[0, t_src], np.r_[0, X[r]])
    return out

sg_ti = interp_rows(t_m, sg_t, t_g)
sg_mi = interp_rows(t_m, sg_m, t_g)
l1_t, linf_t = np.abs(sg_g - sg_ti).mean(axis=0), np.abs(sg_g - sg_ti).max(axis=0)
l1_m, linf_m = np.abs(sg_g - sg_mi).mean(axis=0), np.abs(sg_g - sg_mi).max(axis=0)

print('=== Cellwise |Sg - Sg_GEOS|, interpolated onto GEOS times ===')
print(f'{"t [y]":>7s} {"mean TPFA":>10s} {"mean MPFA":>10s} '
      f'{"max TPFA":>9s} {"max MPFA":>9s}  closer to GEOS')
for i in range(0, len(t_g), max(1, len(t_g)//10)):
    closer = 'MPFA' if l1_m[i] < l1_t[i] else 'TPFA'
    print(f'{t_g[i]:7.1f} {l1_t[i]:10.4f} {l1_m[i]:10.4f} '
          f'{linf_t[i]:9.4f} {linf_m[i]:9.4f}  {closer}')
i_end = -1
print(f'\nAt t=500y: mean|dSg| TPFA-GEOS = {l1_t[i_end]:.4f}, '
      f'MPFA-GEOS = {l1_m[i_end]:.4f}  '
      f'({"MPFA" if l1_m[i_end]<l1_t[i_end] else "TPFA"} closer)')

# ---- CO2 mass in place ----
_tab = np.array([l.split() for l in open(
    '/home/adriano/codes/MRST/punqs3-ref/pvdg_fine.txt') if not l.startswith('#')], float)
pvdg_p, pvdg_b = _tab[1:, 0], 1.0 / _tab[1:, 1]
cr, pref = 5e-5 / 1e5, 1e5
def co2_mass(sg, p):
    rho = 1.8 * np.interp(p, pvdg_p, pvdg_b)
    pv = pv_ref[:, None] * (1 + cr * (p - pref))
    return (sg * pv * rho).sum(axis=0)
mass_t, mass_m, mass_g = co2_mass(sg_t, p_t), co2_mass(sg_m, p_m), mg_g.sum(axis=0)

i10g, i10m = np.argmin(np.abs(t_g - 10)), np.argmin(np.abs(t_m - 10))
print(f'\n=== CO2 mass at end of injection (t=10y) ===')
print(f'GEOS {mass_g[i10g]/1e9:.3f} Mt | TPFA {mass_t[i10m]/1e9:.3f} Mt '
      f'(ratio {mass_t[i10m]/mass_g[i10g]:.4f}) | '
      f'MPFA {mass_m[i10m]/1e9:.3f} Mt (ratio {mass_m[i10m]/mass_g[i10g]:.4f})')

# ---- figure: observation cells ----
OBS = [(13, 18, 1), (7, 21, 1), (11, 11, 1)]
obs_rows = [tk.index(o) for o in OBS]
fig, axs = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
for a, (o, r) in zip(axs, zip(OBS, obs_rows)):
    a.plot(np.r_[0, t_m], np.r_[0, sg_t[r]], '-', color=C_TPFA, lw=2, label='TPFA')
    a.plot(np.r_[0, t_m], np.r_[0, sg_m[r]], '-', color=C_MPFA, lw=2, label='MPFA')
    a.plot(t_g, sg_g[r], '--', color=C_GEOS, lw=2, marker='o', ms=4,
           markevery=4, label='GEOS')
    a.axvline(10, color='0.6', lw=1, ls=':')
    a.set_title(f'cell ({o[0]},{o[1]},{o[2]})')
    a.set_xlabel('time [yr]')
    a.grid(alpha=0.25)
axs[0].set_ylabel('gas saturation $S_g$ [-]')
axs[0].legend(frameon=False)
fig.suptitle('PUNQ-S3, 500y: $S_g$ at observation cells -- TPFA vs MPFA vs GEOS')
fig.tight_layout()
fig.savefig('../figures/cmp_threeway_sg_obs.png', dpi=150)
print('\nsaved figures/cmp_threeway_sg_obs.png')

# ---- figure: mass in place + cellwise agreement with GEOS over time ----
fig, axs = plt.subplots(1, 2, figsize=(11, 4))
axs[0].plot(np.r_[0, t_m], np.r_[0, mass_t]/1e9, '-', color=C_TPFA, lw=2, label='TPFA')
axs[0].plot(np.r_[0, t_m], np.r_[0, mass_m]/1e9, '-', color=C_MPFA, lw=2, label='MPFA')
axs[0].plot(t_g, mass_g/1e9, '--', color=C_GEOS, lw=2, label='GEOS')
axs[0].set_title('CO$_2$ mass in place'); axs[0].set_ylabel('mass [Mt]')
axs[0].set_xlabel('time [yr]'); axs[0].legend(frameon=False); axs[0].grid(alpha=0.25)

axs[1].plot(t_g, l1_t, '-', color=C_TPFA, lw=2, label='mean $|\\Delta S_g|$ TPFA$-$GEOS')
axs[1].plot(t_g, l1_m, '-', color=C_MPFA, lw=2, label='mean $|\\Delta S_g|$ MPFA$-$GEOS')
axs[1].set_title('cellwise agreement with GEOS'); axs[1].set_ylabel('mean $|\\Delta S_g|$')
axs[1].set_xlabel('time [yr]'); axs[1].legend(frameon=False); axs[1].grid(alpha=0.25)
fig.tight_layout()
fig.savefig('../figures/cmp_threeway_inventory.png', dpi=150)
print('saved figures/cmp_threeway_inventory.png')

# ---- figure: layer maps at t=500y ----
nx, ny, nz = 19, 28, 5
full = lambda: np.full(nx*ny*nz, np.nan)
lin = (ijk[:, 0]-1) + (ijk[:, 1]-1)*nx + (ijk[:, 2]-1)*nx*ny
fig, axs = plt.subplots(3, nz, figsize=(2.4*nz, 7.2))
cmap = plt.get_cmap('viridis').copy(); cmap.set_bad('0.85')
for row, (sg, name) in enumerate([(sg_t[:, -1], 'TPFA'), (sg_m[:, -1], 'MPFA'),
                                  (sg_g[:, -1], 'GEOS')]):
    v = full(); v[lin] = sg
    v3 = v.reshape(nz, ny, nx)
    for k_ in range(nz):
        ax = axs[row, k_]
        im = ax.imshow(v3[k_], origin='lower', cmap=cmap, vmin=0, vmax=0.7)
        ax.set_xticks([]); ax.set_yticks([])
        if row == 0:
            ax.set_title(f'layer {k_+1}')
        if k_ == 0:
            ax.set_ylabel(name, fontsize=12)
fig.suptitle('$S_g$ at t = 500 y: TPFA vs MPFA vs GEOS')
fig.colorbar(im, ax=axs, shrink=0.8, label='$S_g$ [-]')
fig.savefig('../figures/cmp_threeway_maps500y.png', dpi=150)
print('saved figures/cmp_threeway_maps500y.png')

# ---- metrics table ----
def metrics(sg):
    return (sg.max(), (sg > 0.05).sum(), (sg > 0.3).sum(), sg[sg > 0.05].mean())

print('\n=== Plume metrics ===')
hdr = f'{"metric":26s} {"TPFA":>10s} {"MPFA":>10s} {"GEOS":>10s}'
for label, im_, ig_ in [('t = 10 y', i10m, i10g), ('t = 500 y', -1, -1)]:
    t_, m_, g_ = metrics(sg_t[:, im_]), metrics(sg_m[:, im_]), metrics(sg_g[:, ig_])
    print(f'\n--- {label} ---\n{hdr}')
    for n_, name in enumerate(['max Sg', 'cells Sg > 0.05', 'cells Sg > 0.3',
                               'mean Sg where > 0.05']):
        fm = '10.3f' if n_ in (0, 3) else '10d'
        print(f'{name:26s} {t_[n_]:{fm}} {m_[n_]:{fm}} {g_[n_]:{fm}}')
