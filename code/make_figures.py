"""make_figures.py -- regenerate the three paper figures from the analysis.
Run after analyze_pockmarks.py (uses the same driver to obtain curves)."""
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from analyze_pockmarks import main, FIELDS

COL = {"Offshore Zannone Island": "#2980B9",
       "northwestern Calabrian margin": "#C0392B",
       "Malta Plateau": "#27AE60"}
SHORT = {"Offshore Zannone Island": "Zannone (isotropic ~random)",
         "northwestern Calabrian margin": "Calabrian (clustered)",
         "Malta Plateau": "Malta Plateau (regular, anisotropic)"}


def run(shp, dbf, outdir):
    curves = main(shp, dbf, outdir)

    # ---- Fig 1: point patterns ----
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    for a, (name, role) in zip(ax, FIELDS):
        p = curves[name]['pts'] / 1000.0
        a.scatter(p[:, 0], p[:, 1], s=10, alpha=0.75, color=COL[name])
        a.set_aspect('equal'); a.set_xlabel('km'); a.set_ylabel('km')
        a.set_title(f"{SHORT[name]}\nN={len(p)}, NNI={curves[name]['NNI']:.2f}", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{outdir}/figures/fig_fields.pdf", bbox_inches="tight")
    fig.savefig(f"{outdir}/figures/fig_fields.png", dpi=110, bbox_inches="tight")

    # ---- Fig 2: order plane + S(k) ----
    fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
    a = ax[0]
    a.axhline(1, color='#888', ls='--', lw=1); a.axvline(1, color='#888', ls='--', lw=1)
    a.scatter([1], [1], marker='s', s=90, color='k', zorder=5)
    a.annotate('CSR (Poisson)', (1, 1), textcoords='offset points', xytext=(8, 6), fontsize=8)
    a.scatter([1.3], [0.7], marker='^', s=90, color='#7D3C98', zorder=5)
    a.annotate('RSA hard-core\n(liquid-like)', (1.3, 0.7), textcoords='offset points', xytext=(6, -22), fontsize=8, color='#7D3C98')
    a.scatter([1.5], [0.03], marker='*', s=220, color='#E67E22', zorder=5)
    a.annotate('disordered\nHYPERUNIFORM', (1.5, 0.03), textcoords='offset points', xytext=(8, 2), fontsize=8, color='#E67E22')
    a.fill_between([0.3, 1.8], [0.02, 0.02], [1, 1], color='#E67E22', alpha=0.06)
    a.text(0.34, 0.04, 'hyperuniform zone (empty)', fontsize=7.5, color='#E67E22')
    for name, _ in FIELDS:
        a.scatter([curves[name]['NNI']], [curves[name]['sk0']], s=150, color=COL[name],
                  zorder=6, edgecolor='k', lw=0.5)
        a.annotate(SHORT[name].split(' (')[0], (curves[name]['NNI'], curves[name]['sk0']),
                   textcoords='offset points', xytext=(8, 4), fontsize=8.5, color=COL[name])
    a.set_yscale('log'); a.set_xlim(0.3, 1.8); a.set_ylim(0.02, 30)
    a.set_xlabel('NNI  (local order:  <1 clustered  →  >1 regular)')
    a.set_ylabel('S(k→0)   (long range:  →0 hyperuniform  |  >1 clustered)')
    a.set_title('(a) order plane: local order varies, but every field is\n'
                'long-range super-Poissonian — none hyperuniform')
    b = ax[1]
    for name, _ in FIELDS:
        cu = curves[name]
        b.plot(cu['kc']/cu['kref'], cu['Sr'], color=COL[name], lw=1.6, label=SHORT[name])
    b.axhline(1, color='#888', ls='--', lw=1, label='Poisson S=1')
    b.set_xlim(0, 4); b.set_ylim(0, None)
    b.set_xlabel('k / k_ρ'); b.set_ylabel('S(k)')
    b.set_title('(b) radial structure factor: low-k rises above 1\n(heterogeneity); no S(k→0)→0 hyperuniform dip')
    b.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{outdir}/figures/fig_orderplane.pdf", bbox_inches="tight")
    fig.savefig(f"{outdir}/figures/fig_orderplane.png", dpi=110, bbox_inches="tight")

    # ---- Fig 3: number variance (long-range: not hyperuniform) ----
    fig, axn = plt.subplots(1, 1, figsize=(7.2, 5.2))
    for name, _ in FIELDS:
        cu = curves[name]; rr = cu['rr']; nv = cu['nv']; ok = (nv > 0)
        axn.plot(rr[ok]/1000, nv[ok], 'o-', color=COL[name], lw=1.4, ms=5,
                 label=f"{SHORT[name].split(' (')[0]} (NNI={cu['NNI']:.2f})")
    cu = curves["northwestern Calabrian margin"]; rr = cu['rr']/1000; nv = cu['nv']
    r0, v0 = rr[0], nv[0]; rline = np.array([rr[0], rr[-1]])
    axn.plot(rline, v0*(rline/r0)**2, 'k--', lw=1, label='slope 2 (area / Poisson)')
    axn.plot(rline, v0*(rline/r0)**1, 'k:', lw=1.2, label='slope 1 (perimeter / hyperuniform)')
    axn.set_xscale('log'); axn.set_yscale('log')
    axn.set_xlabel('window radius R (km)'); axn.set_ylabel(r'number variance $\sigma^2(R)$')
    axn.set_title('Number variance grows ~ area (slope ≈ 2–3) for every field,\n'
                  'not ~ perimeter (slope 1): no hyperuniform suppression')
    axn.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{outdir}/figures/fig_numbervariance.pdf", bbox_inches="tight")
    fig.savefig(f"{outdir}/figures/fig_numbervariance.png", dpi=110, bbox_inches="tight")
    print("figures written to", f"{outdir}/figures/")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shp", default="data/raw/Pockmarks_Med.shp")
    ap.add_argument("--dbf", default="data/raw/Pockmarks_Med.dbf")
    ap.add_argument("--outdir", default="..")
    a = ap.parse_args()
    run(a.shp, a.dbf, a.outdir)
