"""
analyze_pockmarks.py
--------------------
Hyperuniformity-class characterization of single pockmark fields from the
Spatola et al. (2025) Mediterranean geodatabase (Zenodo 10.5281/zenodo.15425907).

For each selected single field we compute, in its projected (UTM, metre) frame:
  * Clark-Evans nearest-neighbour index (NNI) and its z-score  -- local order;
  * pair correlation g(r)                                      -- local exclusion;
  * radially-averaged structure factor S(k)                   -- long-range order;
  * number variance sigma^2(R) and its log-log slope          -- HU vs clustered;
and compare each against three nulls (CSR, RSA hard-core matched to the field's
nearest-neighbour scale, and an inhomogeneous-Poisson heterogeneity null), with a
disordered-hyperuniform perturbed-lattice reference for context.

The raw geodatabase (shapefile) is third-party data (Spatola et al. 2025,
CC-BY-4.0); download it from Zenodo and pass --shp / --dbf. Derived per-field
coordinates and all result tables are written to ../data and ../results.

Usage:
    python analyze_pockmarks.py --shp data/raw/Pockmarks_Med.shp \
                                --dbf data/raw/Pockmarks_Med.dbf --outdir ..
"""
import argparse, csv, os
import numpy as np
from scipy.spatial import cKDTree, ConvexHull
from shp_reader import read_points, read_dbf_column
from ppstats2d import (structure_factor_rect, number_variance_rect,
                       pair_correlation, inhomogeneous_poisson, rsa_hardcore)

# the three representative single fields (Location_1 values in the geodatabase)
FIELDS = [
    ("Offshore Zannone Island",      "isotropic / near-random"),
    ("northwestern Calabrian margin", "clustered (multi-patch)"),
    ("Malta Plateau",                 "regular but anisotropic"),
]
SEED = 0


def nni(pts):
    n = len(pts); A = ConvexHull(pts).volume; rho = n / A
    dnn = cKDTree(pts).query(pts, k=2)[0][:, 1]
    exp = 1.0 / (2*np.sqrt(rho))
    R = dnn.mean() / exp
    z = (dnn.mean() - exp) / (0.26136 / np.sqrt(n*rho))
    return R, z, np.median(dnn), rho


def long_range(pts, box, rho):
    xmin, xmax, ymin, ymax = box; Lx, Ly = xmax-xmin, ymax-ymin
    kref = 2*np.pi*np.sqrt(rho)
    kc, Sr = structure_factor_rect(pts, Lx, Ly, kmax=4*kref, nk=30)
    sk0 = np.nanmean(Sr[kc < 0.4*kref])
    Rmax = 0.22*min(Lx, Ly); radii = np.linspace(0.15*Rmax, Rmax, 7)
    rng = np.random.default_rng(SEED)
    rr, nv = number_variance_rect(pts, box, radii, 500, rng)
    ok = ~np.isnan(nv) & (nv > 0)
    slope = np.polyfit(np.log(rr[ok]), np.log(nv[ok]), 1)[0] if ok.sum() >= 3 else np.nan
    return kc, Sr, sk0, rr, nv, slope, kref


def make_rsa(n, box, dmin, rng, maxtry=400000):
    xmin, xmax, ymin, ymax = box; pts = []; t = 0
    while len(pts) < n and t < maxtry:
        p = np.array([rng.uniform(xmin, xmax), rng.uniform(ymin, ymax)]); t += 1
        if not pts or np.min(np.linalg.norm(np.array(pts)-p, axis=1)) >= dmin:
            pts.append(p)
    return np.array(pts)


def main(shp, dbf, outdir):
    st, P = read_points(shp)
    loc = np.array(read_dbf_column(dbf, "Location_1")[:len(P)])
    rng = np.random.default_rng(SEED)

    os.makedirs(f"{outdir}/data", exist_ok=True)
    os.makedirs(f"{outdir}/results", exist_ok=True)
    os.makedirs(f"{outdir}/figures", exist_ok=True)

    summary, nulls, coords = [], [], []
    curves = {}
    for name, role in FIELDS:
        pts = np.unique(P[loc == name], axis=0).astype(float)
        c = pts.mean(0); pts = pts - c
        box = (pts[:,0].min(), pts[:,0].max(), pts[:,1].min(), pts[:,1].max())
        Lx, Ly = box[1]-box[0], box[3]-box[2]
        R, z, d0, rho = nni(pts)
        kc, Sr, sk0, rr, nv, slope, kref = long_range(pts, box, rho)
        rcent, g = pair_correlation(pts, max(Lx, Ly), rmax=6*d0, dr=0.4*d0)
        summary.append(dict(field=name, role=role, N=len(pts),
                            extent_km=round(max(Lx, Ly)/1000, 1),
                            aspect=round(max(Lx, Ly)/min(Lx, Ly), 2),
                            NNI=round(R, 3), z=round(z, 1), d0_m=int(d0),
                            Sk0=round(sk0, 2), NVslope=round(slope, 2)))
        curves[name] = dict(kc=kc, Sr=Sr, kref=kref, rcent=rcent, g=g, rr=rr, nv=nv,
                            NNI=R, sk0=sk0, pts=pts)
        for x, y in pts:
            coords.append([name, round(x, 1), round(y, 1)])

        # nulls
        def avg_null(gen, m):
            a = []
            for _ in range(m):
                q = gen()
                if len(q) < 8: continue
                Rn, zn, _, rhon = nni(q)
                _, _, sk0n, _, _, sln, _ = long_range(q, box, rhon)
                a.append((Rn, sk0n, sln))
            return np.nanmean(a, axis=0)
        csr = avg_null(lambda: np.column_stack([rng.uniform(box[0], box[1], len(pts)),
                                                rng.uniform(box[2], box[3], len(pts))]), 8)
        rsa = avg_null(lambda: make_rsa(len(pts), box, 0.85*d0, rng), 6)
        ihp = avg_null(lambda: inhomogeneous_poisson(pts, box, 6, rng), 8)
        for lab, vals in [("DATA", (R, sk0, slope)), ("CSR", csr),
                          ("RSA_hardcore", rsa), ("inhomogeneous_Poisson", ihp)]:
            nulls.append(dict(field=name, model=lab, NNI=round(vals[0], 3),
                              Sk0=round(vals[1], 2), NVslope=round(vals[2], 2)))

    # write tables
    with open(f"{outdir}/results/field_summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary[0].keys())); w.writeheader(); w.writerows(summary)
    with open(f"{outdir}/results/null_comparison.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(nulls[0].keys())); w.writeheader(); w.writerows(nulls)
    with open(f"{outdir}/data/selected_field_coordinates.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["field", "x_utm33n_m_centered", "y_utm33n_m_centered"])
        w.writerows(coords)
    for name in curves:
        cu = curves[name]; tag = name.split()[0].lower()
        with open(f"{outdir}/results/Sk_{tag}.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["k_over_krho", "S_k"])
            for kk, ss in zip(cu['kc']/cu['kref'], cu['Sr']):
                if not np.isnan(ss): w.writerow([round(kk, 4), round(ss, 4)])
        with open(f"{outdir}/results/gr_{tag}.csv", "w", newline="") as f:
            w = csv.writer(f); w.writerow(["r_m", "g_r"])
            for rr2, gg in zip(cu['rcent'], cu['g']):
                w.writerow([round(float(rr2), 1), round(float(gg), 4)])
    print("=== field_summary ===")
    for s in summary: print(s)
    print("\n=== null_comparison ===")
    for nrow in nulls: print(nrow)
    return curves


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shp", default="data/raw/Pockmarks_Med.shp")
    ap.add_argument("--dbf", default="data/raw/Pockmarks_Med.dbf")
    ap.add_argument("--outdir", default="..")
    a = ap.parse_args()
    main(a.shp, a.dbf, a.outdir)
