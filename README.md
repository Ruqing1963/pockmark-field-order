# pockmark-field-order

**Pockmark fields are locally ordered but not hyperuniform: a structure-factor analysis of Mediterranean seabed fluid-escape depressions**

Ruqing Chen — GUT Geoservice Inc., Montréal, Québec, Canada

Code, derived data and figures for a hyperuniformity-class analysis of submarine
pockmark fields, using single fields from the open Mediterranean pockmark
geodatabase of Spatola et al. (2025).

---

## Summary

Pockmark fields have been described both as fault-controlled clusters and as
regularly spaced by lateral inhibition. We test a sharper question: are they
**hyperuniform** (long-wavelength density fluctuations suppressed,
S(k→0)→0, as in avian photoreceptor mosaics), or is any regularity purely local?

For three single fields spanning the range of local order we combine the
Clark–Evans nearest-neighbour index (local) with the radial structure factor
S(k) and number variance σ²(R) (long-range), each against CSR, hard-core (RSA)
and inhomogeneous-Poisson nulls.

- **Local order varies and is real:** NW Calabrian margin is genuinely clustered
  (NNI 0.47, below its heterogeneity null); Malta Plateau is genuinely regular
  (NNI 1.26, above a hard-core gas); Offshore Zannone is ≈ random (1.08).
- **No field is hyperuniform:** every field has S(k→0) > 1 and number variance
  growing with window **area** (slope ≈ 2–3), not perimeter; each field's
  long-range fluctuation matches its own **inhomogeneous-Poisson** null, so the
  excess is density heterogeneity, not intrinsic order.
- The one regular field is an **anisotropic lineament**, not an isotropic field.

**Verdict:** pockmark order is local and steric/structural, heterogeneity-limited
at long range — *not* hyperuniform. Nearest-neighbour regularity must not be read
as long-range order without a structure-factor test.

| Field | N | NNI | S(k→0) | NV slope | character |
|---|---:|---:|---:|---:|---|
| Offshore Zannone Island | 112 | 1.08 | 2.5 | 2.8 | isotropic ~random |
| NW Calabrian margin | 418 | 0.47 | 19.4 | 2.9 | clustered (multi-patch) |
| Malta Plateau | 184 | 1.26 | 5.4 | 3.3 | regular but anisotropic |

(Reference: CSR S(k→0)=1, slope 2; hyperuniform S(k→0)→0, slope 1.)

---

## Repository layout

```
pockmark-field-order/
├── README.md
├── LICENSE                       # CC BY 4.0 (raw geodatabase NOT included)
├── requirements.txt              # numpy, scipy, matplotlib
├── code/
│   ├── shp_reader.py             # pure-Python ESRI shapefile (.shp) + .dbf reader
│   ├── ppstats2d.py              # NNI, g(r), S(k), number variance, nulls
│   ├── analyze_pockmarks.py      # field selection + diagnostics + null tables
│   └── make_figures.py           # the three paper figures
├── data/
│   ├── selected_field_coordinates.csv   # DERIVED: the 3 analysed fields' XY (Spatola 2025, CC-BY-4.0)
│   └── raw/README.md             # how to download the full geodatabase
├── results/
│   ├── field_summary.csv         # Table 1
│   ├── null_comparison.csv       # Table 2
│   ├── Sk_*.csv                  # radial structure factor per field
│   └── gr_*.csv                  # pair correlation per field
├── figures/
│   ├── fig_fields.{pdf,png}              # Figure 1 (point patterns)
│   ├── fig_orderplane.{pdf,png}          # Figure 2 (order plane + S(k))
│   └── fig_numbervariance.{pdf,png}      # Figure 3 (number variance)
└── paper/
    ├── pockmark_paper.pdf
    ├── pockmark_paper.tex
    └── (figure pdfs)
```

---

## Data source (not redistributed in full)

D. Spatola, M. Rovere, D. Casalbore & F. L. Chiocci (2025).
*Pockmarks of the Mediterranean region seas: a comprehensive geodatabase for
marine geomorphological analysis.* **Scientific Data 12, 1049.**
Geodatabase on Zenodo: **https://doi.org/10.5281/zenodo.15425907** (CC BY 4.0;
7,516 pockmarks; ESRI point shapefile, WGS84 UTM 33N).

Download the shapefile into `data/raw/` (see `data/raw/README.md`). Only the
derived coordinates of the three analysed fields are included in this repository.

---

## Reproducing the analysis

```bash
pip install -r requirements.txt
# download the geodatabase shapefile into data/raw/ (see data/raw/README.md)
cd code
python analyze_pockmarks.py --shp ../data/raw/Pockmarks_Med.shp \
                            --dbf ../data/raw/Pockmarks_Med.dbf --outdir ..
python make_figures.py      --shp ../data/raw/Pockmarks_Med.shp \
                            --dbf ../data/raw/Pockmarks_Med.dbf --outdir ..
```

This regenerates `results/`, the derived coordinates in `data/`, and the figures.

---

## Method notes

- **Local vs long-range.** Hyperuniformity is a long-range property: it cannot be
  read from nearest-neighbour spacing. A hard-core liquid or a lateral-inhibition
  field has NNI > 1 yet is *not* hyperuniform. The structure factor S(k→0) and
  the number-variance slope are the decisive long-range diagnostics.
- **Compilation → Poisson.** The geodatabase pools many independent fields;
  analysing it as a whole would drive any statistic toward Poisson. We analyse
  single, contiguous, artefact-free fields only.
- **Heterogeneity nulls.** Each field is compared to an inhomogeneous-Poisson
  null reproducing its coarse density variation; matching it means the long-range
  fluctuation is heterogeneity, not intrinsic correlation.

---

## Citation

Please cite the paper (`paper/pockmark_paper.pdf`) and this repository (Zenodo DOI
on release), and the data source, Spatola et al. (2025),
https://doi.org/10.5281/zenodo.15425907.

## License

Code, derived data, figures and paper: **CC BY 4.0** (see `LICENSE`). The raw
geodatabase is © Spatola et al. (2025), CC BY 4.0, obtained separately.
