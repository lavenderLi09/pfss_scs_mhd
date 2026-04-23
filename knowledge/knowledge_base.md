# Knowledge Base

## SCS can be computed from spherical harmonics

Yes. The Schatten Current Sheet (SCS) field can be evaluated from spherical-harmonic expansion coefficients.

In the literature, these coefficients are commonly written as `g'_{nm}`, `h'_{nm}` or `g_l^m`, `h_l^m`. In code, they may also appear as `glm` and `hlm`. These are the same kind of objects: harmonic coefficients used to reconstruct the magnetic field in the SCS region.

### Main takeaway

- PFSS and SCS are both treated as potential, divergence-free magnetic fields in spherical geometry.
- The SCS field outside the cusp / current-sheet interface can be written directly as a spherical-harmonic expansion.
- Therefore, once the SCS coefficients are known, the field can be evaluated directly at arbitrary `(r, theta, phi)` points without requiring an intermediate storage grid in principle.

### Important practical note

In many SCS workflows, the field used to derive the SCS coefficients is first reoriented so that the radial component points outward everywhere at the matching surface. After the harmonic coefficients are obtained, a polarity-restoration step may still be needed when reconstructing the physical field.

### Why this matters for this project

For the PFSS + SCS -> AMRVAC workflow, this supports the idea that SCS does not have to be treated only as an interpolated grid product. Instead, it is reasonable to design the workflow around direct evaluation from harmonic coefficients on the AMRVAC target mesh.

### References

1. Knizhnik, K. (2024), review of PFSS and SCS:
   [Frontiers article](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2024.1476498/full)

2. Nikolic, L. (2017), explicit SCS harmonic formulas:
   [Government of Canada PDF](https://publications.gc.ca/collections/collection_2017/rncan-nrcan/M183-2/M183-2-8007-eng.pdf)

## How global coronal MHD models choose domain and grid

For global coronal simulations extending to roughly `10-25 Rsun`, the literature shows a strong common pattern:

- the radial domain is usually stretched rather than uniform;
- the outer boundary is often placed near `20-25 Rsun` when the model is meant to feed a heliospheric solution;
- the angular mesh is usually designed to avoid polar singularity or excessive pole clustering;
- many modern models do not use a simple uniform `(r, theta, phi)` spherical grid.

### Feng / SIP-CESE style global corona

A 2023 paper by Feng et al. gives a very explicit example of how the radial mesh is designed.

- Domain: `1 Rsun` to `6.7 Rsun`
- Grid type: `Yin-Yang` overlapping spherical grid
- Resolution: for each component, `Nr x Ntheta x Nphi = 85 x 60 x 180`
- Angular spacing: uniform, `Δtheta = Δphi = 1.5 deg`
- Radial spacing: nonuniform / stretched

Their radial spacing is piecewise:

- `Δr = 0.01 Rsun` for `r < 1.1 Rsun`
- a transition rule for `1.1 <= r < 3.5 Rsun`
- `Δr = r * Δtheta` for `r >= 3.5 Rsun`

This is a useful reference because it shows a concrete and physically motivated stretched radial grid, with finer spacing near the lower corona and larger cells farther out.

Reference:
[Feng et al. 2023, MNRAS](https://academic.oup.com/mnras/article/519/4/6297/6967137)

### COCONUT

COCONUT is a strong reference for a modern global coronal model extending to a heliosphere-coupling radius.

- Domain: `1.01 Rsun` to `25 Rsun`
- Grid type: unstructured subdivided geodesic mesh
- Radial grid: stretched
- Radial layers: `73`
- Cells per radial layer: `20,480`
- Total cells: `1,495,040`

This means COCONUT avoids the latitude-longitude pole problem entirely and uses gradual radial stretching instead of a uniform radial mesh.

Related coupling work also uses `21.5 Rsun` as the corona-to-heliosphere handoff radius, motivated by the expectation that the solar wind is already supersonic and super-Alfvenic there.

References:
[Wang et al. 2025, A&A](https://www.aanda.org/articles/aa/full_html/2025/02/aa52279-24/aa52279-24.html)
[COCONUT + EUHFORIA 2025, A&A](https://www.aanda.org/articles/aa/full_html/2025/01/aa51854-24/aa51854-24.html)

### AWSoM / AWSoM-R

AWSoM is an important reference because it represents the mainstream “global corona plus heliosphere coupling” approach.

- Coronal domain: typically `1.0-1.1 Rsun` to `24 Rsun`
- Grid type: spherical block-adaptive grid
- Radial grid: stretched
- Angular refinement: adaptive, with stronger refinement near important structures such as the current sheet

One concrete AWSoM example reports:

- Domain: `1 Rsun` to `24 Rsun`
- Total cells: `4,624,896`
- Finest angular resolution near the inner boundary: about `1.4 deg`
- Mean radial spacing grows strongly with radius:
- about `1e-4 Rsun` at the bottom boundary
- about `0.02 Rsun` at `1.15 Rsun`
- about `0.86 Rsun` at `20 Rsun`

This is another strong indication that a stretched radial grid is standard practice for coronal MHD domains extending to large radius.

References:
[CCMC AWSoM-R model page](https://ccmc.gsfc.nasa.gov/models/SWMF~AWSoM_R~1.0)
[Schad et al. 2024, Science Advances / PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11421591/)

### Cross-model takeaway

Across Feng-style models, COCONUT, and AWSoM, the common practical choices are:

- Use a lower boundary near `1.0-1.1 Rsun`
- Use an outer boundary near `20-25 Rsun` when later heliospheric coupling or extended coronal relaxation is needed
- Use a stretched radial mesh
- Avoid a naive uniform spherical latitude-longitude grid for the full domain

### What this suggests for the current project

For a PFSS + SCS based AMRVAC initial-condition workflow, the literature supports the following high-level direction:

- a radial domain extending to about `20 Rsun` is physically and numerically reasonable;
- a stretched radial mesh is the normal choice, not a uniform radial grid;
- angular discretization should be chosen with pole behavior in mind;
- the final target mesh should be designed around the MHD solver requirements, then PFSS/SCS should be evaluated on that mesh rather than forcing the solver to inherit an intermediate PFSS/SCS storage grid.

## What is commonly used for SCS `lmax`

The literature is much less explicit about SCS truncation order than it is about PFSS truncation order. Many papers use PFSS + SCS operationally, but do not clearly document the exact SCS harmonic cutoff.

### What can be stated with confidence

- `lmax_scs = 10` is a real and commonly used operational baseline.
- It should not be treated as a universal standard, because many papers do not publish the exact value.
- For interface matching problems, `lmax_scs = 10` may be too low even if it is common.

### Explicit example found in the literature

One modern example is Narechania et al. (2021), which states:

- angular grid on the fitting sphere: `181 x 360`
- SCS maximum degree: `N_s = 10`

This corresponds to using **SCS `lmax = 10`** in that implementation.

Reference:
[Narechania et al. 2021, SWSC PDF](https://www.swsc-journal.org/articles/swsc/pdf/2021/01/swsc200078.pdf)

### Practical implication for this project

For this project, the safest interpretation is:

- start from the fact that `lmax_scs = 10` is not unusual;
- do not assume it is automatically sufficient for smooth PFSS-SCS matching;
- if interface mismatch remains significant, test `lmax_scs = 15` or `20` before changing too many other things at once.

### Supporting context

The 2024 SCS review also notes that many publications use SCS without fully specifying implementation details, which is one reason a single universally documented “standard SCS `lmax`” is hard to identify.

Reference:
[Knizhnik 2024 review](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2024.1476498/full)
