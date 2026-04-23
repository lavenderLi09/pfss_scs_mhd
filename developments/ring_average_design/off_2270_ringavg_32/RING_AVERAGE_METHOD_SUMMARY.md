# Ring Average Method Summary

## Goal

Ring Average is designed to remove the severe CFL timestep restriction caused by the clustering of azimuthal cells near the pole or axis in spherical/cylindrical finite-volume MHD.

In spherical 3D, the azimuthal cell size is approximately

\[
\Delta s_\phi \sim r \sin\theta \,\Delta\phi
\]

As \(\theta \to 0\) or \(\pi\), we have \(\sin\theta \to 0\), so the \(\phi\)-direction cells become extremely thin. For an explicit scheme, the CFL condition is roughly

\[
\Delta t \le C_{\mathrm{CFL}} \min\left(\frac{\Delta s}{\lambda_{\max}}\right)
\]

so the smallest \(\Delta s_\phi\) near the pole can dominate the global timestep.

Ring Average addresses this without changing the actual computational mesh. The original solver still advances the solution on the original grid, but after each full update the near-pole solution is post-processed ring by ring so that it behaves as if it lived on a coarser effective azimuthal grid.

## Core Idea

The method does not replace the original mesh. Instead:

1. Keep the original spherical finite-volume grid unchanged.
2. Group thin azimuthal cells in selected polar rings into larger chunks.
3. Average conservative quantities inside each chunk.
4. Reconstruct monotonically back to the original thin cells inside that chunk.
5. Recompute the CFL timestep based on the effective grid implied by the chunking.

So Ring Average is a conservative averaging-reconstruction post-processing algorithm, not just a timestep hack and not a mesh modification.

## Chunk Geometry

For a fixed \((r,\theta)\) ring, split the azimuthal direction into chunks, each containing \(N_C\) original cells.

If the ring has \(N_\phi\) cells total, then the number of chunks is

\[
N_{\mathrm{chunk}}=\frac{N_\phi}{N_C}
\]

Near the pole, \(N_C\) is larger. Away from the pole, \(N_C\) decreases and eventually Ring Average is turned off.

## Fluid Conservative Variables

For fluid variables, Ring Average acts on cell-centered conservative quantities such as

\[
\mathbf{U}=
\begin{pmatrix}
\rho \\
\rho v_r \\
\rho v_\theta \\
\rho v_\phi \\
E
\end{pmatrix}
\]

The operation is one-dimensional along \(\phi\) within each selected ring.

### 1. Chunk Averaging

For a conservative variable \(q\), the chunk mean is

\[
q_{\mathrm{Mean}}=\frac{1}{N_C}\sum_{k_l=1}^{N_C} q(k_l)
\]

where \(k_l=1,\dots,N_C\) is the local cell index within the chunk.

This removes the smallest azimuthal variations inside the chunk.

### 2. Reconstruction

The chunk average is then reconstructed back onto the original thin cells using a monotone reconstruction. The paper discusses:

- PCM: Piecewise Constant Method
- PLM: Piecewise Linear Method
- PPM: Piecewise Parabolic Method

The paper mainly uses PPM. In that case, a parabola is constructed inside each chunk:

\[
q(x)=Ax^2+Bx+C
\]

subject to:

- left interface value \(q_L\)
- right interface value \(q_R\)
- conservation of the chunk-integrated average \(q_{\mathrm{Mean}}\)

The coefficients are

\[
A = 3 \left(q_L-q_R-2q_{\mathrm{Mean}}\right)
\]

\[
B = 2 \left(3q_{\mathrm{Mean}}-q_R-2q_L\right)
\]

\[
C = q_L
\]

### 3. Recover Cell Averages in the Original Thin Cells

For the \(k_l\)-th original cell inside a chunk, the reconstructed cell average is

\[
q_r(k_l)=
\frac{A}{3N_C^2}\left(3k_l^2-3k_l+1\right)
+\frac{B}{2N_C}(2k_l-1)+C
\]

This preserves the total chunk amount:

\[
\sum_{k_l=1}^{N_C} q_r(k_l)=N_C\, q_{\mathrm{Mean}}
\]

so Ring Average is conservative on each chunk.

## PCM, PLM, and PPM in Ring Average

### PCM

Piecewise Constant Method sets every cell in the chunk to the chunk mean:

\[
q_r(k_l)=q_{\mathrm{Mean}}
\]

This is simplest and most robust, but also the most diffusive.

### PLM

Piecewise Linear Method reconstructs a linear profile:

\[
q(x)=ax+b
\]

with slope limiting to avoid oscillations. It preserves more intra-chunk structure than PCM.

### PPM

Piecewise Parabolic Method reconstructs a parabolic profile:

\[
q(x)=Ax^2+Bx+C
\]

This is higher-order and closest to the paper's main implementation.

## Why This Is Not Just a Standard Filter

Ring Average behaves roughly like a localized, non-linear azimuthal filter, but it is not just a standard boxcar or Fourier filter:

- it is only applied in selected polar rings
- it preserves chunk-integrated conservation exactly
- it uses monotone reconstruction
- it can preserve sharp structures better than a simple linear low-pass filter

The paper also gives a Fourier interpretation for PCM averaging:

\[
\tilde f_{\mathrm{avg}}(m)=
\tilde f(m)\,
\frac{1}{L}
\left[
\frac{\sin(L m \delta)}{\sin(m\delta)}
\right]^2
\]

which shows strong damping of high-\(m\) modes.

## Magnetic Field Treatment

The magnetic-field part depends on how the code stores magnetic variables.

### Cell-Centered Magnetic Field

If \(\mathbf{B}\) is cell-centered and the solver uses divergence cleaning, then the magnetic variables can be treated in the same way as fluid conservative variables. After Ring Average, a divergence-cleaning step is applied.

### Staggered / Face-Centered Magnetic Field

This is much harder, because it is necessary to preserve

\[
\nabla\cdot\mathbf{B}=0
\]

If all three magnetic-flux components are independently reconstructed, the local divergence-free condition will in general be violated.

## Paper's Staggered-Field Strategy

Let the curvilinear coordinate directions be \((\hat i,\hat j,\hat k)\), where:

- \(\hat i\): axis direction
- \(\hat j\): transverse direction
- \(\hat k\): azimuthal direction, the direction of Ring Average

The paper applies conservative averaging-reconstruction only to the first two magnetic-flux components:

\[
\Phi_i,\quad \Phi_j
\]

It does not reconstruct \(\Phi_k\) independently. Instead, \(\Phi_k\) is updated through Faraday's law using perturbation electric fields.

### Perturbation Electric Field From \(\Phi_i\)

After reconstructing \(\Phi_i^r(k_l)\), define

\[
\delta E_j(k_l+1)=
\sum_{k_l=1}^{N_C-1}
\frac{\Phi_i^r(k_l)-\Phi_i(k_l)}{\Delta t}
\]

### Perturbation Electric Field From \(\Phi_j\)

Similarly,

\[
\delta E_i(k_l+1)=
-
\sum_{k_l=1}^{N_C-1}
\frac{\Phi_j^r(k_l)-\Phi_j(k_l)}{\Delta t}
\]

### Explicit Choice

\[
\delta E_k=0
\]

Then \(\Phi_k\) is updated consistently through Faraday's law rather than reconstructed independently. This is the key step that preserves local and global

\[
\nabla\cdot\mathbf{B}=0
\]

for the staggered formulation.

## Effective CFL Scale

The larger timestep does not come from arbitrarily forcing a bigger \(\Delta t\). It comes from the fact that the near-pole solution is replaced by a chunk-averaged, reconstructed state that behaves like it lives on an effective coarser azimuthal grid.

The timestep is then recomputed using the effective grid:

\[
\Delta t = \min\left(\frac{\Delta L_{\mathrm{eff}}}{V_{\max}}\right)
\]

with

\[
\Delta L_{\mathrm{eff}}=
\frac{\Delta V}{\min(A_i^\ast,A_j^\ast,A_k)}
\]

where:

- \(\Delta V\) is the cell volume
- \(A_i^\ast, A_j^\ast\) are effective face areas on the chunked grid
- \(A_k\) is the remaining face area contribution

So in the full paper method, the larger timestep is a consequence of the averaging-reconstruction algorithm, not an independent shortcut.

## Choosing Chunk Size

The paper notes that chunk choice is not unique. A practical rule is to choose chunk sizes so that the effective azimuthal scale becomes comparable to the radial or polar scale.

For uniform angular spacing, the paper gives an effective azimuthal length estimate:

\[
\Delta L_{\mathrm{eff}}
=
\left(\frac{2m-1}{2}\right)\left(\frac{2\pi}{N}\right)\Delta R
\]

where:

- \(m\) is the ring index from the pole
- \(N\) is the total number of azimuthal cells
- \(\Delta R\) is the other local grid scale

They propose selecting chunk sizes \(2^k\) such that

\[
\left(\frac{2m-1}{2}\right)\pi\left(\frac{2^k}{N}\right)\approx 1
\]

which keeps the effective aspect ratio near unity.

## Full Algorithmic Sequence

For one timestep, the full Ring Average method is conceptually:

1. Advance the MHD solution with the original finite-volume solver on the original grid.
2. In selected polar rings, group azimuthal cells into chunks.
3. Average fluid conservative variables within each chunk.
4. Reconstruct them monotonically back to the original thin cells.
5. If magnetic field is cell-centered, apply the same idea to magnetic variables and then do divergence cleaning.
6. If magnetic field is staggered, reconstruct \(\Phi_i\) and \(\Phi_j\), construct perturbation electric fields, and update \(\Phi_k\) through Faraday's law.
7. Refresh boundaries and ghost cells.
8. Recompute the timestep using the effective grid.

So Ring Average is a post-processing algorithm inserted after the normal finite-volume update.

## Difference From a CFL-Only Approximation

This is the most important implementation distinction.

### CFL-Only Approximation

A simplified engineering approximation would only change the CFL estimate, for example replacing

\[
\frac{c_{\max}}{\Delta s_\phi}
\quad\rightarrow\quad
\frac{c_{\max}}{\Delta s_{\phi,\mathrm{eff}}}
\]

This only changes the allowed timestep. It does not change the actual near-pole updated state.

### Full Ring Average

The full paper method:

1. changes the near-pole numerical solution after each update by chunk averaging and reconstruction
2. then recomputes the timestep from the resulting effective grid

So the full method changes both:

- the allowed timestep
- the near-pole state itself

That is why a faithful AMRVAC implementation is much deeper than a simple `mod_dt.t` change.

## AMRVAC Interpretation

If Ring Average is implemented in MPI-AMRVAC, the clean conceptual decomposition is:

1. ring metadata:
   - identify which \((r,\theta)\) rings are inside the polar cap
   - choose \(N_C\) for each ring
   - build chunk grouping
2. fluid ring-average:
   - apply averaging-reconstruction to \(\rho\), momentum, and energy
3. magnetic ring-average:
   - cell-centered path: treat magnetic variables similarly, then do divergence cleaning
   - staggered path: use the paper's \(\Phi_i,\Phi_j,\delta E_i,\delta E_j,\Phi_k\) construction
4. post-processing hook:
   - insert after each full update stage
5. CFL update:
   - recompute timestep using the effective grid

## One-Sentence Summary

Ring Average is a conservative post-processing method that groups thin polar azimuthal cells into chunks, averages conservative quantities within each chunk, reconstructs them monotonically back to the original cells, preserves magnetic divergence with a special staggered-field treatment when needed, and then uses the resulting effective grid to permit a larger stable timestep.
