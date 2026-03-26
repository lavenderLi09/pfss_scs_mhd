# Analytic Bipolar 1-20 Rsun

This AMRVAC case is a clean analytical precursor for later PFSS/SCS-to-AMRVAC initialization. It reuses the reviewed `polytropic_bipolar` case structure, but keeps only a Parker solar-wind background plus the analytical charge-pair bipolar magnetic field.

Key setup choices:

- Geometry: `spherical_3D`
- Radial domain: `1.0` to `20.0 Rsun`
- Base mesh: `72 x 72 x 72`
- Radial stretching: geometric, `qstretch_baselevel = 1.005`
- MHD model: polytropic with `mhd_gamma = 1.05`
- Divergence control: constrained transport (`typedivbfix='ct'`)

The magnetic field is the analytical buried-charge bipolar field implemented in `bipolar_field()`. This case does not include the Titov flux rope, electric driving, or PFSS/SCS input files.

Parker background constants in `mod_usr.t`:

- `rc = 3.45d0`
- `Vs = 117.54d0`
- `rhob = 5.0d9/1000.0d0/unit_numberdensity`
- `Tiso = 2.0d0`

User-facing magnetic parameters in `amrvac.par`:

- `f_q`: bipolar field strength scale
- `f_d`: burial depth scale
- `f_L`: half-separation of the two buried charges

Current default geometry uses a more deeply buried charge pair:

- `f_d = 2.5`
- `f_L = 0.5`
