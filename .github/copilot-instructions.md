# PhaseFromInterferograms.jl — Copilot Context

## What this package is

`PhaseFromInterferograms` extracts phase from interferometric fringe patterns using FFT-based analysis, Phase-Shifted Interferometry (PSI), and Phase-Tilted Interferometry (PTI). It also handles tilt extraction and fringe carrier-frequency detection.

## Key types and functions

| Symbol | Purpose |
|---|---|
| `PTIestimate` | State struct for PTI phase reconstruction |
| `RoughTilts` / `FineTilts` | Tilt extraction algorithm variants |
| `diffPSI` / `LSPSI` | Differential PSI and least-squares PSI algorithms |
| `PSIAlg` / `TiltExtractionAlg` | Abstract algorithm types for dispatch |
| `findfirstharmonic` | Detect carrier frequency peak in FFT spectrum |
| `findfirstharmonic2` / `_v2` | Zoom-FFT variants for sub-pixel carrier detection |
| `eraseZerothOrder` / `eraseZerothOrder!` | Remove DC term from fringe spectrum |
| `get_tilt` / `fourier_tilt` | Extract tilt from a single interferogram |
| `get_phase_from_igrams_with_tilts` | Full pipeline: interferograms → phase + tilts |
| `get_phase_from_n_psi` | N-frame PSI phase extraction |
| `get_tilt_dirs` | Determine tilt directions from fringe orientations |
| `get_aperture` | Detect aperture from fringe image |

## Tilt types (re-exported from PhaseUtils)

`Tilt`, `TiltCentered`, `FreeTilt`, `sigma`/`tau` accessors — see `PhaseUtils` context.

## Coordinate descriptors (re-exported from PhaseUtils)

`ArrayAxes`, `FourierAxes`, `DataAxes`, `DataAxesCentered` — describe physical vs. Fourier-space grids.

## Zoom-FFT extension

`ZoomFFTVisualExt.jl` (in `ext/`) provides visualization for zoom-FFT intermediate steps; loaded when CairoMakie is available.

## Relationships

- Depends on: `PhaseUtils`, `FFTW`, `FFTViews`, `StatsBase`
- Used by: analysis scripts and notebooks in `Feedback14AMI`

## Source layout

```
src/
    PhaseFromInterferograms.jl  ← module entry, exports
    FindHarmonics.jl            ← carrier frequency detection (submodule)
    algorithms.jl               ← PSI, tilt extraction algorithms
    methods.jl                  ← get_tilt, fourier_tilt
    PTI.jl / PTI_tensor.jl      ← PTI phase reconstruction
    zoomFFTmodule.jl            ← zoom-FFT submodule
ext/
    ZoomFFTVisualExt.jl         ← visualization extension
```
