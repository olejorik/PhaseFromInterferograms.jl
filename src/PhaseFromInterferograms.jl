module PhaseFromInterferograms
using FFTW
using FFTViews
using PhaseUtils: phwrap, dotproduct
using StatsBase

# Overview of the exported symbols
include("utils.jl")
export get_aperture

include("PTI_tensor.jl")

include("Windowing.jl")
using .Windowing
export GaussianWindow

include("FindHarmonics.jl")
using .FindHarmonics
export eraseZerothOrder, eraseZerothOrder!, findfirstharmonic

include("zoomFFTmodule.jl")
using .zoomFFT2D
export findfirstharmonic2

include("methods.jl")
export get_tilt, fourier_tilt

include("algorithms.jl")
export RoughTilts, FineTilts, diffPSI, LSPSI, PSIAlg, TiltExtractionAlg
export get_phase_from_igrams_with_tilts, get_tilt_dirs, get_phase_from_n_psi

end
