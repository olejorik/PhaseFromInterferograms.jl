
# # Introduction
#
# The Zoomed FFT (based on the chirp-z transform) is a powerful technique for achieving subpixel precision in frequency detection. Instead of zero-padding the entire signal (which is computationally expensive), the zoomed FFT computes the FFT only in a small region of interest around a target frequency.
#
# This tutorial demonstrates the `FFT2Zoom` implementation and the `findfirstharmonic2_v2` function for finding the dominant frequency in a signal with high precision.

using PhaseFromInterferograms
using PhaseFromInterferograms.zoomFFT2D: FFT2Zoom, findfirstharmonic2_v2
using FFTW
using CairoMakie

# # Part 1: Understanding Zero-Padding vs Zoomed FFT
#
# ## The Zero-Padding Principle
#
# Zero-padding in the spatial domain is equivalent to interpolation in the frequency domain. If we have a signal of size M×N and pad it to size pM×pN, the FFT will have p times more frequency samples.

## Create a test signal
M, N = 16, 20
signal = randn(ComplexF64, M, N)

## Method 1: Zero-pad and compute FFT
padfactor = 4
padM, padN = padfactor * M, padfactor * N
padded_signal = zeros(ComplexF64, padM, padN)
padded_signal[1:M, 1:N] .= signal
reference_fft = fft(padded_signal)

println("Original signal size: $(M)×$(N)")
println("Padded signal size: $(padM)×$(padN)")
println("FFT size: $(size(reference_fft))")

# ## Using FFT2Zoom
#
# Instead of computing the entire padded FFT, we can use `FFT2Zoom` to compute only the frequencies we need. The key is to match the frequency ranges correctly using `fftshift`.

## Method 2: Use FFT2Zoom to compute the interpolated FFT
Mset = 0:(M - 1)
Nset = 0:(N - 1)

## Create frequency ranges that match the padded FFT
freqs_x = fftshift(fftfreq(padM))
freqs_y = fftshift(fftfreq(padN))
Rset = freqs_x[1]:step(freqs_x):freqs_x[end]
Sset = freqs_y[1]:step(freqs_y):freqs_y[end]

fft_zoom = FFT2Zoom(Mset, Nset, Rset, Sset)
zoomed_fft = fft_zoom(signal)

## Compare: the zoomed FFT should match the fftshifted reference
max_error = maximum(abs.(zoomed_fft .- fftshift(reference_fft)))
println("Maximum error between zero-padded and zoomed FFT: $(max_error)")

# # Part 2: Subregion Extraction
#
# The real power of zoomed FFT is that we can compute just a small region around a frequency of interest.

## Extract a quarter region from the center
center_M, center_N = padM ÷ 2, padN ÷ 2
halfM, halfN = padM ÷ 8, padN ÷ 8
M_range = (center_M - halfM + 1):(center_M + halfM)
N_range = (center_N - halfN + 1):(center_N + halfN)

## Reference: extract from zero-padded FFT
ref_region = fftshift(reference_fft)[M_range, N_range]

## Zoomed FFT: compute only the region of interest
Rset_region = freqs_x[M_range]
Sset_region = freqs_y[N_range]

fft_zoom_region = FFT2Zoom(Mset, Nset, Rset_region, Sset_region)
zoomed_region = fft_zoom_region(signal)

## Compare
max_error_region = maximum(abs.(zoomed_region .- ref_region))
println("Maximum error for subregion: $(max_error_region)")
println(
    "Computed only $(length(M_range))×$(length(N_range)) instead of $(padM)×$(padN) samples!",
)

# # Part 3: Finding Tilts with Subpixel Precision
#
# ## Creating a Test Signal
#
# Let's create a signal with a known linear phase (tilt):

arrsize = (25, 31)  # Use odd dimensions
f1, f2 = 0.18, -0.15  # Frequency components
offset = π / 5

tilt = fourier_tilt(2π .* (f1, f2), offset, arrsize)
sig = cis.(tilt)

## Visualize the phase
fig = Figure(; size=(800, 400))
ax1 = Axis(fig[1, 1]; title="Real part", aspect=DataAspect())
ax2 = Axis(fig[1, 2]; title="Phase", aspect=DataAspect())

heatmap!(ax1, real.(sig))
heatmap!(ax2, angle.(sig))

fig

# ## Coarse Frequency Detection
#
# First, let's see what a standard FFT gives us:

spectrum = fft(sig)
abs_spectrum = abs.(spectrum)

## Find peak
peak_idx = argmax(abs_spectrum)
freqs_x_coarse = fftfreq(arrsize[1])
freqs_y_coarse = fftfreq(arrsize[2])
f1_coarse = freqs_x_coarse[peak_idx[1]]
f2_coarse = freqs_y_coarse[peak_idx[2]]

println("True frequencies: f1 = $(f1), f2 = $(f2)")
println("Coarse detection: f1 = $(f1_coarse), f2 = $(f2_coarse)")
println("Errors: Δf1 = $(abs(f1 - f1_coarse)), Δf2 = $(abs(f2 - f2_coarse))")

# ## Zoomed FFT Refinement
#
# Now let's use zoomed FFT to refine the estimate:

zoom_factor = 32
padsize = zoom_factor .* arrsize

## Create zoomed frequency grid
freqs_x_zoom = fftshift(fftfreq(padsize[1]))
freqs_y_zoom = fftshift(fftfreq(padsize[2]))

## Find indices closest to our coarse estimate
idx_f1 = argmin(abs.(freqs_x_zoom .- f1_coarse))
idx_f2 = argmin(abs.(freqs_y_zoom .- f2_coarse))

## Extract a region around the peak
half_width = 16
M_range = max(1, idx_f1 - half_width):min(padsize[1], idx_f1 + half_width)
N_range = max(1, idx_f2 - half_width):min(padsize[2], idx_f2 + half_width)

## Compute zoomed FFT
Mset = 0:(arrsize[1] - 1)
Nset = 0:(arrsize[2] - 1)
Rset = freqs_x_zoom[M_range]
Sset = freqs_y_zoom[N_range]

fft_zoom = FFT2Zoom(Mset, Nset, Rset, Sset)
zoomed_spectrum = fft_zoom(sig)

## Find peak in zoomed region
zoom_peak_idx = argmax(abs.(zoomed_spectrum))
f1_refined = Rset[zoom_peak_idx[1]]
f2_refined = Sset[zoom_peak_idx[2]]
phase_refined = angle(zoomed_spectrum[zoom_peak_idx])

println("\nRefined detection (zoom=$(zoom_factor)):")
println("f1 = $(f1_refined), f2 = $(f2_refined)")
println("Errors: Δf1 = $(abs(f1 - f1_refined)), Δf2 = $(abs(f2 - f2_refined))")

# Visualize the zoomed spectrum
fig = Figure(; size=(800, 400))
ax1 = Axis(fig[1, 1]; title="Coarse spectrum (full)", aspect=DataAspect())
ax2 = Axis(fig[1, 2]; title="Zoomed spectrum ($(zoom_factor)x)", aspect=DataAspect())

heatmap!(ax1, fftshift(abs_spectrum))
heatmap!(ax2, abs.(zoomed_spectrum))

fig


# # Part 4: Reconstruction and Error Analysis
#
# Let's reconstruct the tilt from the detected parameters and measure the error:

reconstructed_tilt = fourier_tilt(2π .* (f1_refined, f2_refined), phase_refined, arrsize)

## Calculate RMS error
phase_diff = angle.(exp.(im .* (tilt .- reconstructed_tilt)))
rms_error = sqrt(sum(phase_diff .^ 2) / length(phase_diff))

println("\nReconstruction RMS error: $(rms_error) radians")

## Visualize
fig = Figure(; size=(1200, 400))
ax1 = Axis(fig[1, 1]; title="Original tilt", aspect=DataAspect())
ax2 = Axis(fig[1, 2]; title="Reconstructed tilt", aspect=DataAspect())
ax3 = Axis(fig[1, 3]; title="Error", aspect=DataAspect())

heatmap!(ax1, tilt; colormap=:twilight)
heatmap!(ax2, reconstructed_tilt; colormap=:twilight)
hm = heatmap!(ax3, phase_diff; colormap=:RdBu)
Colorbar(fig[1, 4], hm; label="Phase error (rad)")

fig

# # Part 5: Automatic Zoom Levels with findfirstharmonic2_v2
#
# The `findfirstharmonic2_v2` function automates the iterative refinement process:

## Test with increasing zoom levels
zoom_factors = [4, 8, 16, 32]
errors = Float64[]

for zf in zoom_factors
    (fhat, phase_est), _, _ = findfirstharmonic2_v2(sig; zoomlevels=[1, zf], erasesize=0)

    ## Reconstruct and calculate error
    reconstructed = fourier_tilt(2π .* (fhat[1], fhat[2]), phase_est, arrsize)
    phase_diff = angle.(exp.(im .* (tilt .- reconstructed)))
    rms = sqrt(sum(phase_diff .^ 2) / length(phase_diff))

    push!(errors, rms)
    println("Zoom factor $(zf): RMS error = $(round(rms, digits=6)) rad")
end

## Plot convergence
fig = Figure(; size=(600, 400))
ax = Axis(
    fig[1, 1];
    xlabel="Zoom factor",
    ylabel="RMS error (radians)",
    title="Reconstruction accuracy vs zoom level",
    xscale=log10,
    yscale=log10,
)

scatterlines!(ax, zoom_factors, errors; marker=:circle, markersize=12)
fig

# # Part 6: Handling Real Signals with DC Components
#
# Real non-negative signals (like squared interferogram differences) have strong DC components that must be excluded.

## Create a realistic interferogram-like signal
slow_x = range(0, 2π; length=arrsize[1])
slow_y = range(0, 2π; length=arrsize[2])
slow_background = 1.0 .+ 0.3 .* cos.(slow_x) .+ 0.3 .* cos.(slow_y)'

contrast = 0.8
real_signal = slow_background .* (1.0 .+ contrast .* cos.(tilt))

## Visualize
fig = Figure(; size=(800, 400))
ax1 = Axis(fig[1, 1]; title="Real signal", aspect=DataAspect())
ax2 = Axis(fig[1, 2]; title="FFT magnitude (log)", aspect=DataAspect())

heatmap!(ax1, real_signal)

## Show FFT with DC
spec_real = fft(real_signal)
heatmap!(ax2, log10.(1 .+ abs.(fftshift(spec_real))))

fig

# ## Finding the Sidelobe (not DC)
#
# With `erasesize > 0`, the function removes DC before searching:

## Without DC removal (finds DC)
(fhat_dc, _), _, _ = findfirstharmonic2_v2(real_signal; zoomlevels=[1], erasesize=0)
println("Without DC removal: f1 = $(fhat_dc[1]), f2 = $(fhat_dc[2])")

## With DC removal (finds sidelobe)
(fhat_sidelobe, phase_est), _, _ = findfirstharmonic2_v2(
    real_signal; zoomlevels=[1, 8], erasesize=5
)
println("With DC removal: f1 = $(fhat_sidelobe[1]), f2 = $(fhat_sidelobe[2])")
println("True frequencies: f1 = $(f1), f2 = $(f2)")

# # Summary
#
# The zoomed FFT technique provides:
#
# 1. **Computational efficiency**: Only computes frequencies in regions of interest
# 2. **Subpixel precision**: Achieves much finer frequency resolution than standard FFT
# 3. **Iterative refinement**: Can progressively increase zoom for better accuracy
# 4. **DC handling**: Can exclude DC regions for real non-negative signals
#
# Key parameters:
#
# - `zoomlevels`: Array of zoom factors to apply sequentially (default: `[1, 2, 4, 8, 16, ...]`)
# - `half_width`: Half-width of search window around peak (default: 16 samples)
# - `erasesize`: Size of DC region to exclude (default: 5, use 0 for complex signals)
#
# The `findfirstharmonic2_v2` function automates this process, making it easy to find dominant frequencies with high precision in both complex and real signals.
