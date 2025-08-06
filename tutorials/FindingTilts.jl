# # Finding Tilts from the Autocorrelation Function
#
# ## Introduction
# This tutorial demonstrates how to detect and reconstruct linear phase components (tilts) in signals using Fourier analysis and autocorrelation. The algorithm uses the maxima of the autocorrelation spectra to deduce the tilt parameters. For this purpose, a simple function [`fourier_tilt`](@ref) is provided.
#
# ## 1. Creating a Basic Tilt
#
# We first construct a tilt (a linear function) with periodic frequencies. This section demonstrates how a tilt is represented and visualized.

using CairoMakie
using PhasePlots, PhaseUtils

using FFTW, FFTViews
using PhaseFromInterferograms
using PhaseFromInterferograms: fourier_tilt, getslopes
using LaTeXStrings
CairoMakie.activate!()

arrsize = (50, 35)
i, j = 3, 8
f1, f2 = getindex.(fftfreq.(arrsize), [i, j])
offset = 2
## Generate a linear tilt with specified frequencies and offset
## The tilt is a linear function: t(x) = 2π(f₁x₁ + f₂x₂) + offset
## This is useful for simulating phase ramps in interferometric data

tilt = fourier_tilt(2π .* (f1, f2), offset, arrsize)
## Figure 1: Linear function with periodic frequencies
showarray(
    tilt;
    axis=(; title=L"Linear function $t(x)$ with frequencies $f_1=3/50$, $f_2=8/35$",),
    rot=0,
)

# ## 2. Visualizing the Wrapped Phase
#
# The exponentiation creates a linear phase term. The wrapped phase is visualized below.

sig = cis.(tilt)
## Figure 2: Wrapped phase pattern
showarray(angle.(sig); axis=(; title=L"Wrapped phase $\phi(x)=w(t(x))$",), rot=0)

# ## 3. Fourier Spectrum of the Tilt
#
# The Fourier spectrum of the complex exponential of the tilt shows a single peak at the tilt frequency.

spec = fft(sig)
## Figure 3: Power spectrum with single peak at tilt frequency
showarray(
    abs.(spec);
    axis=(; title="Power spectrum showing single peak at tilt frequency (3,8)"),
    rot=0,
)

## Highlight the peak location in the spectrum
scatter!(i, j; color=:red)
current_figure()

## The argument of the peak gives the original offset
angle(spec[i, j]) ≈ offset

# ## 4. Introducing Subpixel Coordinates
#
# In real-world applications, the signal may have subpixel shifts. This section explores the impact of subpixel coordinates on tilt reconstruction.

subpixelshifts = [0.1, 0.52]
## subpixelshifts = rand(2) # or use rand(2)
f1s, f2s = [f1, f2] .+ subpixelshifts ./ arrsize

# Construct the tilt with subpixel frequencies.
tilts = fourier_tilt(2π .* (f1s, f2s), offset, arrsize)
## Figure 4: Linear function with subpixel frequencies
showarray(tilts; axis=(title=L"Linear function $t(x)$ ",), rot=0)

# It doesn't produce a periodic signal anymore!
sigs = cis.(tilts)
## Figure 5: Wrapped phase of the signal with subpixel frequencies
showarray(angle.(sigs); axis=(title=L"Linear function $2\pi t(x)$ , wrapped",), rot=0)

# And that's why the Fourier spectrum demonstrates aliasing
specs = fft(sigs)
## Figure 6: Power spectrum of the complex exponent of linear function with subpixel frequencies
showarray(
    abs.(specs);
    axis=(title="Power spectrum of the complex exponent of linear function",),
    rot=0,
)

# We can detect frequencies both in periodic and non-periodic signals.
# We start with the periodic complex signal

# ## 5. Frequency Detection in Periodic Signals
#
# This section demonstrates the frequency detection capabilities of the algorithm on signals with known periodicity.

# Our algorithm  for the complex signal should behave the same at zoom level 1 as making the Fourier transform and taking the component with the max coordinate
(fhat, sigma), amp_hist, freqs_hist = findfirstharmonic2(sig; zoomlevels=[1])
@show fhat
@show sigma
@show all(fhat .≈ (f1, f2))
@show sigma .≈ offset

# The results are, of course, the same for other zoom levels
for zl in [[1], [1, 2], [1, 8], [1, 2, 16], nothing]
    fhat, sigma = findfirstharmonic2(sig; zoomlevels=zl)[1]
    fhat = flipsign.(fhat, fhat[1])
    ## @test all(fhat .≈ [f1, f2])
    @show fhat
    @show sigma
end

# ## 6. Frequency Detection in Real Signals
#
# The following examples show the application of the algorithm to real signals, emphasizing the robustness and accuracy of frequency detection.

# Check it on the real signal
for zl in [[1], [1, 2], [1, 8], [1, 2, 16], nothing]
    fhat, sigma = findfirstharmonic2(real.(sig); zoomlevels=zl)[1]
    sigma = flipsign(sigma, fhat[1])
    fhat = flipsign.(fhat, fhat[1])
    ## @test all(fhat .≈ [f1, f2])
    @show fhat
    @show sigma
end

# ## 7. Frequency Detection in Non-Periodic Signals
#
# Non-periodic signals present additional challenges for frequency detection. This section explores the performance of the algorithm on non-periodic signals.

#  And now check on the real non-periodic signal
scales = Int32[]
relerrsX = Float64[]
relerrsY = Float64[]
sigmas = Float32[]
for zl in [[1], [1, 2], [1, 4], [1, 8], [1, 4, 16], nothing]
    fhat, sigma = findfirstharmonic2(real.(sigs); zoomlevels=zl)[1]
    sigma = flipsign(sigma, fhat[1])
    fhat = flipsign.(fhat, fhat[1])
    scale = isnothing(zl) ? minimum(arrsize) : last(zl)
    ## @test all(abs.(fhat .- [f1s, f2s]) .* arrsize .* scale .< 0.50001) # approximately 0.5
    relerr = abs.(fhat .- [f1s, f2s]) .* arrsize
    @show scale
    @show fhat
    @show sigma
    push!(scales, scale)
    push!(relerrsX, relerr[1])
    push!(relerrsY, relerr[2])
    push!(sigmas, sigma)
end

# We see that the error is decreasing with scale
fig, ax, p = lines(scales, 0.5 ./ scales; label="0.5/scale", linestyle=:dot);
scatter!(scales, relerrsX; label="X")
scatter!(scales, relerrsY; label="Y")
lines!(scales, sqrt.(relerrsX .^ 2 .+ relerrsY .^ 2); label="joint x and y")
ax.title = "Relative error in frequency detection"
axislegend()
fig

# The offset error also decreases
fig, ax, l = lines(scales, sigmas; label="restored")
hlines!(offset; label="GT", color=:orange)
axislegend()
ax.title = "Offset detection"
fig

# ## 8. Tilt Reconstruction
#
# Once the frequencies are detected, the tilt can be reconstructed. This section demonstrates the reconstruction of the tilt from the detected frequencies.

# Thus, only from the real signal we have restored the parameters of its main harmonics (we have used however the _a priory_ knowledge about the sign of the tilt).
# Finally, we can reconstruct the tilt from the found frequencies using the same function
fhat, sigma = findfirstharmonic2(real.(sigs))[1]
sigma = flipsign(sigma, fhat[1])
fhat = flipsign.(fhat, fhat[1])
restored_tilt = fourier_tilt(2π * fhat, sigma, arrsize)
fig, ax, hm = showarray(restored_tilt; axis=(title=L"Restored function $t(x)$ ",), rot=0);
Colorbar(fig[1, 2], hm)
fig

# ## 9. Error Analysis in Tilt Reconstruction
#
# Analyzing the error in tilt reconstruction is crucial for evaluating the performance of the detection algorithm. The following figures illustrate the error between the original and reconstructed tilts.

# And we check the error in the restoration
err_tilt = restored_tilt .- tilts
fig, ax, hm = showarray(err_tilt; axis=(title=L"Error $t(x) - \hat{t}(x)$ ",), rot=0);
Colorbar(fig[1, 2], hm)
fig

#  We see that the error is quite small compared with the size of the tilt itself.
