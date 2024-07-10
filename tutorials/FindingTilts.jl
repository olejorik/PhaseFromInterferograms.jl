# # Finding tilts from the autocorrelation function
#
# The algorithm uses the maxima of the autocorrelation spectra to deduce the tilt parameters.
# For this purpoes, a simple function [`fourier_tilt`](@ref) is provided.
# Here we demonstrate its work.
#
# First, construct a tilt (a linear function) with this function

using PhasePlots, PhaseUtils
using CairoMakie
using FFTW, FFTViews
using PhaseFromInterferograms
using PhaseFromInterferograms: fourier_tilt, getslopes
using LaTeXStrings

arrsize = (50, 35)
i, j = 3, 8
f1, f2 = getindex.(fftfreq.(arrsize), [i, j])
offset = 2
tilt = fourier_tilt(2π .* (f1, f2), offset, arrsize)
showarray(tilt; axis=(title=L"Linear function $t(x)$ ",), rot=0)

# and we
# can see that in the exponent it will create the linear phase term.
sig = cis.(tilt)
showarray(angle.(sig); axis=(title=L"Linear function $t(x)$ , wrapped",), rot=0)

# Now examine its Fourier spectrum
spec = fft(sig)
showarray(
    abs.(spec);
    axis=(title="Power spectrum of the complex exponent of linear function",),
    rot=0,
)

# Note that the only non-zero component has index (i,j)
scatter!(i, j; color=:red)
current_figure()

# And its argument produces the original offset
angle(spec[i, j]) ≈ offset

# Repeat the same example with subpixel coordinates

subpixelshifts = [0.1, 0.52]
## subpixelshifts = rand(2) # or use rand(2)
f1s, f2s = [f1, f2] .+ subpixelshifts ./ arrsize

# Construct the tilt with subpixel frequencies.
tilts = fourier_tilt(2π .* (f1s, f2s), offset, arrsize)
showarray(tilts; axis=(title=L"Linear function $t(x)$ ",), rot=0)

# It doesn't produce a periodic signal anymore!
sigs = cis.(tilts)
showarray(angle.(sigs); axis=(title=L"Linear function $2\pi t(x)$ , wrapped",), rot=0)

# And that's why the Fourier spectrum demonstrates aliasing
specs = fft(sigs)
showarray(
    abs.(specs);
    axis=(title="Power spectrum of the complex exponent of linear function",),
    rot=0,
)




# We can detect frequencies both in periodic and non-periodic signals.
# We start with the periodic complex signal

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

# Check it on the real signal
for zl in [[1], [1, 2], [1, 8], [1, 2, 16], nothing]
    fhat, sigma = findfirstharmonic2(real.(sig); zoomlevels=zl)[1]
    sigma = flipsign(sigma, fhat[1])
    fhat = flipsign.(fhat, fhat[1])
    ## @test all(fhat .≈ [f1, f2])
    @show fhat
    @show sigma
end

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

# Thus, only from the real signal we have restored the parameters of its main harmonics (we have used however the _a priory_ knowledge about the sign of the tilt).
# Finally, we can reconstruct the tilt from the found frequencies using the same function
fhat, sigma = findfirstharmonic2(real.(sigs))[1]
sigma = flipsign(sigma, fhat[1])
fhat = flipsign.(fhat, fhat[1])
restored_tilt = fourier_tilt(2π * fhat, sigma, arrsize)
fig, ax, hm = showarray(restored_tilt; axis=(title=L"Restored function $t(x)$ ",), rot=0);
Colorbar(fig[1, 2], hm)
fig

# And we check the error in the restoration
err_tilt = restored_tilt .- tilts
fig, ax, hm = showarray(err_tilt; axis=(title=L"Error $t(x) - \hat{t}(x)$ ",), rot=0);
Colorbar(fig[1, 2], hm)
fig

#  We see that the errror is quite small compared with the size of the tilt itself.
