```@meta
EditURL = "../../../tutorials/FindingTilts.jl"
```

# Fnding tilts from the autocorrelation function

The algorithm uses the maxima of the autocorrelation spectra to deduce the tilt parameters.
For this purpoes, a simple function [`fourier_tilt`](@ref) is provided.
Here we demonstrate its work.

First, construct a tilt (a linear function) with this function

````@example FindingTilts
using PhasePlots, PhaseUtils
using CairoMakie
using FFTW, FFTViews
using PhaseFromInterferograms
using PhaseFromInterferograms: fourier_tilt, getslopes
using LaTeXStrings

arrsize = (42, 30)
i, j = 3, 7
f1, f2 = getindex.(fftfreq.(arrsize), [i, j])
offset = π / 3
tilt = fourier_tilt(2π .* (f1, f2), offset, arrsize)
showarray(tilt; axis=(title=L"Linear function $t(x)$ ",), rot=0)
````

and we
can see that in the exponent it will create the linear phase term.

````@example FindingTilts
sig = cis.(tilt)
showarray(angle.(sig); axis=(title=L"Linear function $t(x)$ , wrapped",), rot=0)
````

Now examine its Fourier spectrum

````@example FindingTilts
spec = fft(sig)
showarray(
    abs.(spec);
    axis=(title="Power spectrum of the complex exponent of linear function",),
    rot=0,
)
````

Note that the only non-zero component has index (i,j)

````@example FindingTilts
scatter!(i, j; color=:red)
current_figure()
````

And its argument produces the original offset

````@example FindingTilts
angle(spec[i, j]) ≈ offset
````

Repeat the same example with subpixel coordinates

````@example FindingTilts
subpixelshifts = [0.1, 0.52]
# subpixelshifts = rand(2) # or use rand(2)
f1s, f2s = [f1, f2] .+ subpixelshifts ./ arrsize
````

Construct the tilt with subpixel frequencies.

````@example FindingTilts
tilts = fourier_tilt(2π .* (f1s, f2s), offset, arrsize)
showarray(tilts; axis=(title=L"Linear function $t(x)$ ",), rot=0)
````

It doesn't produce a periodic signal anymore!

````@example FindingTilts
sigs = cis.(tilts)
showarray(angle.(sigs); axis=(title=L"Linear function $2\pi t(x)$ , wrapped",), rot=0)
````

And that's why the Fourier spectrum demonstrates aliasing

````@example FindingTilts
specs = fft(sigs)
showarray(
    abs.(specs);
    axis=(title="Power spectrum of the complex exponent of linear function",),
    rot=0,
)
````

We can detect frequencies both in periodic and non-periodic signals.
We start with the periodic complex signal

````@example FindingTilts
for zl in [[1], [1, 2], [1, 8], [1, 2, 16], nothing]
    fhat = findfirstharmonic2(sig; zoomlevels=zl)[1]
    fhat = flipsign.(fhat, fhat[1])
    # @test all(fhat .≈ [f1, f2])
    @show fhat
end
````

Check it on the real signal

````@example FindingTilts
for zl in [[1], [1, 2], [1, 8], [1, 2, 16], nothing]
    fhat = findfirstharmonic2(real.(sig); zoomlevels=zl)[1]
    fhat = flipsign.(fhat, fhat[1])
    # @test all(fhat .≈ [f1, f2])
    @show fhat
end
````

 And now check on the real non-periodic signal

````@example FindingTilts
scales = Int32[]
relerrsX = Float64[]
relerrsY = Float64[]
for zl in [[1], [1, 2], [1, 4], [1, 8], [1, 4, 16], nothing]
    fhat = findfirstharmonic2(real.(sigs); zoomlevels=zl)[1]
    fhat = flipsign.(fhat, fhat[1])
    scale = isnothing(zl) ? minimum(arrsize) : last(zl)
    # @test all(abs.(fhat .- [f1s, f2s]) .* arrsize .* scale .< 0.50001) # approximately 0.5
    relerr = abs.(fhat .- [f1s, f2s]) .* arrsize
    @show scale
    @show fhat
    push!(scales, scale)
    push!(relerrsX, relerr[1])
    push!(relerrsY, relerr[2])
end
````

We see that the error is decreasing with scale

````@example FindingTilts
fig, ax, p = lines(scales, 0.5 ./ scales; label="0.5/scale", linestyle=:dot);
scatter!(scales, relerrsX; label="X")
scatter!(scales, relerrsY; label="Y")
lines!(scales, sqrt.(relerrsX .^ 2 .+ relerrsY .^ 2); label="joint x and y")
ax.title = "Relative error in frequency detection"
axislegend()
fig
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

