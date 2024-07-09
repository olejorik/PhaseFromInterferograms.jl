using PhaseFromInterferograms
using Test

@testset "PhaseFromInterferograms.jl" begin
    # Write your tests here.
end

@testset "tilts and zoomFFT2D" begin
    arrsize = (20, 25)
    i, j = 3, 7
    f1, f2 = getindex.(fftfreq.(arrsize), [i, j])
    offset = π / 3
    @testset "pixel-level detection" begin
        tilt = fourier_tilt(2π .* (f1, f2), offset, arrsize)
        sig = cis.(tilt)
        spec = fft(sig)
        @test angle(spec[i, j]) ≈ offset
        spec[i, j] = 0
        @test all(1 .+ abs.(spec) .≈ 1)

        for zl in [[1], [1, 2], [1, 8], [1, 2, 16], nothing]
            fhat, sigma = findfirstharmonic2(sig; zoomlevels=zl)[1]
            fhat = flipsign.(fhat, fhat[1])
            @test all(fhat .≈ [f1, f2])
            # @show fhat
        end
    end

    @testset "subpixel detection" begin
        f1s, f2s = [f1, f2] .+ [0.1, 0.77] ./ arrsize
        # @show f1s, f2s
        tilts = fourier_tilt(2π .* (f1s, f2s), offset, arrsize)
        sigs = cis.(tilts)
        specs = fft(sigs)
        for zl in [[1], [1, 2], [1, 8], [1, 2, 16], nothing]
            fhat, sigma = findfirstharmonic2(sigs; zoomlevels=zl)[1]
            fhat = flipsign.(fhat, fhat[1])
            scale = isnothing(zl) ? minimum(arrsize) : last(zl)
            @test all(abs.(fhat .- [f1s, f2s]) .* arrsize .* scale .< 0.50001) # approximately 0.5
            # @show fhat
            # @show scale
        end
    end
end
