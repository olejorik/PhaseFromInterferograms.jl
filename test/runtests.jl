using PhaseFromInterferograms
using Test
using FFTW
using PhaseFromInterferograms.zoomFFT2D: FFT2Zoom, findfirstharmonic2_v2

@testset "PhaseFromInterferograms.jl" begin
    # Write your tests here.
end

@testset "zoomFFT2D consistency" begin
    ## Test that zoomed FFT matches zero-padded FFT
    ## Based on the principle that zero-padding in spatial domain = interpolation in frequency domain

    @testset "Full FFT interpolation" begin
        ## Create test signals with various sizes
        for (M, N) in [(16, 20), (12, 15), (20, 24)]
            ## Create a simple complex signal
            signal = randn(ComplexF64, M, N)

            ## Method 1: Zero-pad and compute FFT (reference method)
            padfactor = 4
            padM, padN = padfactor * M, padfactor * N
            padded_signal = zeros(ComplexF64, padM, padN)
            padded_signal[1:M, 1:N] .= signal
            reference_fft = fft(padded_signal)

            ## Method 2: Use FFT2Zoom to compute the interpolated FFT
            ## The key insight: FFT2Zoom with fftshift-ed frequency ranges
            ## should match the fftshift-ed zero-padded FFT
            Mset = 0:(M - 1)
            Nset = 0:(N - 1)
            freqs_x = fftshift(fftfreq(padM))
            freqs_y = fftshift(fftfreq(padN))
            Rset = freqs_x[1]:step(freqs_x):freqs_x[end]
            Sset = freqs_y[1]:step(freqs_y):freqs_y[end]

            fft_zoom = FFT2Zoom(Mset, Nset, Rset, Sset)
            zoomed_fft = fft_zoom(signal)

            ## Compare with fftshifted reference
            @test zoomed_fft ≈ fftshift(reference_fft) rtol = 1e-4
        end
    end

    @testset "Subregion zoom" begin
        ## Test zooming into a specific frequency subregion
        M, N = 20, 24
        signal = randn(ComplexF64, M, N)

        ## Create high-resolution reference via zero-padding
        padM, padN = 8M, 8N
        padded = zeros(ComplexF64, padM, padN)
        padded[1:M, 1:N] .= signal
        ref_fft = fftshift(fft(padded))

        ## Extract a centered quarter region
        center_M, center_N = padM ÷ 2, padN ÷ 2
        halfM, halfN = padM ÷ 8, padN ÷ 8
        M_range = (center_M - halfM + 1):(center_M + halfM)
        N_range = (center_N - halfN + 1):(center_N + halfN)
        ref_region = ref_fft[M_range, N_range]

        ## Compute same region with FFT2Zoom
        Mset = 0:(M - 1)
        Nset = 0:(N - 1)
        freqs_x = fftshift(fftfreq(padM))
        freqs_y = fftshift(fftfreq(padN))
        Rset = freqs_x[M_range]
        Sset = freqs_y[N_range]

        fft_zoom = FFT2Zoom(Mset, Nset, Rset, Sset)
        zoomed_region = fft_zoom(signal)

        @test zoomed_region ≈ ref_region rtol = 1e-4
    end

    @testset "Subregion zoom with linear phase" begin
        ## Test zooming around a known frequency peak from a linear phase signal
        ## and validate reconstruction quality at different zoom levels
        arrsize = (25, 31)  # Use odd dimensions to test fftshift/ifftshift behavior
        ## Choose frequencies that will create a clear peak
        f1, f2 = 0.18, -0.15
        offset = π / 5

        ## Create linear phase signal
        tilt = fourier_tilt(2π .* (f1, f2), offset, arrsize)
        signal = cis.(tilt)

        ## Test with increasing zoom levels
        zoom_factors = [4, 8, 16, 32]

        for padfactor in zoom_factors
            ## Create frequency grid for the zoomed region
            padsize = padfactor .* arrsize
            freqs_x = fftshift(fftfreq(padsize[1]))
            freqs_y = fftshift(fftfreq(padsize[2]))

            ## Find indices closest to our target frequencies
            idx_f1 = argmin(abs.(freqs_x .- f1))
            idx_f2 = argmin(abs.(freqs_y .- f2))

            ## Extract a region around the peak (±16 samples)
            half_width = 16
            M_range = max(1, idx_f1 - half_width):min(padsize[1], idx_f1 + half_width)
            N_range = max(1, idx_f2 - half_width):min(padsize[2], idx_f2 + half_width)

            ## Compute zoomed FFT in the region
            Mset = 0:(arrsize[1] - 1)
            Nset = 0:(arrsize[2] - 1)
            Rset = freqs_x[M_range]
            Sset = freqs_y[N_range]

            fft_zoom = FFT2Zoom(Mset, Nset, Rset, Sset)
            zoomed_region = fft_zoom(signal)

            ## Find peak and extract parameters
            zoom_peak_idx = argmax(abs.(zoomed_region))
            peak_freq_x = Rset[zoom_peak_idx[1]]
            peak_freq_y = Sset[zoom_peak_idx[2]]
            peak_phase = angle(zoomed_region[zoom_peak_idx])

            ## Reconstruct the tilt from detected parameters
            reconstructed_tilt = fourier_tilt(
                2π .* (peak_freq_x, peak_freq_y), peak_phase, arrsize
            )

            ## Calculate RMS error between original and reconstructed tilt
            phase_diff = angle.(exp.(im .* (tilt .- reconstructed_tilt)))
            rms_error = sqrt(sum(phase_diff .^ 2) / length(phase_diff))

            @info "Zoom factor: $padfactor, RMS error: $(round(rms_error, digits=6)) rad"
        end

        ## Final validation: highest zoom should give very good reconstruction
        padfactor = 32
        padsize = padfactor .* arrsize
        freqs_x = fftshift(fftfreq(padsize[1]))
        freqs_y = fftshift(fftfreq(padsize[2]))
        idx_f1 = argmin(abs.(freqs_x .- f1))
        idx_f2 = argmin(abs.(freqs_y .- f2))
        half_width = 16
        M_range = max(1, idx_f1 - half_width):min(padsize[1], idx_f1 + half_width)
        N_range = max(1, idx_f2 - half_width):min(padsize[2], idx_f2 + half_width)

        Mset = 0:(arrsize[1] - 1)
        Nset = 0:(arrsize[2] - 1)
        Rset = freqs_x[M_range]
        Sset = freqs_y[N_range]

        fft_zoom = FFT2Zoom(Mset, Nset, Rset, Sset)
        zoomed_region = fft_zoom(signal)
        zoom_peak_idx = argmax(abs.(zoomed_region))
        peak_freq_x = Rset[zoom_peak_idx[1]]
        peak_freq_y = Sset[zoom_peak_idx[2]]
        peak_phase = angle(zoomed_region[zoom_peak_idx])
        reconstructed_tilt = fourier_tilt(
            2π .* (peak_freq_x, peak_freq_y), peak_phase, arrsize
        )
        phase_diff = angle.(exp.(im .* (tilt .- reconstructed_tilt)))
        rms_error = sqrt(sum(phase_diff .^ 2) / length(phase_diff))

        ## Verify high accuracy with sufficient zoom
        @test peak_freq_x ≈ f1 atol = 0.005
        @test peak_freq_y ≈ f2 atol = 0.005
        @test rms_error < 0.05  # Less than 0.05 radians RMS error
    end
end

# @testset "tilts and zoomFFT2D" begin
#     arrsize = (20, 25)
#     i, j = 3, 7
#     f1, f2 = getindex.(fftfreq.(arrsize), [i, j])
#     offset = π / 3
#     @testset "pixel-level detection" begin
#         tilt = fourier_tilt(2π .* (f1, f2), offset, arrsize)
#         sig = cis.(tilt)
#         spec = fft(sig)
#         @test angle(spec[i, j]) ≈ offset
#         spec[i, j] = 0
#         @test all(1 .+ abs.(spec) .≈ 1)

#         for zl in [[1], [1, 2], [1, 8], [1, 2, 16], nothing]
#             fhat, sigma = findfirstharmonic2(sig; zoomlevels=zl)[1]
#             fhat = flipsign.(fhat, fhat[1])
#             @test all(fhat .≈ [f1, f2])
#             # @show fhat
#         end
#     end

#     @testset "subpixel detection" begin
#         f1s, f2s = [f1, f2] .+ [0.1, 0.77] ./ arrsize
#         # @show f1s, f2s
#         tilts = fourier_tilt(2π .* (f1s, f2s), offset, arrsize)
#         sigs = cis.(tilts)
#         specs = fft(sigs)
#         for zl in [[1], [1, 2], [1, 8], [1, 2, 16], nothing]
#             fhat, sigma = findfirstharmonic2(sigs; zoomlevels=zl)[1]
#             fhat = flipsign.(fhat, fhat[1])
#             scale = isnothing(zl) ? minimum(arrsize) : last(zl)
#             @test all(abs.(fhat .- [f1s, f2s]) .* arrsize .* scale .< 0.50001) # approximately 0.5
#             # @show fhat
#             # @show scale
#         end
#     end
# end

@testset "findfirstharmonic2_v2" begin
    arrsize = (25, 31)
    ## Test with exact pixel frequencies
    i, j = 5, 8
    f1, f2 = fftfreq(arrsize[1])[i], fftfreq(arrsize[2])[j]
    offset = π / 4

    @testset "Pixel-level detection v2" begin
        tilt = fourier_tilt(2π .* (f1, f2), offset, arrsize)
        sig = cis.(tilt)

        ## Test with single zoom level (should match FFT result)
        (fhat, phase), _, _ = findfirstharmonic2_v2(sig; zoomlevels=[1])
        @test fhat[1] ≈ f1 rtol = 1e-10
        @test fhat[2] ≈ f2 rtol = 1e-10
        @test phase ≈ offset rtol = 1e-10
    end

    @testset "Subpixel detection v2" begin
        ## Add fractional frequency offsets
        f1_sub = f1 + 0.12 / arrsize[1]
        f2_sub = f2 + 0.37 / arrsize[2]
        offset_sub = π / 5

        tilt = fourier_tilt(2π .* (f1_sub, f2_sub), offset_sub, arrsize)
        sig = cis.(tilt)

        ## Test with multiple zoom levels
        for zl in [[1, 4], [1, 8], [1, 4, 16], [1, 8, 32]]
            (fhat, phase), amp_hist, freqs_hist = findfirstharmonic2_v2(sig; zoomlevels=zl)

            ## Check frequency detection accuracy improves with zoom
            freq_error = sqrt(sum((fhat .- [f1_sub, f2_sub]) .^ 2))
            max_zoom = last(zl)
            ## Use a more realistic accuracy bound (factor of 2 for discretization effects)
            expected_accuracy = 2.0 / (2 * max_zoom * minimum(arrsize))

            @test freq_error < expected_accuracy

            ## Verify history length matches zoom levels
            @test length(freqs_hist) == length(zl)
            @test length(amp_hist) == length(zl)
        end
    end

    @testset "Reconstruction accuracy v2" begin
        ## Use ground truth parameters for full reconstruction test
        f1_gt, f2_gt = 0.18, -0.15
        offset_gt = π / 5

        tilt = fourier_tilt(2π .* (f1_gt, f2_gt), offset_gt, arrsize)
        sig = cis.(tilt)

        ## High zoom for accurate reconstruction
        (fhat, phase), _, _ = findfirstharmonic2_v2(sig; zoomlevels=[1, 8, 32], erasesize=0)

        ## Reconstruct and check RMS error
        reconstructed_tilt = fourier_tilt(2π .* (fhat[1], fhat[2]), phase, arrsize)
        phase_diff = angle.(exp.(im .* (tilt .- reconstructed_tilt)))
        rms_error = sqrt(sum(phase_diff .^ 2) / length(phase_diff))

        @test fhat[1] ≈ f1_gt atol = 0.005
        @test fhat[2] ≈ f2_gt atol = 0.005
        @test rms_error < 0.05  # Less than 0.05 radians RMS
    end

    @testset "DC exclusion for real signals v2" begin
        ## Test with intensity modulation pattern (like interferogram squared)
        ## Create a realistic signal: slowly varying envelope with oscillating tilt
        f1_gt, f2_gt = 0.18, -0.15
        offset_gt = π / 5

        tilt = fourier_tilt(2π .* (f1_gt, f2_gt), offset_gt, arrsize)

        ## Create slowly varying background (low frequency component)
        slow_x = range(0, 2π; length=arrsize[1])
        slow_y = range(0, 2π; length=arrsize[2])
        slow_background = 1.0 .+ 0.3 .* cos.(slow_x) .+ 0.3 .* cos.(slow_y)'

        ## Create interferogram-like signal with modulation: I = background * (1 + contrast*cos(tilt))
        ## This creates a real signal with DC + sidelobes
        contrast = 0.8
        real_signal = slow_background .* (1.0 .+ contrast .* cos.(tilt))
        ## The FFT of real_signal will have large DC component
        ## With erasesize=2 (default), should find the fundamental frequency
        (fhat, phase), _, _ = findfirstharmonic2_v2(
            real_signal; zoomlevels=[1, 8], erasesize=2
        )

        ## For squared interferogram, dominant non-DC peak is at 2*f
        ## Actually for (1 + 0.5*cos(tilt))^2 = 1.25 + cos(tilt) + 0.125*cos(2*tilt)
        ## The fundamental at f should be strongest non-DC component
        @test abs(fhat[1]) ≈ abs(f1_gt) atol = 0.02
        @test abs(fhat[2]) ≈ abs(f2_gt) atol = 0.02

        ## Test with erasesize=0 should find DC (frequency close to zero)
        (fhat_dc, _), _, _ = findfirstharmonic2_v2(real_signal; zoomlevels=[1], erasesize=0)
        @test abs(fhat_dc[1]) < 0.05  # Should be near zero (DC)
        @test abs(fhat_dc[2]) < 0.05
    end
end

@testset "PTIestimate interface" begin
    # These are re-exported by PhaseFromInterferograms
    using PhaseFromInterferograms: FreeTilt, sigma, tau
    # Import internal accessors
    import PhaseFromInterferograms: tilts, gettilts

    @testset "Basic construction and dimensions" begin
        ## Create PTIestimate with 2D frames and 1D set
        framesize = (32, 40)
        setsize = (5,)
        pti = PTIestimate(framesize, setsize)

        @test PhaseFromInterferograms.framesize(pti) == framesize
        @test PhaseFromInterferograms.setsize(pti) == setsize
        @test pti.framesize == 2
        @test pti.setsize == 1
        @test pti.fullsize == (32, 40, 5)
    end

    @testset "tilts() with view keyword" begin
        framesize = (20, 25)
        setsize = (3, 4)
        pti = PTIestimate(framesize, setsize)

        ## Test :raw view (default, broadcasting shape)
        t_raw = tilts(pti)
        @test size(t_raw) == (1, 1, 3, 4)
        @test t_raw === pti.tilts  # Should return the field directly

        ## Test :raw view (explicit)
        t_raw_explicit = tilts(pti; view=:raw)
        @test t_raw_explicit === t_raw

        ## Test :set view (set-only dimensions)
        t_set = tilts(pti; view=:set)
        @test size(t_set) == (3, 4)
        @test ndims(t_set) == 2

        ## Verify it's a reshaped view (not copying data)
        @test t_set isa AbstractArray

        ## Test that invalid view throws error
        @test_throws ArgumentError tilts(pti; view=:invalid)
    end

    @testset "tilts() iteration with :set view" begin
        framesize = (15, 20)
        setsize = (2, 3)
        pti = PTIestimate(framesize, setsize)

        ## Set some test tilts with known coefficients
        for (i, idx) in enumerate(CartesianIndices(setsize))
            pti.tilts[1, 1, idx] = FreeTilt([i * 0.1, i * 0.2, i * 0.3])
        end

        ## Iterate with :set view and collect results
        t_set = tilts(pti; view=:set)
        indices_and_tilts = collect(pairs(t_set))

        ## Check we got the right number of elements
        @test length(indices_and_tilts) == prod(setsize)

        ## Check one example to verify structure
        idx, tilt = first(indices_and_tilts)
        @test idx isa CartesianIndex{2}
        @test all(1 .<= Tuple(idx) .<= setsize)
        @test sigma(tilt) isa Real
        @test tau(tilt) isa AbstractVector
    end

    @testset "gettilts() with view keyword" begin
        framesize = (10, 12)
        setsize = (2, 3)
        pti = PTIestimate(framesize, setsize)

        ## Set different tilts for each set index
        for idx in CartesianIndices(setsize)
            i = LinearIndices(setsize)[idx]
            pti.tilts[1, 1, idx] = FreeTilt([i * π / 6, i * 0.1, i * 0.05])
        end

        ## Test :raw view (default, full tensor)
        t_eval = gettilts(pti)
        @test size(t_eval) == (framesize..., setsize...)
        @test t_eval isa Array{Float64}

        ## Test :raw view (explicit)
        t_eval_explicit = gettilts(pti; view=:raw)
        @test t_eval_explicit == t_eval

        ## Test :set view (slices by set dimensions)
        t_slices = gettilts(pti; view=:set)
        @test length(t_slices) == prod(setsize)

        ## Check one slice to verify structure
        idx, slice = first(pairs(t_slices))
        @test size(slice) == framesize
        @test slice isa AbstractArray{Float64}

        ## Test that invalid view throws error
        @test_throws ArgumentError gettilts(pti; view=:invalid)
    end

    @testset "gettilts() evaluation correctness" begin
        framesize = (8, 10)
        setsize = (2,)
        pti = PTIestimate(framesize, setsize)

        ## Set known tilt coefficients
        σ1, τ1x, τ1y = π / 4, 0.1, 0.2
        pti.tilts[1, 1, 1] = FreeTilt([σ1, τ1x, τ1y])

        ## Get evaluated tilts
        t_eval = gettilts(pti; view=:raw)

        ## Manually compute expected values for first set index
        coords_x = pti.frameaxes[1]
        coords_y = pti.frameaxes[2]
        expected = [σ1 + τ1x * x + τ1y * y for x in coords_x, y in coords_y]

        ## Compare
        @test t_eval[:, :, 1] ≈ expected
    end

    @testset "Frame-only accessors" begin
        framesize = (16, 20)
        setsize = (3,)
        pti = PTIestimate(framesize, setsize)

        ## Test getbackground
        bg = PhaseFromInterferograms.getbackground(pti)
        @test size(bg) == framesize
        @test bg isa Array{Float64}

        ## Test getphase
        phase = PhaseFromInterferograms.getphase(pti)
        @test size(phase) == framesize
        @test phase isa Array{Float64}

        ## Test getapodization
        apod = PhaseFromInterferograms.getapodization(pti)
        @test size(apod) == framesize
        @test apod isa Array{Float64}

        ## Verify values are sensible
        @test all(apod .>= 0)  # Apodization should be non-negative
    end

    @testset "Backward compatibility" begin
        ## Ensure default behavior matches old interface
        framesize = (12, 15)
        setsize = (4,)
        pti = PTIestimate(framesize, setsize)

        ## tilts() with no keyword should return raw field
        @test tilts(pti) === pti.tilts

        ## gettilts() with no keyword should return evaluated full tensor
        t_eval_default = gettilts(pti)
        t_eval_raw = gettilts(pti; view=:raw)
        @test t_eval_default == t_eval_raw
    end

    @testset "Multi-dimensional sets" begin
        ## Test with 2D set dimensions
        framesize = (10, 12)
        setsize = (3, 4)
        pti = PTIestimate(framesize, setsize)

        ## Test tilts :set view
        t_set = tilts(pti; view=:set)
        @test size(t_set) == setsize

        ## Test gettilts :set view
        t_eval_set = gettilts(pti; view=:set)
        @test length(t_eval_set) == prod(setsize)

        ## Check one example
        idx, slice = first(pairs(t_eval_set))
        @test size(slice) == framesize
        @test idx isa CartesianIndex{2}
    end
end
