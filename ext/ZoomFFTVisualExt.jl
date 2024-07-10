module ZoomFFTVisualExt

using PhaseFromInterferograms, CairoMakie
using FFTW, FFTViews
using PhaseFromInterferograms.zoomFFT2D: removeDC, FFT2Zoom
import PhaseFromInterferograms.zoomFFT2D: findfirstharmonic2
function findfirstharmonic2visual(
    a; zoomlevels=nothing, visualdebug=false, erasesize=2, cropsize=2
)
    arrsize = size(a)

    if isnothing(zoomlevels)
        zoomlevels = [
            1,
            2,
            4,
            8,
            16,
            floor(maximum(arrsize) / 4),
            floor(maximum(arrsize) / 2),
            minimum(arrsize),
        ]
    end

    freqrange = fftfreq.(arrsize)

    if visualdebug
        display(
            heatmap(
                fftshift(freqrange[1]),
                fftshift(freqrange[2]),
                abs.(fftshift(fft(a)));
                axis=(
                    aspect=DataAspect(), title="Abs of the Fourier transform of the signal"
                ),
            ),
        )
    end

    signal = removeDC(a, erasesize)
    # signal = a
    spectrum = fft(signal)
    if visualdebug
        display(
            heatmap(
                fftshift(freqrange[1]),
                fftshift(freqrange[2]),
                abs.(fftshift(spectrum));
                axis=(aspect=DataAspect(), title="DC removed from the signal"),
            ),
        )
    end

    aaa = FFTView(spectrum)
    freqrange0 = FFTView.(freqrange)
    iii = argmax(abs.(aaa[0:(arrsize[1] ÷ 2), :]))
    rough_freqs = [freqrange0[i][iii[i]] for i in 1:length(iii)]

    # cropsize = 2
    indexcrop = CartesianIndex(cropsize, cropsize)
    xfreqs = freqrange0[1][(iii[1] - cropsize):(iii[1] + cropsize)]
    yfreqs = freqrange0[2][(iii[2] - cropsize):(iii[2] + cropsize)]

    if visualdebug
        hm = heatmap(
            xfreqs,
            yfreqs,
            abs.(spectrum[(iii - indexcrop):(iii + indexcrop)]);
            axis=(aspect=DataAspect(), title="Spectrum around the maximum"),
        )
        # scatter!(gt_freqs..., label = "True frequency")
        # axislegend()
        display(hm)
    end

    Mset = 1:arrsize[1]
    Nset = 1:arrsize[2]

    last_freqs = copy(rough_freqs)
    freqs_hyst = []
    amp_hyst = []
    for zoom in zoomlevels
        # zoom = zoomlevels[2]

        Rset = (fftshift(freqrange[1]) .- last_freqs[1]) / zoom .+ last_freqs[1]
        Sset = (fftshift(freqrange[2]) .- last_freqs[2]) / zoom .+ last_freqs[2]

        fftX2 = FFT2Zoom(Mset, Nset, Rset, Sset)
        zoomedspectrum = fftX2(signal)
        jjj = argmax(abs.(zoomedspectrum))
        amp = zoomedspectrum[jjj]
        fine_freqs = [Rset[jjj[1]], Sset[jjj[2]]]
        # select the frequncy in the positive half
        # if fine_freqs[1] ≈ 0
        #     if sign(fine_freqs[2]) == -1
        #         fine_freqs = fine_freqs .* sign(fine_freqs[2])
        #         amp = conj(amp)
        #     end
        # else
        #     if sign(fine_freqs[1]) == -1
        #         fine_freqs = fine_freqs .* sign(fine_freqs[1])
        #         amp = conj(amp)
        #     end
        # end
        # if fine_freqs[1] * n[1] + fine_freqs[2] * n[2] < 0
        #     amp = conj(amp)
        #     fine_freqs .*= -1
        # end

        if visualdebug
            fig, ax, hm = heatmap(
                Rset,
                Sset,
                abs.(zoomedspectrum);
                axis=(aspect=DataAspect(), title="zoomed x$zoom part of the spectrum"),
            )
            vlines!(last_freqs[1]; color=:gray)
            hlines!(last_freqs[2]; color=:gray, label="Previous fine frequency")
            # scatter!(gt_freqs..., color=:red, label = "True frequency")
            scatter!(fine_freqs...; color=:blue, label="New fine frequency")
            axislegend()
            display(fig)
        end

        last_freqs .= fine_freqs
        push!(freqs_hyst, fine_freqs)
        push!(amp_hyst, amp)
    end
    return last_freqs, amp_hyst, freqs_hyst
end  # function findfirstharmonic

end
