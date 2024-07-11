module zoomFFT2D
# Based on Jurling article "D:\Documents\articles unsorted\josaa-35-11-1784.pdf"

using FFTW
using FFTViews
# using CairoMakie # for debug only

struct FFT2Zoom
    Mset::Vector{Real}
    Nset::Vector{Real}
    Rset::Vector{Real}
    Sset::Vector{Real}
    αx::Real
    αy::Real
    Kprime::Int
    Lprime::Int
    a::Array{Complex,2}
    bhat::Array{Complex,2}
    Hhat::Array{Complex,2}
end

function nextfftlength(n)
    return nextprod((2, 3, 5), n)
end

"""
    expsquared(x, coef)

Calculate exp(coef * x^2)
"""
function expsquared(x, coef)
    return exp.(coef .* (x .^ 2))
end  # function expsquared

"""
    FFT2Zoom(M,N,R,S,αx, αy)

Structure containing necessary arrays for calculation of zoomed 2D FFT
"""
function FFT2Zoom(Mset, Nset, Rset, Sset, αx, αy)
    M, N, R, S = length.((Mset, Nset, Rset, Sset))
    Kprime = nextfftlength(M + R - 1)
    Lprime = nextfftlength(N + S - 1)

    bhat = zeros(ComplexF32, (Kprime, Lprime))
    # for m in 0:M-1, n in 0:N-1
    #     bhat[m+1,n+1] = exp(-im*π*(m^2*αx+n^2*αy))
    # end
    bhat[1:M, 1:N] .=
        expsquared(Mset, -im * π * αx) .* transpose(expsquared(Nset, -im * π * αy))

    a = zeros(ComplexF32, (R, S))
    # for r in 0:R-1, s in 0:S-1
    #     a[r+1,s+1] = exp(-im*π*(r^2*αx+s^2*αy))
    # end
    a[1:R, 1:S] .=
        expsquared(Rset, -im * π * αx) .* transpose(expsquared(Sset, -im * π * αy))

    hhat = zeros(ComplexF32, (Kprime, Lprime))
    hhatx = zeros(ComplexF32, Kprime)
    hhaty = zeros(ComplexF32, Lprime)
    for m in (0:(R - 1))
        hhatx[m + 1] = exp(im * π * (m + first(Rset) - first(Mset))^2 * αx)
    end
    for m in (Kprime - M + 1):(Kprime - 1)
        hhatx[m + 1] = exp(im * π * (m - Kprime + first(Rset) - first(Mset))^2 * αx)
    end
    for n in (0:(S - 1))
        hhaty[n + 1] = exp(im * π * (n + first(Sset) - first(Nset))^2 * αy)
    end
    for n in (Lprime - N + 1):(Lprime - 1)
        hhaty[n + 1] = exp(im * π * (n - Lprime + first(Sset) - first(Nset))^2 * αy)
    end

    Hhat = fft(hhatx .* transpose(hhaty))
    # for m in 1:Kprime, n in 1:Lprime
    #     hhat[m,n] = hhatx[m] * hhaty[n]
    # end
    # Hhat = fft(hhat)

    return FFT2Zoom(Mset, Nset, Rset, Sset, αx, αy, Kprime, Lprime, a, bhat, Hhat)
end  # function FFT2ZoomM,N,R,S,αx, αy

function FFT2Zoom(
    Mset::AbstractRange, Nset::AbstractRange, Rset::AbstractRange, Sset::AbstractRange
)
    dm, dn, dr, ds = step.((Mset, Nset, Rset, Sset))
    Mset /= dm
    Nset /= dn
    Rset /= dr
    Sset /= ds
    αx = dm * dr
    αy = dn * ds
    M, N, R, S = length.((Mset, Nset, Rset, Sset))
    Kprime = nextfftlength(M + R - 1)
    Lprime = nextfftlength(N + S - 1)

    bhat = zeros(ComplexF32, (Kprime, Lprime))
    # for m in 0:M-1, n in 0:N-1
    #     bhat[m+1,n+1] = exp(-im*π*(m^2*αx+n^2*αy))
    # end
    bhat[1:M, 1:N] .=
        expsquared(Mset, -im * π * αx) .* transpose(expsquared(Nset, -im * π * αy))

    a = zeros(ComplexF32, (R, S))
    # for r in 0:R-1, s in 0:S-1
    #     a[r+1,s+1] = exp(-im*π*(r^2*αx+s^2*αy))
    # end
    a[1:R, 1:S] .=
        expsquared(Rset, -im * π * αx) .* transpose(expsquared(Sset, -im * π * αy))

    hhat = zeros(ComplexF32, (Kprime, Lprime))
    hhatx = zeros(ComplexF32, Kprime)
    hhaty = zeros(ComplexF32, Lprime)
    for m in (0:(R - 1))
        hhatx[m + 1] = exp(im * π * (m + first(Rset) - first(Mset))^2 * αx)
    end
    for m in (Kprime - M + 1):(Kprime - 1)
        hhatx[m + 1] = exp(im * π * (m - Kprime + first(Rset) - first(Mset))^2 * αx)
    end
    for n in (0:(S - 1))
        hhaty[n + 1] = exp(im * π * (n + first(Sset) - first(Nset))^2 * αy)
    end
    for n in (Lprime - N + 1):(Lprime - 1)
        hhaty[n + 1] = exp(im * π * (n - Lprime + first(Sset) - first(Nset))^2 * αy)
    end

    Hhat = fft(hhatx .* transpose(hhaty))
    # for m in 1:Kprime, n in 1:Lprime
    #     hhat[m,n] = hhatx[m] * hhaty[n]
    # end
    # Hhat = fft(hhat)

    return FFT2Zoom(Mset, Nset, Rset, Sset, αx, αy, Kprime, Lprime, a, bhat, Hhat)
end  # function FFT2ZoomM,N,R,S

"""
    (f::FFT2Zoom)(x)

Calculate zoomed FFT of x
"""
function (f::FFT2Zoom)(x)
    ghat = zeros(ComplexF32, (f.Kprime, f.Lprime))
    M, N, R, S = length.((f.Mset, f.Nset, f.Rset, f.Sset))
    for m in 0:(M - 1), n in 0:(N - 1)
        ghat[m + 1, n + 1] = x[m + 1, n + 1]
    end
    return f.a .* (ifft(fft(f.bhat .* ghat) .* f.Hhat))[1:R, 1:S]
end

"""
    findfirstharmonic2(a, zoomlevels = nothing) ->(freqs, amp), amp_hist, freqs_hist

Find first harmonic in a 2D array `a` with subpixel precision, using sequential FFT2Zoom with `zoomlevels`
"""
function findfirstharmonic2(
    a; zoomlevels=nothing, visualdebug=false, erasesize=2, cropsize=2
)
    arrsize = size(a)

    if isnothing(zoomlevels)
        zl = [1, 2, 4, 8, 16, minimum(arrsize) ÷ 4, minimum(arrsize) ÷ 2, minimum(arrsize)]
    else
        zl = zoomlevels
    end

    freqrange = fftfreq.(arrsize)



    signal = removeDC(a, erasesize)
    # signal = a
    spectrum = fft(signal)


    aaa = FFTView(spectrum)
    freqrange0 = FFTView.(freqrange)
    iii = argmax(abs.(aaa))
    rough_freqs = [freqrange0[i][iii[i]] for i in 1:length(iii)]

    # cropsize = 2
    indexcrop = CartesianIndex(cropsize, cropsize)
    xfreqs = freqrange0[1][(iii[1] - cropsize):(iii[1] + cropsize)]
    yfreqs = freqrange0[2][(iii[2] - cropsize):(iii[2] + cropsize)]



    Mset = (1:arrsize[1]) .- 1
    Nset = (1:arrsize[2]) .- 1

    last_freqs = copy(rough_freqs)
    freqs_hist = []
    amp_hist = []
    for zoom in zl
        # zoom = zl[2]

        Rset = (fftshift(freqrange[1]) .- last_freqs[1]) / zoom .+ last_freqs[1]
        Sset = (fftshift(freqrange[2]) .- last_freqs[2]) / zoom .+ last_freqs[2]

        fftX2 = FFT2Zoom(Mset, Nset, Rset, Sset)
        zoomedspectrum = fftX2(signal)
        jjj = argmax(abs.(zoomedspectrum))
        amp = zoomedspectrum[jjj]
        fine_freqs = [Rset[jjj[1]], Sset[jjj[2]]]

        last_freqs .= fine_freqs
        # @show last_freqs
        # @show fine_freqs
        push!(freqs_hist, fine_freqs)
        push!(amp_hist, amp)
    end
    return (last_freqs, angle(last(amp_hist))), amp_hist, freqs_hist
end  # function findfirstharmonic

"""
    removeDC(a, size = 2)

Remove the zeroth lobe of width 2*`size`+1 from the signal
"""
function removeDC(a, erasesize=0)
    spectrum = fft(a)
    aaa = FFTView(spectrum)
    I1 = one(CartesianIndex{ndims(aaa)})
    aaa[(-erasesize * I1):(erasesize * I1)] .= 0
    return ifft(spectrum)
end  # function removeDC

# to make it compatible with Quarto, not needed in the final version of the module
# dpng(x) = display("image/png", x)

dpng(x) = nothing

export findfirstharmonic2

end
