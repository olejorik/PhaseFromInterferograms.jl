# Methods usded for the algorithms

######################
# Tilt from parameters
######################

"""
    fourier_tilt(τ, σ, arrsize, fftshift = false)

Construct an array of values of linear function `t(x) = τ⋅x +σ` compatible with Fourier transform coordinates, that is `T(ξ)=F(exp(i t(x)))` has maximum at `ξ=τ/(2π)`. Coordinates: `ξ` is defined by `fftfreq`, `x` defined with the origin at `arrsize÷2+1` for `fftshift = true` and with the origin at at the first element of the array if `fftshift = false`.
"""
function fourier_tilt(τ, σ, arrsize, fftshift=false)
    @assert length(τ) == length(arrsize)
    if fftshift
        coord = [fftshift(fftfreq(d, d)) for d in arrsize]
    else
        coord = [0:(d - 1) for d in arrsize]
    end
    tilt = σ .+ reshape([[τ...]' * [i...] for i in Iterators.product(coord...)], arrsize)

    return tilt
end

# ###############
# Tilt extraction
# ###############

"""
    get_tilt(idiff, erasesize=3, selectsize=3, debug=false)

Get the approximate tilt from the interferogram difference through the inverse Fourier of the part of the Fourier spectrum of `idiff` loacated around the position of the maximal side lobe.

Parameters:
 - `erasesize` -- minimal distance from the sidelobe to the DC component
 - `selectsize` -- number of the neigbour Fourier components to include in the tilt reconstruction.
"""
function get_tilt(idiff, erasesize=3, selectsize=3, debug=false)
    spectrum = fft(idiff .^ 2)
    eraseZerothOrder!(spectrum; erasesize=erasesize)

    # Calculate frequencies
    iii = argmax(abs.(spectrum))
    freqs = [fftfreq(size(idiff)[d])[iii[d]] for d in eachindex(Tuple(iii))]

    if debug
        println("max position $iii")
        println("max value is $(spectrum[iii])\n")
    end

    # select

    aaa = FFTView(spectrum)
    bbb = zeros(ComplexF32, size(spectrum))
    bbbFV = FFTView(bbb)

    bbox = selectsize
    ind1 = one(iii)
    bbbFV[(iii - bbox * ind1):(iii + bbox * ind1)] .= aaa[(iii - bbox * ind1):(iii + bbox * ind1)]

    # delta = angle.(ifft(bbb))
    delta = angle.(ifft(-bbb)) # add π
    return delta, freqs
end  # function get_tilt


"""
    getfinetilt(idiff; n=[1, 0], zoomlevels=nothing, visualdebug=false, erasesize=2, cropsize=2) -> tilt, τ, σ

TBW
"""
function getfinetilt(
    idiff; n=[1, 1], zoomlevels=nothing, visualdebug=false, erasesize=2, cropsize=2
)
    # process default normal
    if isnothing(n)
        n = [1, 1]
    end

    arrsize = size(idiff)
    (f, σ) = first(
        findfirstharmonic2(
            idiff .^ 2;
            zoomlevels=zoomlevels,
            visualdebug=visualdebug,
            erasesize=erasesize,
            cropsize=cropsize,
        ),
    )
    # @show τ, σ
    if dotproduct(f, n) < 0
        f .*= -1
        σ *= -1
    end
    τ = 2π .* f
    σ = phwrap(σ + π)
    # tilt = σ .+ [2π * (i * τ[1] + j * τ[2]) for i in 1:arrsize[1], j in 1:arrsize[2]]
    tilt = fourier_tilt(2π * f, σ, arrsize)
    return tilt, τ, σ
end

# ################
# Phase extraction
# ################


"""
    get_diff_phase_from_n_psi(images, deltas)

Given n+1 interferometric `images` with n known phase shifts `deltas` between them,
reconstruct the phase using the phase-difference-based approach.
"""
function get_diff_phase_from_n_psi(images, deltas)
    n = length(deltas)
    @assert length(images) == n + 1 "Not compatible number of images and phase shifts!"

    cosdelta = [cos.(δ) for δ in deltas]
    sindelta = [sin.(δ) for δ in deltas]

    uC = sum(cosdelta)
    vC = sum([cos.(2 * δ) for δ in deltas])
    uS = sum(sindelta)
    vS = sum([sin.(2 * δ) for δ in deltas])

    Idiff = -diffirst(images)

    a2 = dotproduct(sindelta, Idiff)
    a1 = sum(Idiff) - dotproduct(cosdelta, Idiff)

    bi11 = n .- vC
    bi12 = -2 * uS + vS
    bi22 = 3n .- 4 * uC + vC
    x = a1 .* bi11 + a2 .* bi12
    y = a1 .* bi12 + a2 .* bi22

    return atan.(y, x)
end

"""
    get_LS_ phase_from_n_psi(images, deltas)

Given n interferometric `images` with known phase shifts `deltas` between them,
reconstruct the phase using the least-squares diversity-based method.
"""
function get_LS_phase_from_n_psi(images, deltas)
    n = length(deltas)
    @assert (length(images) == n) "Not compatible number of images and phase shifts!"

    #construct the diversity factors
    ds = [exp.(im * δ) for δ in deltas]

    s1 = sum(ds)
    s2 = sum(x -> x .^ 2, ds)

    ai = [dotproduct(conj.(ds), images), dotproduct(ds, images), sum(images)]

    den = @. 2real(s1^2 * conj(s2)) - n * (2 * abs2(s1) + abs2(s2)) + n^3
    m11 = @. n^2 - abs2(s1)
    m12 = @. conj(s1)^2 - n * conj(s2)
    m13 = @. s1 * conj(s2) - n * conj(s1)

    m31 = conj(m13)
    m32 = @. s1 * conj(s2) - n * conj(s1)
    m33 = @. n^2 - abs2(s2)

    x = dotproduct([m11, m12, m13], ai) ./ den
    b = real(dotproduct([m31, m32, m33], ai) ./ den)

    return angle.(x), x, b
end
