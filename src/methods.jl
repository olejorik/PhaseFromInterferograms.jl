# Methods usded for the algorithms

# ###############
# Tilt extraction
# ###############

"""
    get_tilt(idiff)

Get the tilt from the interferogram difference
"""
function get_tilt(idiff, erasesize=3, selectsize=3, debug=false)
    spectrum = fft(idiff .^ 2)
    eraseZerothOrder!(spectrum)
    aaa = FFTView(spectrum)
    bbb = zeros(ComplexF32, size(spectrum))
    bbbFV = FFTView(bbb)

    iii = argmax(abs.(spectrum))
    if debug
        println("max position $iii")
        println("max value is $(spectrum[iii])\n")
    end

    # Calculate frequencies
    freqs = [fftfreq(size(idiff)[d])[iii[d]] for d in eachindex(Tuple(iii))]

    # select

    bbox = selectsize
    bbbFV[(iii - CartesianIndex(bbox, bbox)):(iii + CartesianIndex(bbox, bbox))] .= aaa[(iii - CartesianIndex(
        bbox, bbox
    )):(iii + CartesianIndex(bbox, bbox))]

    delta = angle.(ifft(-bbb)) # add π
    return delta, freqs
end  # function get_tilt


"""
    getfinetilt(idiff, n=[1, -1]) -> tilt, τ, σ

TBW
"""
function getfinetilt(idiff, n=[1, 1])
    arrsize = size(idiff)
    res = findfirstharmonic2(idiff .^ 2; visualdebug=false)
    τ = res[1]
    σ = angle(res[2][end]) - π
    if dotproduct(τ, n) < 0
        τ .*= -1
        σ *= -1
    end
    tilt = σ .+ [2π * (i * τ[1] + j * τ[2]) for i in 1:arrsize[1], j in 1:arrsize[2]]
    return tilt, τ, σ
end

# ################
# Phase extraction
# ################


"""
    get_phase_from_n_psi(images, deltas)

Given n interferometric `images` with known phase shifts `deltas` between them,
reconstruct the phase.
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
reconstruct the phase using the least-squares diversity -based method.
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
    b = dotproduct([m31, m32, m33], ai) ./ den

    return angle.(x), x, b
end
