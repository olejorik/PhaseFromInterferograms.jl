module PhaseFromInterferograms
using FFTW
using FFTViews
using PhaseUtils: phwrap, dotproduct

# Write your package code here.
include("utils.jl")

include("Windowing.jl")
using .Windowing
export GaussianWindow

include("FindHarmonics.jl")
using .FindHarmonics
export eraseZerothOrder, eraseZerothOrder!, findfirstharmonic

include("zoomFFTmodule.jl")
using .zoomFFT2D
export findfirstharmonic2

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

    # select

    bbox = selectsize
    bbbFV[(iii - CartesianIndex(bbox, bbox)):(iii + CartesianIndex(bbox, bbox))] .= aaa[(iii - CartesianIndex(
        bbox, bbox
    )):(iii + CartesianIndex(bbox, bbox))]

    delta = angle.(ifft(-bbb)) # add π
    return delta
end  # function get_tilt

export get_tilt

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

function getfinetilt(idiff)
    arrsize = size(idiff)
    res = findfirstharmonic2(idiff .^ 2; visualdebug=false)
    return that =
        angle(res[2][end]) .+
        [2π * (i * res[1][1] + j * res[1][2]) .- π for i in 1:arrsize[1], j in 1:arrsize[2]]
end

abstract type AbstractAlg end
# Algorithms for tilt extraction
abstract type TiltExtractionAlg <: AbstractAlg end
struct RoughTilts <: TiltExtractionAlg
    erasesize::Int
    selectsize::Int
end
RoughTilts() = RoughTilts(2, 2)
struct FineTilts <: TiltExtractionAlg end

function (alg::RoughTilts)(idiffs)
    tilts = [get_tilt(id, alg.erasesize, alg.selectsize) for id in idiffs]
    # get signs of the restored tilts
    s = sign.([phwrap.(diff(dd; dims=2))[1] for dd in tilts])
    return tilts .*= s
end
(::FineTilts)(idiffs) = [getfinetilt(id) for id in idiffs]

# Algorithms for Phase-Shifting Interferometry (PSI)
abstract type PSIAlg <: AbstractAlg end
struct diffPSI <: PSIAlg end
struct LSPSI <: PSIAlg end

(::diffPSI)(images, deltas) = get_diff_phase_from_n_psi(images, deltas)
(::LSPSI)(images, deltas) =
    if length(images) == length(deltas)
        get_LS_phase_from_n_psi(images, deltas)[1]
    else
        @info "Assuming the first delta is zero"
        fulldeltas = [[zero(deltas[1])]; deltas]
        get_LS_phase_from_n_psi(images, fulldeltas)[1]
    end

export RoughTilts, FineTilts, diffPSI, LSPSI

function get_phase_from_igrams_with_tilts(
    igramsF, tiltsmethod::TiltExtractionAlg, psimethod::PSIAlg
)
    idiffs = diffirst(igramsF)
    # Extract  tilts with tilts method
    tilts = tiltsmethod(idiffs)
    # get signs of the restored tilts
    # s = sign.([phwrap.(diff(dd; dims=2))[1] for dd in tilts])
    # tilts .*= s

    phase = psimethod(igramsF, tilts)
    return phase
end

export get_phase_from_igrams_with_tilts

end
