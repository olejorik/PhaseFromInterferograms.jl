module PhaseFromInterferograms
using FFTW
using FFTViews
using PhaseUtils: phwrap, dotproduct
using StatsBase

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
    tiltsandfreqs = [get_tilt(id, alg.erasesize, alg.selectsize) for id in idiffs]
    tilts = [tf[1] for tf in tiltsandfreqs]
    freqs = [tf[2] for tf in tiltsandfreqs]
    # get signs of the restored tilts
    s = [sign(dotproduct(f, [1, 1])) for f in freqs]
    tilts .*= s
    return tilts
end
(::FineTilts)(idiffs) = [getfinetilt(id)[1] for id in idiffs]

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

function reorder(images, deltas, ref)
    if length(images) == length(deltas)
        fulldeltas = deltas
    else
        @info "Assuming the first delta is zero"
        fulldeltas = [[zero(deltas[1])]; deltas]
        removefirst = true
    end
    n = length(images)
    @assert ref <= n
    perm = collect(1:n)
    perm[1] = ref
    perm[ref] = 1
    imnew = images[perm]
    deltasnew = fulldeltas[perm] .- [fulldeltas[ref]]
    if removefirst
        return imnew, deltasnew[2:end]
    else
        return imnew, deltasnew
    end
end

export RoughTilts, FineTilts, diffPSI, LSPSI

function get_phase_from_igrams_with_tilts(
    igramsF, tiltsmethod::TiltExtractionAlg, psimethod::PSIAlg
)
    idiffs = diffirst(igramsF)
    # Extract  tilts with tilts method
    tilts = tiltsmethod(idiffs)
    # get signs of the restored tilts
    s = sign.([phwrap.(diff(dd; dims=2))[1] for dd in tilts])
    tilts .*= s

    phase = psimethod(igramsF, tilts)
    return phase
end

function get_phase_from_igrams_with_tilts(
    igramsF, dirs::Vector{String}, tiltsmethod::TiltExtractionAlg, psimethod::PSIAlg
)
    idiffs = diffirst(igramsF)
    # Extract  tilts with tilts method
    tilts = tiltsmethod(idiffs)
    # get signs of the restored tilts
    s = getsign.(tilts, dirs[2:end])
    tilts .*= s
    phase = psimethod(igramsF, tilts)
    return phase
end

function get_phase_from_igrams_with_tilts(
    igramsF,
    ref::Integer,
    dirs::Vector{String},
    tiltsmethod::TiltExtractionAlg,
    psimethod::PSIAlg,
)
    idiffs = diffirst(igramsF)
    # Extract  tilts with tilts method
    tilts = tiltsmethod(idiffs)

    # Change the order of the interferograms

    phase = psimethod(reorder(igramsF, tilts, ref)...)
    return phase
end

export get_phase_from_igrams_with_tilts, get_tilt_dirs
export get_aperture

end
