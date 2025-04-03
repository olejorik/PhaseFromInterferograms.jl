
using StaticArrays
using LinearAlgebra: dot
using FFTW
using PhaseUtils: diffirst

import Base.zero

"""
    PTIestimate

Model for the interferograms obtained with in phase-tilted interferometry.
Each interferogram `Iₙ` is modeled as `Iₙ = a + 2 Re (c dₙ) a + c dₙ +̄c̄ ̄dₙ`, where `dₙ = exp(i δₙ)`.
`a` is the background, `c` is the complex amplitude, and `dₙ` is the linear in the image coordinates phase shift, `dₙ = exp(i δₙ)`, where `δₙ` is the phase tilt.

mask is the same for all the interferograms.

"""
struct PTIestimate{M,TT,TFA,TSA}
    framesize::Int
    setsize::Int
    fullsize::NTuple{M,Int}
    background::Array{Float64,M}
    complexamplitude::Array{ComplexF64,M}
    mask::BitArray{M}
    tilts::Array{TT,M}
    frameaxes::TFA
    setaxes::TSA
    igrams::Array{Float64,M}
    insync::Vector{Bool}
end

function PTIestimate(
    framesize::NTuple, setsize::NTuple; frameaxes=DataAxesCentered(), setaxes=DataAxes()
)
    K = length(framesize)
    M = length(setsize)
    fullsize = (framesize..., setsize...)
    frameax = frameaxes(framesize)
    setax = setaxes(setsize)
    return PTIestimate(
        K,
        M,
        fullsize,
        ones(Float64, framesize..., fill(1, M)...),
        0.5 * ones(ComplexF64, framesize..., fill(1, M)...),
        trues(framesize..., fill(1, M)...),
        reshape(
            [FreeTilt(zeros(Float64, K + 1)) for _ in CartesianIndices(setsize)],
            fill(1, K)...,
            setsize...,
        ),
        frameax,
        setax,
        2 * ones(Float64, framesize..., setsize...),
        [true],
    )
end

PTIestimate(igrams::Union{Array{T} where {T<:Array},Slices}; axes...) =
    PTIestimate(size(igrams[1]), size(igrams); axes...)

#  Interfaces

background(p::PTIestimate) = p.background
complexamplitude(p::PTIestimate) = p.complexamplitude
mask(p::PTIestimate) = p.mask
tilts(p::PTIestimate) = p.tilts
framesize(p::PTIestimate) = p.fullsize[1:(p.framesize)]
setsize(p::PTIestimate) = p.fullsize[(p.framesize + 1):(p.setsize + p.framesize)]

getphase(p::PTIestimate) = reshape(angle.(complexamplitude(p)), framesize(p))
getigrams(p::PTIestimate) = p.insync[1] ? p.igrams : materialize!(p).igrams
getigramssliced(p::PTIestimate) = eachslice(getigrams(p); dims=setdims(p))


setbackground!(p::PTIestimate, b) = (p.insync .= false; p.background .= b)
setcomplexamplitude!(p::PTIestimate, c) = (p.insync .= false; p.complexamplitude .= c)
setmask!(p::PTIestimate, m) = (p.insync .= false; p.mask .= m)
settilts!(p::PTIestimate, t) = (p.insync .= false; p.tilts .= t)
setphase!(p, φ) = (setcomplexamplitude!(p, abs.(complexamplitude(p)) .* cis.(φ)))


function materialize!(p::PTIestimate)
    coords = Iterators.product(p.frameaxes...)

    p.igrams .= background(p) .+ real.(complexamplitude(p) .* cis.(apply.(p.tilts, coords)))
    p.igrams .*= mask(p)
    p.insync .= true
    return p
end

"""
    get_diversed_complex_amplitude(p::PTIestimate, dims...)

Get the complex amplitude of the interferograms with the tilts (component speciefed by `dims` applied).
"""
function get_diversed_complex_amplitude(p::PTIestimate, dims...)
    coords = Iterators.product(p.frameaxes...)
    return complexamplitude(p) .* cis.(apply.(p.tilts, coords, (dims...,)))
end

function get_single_diversed_complex_amplitude!(arr, p::PTIestimate, i, dims...)
    coords = Iterators.product(p.frameaxes...)
    aaa = apply.((p.tilts[i],), coords, (dims...,))
    @show size(aaa)
    @show size(complexamplitude(p))
    arr .= complexamplitude(p) .* cis.(apply.((p.tilts[i],), coords, (dims...,)))
    return arr
end

function update_igrams!(p::PTIestimate)
    return p.igrams .= materialize.(p.tilts, (p.axes[1:(p.framesize)],))
end

framedims(p::PTIestimate) = Tuple(i for i in 1:(p.framesize))
setdims(p::PTIestimate) = Tuple((i + p.framesize) for i in 1:(p.setsize))


# Main functions
function initialize!(p::PTIestimate, data, alg=FFTcrop1(); refframe=1)

    for (i, igram) in pairs(IndexCartesian(), data)
        if i == CartesianIndex(refframe)
            tiltguess = FreeTilt([0.0, 0.0, 0.0])
        else
            idiffsq = (igram - data[refframe]) .^ 2
            pos, freq, amp = get_side_lobe_freq(idiffsq, alg)
            tiltguess = FreeTilt([-π - angle(amp), (2π * freq)...])
        end
        p.tilts[i] = tiltguess
    end
    p.insync .= false
    return p
end

function set_tilt_signs!(p::PTIestimate, normals)
    for (tp, n) in zip(p.tilts, normals)
        if dot(tau(tp), n) < 0
            tp.coefs .*= -1
        end
    end
end


function update_background_amplitude!(p::PTIestimate, data, alg)
    a, c = alg(data, p.tilts)
    setbackground!(p, a)
    return setcomplexamplitude!(p, c)
end

function update_tilts!(p::PTIestimate, igrams, alg)
    tiltguess = (alg)(igrams, background(p), complexamplitude(p), p.frameaxes)
    for (tp, tg) in zip(p.tilts, tiltguess)
        setall!(tp, tg)
    end

end



abstract type Tilt end


sigma(t::Tilt) = t.coefs[1]
tau(t::Tilt) = t.coefs[2:end]
tau(t::Tilt, j) = t.coefs[1 + j]
setsigma!(t::Tilt, s) = (t.coefs[1] = s)
settau!(t::Tilt, τ) = (t.coefs[2:end] .= τ)
setall!(t::Tilt, v) = (t.coefs .= v)
materialize(t::Tilt, dims) =
    [sigma(t) + dot(tau(t), x) for x in Iterators.product(fftshift.(fftfreq.(dims))...)]

materialize(t::Tilt, axes::Vector{T} where {T<:AbstractVector}) =
    [sigma(t) + dot(tau(t), x) for x in Iterators.product(axes...)]



"""
    TiltCentered


"""
struct TiltCentered{N} <: Tilt
    coefs::MVector{N,Float64}
end

TiltCentered(coefs::AbstractVector) = TiltCentered(MVector(coefs...))

struct FreeTilt{N} <: Tilt
    coefs::MVector{N,Float64}
end

FreeTilt(coefs::AbstractVector) = FreeTilt(MVector(coefs...))

using LinearAlgebra: dot
apply(t::Tilt, x) = sigma(t) + dot(tau(t), x)
apply(t::Tilt, x, dims) = sum(t.coefs[i + 1] * get(x, i, 1) for i in dims)


## Fourier transform coordinates
abstract type ArrayAxes end
struct FourierAxes <: ArrayAxes end
struct DataAxes <: ArrayAxes end
struct DataAxesCentered <: ArrayAxes end


(alg::FourierAxes)(dims::NTuple) = [fftshift(fftfreq(d)) for d in dims]
(alg::DataAxes)(dims::NTuple) = [1:d for d in dims]
(alg::DataAxesCentered)(dims::NTuple) = [fftshift(fftfreq(d, d)) for d in dims]

# end # module PTI
