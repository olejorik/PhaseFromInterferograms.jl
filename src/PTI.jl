module PTI
using StaticArrays
using LinearAlgebra: dot
using FFTW

import Base.zero

"""
    PTIestimate

Model for the interferograms obtained with in phase-tilted interferometry.
Each interferogram `Iₙ` is modeled as `Iₙ = a + 2 Re (c dₙ) a + c dₙ +̄c̄ ̄dₙ`, where `dₙ = exp(i δₙ)`.
`a` is the background, `c` is the complex amplitude, and `dₙ` is the linear in the image coordinates phase shift, `dₙ = exp(i δₙ)`, where `δₙ` is the phase tilt.

"""
struct PTIestimate{M,N,TT}
    framesize::NTuple{M,Int}
    setsize::NTuple{N,Int}
    background::Array{Float64,M}
    complexamplitude::Array{ComplexF64,M}
    mask::BitArray{M}
    tilts::Array{TT,N}
end

PTIestimate(igrams) = PTIestimate(
    size(igrams[1]),
    size(igrams),
    ones(Float64, size(igrams[1])),
    0.5 * ones(ComplexF64, size(igrams[1])),
    trues(size(igrams[1])),
    [TiltCentered(zeros(Float64, length(size(igrams[1])) + 1)) for _ in igrams],
)

background(p::PTIestimate) = p.background
complexamplitude(p::PTIestimate) = p.complexamplitude
mask(p::PTIestimate) = p.mask
tilts(p::PTIestimate) = p.tilts
framesize(p::PTIestimate) = p.framesize
setsize(p::PTIestimate) = p.setsize

setbackground!(p::PTIestimate, b) = (p.background .= b)
setcomplexamplitude!(p::PTIestimate, c) = (p.complexamplitude .= c)
setmask!(p::PTIestimate, m) = (p.mask .= m)
settilts!(p::PTIestimate, t) = (p.tilts .= t)

materialize(p::PTIestimate) = [
    background(p) .+
    2real(complexamplitude(p) .* cis.(PTI.materialize(tilt, framesize(p)))) for
    tilt in tilts(p)
]

"""
    TiltCentered


"""
struct TiltCentered{N}
    coefs::MVector{N,Float64}
end

TiltCentered(coefs::AbstractVector) = TiltCentered(MVector(coefs...))

sigma(t::TiltCentered) = t.coefs[1]
tau(t::TiltCentered) = t.coefs[2:end]
setsigma!(t::TiltCentered, s) = (t.coefs[1] = s)
settau!(t::TiltCentered, τ) = (t.coefs[2:end] .= τ)
setall!(t::TiltCentered, v) = (t.coefs .= v)
materialize(t::TiltCentered, dims) =
    [sigma(t) + dot(tau(t), x) for x in Iterators.product(fftshift.(fftfreq.(dims))...)]



end # module PTI
