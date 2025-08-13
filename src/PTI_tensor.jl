using StaticArrays
using LinearAlgebra: dot
using FFTW
using PhaseUtils:
    Tilt,
    TiltCentered,
    FreeTilt,
    sigma,
    tau,
    setsigma!,
    settau!,
    setall!,
    materialize,
    apply
using PhaseUtils: ArrayAxes, FourierAxes, DataAxes, DataAxesCentered

import Base.zero

"""
    PTIestimate{M,TT,TFA,TSA}

Comprehensive model and state container for Phase-Tilted Interferometry (PTI) analysis.

# Physical Model
Each interferogram Iₙ is modeled as:
Iₙ = a + 2Re(c·dₙ) = a + c·dₙ + c̄·d̄ₙ
where:
- a is the background intensity (shared across all interferograms)
- c is the complex amplitude (contains the phase to be estimated)
- dₙ = exp(iδₙ) is the tilt factor for interferogram n
- δₙ(x, y) = σₙ + τₙ¹ x + τₙ² y is a linear phase function

# Structure Fields
- `framesize::Int`: Number of spatial dimensions
- `setsize::Int`: Number of interferogram set dimensions
- `fullsize::NTuple{M,Int}`: Complete dimensions (frame + set)
- `background::Array{Float64,M}`: Background intensity component
- `complexamplitude::Array{ComplexF64,M}`: Complex amplitude containing phase
- `mask::BitArray{M}`: Binary mask for valid pixels (shared across interferograms)
- `tilts::Array{TT,M}`: Array of tilt objects for each interferogram
- `frameaxes::TFA`: Coordinate system for frame dimensions
- `setaxes::TSA`: Coordinate system for set dimensions
- `igrams::Array{Float64,M}`: Forward model output (synthesized interferograms)
- `data::Union{Nothing,Array{Float64,M}}`: Measured interferograms (optional)
- `insync::Vector{Bool}`: Flag indicating if `igrams` is synchronized with current parameters

# Usage
PTIestimate serves as:
1. A problem formulation for PTI analysis
2. A state container for iterative algorithms
3. A unified interface for both iterative and non-iterative algorithms
4. A forward model for generating synthetic interferograms

# Example
```julia
## Create a PTIestimate for a 2D problem with 4 interferograms (forward modeling)
pti = PTIestimate((256, 256), (4,))  # data = nothing

## Create a PTIestimate from measured data (inverse problem)
pti = PTIestimate(measured_interferograms)  # data = measured_interferograms

## Run parameter estimation (requires data)
update_background_amplitude!(pti, LSPhaseAlg())
update_tilts!(pti, GradientDescentTilt())

## Extract the estimated phase
phase_estimate = getphase(pti)
```

See also: `initialize!`, `update_background_amplitude!`, `update_tilts!`, `materialize!`
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
    data::Union{Nothing,Array{Float64,M}}
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
        nothing,  # No measured data for forward modeling
        [true],
    )
end

function PTIestimate(
    igrams::Union{Array{T} where {T<:Array},Slices};
    frameaxes=DataAxesCentered(),
    setaxes=DataAxes(),
)
    # Create the basic structure
    framesize = size(igrams[1])
    setsize = size(igrams)
    K = length(framesize)
    M = length(setsize)
    fullsize = (framesize..., setsize...)
    frameax = frameaxes(framesize)
    setax = setaxes(setsize)

    # Convert igrams to array format for data storage
    data_array = Array{Float64}(undef, fullsize...)
    for (i, igram) in pairs(IndexCartesian(), igrams)
        # Use linear indexing for the last dimensions
        indices = (Colon() for _ in 1:K)..., Tuple(i)...
        data_array[indices...] = igram
    end

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
        copy(data_array),  # Initialize igrams with copy of data
        copy(data_array),  # Store measured data
        [false],  # igrams need to be updated
    )
end

#  Interfaces

background(p::PTIestimate) = p.background
complexamplitude(p::PTIestimate) = p.complexamplitude
mask(p::PTIestimate) = p.mask
tilts(p::PTIestimate) = p.tilts
data(p::PTIestimate) = p.data
framesize(p::PTIestimate) = p.fullsize[1:(p.framesize)]
setsize(p::PTIestimate) = p.fullsize[(p.framesize + 1):(p.setsize + p.framesize)]

hasdata(p::PTIestimate) = p.data !== nothing
function requiredata(p::PTIestimate)
    hasdata(p) ||
        error("This operation requires measured data. Use setdata!(p, data) first.")
    return data(p)
end
function requiredatasliced(p::PTIestimate)
    hasdata(p) ||
        error("This operation requires measured data. Use setdata!(p, data) first.")
    return eachslice(data(p); dims=1:(p.framesize))
end

getphase(p::PTIestimate) = reshape(angle.(complexamplitude(p)), framesize(p))
getigrams(p::PTIestimate) = p.insync[1] ? p.igrams : materialize!(p).igrams
getigramssliced(p::PTIestimate) = eachslice(getigrams(p); dims=setdims(p))


setbackground!(p::PTIestimate, b) = (p.insync .= false; p.background .= b)
setcomplexamplitude!(p::PTIestimate, c) = (p.insync .= false; p.complexamplitude .= c)
setmask!(p::PTIestimate, m) = (p.insync .= false; p.mask .= m)
settilts!(p::PTIestimate, t) = (p.insync .= false; p.tilts .= t)
setdata!(p::PTIestimate, d) = (p.data = copy(d))
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

"""
    update_igrams!(p::PTIestimate)

TBW
"""
function update_igrams!(p::PTIestimate)
    return p.igrams .= materialize.(p.tilts, (p.frameaxes))
end

framedims(p::PTIestimate) = Tuple(i for i in 1:(p.framesize))
setdims(p::PTIestimate) = Tuple((i + p.framesize) for i in 1:(p.setsize))


# Main functions
function initialize!(p::PTIestimate, data, alg=FFTcrop1(); refframe=1)
    data = eachslice(data; dims=3)
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

# Convenience method that uses internal data - more specific signature
function initialize!(p::PTIestimate; refframe::Int=1)
    return initialize!(p, requiredata(p), FFTcrop1(); refframe=refframe)
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

# Convenience method that uses internal data - more specific signature
function update_background_amplitude!(p::PTIestimate)
    return update_background_amplitude!(p, requiredata(p), LSPhaseAlg())
end

function update_tilts!(p::PTIestimate, igrams, alg)
    tiltguess = (alg)(igrams, background(p), complexamplitude(p), p.frameaxes)
    for (tp, tg) in zip(p.tilts, tiltguess)
        setall!(tp, tg)
    end
end

# Convenience method that uses internal data - more specific signature
function update_tilts!(p::PTIestimate)
    return update_tilts!(p, requiredata(p), SymmetricLS())
end



## Tilt and axes types now live in PhaseUtils; keep using them via imports above

# end # module PTI
