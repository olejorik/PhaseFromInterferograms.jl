# Abstract and concrete types of the algorithms used in the package
# #####################



abstract type AbstractAlg end


# ##############################
# Algorithms for tilt extraction
# ##############################

"""
    TiltExtractionAlg <: AbstractAlg

Abstract algorithm for retrieval tilts from interferogram differences.

Concrete instances of the type are callable:

    `(alg::TiltExtractionAlg)(idiffs)` -> tilts

or can be used in [`gettilts`](@ref) function.

See also [`gettilts`](@ref).
"""
abstract type TiltExtractionAlg <: AbstractAlg end

"""
    gettilts(idiffs, alg::TiltExtractionAlg)

Extract tilts from the interferogram differences `idiffs` using `alg` method. Return tilts and additional information depending on the method used.
"""
gettilts(idiffs, alg::TiltExtractionAlg) = error("Method $(typeof(alg)) is not implemented")


(alg::TiltExtractionAlg)(idiffs) = gettilts(idiffs, alg)

"""
     RoughTilts <: TiltExtractionAlg

Simple method of extracting approximate tilts from an array of the interferogram differences based on the Fourier filtering.

Paremeters:
 - `erasesize` = 2 -- half-diameter of the central lobe to be excluded: block with size `2 erasesize + 1` centered around the origin
 - `selectsize` = 2 -- half diameter of the side lobe used to reconstruct the tilt: block with size `2 selectsize + 1` centerd around the maximum of the spectrum.
"""
@kwdef struct RoughTilts <: TiltExtractionAlg
    erasesize::Int = 2
    selectsize::Int = 2
end

"""
     FineTilts <: TiltExtractionAlg

Method of extracting tilts from an array of the interferogram differences based on locating maxima in the spectra with subpixel accuracy using zoomFFT (aka CZT or Bluestein algorithm).

Paremeters:
 - `erasesize` = 2 -- half-diameter of the central lobe to be excluded: block with size `2 erasesize + 1` centered around the origin
 - `cropsize` = 2 -- half diameter of the side lobe used tolocate the maximum: block with size `2 cropsize + 1` centerd around the maximum of the spectrum.
 - `zoomlevels` = nothing -- array of zoomlevels used for iterational subpixel approximation. E.g. `zoomlevels = [1,4]` will first find the side lobe maximum in the original Fourier spectrum , and then in the spectrum sampled with 4 times higher rate. If set to `nothing` uses automatic sequence of the zoom levels.
 - `normals = [nothing]` -- array of the normals defining the half plane where the maximum is located.
"""
@kwdef struct FineTilts <: TiltExtractionAlg
    zoomlevels = nothing
    erasesize = 2
    cropsize = 2
    normals = [nothing]
end



function gettilts(idiffs, alg::RoughTilts)
    tiltsandfreqs = [get_tilt(id, alg.erasesize, alg.selectsize) for id in idiffs]
    tilts = [tf[1] for tf in tiltsandfreqs]
    freqs = [tf[2] for tf in tiltsandfreqs]
    # get signs of the restored tilts
    s = [sign(dotproduct(f, [1, 1])) for f in freqs]
    tilts .*= s
    return tilts
end



function gettilts(idiffs, alg::FineTilts)
    tilts = [similar(id) for id in idiffs]
    taus = [zeros(2) for i in 1:length(idiffs)]
    sigmas = zeros(length(idiffs))
    arrsize = size(idiffs[1])
    if length(alg.normals) == 1
        normals = fill(alg.normals[1], length(idiffs))
    else
        normals = alg.normals
    end
    for (i, id) in enumerate(idiffs)
        tilt, τ, σ = getfinetilt(
            id;
            n=normals[i],
            zoomlevels=alg.zoomlevels,
            erasesize=alg.erasesize,
            cropsize=alg.cropsize,
            visualdebug=false,
        )
        tilts[i] .= tilt
        taus[i] .= τ
        sigmas[i] = σ
    end
    return tilts, taus, sigmas
end



# ##################################################
# Algorithms for Phase-Shifting Interferometry (PSI)
# ##################################################

"""
    PSIAlg <: AbstractAlg

An abstract algorithm for Phase-shifting Interferometry (PSI).

Concrete instances are callable or can be used in [`get_phase_from_n_psi`](@ref) functions:

    `(alg::PSIAlg)(measurements, deltas)`
extracts the phase `ϕ` from several `measurements` Iₖ = b + a cos(ϕ + δₖ).

`deltas` is an array of δₖ and it should be the same length as the number of `measurements` or one less. In this case, δ₁ is asuumed to be zero.

See also `get_phase_from_n_psi`](@ref).
"""
abstract type PSIAlg <: AbstractAlg end


"""
    get_phase_from_n_psi(images, deltas, alg::PSIAlg; kwargs...)

Retreive phase from several interferograms obtained with phase-shifting interferometry.
"""
get_phase_from_n_psi(images, deltas, alg::PSIAlg; kwargs...) =
    error("$(typeof(alg)) is not implemented")

function (alg::PSIAlg)(images, deltas)
    if length(images) == length(deltas)
        ret = get_phase_from_n_psi(images, deltas, alg)
    else
        @info "Assuming the first delta is zero"
        fulldeltas = [[zero(deltas[1])]; deltas]
        ret = get_phase_from_n_psi(images, fulldeltas, alg)
    end
    return ret
end

struct diffPSI <: PSIAlg end
struct LSPSI <: PSIAlg end

(::diffPSI)(images, deltas) = get_diff_phase_from_n_psi(images, deltas)

function (::LSPSI)(images, deltas; full=false)
    if length(images) == length(deltas)
        ret = get_LS_phase_from_n_psi(images, deltas)
    else
        @info "Assuming the first delta is zero"
        fulldeltas = [[zero(deltas[1])]; deltas]
        ret = get_LS_phase_from_n_psi(images, fulldeltas)
    end
    if full
        return ret
    else
        return ret[1]
    end
end




# ############################################
# Algorithms for phase extraction from several interferograms with unknown tilts
# ############################################


function get_phase_from_igrams_with_tilts(
    igramsF, tiltsmethod::TiltExtractionAlg, psimethod::PSIAlg
)
    idiffs = diffirst(igramsF)
    # Extract  tilts with tilts method
    tilthats, tauhats, sigmahats = gettilts(idiffs, tiltsmethod)
    phase = psimethod(igramsF, tilthats)
    return phase
end

function get_phase_from_igrams_with_tilts(
    igramsF, dirs::Vector{String}, tiltsmethod::TiltExtractionAlg, psimethod::PSIAlg
)
    idiffs = diffirst(igramsF)
    # Extract  tilts with tilts method
    tilts = first(tiltsmethod(idiffs))
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
    psimethod::PSIAlg;
    kwargs...,
)
    idiffs = diffirst(igramsF, ref)
    # Extract  tilts with tilts method
    tilts = tiltsmethod(idiffs)
    # get signs of the restored tilts
    s = getsign.(tilts, deleteat!(copy(dirs), ref); kwargs...)
    tilts .*= s

    # Change the order of the interferograms
    # phase = psimethod(reorder(igramsF, tilts, ref)...)
    # n = length(igramsF)
    # @assert ref <= n
    # perm = collect(1:n)
    # perm[1] = ref
    # perm[ref] = 1
    #
    # No, just insert zero tilt in the needed position
    phase = psimethod(igramsF, insert!(tilts, ref, zero(tilts[1])); kwargs...)
    #
    return phase
end
