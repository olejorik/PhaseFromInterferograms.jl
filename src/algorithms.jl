# Abstract and concrete types of the algorithms used in the package
# #####################



abstract type AbstractAlg end


# ##############################
# Algorithms for tilt extraction
# ##############################

abstract type TiltExtractionAlg <: AbstractAlg end

@kwdef struct RoughTilts <: TiltExtractionAlg
    erasesize::Int = 2
    selectsize::Int = 2
end

@kwdef struct FineTilts <: TiltExtractionAlg
    zoomlevels = nothing
    erasesize = 2
    cropsize = 2
    normals = [nothing]
end



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
            n=alg.normals[i],
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

`(alg::PSIAlg)(measurements, deltas)` extracts the phase `ϕ` from several `measurements` Iₖ = b + a cos(ϕ + δₖ).

`deltas` is an array of δₖ and it should be the same length as the number of `measurements` or one less. In this case, δ₁ is asuumed to be zero.

"""
abstract type PSIAlg <: AbstractAlg end

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

get_phase_from_n_psi(images, deltas, alg::PSIAlg; kwargs...) =
    error("$(typeof(alg)) is not implemented")

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
