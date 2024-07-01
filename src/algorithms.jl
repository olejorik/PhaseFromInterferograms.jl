# Abstract and concrete types of the algorithms used in the package
# #####################



abstract type AbstractAlg end


# Algorithms for tilt extraction
# ##############################

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
# ##################################################


abstract type PSIAlg <: AbstractAlg end
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
