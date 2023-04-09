module Windowing

abstract type Window{T<:Union{Real,Tuple}} end

(w::Window)(dims::Tuple{Vararg{Int}}) = error("Windowing for type $(typeof(w)) is not defined")
    
"""
    GaussianWindow(w) is an object corresponding to a centerd gaussian window  `w` pixels wide.

    `w` can be a single number or a tuple, with specified width for each dimension.

"""
struct GaussianWindow{T} <: Window{T}
    width::T
end

_grange(len,wid) = range(-(len-1)/wid,(len-1)/wid,step = 2/wid)
_gwidth(w::GaussianWindow{<:Real}, dims) = fill(w.width, length(dims)) 
_gwidth(w::GaussianWindow{<:Tuple}, dims) = length(w.width) == length(dims) ? w.width : error("Incompatible dimensions of the window and the array")

function (w::GaussianWindow)(dims::Tuple{Vararg{Int}})
    a = zeros(dims) .+ 1
    widths = _gwidth(w, dims)
    for (id,d) in enumerate(dims)
        wd = _grange(d, widths[id])
        for  (is,s) in enumerate(eachslice(a, dims = id))
            s .*= exp(- wd[is]^2/2)
        end
    end
    return a
end

export GaussianWindow

end