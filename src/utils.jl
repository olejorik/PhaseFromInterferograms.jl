

diffirst(v) = [t .- v[1] for t in v[2:end]]
diffirst(v, i) = [t .- v[i] for t in v[vcat(1:(i - 1), (i + 1):end)]]

function get_direction_from_filename(s::String, base::String)
    bname = splitext(s)[1]
    if bname == base
        return "o"
    else
        return "$(bname[length(base) + 1])"
    end
end

function longest_common_subsequence(a::Vector)
    l = 0
    n = minimum(length, a)
    for i in 1:n
        allequal(getindex.(a, [1:i])) || break
        l = i
    end
    return a[1][1:l]
end

function get_tilt_dirs(fn::Vector{String})
    base = longest_common_subsequence(fn)
    return get_direction_from_filename.(fn, base)
end

dirtov = (o=[0, 0], v=[0, 1], h=[1, 0], d=[1, 1])
dirtov = (o=[0, 0], l=[0, -1], r=[0, 1], u=[1, 0], d=[-1, 0])


function getsign(tilt, n::Vector)
    # return 2Int(dotproduct([tilt[1, 2] - tilt[1, 1], tilt[2, 1] - tilt[1, 1]], n) >= 0) - 1
    s = getslopes(tilt)
    return 2Int((n' * [s...]) >= 0) - 1
end

getsign(tilt, dir::String) = getsign(tilt, dirtov[Symbol(dir)])

function get_aperture(igramsF, relative_threshold=1.3)
    i3d = stack(igramsF)
    variancei3d = var(i3d; dims=(3))[:, :, 1]
    h = fit(Histogram, variancei3d[:]; nbins=500)
    w = h.weights
    edges = h.edges[1]
    pos = find_first_minimum(w)
    tr = (edges[pos] + edges[pos + 1]) / 2
    tr = tr * relative_threshold # empirical
    return m = threshold!(ones(Bool, size(variancei3d)), variancei3d, tr)
end

function find_first_minimum(a)
    i = 1
    while a[i + 3] + a[i + 2] - a[i + 1] - a[i] <= 0 # derivative of a smoothed array
        i += 1
    end
    return i + 1
end

threshold!(mask, arr, t) = (mask[arr .< t] .= 0; mask)
threshold(arr, t) = threshold!(ones(size(arr)), arr, t)

getslopes(tilt, x, y) =
    (phwrap(tilt[x + 1, y] - tilt[x, y]), phwrap(tilt[x, y + 1] - tilt[x, y]))
function getslopes(tilt)
    x, y = Int.(round.(size(tilt) ./ 2))
    return getslopes(tilt, x, y)
end
