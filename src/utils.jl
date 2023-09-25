
diffirst(v) = [t .- v[1] for t in v[2:end]]

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

dirtov = (o=[0, 0], h=[1, 0], v=[0, 1], d=[1, 1])

function getsign(tilt, n::Vector)
    return 2Int(dotproduct([tilt[1, 2] - tilt[1, 1], tilt[2, 1] - tilt[1, 1]], n) >= 0) - 1
end

getsign(tilt, dir::String) = getsign(tilt, dirtov[Symbol(dir)])
