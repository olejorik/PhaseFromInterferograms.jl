module FindHarmonics
# include("windowing.jl")
using ..Windowing
using FFTViews
using FFTW
using LinearAlgebra

"""
    eraseZerothOrder(spectrum, method = "block"; erasesize = 3)

Erase the zeroth order in the spectrum.
"""
function eraseZerothOrder!(spectrum; method = "block", erasesize = 3)
    if method == "block"
        aaa = FFTView(spectrum)
        I1 = one(CartesianIndex{ndims(aaa)})
        aaa[-erasesize* I1: erasesize *I1] .= 0
        return spectrum
    else
            error("method $method is not implemented")
    end
    
end  # function eraseZerothOrder

eraseZerothOrder(spectrum; args...) = eraseZerothOrder!(copy(spectrum);args...)


function findfirstharmonic(array, params)
    hr = findroughharmonic(array; params)
    hf = findfineharmonic(array, h, params)
    return hf
end

"""
    findroughharmonic(array, params)

Document this function
"""
function findroughharmonic(array; window_scale = 1/3, erasesize = 3, halfplane = "vertical")
    if halfplane == "vertical" halfplane = zeros(length(size(array))); halfplane[1] =1 end
    if halfplane == "horisontal" halfplane = zeros(length(size(array))); halfplane[2] =1 end
    w = GaussianWindow(Tuple(collect(size(array)) .* window_scale))
    s = fft(w(size(array)) .* array)
    eraseZerothOrder!(s, erasesize = erasesize)
    posmax = argmax(FFTView(abs.(s)))
    @show posmax
    rough_freq = [FFTView(fftfreq(size(array)[i]))[posmax[i]] for i in 1:length(size(array))]
    signature = Int(sign(dot(halfplane, rough_freq))) # Select "positive tilt")
    posmax *= signature
    rough_freq *= signature
    @show posmax

    complex_amplitude = FFTView(s)[posmax]/sum(w(size(array)))
    
    return posmax, rough_freq, complex_amplitude  
    
    
end  # function findroughharmonic
    
"""
    findfineharmonic(array, h, params)

Document this function
"""
function findfineharmonic(array, h, params)
    r = setrangesCST(h, array, params)
    fftX = FFTZoom(r)
    
end  # function findfineharmonic

export eraseZerothOrder, eraseZerothOrder!, findfirstharmonic
end # module