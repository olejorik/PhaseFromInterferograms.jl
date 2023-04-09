module PhaseFromInterferograms
using FFTW
using FFTViews

# Write your package code here.
include("Windowing.jl")
using .Windowing
export GaussianWindow

include("FindHarmonics.jl")
using .FindHarmonics
export eraseZerothOrder, eraseZerothOrder!, findfirstharmonic


"""
    get_tilt(idiff)

Get the tilt from the interferogram difference
"""
function get_tilt(idiff, erasesize=3, selectsize=3, debug = false)
    spectrum = fft(idiff .^ 2)
    eraseZerothOrder!(spectrum) 
    aaa = FFTView(spectrum)
    bbb= zeros(ComplexF32, size(spectrum))
    bbbFV = FFTView(bbb)

    iii = argmax(abs.(spectrum))
    if debug
        println("max position $iii")
        println("max value is $(spectrum[iii])\n")
    end

    # select 

    bbox=selectsize
    bbbFV[iii - CartesianIndex(bbox,bbox) : iii + CartesianIndex(bbox,bbox)] .= aaa[iii - CartesianIndex(bbox,bbox) : iii + CartesianIndex(bbox,bbox)]

    delta = angle.(ifft(- bbb)) # add π
    return delta
    
end  # function get_tilt

export get_tilt


end
