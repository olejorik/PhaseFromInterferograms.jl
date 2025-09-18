module FindHarmonics
using PhaseUtils: GaussianWindow
using FFTViews
using FFTW
using LinearAlgebra

"""
    eraseZerothOrder(spectrum, method = "block"; erasesize = 3)

Erase the zeroth order in the spectrum. Currently impemented methods:
 - "block" -- set to zero a cube of width 2*`erasesize` +1 around the DC frequency.
"""
function eraseZerothOrder!(spectrum; method="block", erasesize=3)
    if method == "block"
        aaa = FFTView(spectrum)
        I1 = one(CartesianIndex{ndims(aaa)})
        aaa[(-erasesize * I1):(erasesize * I1)] .= 0
        return spectrum
    else
        error("method $method is not implemented")
    end

end  # function eraseZerothOrder

eraseZerothOrder(spectrum; args...) = eraseZerothOrder!(copy(spectrum); args...)


function findfirstharmonic(array, params)
    hr = findroughharmonic(array; params)
    hf = findfineharmonic(array, hr, params)
    return hf
end

"""
    findroughharmonic(array; window_scale=1/3, erasesize=3, halfplane="vertical")

Find the rough position and frequency of the first harmonic in interferometric data.

This function performs a coarse search for the dominant spatial frequency component
in an interferogram by analyzing its Fourier spectrum. It's typically used as the
first step in extracting phase information from interference patterns.

# Arguments
- `array`: Input interferogram data (2D or N-D array)
- `window_scale=1/3`: Scale factor for the Gaussian window size relative to array dimensions.
  Smaller values create sharper windows with better frequency localization.
- `erasesize=3`: Half-width of the region around DC frequency to zero out (removes
  background illumination effects)
- `halfplane="vertical"`: Constraint on the search direction:
  - `"vertical"`: Force frequency to have positive component in first dimension
  - `"horisontal"`: Force frequency to have positive component in second dimension
  - `"none"`: No directional constraint

# Returns
A tuple `(posmax, rough_freq, complex_amplitude)` where:
- `posmax`: CartesianIndex of the peak position in the FFT spectrum
- `rough_freq`: Vector of spatial frequencies [fx, fy, ...] in cycles per pixel
- `complex_amplitude`: Complex amplitude of the harmonic (normalized by window sum)

# Algorithm
1. Apply Gaussian windowing to reduce spectral leakage
2. Compute FFT of windowed data
3. Remove DC component to avoid interference from background illumination
4. Find the maximum peak in the magnitude spectrum
5. Convert peak position to spatial frequency coordinates
6. Apply halfplane constraint if specified to resolve directional ambiguity
7. Extract and normalize the complex amplitude

# Notes
The halfplane constraint is useful when the fringe orientation is known a priori,
helping to avoid sign ambiguity in the detected frequency. This is common in
interferometric setups where the tilt direction is controlled.

# Example
```julia
# Find rough harmonic in a 2D interferogram
pos, freq, amp = findroughharmonic(interferogram; window_scale=0.25, halfplane="vertical")
println("Detected frequency: ", freq, " cycles/pixel")
```
"""
function findroughharmonic(array; window_scale=1 / 3, erasesize=3, halfplane="vertical")
    if halfplane == "vertical"
        halfplane = zeros(length(size(array)))
        halfplane[1] = 1
    elseif halfplane == "horisontal"
        halfplane = zeros(length(size(array)))
        halfplane[2] = 1
    end
    w = GaussianWindow(Tuple(collect(size(array)) .* window_scale))
    s = fft(w(size(array)) .* array)
    eraseZerothOrder!(s; erasesize=erasesize)
    posmax = argmax(FFTView(abs.(s)))
    # @show posmax
    rough_freq = [
        FFTView(fftfreq(size(array)[i]))[posmax[i]] for i in 1:length(size(array))
    ]
    if halfplane != "none"
        signature = Int(sign(dot(halfplane, rough_freq))) # Select "positive tilt")
        posmax *= signature
        rough_freq *= signature
    end
    # @show posmax

    complex_amplitude = FFTView(s)[posmax] / sum(w(size(array)))

    return posmax, rough_freq, complex_amplitude


end  # function findroughharmonic

"""
    findfineharmonic(array, h, params)

Document this function
"""
function findfineharmonic(array, h, params)
    r = setrangesCST(h, array, params)
    return fftX = FFTZoom(r)

end  # function findfineharmonic

export eraseZerothOrder, eraseZerothOrder!, findfirstharmonic
end # module
