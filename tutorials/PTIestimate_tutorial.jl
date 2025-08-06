# # PTIestimate Structure Tutorial
#
# ## Introduction
#
# This tutorial demonstrates the use of the [`PTIestimate`](@ref) structure for Phase-Tilted Interferometry (PTI) analysis. The PTIestimate structure serves as a comprehensive model representation that can be used for:
# 1. Problem formulation for PTI analysis
# 2. State container for iterative algorithms
# 3. Forward model for generating synthetic interferograms
# 4. Unified interface for both iterative and non-iterative algorithms

using PhaseFromInterferograms
const PTI = PhaseFromInterferograms
using Statistics
using CairoMakie
using PhasePlots
using PhaseUtils
using ImageFiltering

CairoMakie.activate!(; type="png")

# Define colormaps for consistent visualization
coligram = :linear_blue_95_50_c20_n256
coligramdiff = :diverging_gwr_55_95_c38_n256;

# ## 1. Creating a PTIestimate from Interferograms
#
# We start by creating a set of synthetic interferograms and then constructing a PTIestimate object.

# Create synthetic interferogram data
igrams = [rand(500, 800) for _ in 1:15];

# Create the PTIestimate object with Fourier coordinate system
pti_est = PTIestimate(igrams; frameaxes=PTI.FourierAxes())
@show typeof(pti_est)

# Examine the coordinate system
@show pti_est.frameaxes
coords = Iterators.product(pti_est.frameaxes...)
@show size(coords)

# ## 2. Setting the Aperture Mask
#
# The mask defines the valid region for phase analysis. Here we create a circular aperture in Fourier space.


# Set circular mask in Fourier space
r = 0.45
PTI.setmask!(pti_est, map(x -> x[1]^2 + x[2]^2 <= r^2, coords));

# Visualize the mask across all frames
plot_heatmaps_table(eachslice(PTI.mask(pti_est); dims=3); aspect=AxisAspect(1))

# ## 3. Understanding the Forward Model
#
# Initially, the phase is zero, so all interferograms are identical to the background.

# Display initial interferograms (all identical to background)
plot_heatmaps_table(
    eachslice(PTI.getigrams(pti_est); dims=3); colormap=coligram, aspect=AxisAspect(1)
)

# Show the initial (zero) phase
showphase(PTI.getphase(pti_est))[1]

# ## 4. Setting a Complex Phase Pattern
#
# We create a phase composed of astigmatism and coma terms to demonstrate realistic wavefront aberrations.

# Create a complex phase pattern: astigmatism + coma
complex_phase = 4π * map(x -> -prod(x) / r + 3x[1]^2 * x[2] - x[2]^3, coords) ./ r^3
PTI.setphase!(pti_est, complex_phase);

# Display the phase pattern
showphase(PTI.getphase(pti_est))[1]

# Show how interferograms change with the new phase
plot_heatmaps_table(
    eachslice(PTI.getigrams(pti_est); dims=3); colormap=coligram, aspect=AxisAspect(1)
)

# ## 5. Adding Random Tilts
#
# Tilts represent linear phase ramps that are different for each interferogram. This simulates realistic experimental conditions.

# Add random tilts to each interferogram
PTI.settilts!(pti_est, PTI.FreeTilt.([6π * (randn(3) .- 0.5) for _ in PTI.tilts(pti_est)]));

# Display interferograms with tilts applied
plot_heatmaps_table(
    eachslice(PTI.getigrams(pti_est); dims=3); colormap=coligram, aspect=AxisAspect(1)
)

# Update tilts to demonstrate variability
PTI.settilts!(pti_est, PTI.FreeTilt.([10π * (randn(3) .- 0.5) for _ in PTI.tilts(pti_est)]))
plot_heatmaps_table(
    eachslice(PTI.getigrams(pti_est); dims=3); colormap=coligram, aspect=AxisAspect(1)
) |> display

print("Tilts: ", pti_est.tilts)

# ## 6. Analyzing Interferogram Differences
#
# Interferogram differences are crucial for tilt estimation. We examine the structure of these differences.

# Extract interferograms for analysis
igrams_analysis = PTI.getigrams(pti_est);

# Helper function for plotting heatmap tables
pht(arr; kwargs...) =
    plot_heatmaps_table(eachslice(arr; dims=3); aspect=AxisAspect(1), kwargs...)

# Display the interferograms
pht(igrams_analysis; colormap=:grays)

# Compute and display interferogram differences
deltas = diff(igrams_analysis; dims=3)
pht(deltas; colormap=coligramdiff)

# ## 7. Background Estimation
#
# The background can be estimated by averaging all interferograms, assuming random tilts.

# Simple average for background estimation
aest = sum(igrams_analysis; dims=3) / prod(PTI.setsize(pti_est))
pht(aest; colormap=:grays)

# Improved background estimation with Gaussian filtering
aestf = [imfilter(aest[:, :, 1], Kernel.gaussian(s)) for s in 1:16]
plot_heatmaps_table(
    aestf; colormap=:grays, titles=["σ = $s" for s in 1:16], aspect=AxisAspect(1)
)

# ## 8. Working with Diversed Complex Amplitudes
#
# The diversed complex amplitude represents how the complex amplitude appears under specific tilt conditions.

# Get complex amplitude with no tilt components
ttt = PTI.get_diversed_complex_amplitude(pti_est, (0,))
pht(real.(ttt); colormap=coligram)

# Get complex amplitude with x and y tilt components
ttt = PTI.get_diversed_complex_amplitude(pti_est, (1, 2))
pht(real.(ttt); colormap=coligram)

# ## 9. Phase Reconstruction from Individual Interferograms
#
# This section demonstrates how to reconstruct phase information from individual interferograms using linear algebra.

# Set up linear system for the first interferogram
A1 = [vec(ttt[:, :, 1]) conj.(vec(ttt[:, :, 1]))]
b1 = vec(igrams_analysis[:, :, 1]) - vec(PTI.background(pti_est)[:, :, 1])
d1 = A1 \ b1;

# Compare estimated and ground truth phase offset
estimated_sigma = angle(d1[1])
ground_truth_sigma = PTI.sigma(pti_est.tilts[1])

@show estimated_sigma
@show ground_truth_sigma
@show phwrap(estimated_sigma - ground_truth_sigma)  # Error should be small

# ## 10. Improved Reconstruction Using Mask
#
# Limiting the reconstruction to valid pixels (within the mask) improves accuracy significantly.

# Use only masked pixels for reconstruction
ttt1 = ttt[:, :, 1][pti_est.mask]
A1_masked = [ttt1 conj(ttt1)]
b1_masked =
    igrams_analysis[:, :, 1][pti_est.mask] - PTI.background(pti_est)[:, :, 1][pti_est.mask]
d1_masked = A1_masked \ b1_masked;

# Check improved accuracy
improved_sigma = angle(d1_masked[1])
@show improved_sigma
@show phwrap(improved_sigma - ground_truth_sigma)  # Error should be near machine precision

# ## 11. PTI Workflow: Initialization and Iterative Refinement
#
# This section demonstrates the complete PTI workflow using the PTIestimate structure.

# Create a new estimate for testing the workflow
test_data = eachslice(igrams_analysis; dims=3)
est = PTIestimate(test_data);

# Initialize with rough tilt estimates
refframe = 2
PTI.initialize!(est, test_data; refframe=refframe);

# Display initial tilt estimates
@show est.tilts;

# Show interferograms after initialization
plot_heatmaps_table(
    PTI.getigramssliced(est); colormap=coligram, title="After initialization"
)

# Adjust tilt signs based on ground truth (this would normally be done using other methods)
# For this tutorial, we'll use the ground truth to set proper signs for convergence demonstration
normals = [
    PTI.tau(t) - PTI.tau(pti_est.tilts[1, 1, refframe]) for t in pti_est.tilts[1, 1, :]
]
PTI.set_tilt_signs!(est, normals)

# ## 12. Phase and Background Estimation
#
# Use phase-shifting interferometry to estimate phase and background.

# Apply least-squares PSI algorithm
psialg = PTI.LSPSI()
coords_analysis = Iterators.product(est.frameaxes...)
deltas_est = eachslice(PTI.apply.(est.tilts, coords_analysis); dims=PTI.setdims(est))
PhaseBgAmp = psialg(test_data, deltas_est; full=true);

# Display results
fig = Figure(; size=(1200, 400))
fig[1, 1] = Axis(fig; title="Estimated Phase")
fig[1, 2] = Axis(fig; title="Estimated Background")
fig[1, 3] = Axis(fig; title="Estimated Amplitude")
showphase!(fig[1, 1], PhaseBgAmp[1])
showarray!(fig[1, 2], PhaseBgAmp[3])
showarray!(fig[1, 3], abs.(PhaseBgAmp[2]))
fig

# Update the estimate with new values
PTI.setcomplexamplitude!(est, PhaseBgAmp[2])
PTI.setbackground!(est, PhaseBgAmp[3]);

# ## 13. Iterative Refinement
#
# Demonstrate the iterative refinement process that alternates between phase/background estimation and tilt refinement.

# Exact tilt estimation functions from the original PTIestimate_test.jl
function get_taux(qqq, n)
    ## n is the index of the tilt
    ## we fix first coordinate of the frame, and compose the matrix of the system of equations by iterating by the second coordinate
    A1 = zeros(ComplexF64, PTI.framesize(qqq)[2], 2)
    b1 = zeros(ComplexF64, PTI.framesize(qqq)[2])
    alltau = [get_taux(qqq, n, mx, A1, b1) for mx in 1:PTI.framesize(qqq)[1]]
    ## function gettau(qqq, n, mx) returns the x tau component of the tilt n at the point mx or NaN if the tilt is not defined at the point
    ## Now we extract the slope taux from alltau
    t1est = phwrap(diff(alltau[(!isnan).(alltau)]))
    return mean(t1est) / step(qqq.frameaxes[1])
end

function get_taux(qqq, n, mx, A1, b1)
    ## n is the index of the tilt
    ## mx is the index of the point in the first coordinate of the frame
    igrams = PTI.getigrams(qqq)
    for my in 1:PTI.framesize(qqq)[2]
        A1[my, 1] = PTI.complexamplitude(qqq)[mx, my] * qqq.mask[mx, my]
        A1[my, 2] = conj(A1[my, 1])
        b1[my] = (igrams[mx, my, n] - PTI.background(qqq)[mx, my]) * qqq.mask[mx, my]
    end
    if sum(b1 .!= 0) > 2
        d1 = A1 \ b1
        return angle(d1[1])
    else
        return NaN
    end
end

function get_tauy(qqq, n)
    ## n is the index of the tilt
    ## we fix second coordinate of the frame, and compose the matrix of the system of equations by iterating by the first coordinate
    A1 = zeros(ComplexF64, PTI.framesize(qqq)[1], 2)
    b1 = zeros(ComplexF64, PTI.framesize(qqq)[1])
    alltau = [get_tauy(qqq, n, my, A1, b1) for my in 1:PTI.framesize(qqq)[2]]
    ## function gettau(qqq, n, mx) returns the x tau component of the tilt n at the point mx or NaN if the tilt is not defined at the point
    ## Now we extract the slope taux from alltau
    t1est = phwrap(diff(alltau[(!isnan).(alltau)]))
    return mean(t1est) / step(qqq.frameaxes[2])
end

function get_tauy(qqq, n, my, A1, b1)
    ## n is the index of the tilt
    ## my is the index of the point in the second coordinate of the frame
    igrams = PTI.getigrams(qqq)
    for mx in 1:PTI.framesize(qqq)[1]
        A1[mx, 1] = PTI.complexamplitude(qqq)[mx, my] * qqq.mask[mx, my]
        A1[mx, 2] = conj(A1[mx, 1])
        b1[mx] = (igrams[mx, my, n] - PTI.background(qqq)[mx, my]) * qqq.mask[mx, my]
    end
    if sum(b1 .!= 0) > 2
        d1 = A1 \ b1
        return angle(d1[1])
    else
        return NaN
    end
end

function get_sigma(qqq, n)
    ## n is the index of the tilt
    ## we fix first coordinate of the frame, and compose the matrix of the system of equations by iterating by the second coordinate
    A1 = zeros(ComplexF64, prod(PTI.framesize(qqq)), 2)
    b1 = zeros(ComplexF64, prod(PTI.framesize(qqq)))
    return get_sigma(qqq, n, A1, b1)
end

function get_sigma(qqq, n, A1, b1)
    ## n is the index of the tilt
    ## mx is the index of the point in the first coordinate of the frame
    coords = Iterators.product(qqq.frameaxes...)
    igrams = PTI.getigrams(qqq)
    for (i, x) in enumerate(coords)
        A1[i, 1] = PTI.complexamplitude(qqq)[i] * qqq.mask[i]
        A1[i, 2] = conj(A1[i, 1])
        b1[i] =
            (igrams[i + (n - 1) * length(coords)] - PTI.background(qqq)[i]) * qqq.mask[i]
    end
    if sum(b1 .!= 0) > 2
        d1 = A1 \ b1
        return angle(d1[1])
    else
        return NaN
    end
end

# Perform iterations of refinement using exact tilt estimation
for k in 1:10
    @info "Iteration $k"

    ## Phase and background estimation
    deltas_est = eachslice(PTI.apply.(est.tilts, coords_analysis); dims=PTI.setdims(est))
    PhaseBgAmp = psialg(test_data, deltas_est; full=true)

    ## Update estimates
    PTI.setcomplexamplitude!(est, PhaseBgAmp[2])
    PTI.setbackground!(est, PhaseBgAmp[3])

    ## Exact tilt refinement using the functions from PTIestimate_test.jl
    all_taux = [get_taux(est, tiltind) for tiltind in est.setaxes[1]]
    all_tauy = [get_tauy(est, tiltind) for tiltind in est.setaxes[1]]
    all_sigma = [get_sigma(est, tiltind) for tiltind in est.setaxes[1]]

    ## Update tilts
    newtilts = reshape(
        [PTI.FreeTilt([c]) for c in zip(all_sigma, all_taux, all_tauy)], size(est.tilts)
    )
    est.tilts .= newtilts

    ## Display progress
    if k <= 2  # Show first two iterations
        plot_heatmaps_table(
            PTI.getigramssliced(est); colormap=coligram, title="Iteration $k"
        ) |> display
    end
end

# ## 14. Final Results and Convergence Analysis
#
# Display the final estimated phase and compare with the ground truth.

## Extract final phase estimate
final_phase = PTI.getphase(est)

## Display final results
fig = Figure(; size=(800, 400))
fig[1, 1] = Axis(fig; title="Final Phase Estimate")
fig[1, 2] = Axis(fig; title="Ground Truth Phase")
showphase!(fig[1, 1], final_phase)
showphase!(fig[1, 2], PTI.getphase(pti_est))
fig

# Compare tilt estimation accuracy
gt_taux = [PTI.tau(tilt)[1] for tilt in PTI.tilts(pti_est)[1, 1, :]]
gt_tauy = [PTI.tau(tilt)[2] for tilt in PTI.tilts(pti_est)[1, 1, :]]
gt_sigma = [PTI.sigma(tilt) for tilt in PTI.tilts(pti_est)[1, 1, :]]

est_taux = [PTI.tau(tilt)[1] for tilt in est.tilts[1, 1, :]]
est_tauy = [PTI.tau(tilt)[2] for tilt in est.tilts[1, 1, :]]
est_sigma = [PTI.sigma(tilt) for tilt in est.tilts[1, 1, :]];

# Display tilt estimation errors
fig_tilts = Figure(; size=(1200, 300))
scatter(fig_tilts[1, 1], est_taux - gt_taux; axis=(; title="Error in τₓ"))
scatter(fig_tilts[1, 2], est_tauy - gt_tauy; axis=(; title="Error in τᵧ"))
scatter(
    fig_tilts[1, 3], phwrap.(est_sigma - gt_sigma); axis=(; title="Error in σ (wrapped)")
)
fig_tilts

# ## Conclusions
#
# This tutorial has demonstrated:
#
# 1. **PTIestimate Construction**: How to create PTIestimate objects from interferogram data
# 2. **Forward Modeling**: Using PTIestimate to generate synthetic interferograms with controlled phase and tilt patterns
# 3. **Mask Application**: Setting up aperture masks for realistic experimental conditions
# 4. **Phase Reconstruction**: Extracting phase information using both direct linear algebra and iterative algorithms
# 5. **Iterative Refinement**: The alternating optimization approach for simultaneously estimating phase, background, and tilts
#
# The PTIestimate structure provides a unified framework for PTI analysis that can handle both forward modeling for simulation and inverse problems for phase reconstruction from experimental data.
#
# ## Next Steps
#
# - Explore more sophisticated tilt estimation algorithms
# - Apply the framework to real experimental data
# - Investigate convergence properties of iterative algorithms
# - Extend to more complex aperture geometries and phase patterns
