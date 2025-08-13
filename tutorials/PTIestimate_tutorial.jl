# # PTIestimate Structure Tutorial
#
# ## Introduction
#
# This tutorial demonstrates the use of the [`PTIestimate`](@ref) structure for Phase-Tilted Interferometry (PTI) analysis. The PTIestimate structure serves as a comprehensive framework that supports two primary use cases:
#
# 1. **Forward Modeling**: Generate synthetic interferograms from known phase and tilt parameters
# 2. **Inverse Problems**: Estimate phase and tilts from measured interferogram data
#
# The structure provides a unified interface for both scenarios, with the key distinction being whether measured data is provided or not.

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

# # Part I: Forward Modeling with PTIestimate
#
# In forward modeling, we create a PTIestimate structure to generate synthetic interferograms from known parameters. This is useful for algorithm development, testing, and educational purposes.

# ## 1. Creating PTIestimate for Forward Modeling
#
# We start by creating a PTIestimate structure from dimensions only, without any measured data.

# Create PTIestimate for forward modeling: 500×800 pixels, 15 interferograms
pti_forward = PTIestimate((500, 800), (15,); frameaxes=PTI.FourierAxes())
@show typeof(pti_forward)

# Check that no measured data is present
@show PTI.hasdata(pti_forward)  # Should be false

# Examine the coordinate system and structure
@show pti_forward.frameaxes
@show PTI.framesize(pti_forward)
@show PTI.setsize(pti_forward)

# Create coordinate iterator for later use
coords = Iterators.product(pti_forward.frameaxes...)
@show size(coords)

# ## 2. Setting Up the Forward Model
#
# We now configure the PTIestimate structure with the parameters needed for forward modeling.

# Set circular aperture mask in Fourier space
r = 0.45
PTI.setmask!(pti_forward, map(x -> x[1]^2 + x[2]^2 <= r^2, coords));

# ## 3. Initial State: Understanding the Default Parameters
#
# Let's examine what the PTIestimate structure contains initially and understand the forward model.

# Initially, all parameters are at default values
# Visualize the aperture mask, initial phase, and background together
plot_heatmaps_table(
    [
        PTI.mask(pti_forward)[:, :, 1],
        PTI.getphase(pti_forward),
        PTI.background(pti_forward)[:, :, 1],
    ];
    ncols=3,
    aspect=AxisAspect(1),
    titles=["Aperture Mask", "Initial Phase", "Initial Background"],
)

# Display initial interferograms (all identical due to zero phase and tilts)
plot_heatmaps_table(
    eachslice(PTI.getigrams(pti_forward); dims=3)[1:4];
    colormap=coligram,
    aspect=AxisAspect(1),
    titles=["Frame $i" for i in 1:4],
)

# Show the initial (zero) phase
showphase(PTI.getphase(pti_forward))[1]

# ## 4. Setting a Complex Phase Pattern
#
# Now we'll set a realistic complex phase pattern using Zernike polynomials with Born & Wolf ordering.

phase_pattern = map(coords) do x
    ## Extract coordinates
    x_coord, y_coord = x[1] / r, x[2] / r

    ## Define Zernike coefficients (Born & Wolf ordering, j=1 to 15)
    coefn = [
        0.0,                      # j=1: Piston
        0.04487997792421934,      # j=2: x-tilt
        -0.061213166549227316,    # j=3: y-tilt
        0.0372359177320852,       # j=4: 0° Primary astigmatism
        -0.025807672973317094,    # j=5: Defocus
        -0.10907698627242028,     # j=6: 45° Primary astigmatism
        -0.0025831843312457176,   # j=7: trefoil 1
        -1.2854500898834675,      # j=8: Primary x-coma
        0.5100406297505503,       # j=9: Primary y-coma
        -0.02063276561330821,     # j=10: trefoil 2
        -0.14316991684722485,     # j=11:
        -0.025345357367834,       # j=12: 0° Secondary astigmatism
        0.1393506816540738,       # j=13:  spherical aberration
        -0.19999956756592352,     # j=14: 45° Secondary astigmatism
        0.13278232836541945,       # j=15:
    ]

    ## Calculate polynomial value using weighted Zernike polynomials (compact Cartesian form)
    phase_val =
        coefn[1] * 1 +
        coefn[2] * x_coord +
        coefn[3] * y_coord +
        coefn[4] * (x_coord^2 - y_coord^2) +
        coefn[5] * (2 * (x_coord^2 + y_coord^2) - 1) +
        coefn[6] * (2 * x_coord * y_coord) +
        coefn[7] * (x_coord^3 - 3 * x_coord * y_coord^2) +
        coefn[8] * (3 * x_coord * (x_coord^2 + y_coord^2) - 2 * x_coord) +
        coefn[9] * (3 * y_coord * (x_coord^2 + y_coord^2) - 2 * y_coord) +
        coefn[10] * (3 * x_coord^2 * y_coord - y_coord^3) +
        coefn[11] * (x_coord^4 - 6 * x_coord^2 * y_coord^2 + y_coord^4) +
        coefn[12] * (
            4 * (x_coord^2 + y_coord^2) * (x_coord^2 - y_coord^2) -
            3 * (x_coord^2 - y_coord^2)
        ) +
        coefn[13] * (6 * (x_coord^2 + y_coord^2)^2 - 6 * (x_coord^2 + y_coord^2) + 1) +
        coefn[14] *
        (8 * x_coord * y_coord * (x_coord^2 + y_coord^2) - 6 * x_coord * y_coord) +
        coefn[15] * (4 * x_coord^3 * y_coord - 4 * x_coord * y_coord^3)

    return phase_val
end

# Generate the Zernike phase pattern

PTI.setphase!(pti_forward, 10π * phase_pattern);
showphase(PTI.getphase(pti_forward))[1]

# Display the phase pattern and show interferograms with new phase (still identical due to zero tilts)
plot_heatmaps_table(
    eachslice(PTI.getigrams(pti_forward); dims=3)[1:4];
    colormap=coligram,
    aspect=AxisAspect(1),
    titles=["Frame 1", "Frame 2", "Frame 3", "Frame 4"],
)

# ## 5. Adding Tilts: Creating Interferogram Diversity
#
# Tilts create the diversity between interferograms that makes PTI analysis possible.

# Add random tilts to each interferogram
PTI.settilts!(
    pti_forward, PTI.FreeTilt.([10π * (randn(3) .- 0.5) for _ in PTI.tilts(pti_forward)])
);

# Now interferograms show clear differences due to tilts
plot_heatmaps_table(
    eachslice(PTI.getigrams(pti_forward); dims=3);
    colormap=coligram,
    aspect=AxisAspect(1),
    titles=["Frame $i" for i in 1:15],
)

# Update tilts to demonstrate different tilt patterns
PTI.settilts!(
    pti_forward, PTI.FreeTilt.([10π * (randn(3) .- 0.5) for _ in PTI.tilts(pti_forward)])
)
plot_heatmaps_table(
    eachslice(PTI.getigrams(pti_forward); dims=3);
    colormap=coligram,
    aspect=AxisAspect(1),
    titles=["Frame $i" for i in 1:15],
)

print("Some current tilts: ", pti_forward.tilts[1, 1, 1:3])

# ## 6. Forward Model Applications
#
# The forward model can be used for various analysis tasks.

# Extract the synthetic interferograms for further analysis
synthetic_igrams = PTI.getigrams(pti_forward);

# Helper function for plotting
pht(arr; kwargs...) =
    plot_heatmaps_table(eachslice(arr; dims=3); aspect=AxisAspect(1), kwargs...)

# Analyze interferogram differences (important for tilt estimation)
deltas = diff(synthetic_igrams; dims=3)
pht(deltas; colormap=coligramdiff, titles=["Δ$i" for i in 1:size(deltas, 3)])

# ## 7. Understanding the Physical Model
#
# Each interferogram follows the model: I = background + Re(complex_amplitude × exp(i×tilt))

# Examine different components
@show size(PTI.background(pti_forward))
@show size(PTI.complexamplitude(pti_forward))
@show size(PTI.tilts(pti_forward))

# Store our forward model results for later comparison
ground_truth_phase = PTI.getphase(pti_forward)
ground_truth_tilts = deepcopy(PTI.tilts(pti_forward))
ground_truth_igrams = copy(synthetic_igrams);

# # Part II: Inverse Problems with PTIestimate
#
# Now we'll demonstrate how to use PTIestimate for inverse problems - estimating parameters from measured data.

# ## 8. Creating PTIestimate from Measured Data
#
# We'll use our synthetic data as "measured" interferograms to test the inverse algorithms.

# Create PTIestimate from interferogram data (simulating measured data)
measured_data = eachslice(ground_truth_igrams; dims=3)
pti_inverse = PTIestimate(measured_data; frameaxes=PTI.FourierAxes())

# Check that measured data is now stored in the structure
@show PTI.hasdata(pti_inverse)  # Should be true
@show size(PTI.data(pti_inverse))

# The structure now contains both model parameters and measured data
@show typeof(pti_inverse)
@show PTI.framesize(pti_inverse)
@show PTI.setsize(pti_inverse)

# ## 9. Initial Parameter Estimation
#
# The first step in inverse PTI is to obtain rough estimates of the tilts.

# Set the same aperture mask as used in forward modeling
PTI.setmask!(pti_inverse, map(x -> x[1]^2 + x[2]^2 <= r^2, coords));

# Initialize with rough tilt estimates using the new convenient API
refframe = 2
PTI.initialize!(pti_inverse; refframe=refframe)  # Uses internal data automatically

# Display initial tilt estimates
@show pti_inverse.tilts[1, 1, 1:3]

# Show what the interferograms look like with initial estimates
plot_heatmaps_table(
    PTI.getigramssliced(pti_inverse)[1:6];
    colormap=coligram,
    aspect=AxisAspect(1),
    titles=["Initial Est. $i" for i in 1:6],
)

# Adjust tilt signs for proper convergence (using ground truth for demonstration)
normals = [
    PTI.tau(t) - PTI.tau(ground_truth_tilts[1, 1, refframe]) for
    t in ground_truth_tilts[1, 1, :]
]
PTI.set_tilt_signs!(pti_inverse, normals)

# ## 10. Phase and Background Estimation
#
# Use phase-shifting interferometry to estimate phase and background from the current tilt estimates.

# Apply least-squares PSI algorithm using the convenient API
psialg = PTI.LSPSI()
coords_inv = Iterators.product(pti_inverse.frameaxes...)
deltas_est = eachslice(
    PTI.apply.(pti_inverse.tilts, coords_inv); dims=PTI.setdims(pti_inverse)
)
PhaseBgAmp = psialg(eachslice(PTI.data(pti_inverse); dims=3), deltas_est; full=true);

# Display initial estimation results
fig = Figure(; size=(1200, 400))
fig[1, 1] = Axis(fig; title="Estimated Phase")
fig[1, 2] = Axis(fig; title="Estimated Background")
fig[1, 3] = Axis(fig; title="Estimated Amplitude")
showphase!(fig[1, 1], PhaseBgAmp[1])
showarray!(fig[1, 2], PhaseBgAmp[3])
showarray!(fig[1, 3], abs.(PhaseBgAmp[2]))
fig

# Update the estimate with new values
PTI.setcomplexamplitude!(pti_inverse, PhaseBgAmp[2])
PTI.setbackground!(pti_inverse, PhaseBgAmp[3]);

# ## 11. Iterative Refinement
#
# The key to accurate PTI analysis is iterative refinement of both phase/background and tilts.

# Exact tilt estimation functions (from PTIestimate_test.jl)
function get_taux(qqq, n)
    ## n is the index of the tilt
    A1 = zeros(ComplexF64, PTI.framesize(qqq)[2], 2)
    b1 = zeros(ComplexF64, PTI.framesize(qqq)[2])
    alltau = [get_taux(qqq, n, mx, A1, b1) for mx in 1:PTI.framesize(qqq)[1]]
    t1est = phwrap(diff(alltau[(!isnan).(alltau)]))
    return mean(t1est) / step(qqq.frameaxes[1])
end

function get_taux(qqq, n, mx, A1, b1)
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
    A1 = zeros(ComplexF64, PTI.framesize(qqq)[1], 2)
    b1 = zeros(ComplexF64, PTI.framesize(qqq)[1])
    alltau = [get_tauy(qqq, n, my, A1, b1) for my in 1:PTI.framesize(qqq)[2]]
    t1est = phwrap(diff(alltau[(!isnan).(alltau)]))
    return mean(t1est) / step(qqq.frameaxes[2])
end

function get_tauy(qqq, n, my, A1, b1)
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
    A1 = zeros(ComplexF64, prod(PTI.framesize(qqq)), 2)
    b1 = zeros(ComplexF64, prod(PTI.framesize(qqq)))
    return get_sigma(qqq, n, A1, b1)
end

function get_sigma(qqq, n, A1, b1)
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

# Perform iterative refinement
for k in 1:10
    @info "Iteration $k"

    ## Phase and background estimation
    deltas_est = eachslice(
        PTI.apply.(pti_inverse.tilts, coords_inv); dims=PTI.setdims(pti_inverse)
    )
    PhaseBgAmp = psialg(PTI.data(pti_inverse), deltas_est; full=true)

    ## Update estimates
    PTI.setcomplexamplitude!(pti_inverse, PhaseBgAmp[2])
    PTI.setbackground!(pti_inverse, PhaseBgAmp[3])

    ## Exact tilt refinement
    all_taux = [get_taux(pti_inverse, tiltind) for tiltind in pti_inverse.setaxes[1]]
    all_tauy = [get_tauy(pti_inverse, tiltind) for tiltind in pti_inverse.setaxes[1]]
    all_sigma = [get_sigma(pti_inverse, tiltind) for tiltind in pti_inverse.setaxes[1]]

    ## Update tilts
    newtilts = reshape(
        [PTI.FreeTilt([c]) for c in zip(all_sigma, all_taux, all_tauy)],
        size(pti_inverse.tilts),
    )
    pti_inverse.tilts .= newtilts

    ## Display progress for first few iterations
    if k <= 2
        plot_heatmaps_table(
            PTI.getigramssliced(pti_inverse)[1:4];
            colormap=coligram,
            aspect=AxisAspect(1),
            titles=["Iter $k: Frame $i" for i in 1:4],
        )
    end
end

# ## 12. Validation Against Ground Truth
#
# Now we can compare our estimated parameters with the known ground truth from the forward model.

# Extract final estimated parameters
final_estimated_phase = PTI.getphase(pti_inverse)

# Compare phases
fig = Figure(; size=(1200, 400))
fig[1, 1] = Axis(fig; title="Estimated Phase")
fig[1, 2] = Axis(fig; title="Ground Truth Phase")
fig[1, 3] = Axis(fig; title="Phase Difference")
showphase!(fig[1, 1], final_estimated_phase)
showphase!(fig[1, 2], ground_truth_phase)
showphase!(fig[1, 3], phwrap.(final_estimated_phase - ground_truth_phase))
fig

# Compare tilt estimation accuracy
gt_taux = [PTI.tau(tilt)[1] for tilt in ground_truth_tilts[1, 1, :]]
gt_tauy = [PTI.tau(tilt)[2] for tilt in ground_truth_tilts[1, 1, :]]
gt_sigma = [PTI.sigma(tilt) for tilt in ground_truth_tilts[1, 1, :]]

est_taux = [PTI.tau(tilt)[1] for tilt in pti_inverse.tilts[1, 1, :]]
est_tauy = [PTI.tau(tilt)[2] for tilt in pti_inverse.tilts[1, 1, :]]
est_sigma = [PTI.sigma(tilt) for tilt in pti_inverse.tilts[1, 1, :]]

# Display tilt estimation errors
fig_tilts = Figure(; size=(1200, 300))
scatter(fig_tilts[1, 1], est_taux - gt_taux; axis=(; title="Error in τₓ"))
scatter(fig_tilts[1, 2], est_tauy - gt_tauy; axis=(; title="Error in τᵧ"))
scatter(
    fig_tilts[1, 3], phwrap.(est_sigma - gt_sigma); axis=(; title="Error in σ (wrapped)")
)
fig_tilts

# Print summary statistics
println(
    "Phase reconstruction RMS error: ",
    sqrt(mean((phwrap.(final_estimated_phase - ground_truth_phase)) .^ 2)),
)
println("Tilt τₓ RMS error: ", sqrt(mean((est_taux - gt_taux) .^ 2)))
println("Tilt τᵧ RMS error: ", sqrt(mean((est_tauy - gt_tauy) .^ 2)))
println("Tilt σ RMS error: ", sqrt(mean((phwrap.(est_sigma - gt_sigma)) .^ 2)))

# # Part III: Key Features and API Summary

# ## 13. PTIestimate API Summary
#
# The tutorial has demonstrated the two primary usage patterns:

# **Forward Modeling (no measured data):**
# ```julia
# pti = PTIestimate((nx, ny), (nframes,))  # hasdata(pti) == false
# PTI.setphase!(pti, phase_pattern)
# PTI.settilts!(pti, tilt_array)
# synthetic_data = PTI.getigrams(pti)
# ```

# **Inverse Problems (with measured data):**
# ```julia
# pti = PTIestimate(measured_interferograms)  # hasdata(pti) == true
# PTI.initialize!(pti)  # Uses internal data
# # Iterative refinement...
# estimated_phase = PTI.getphase(pti)
# ```

# ## Conclusions
#
# This tutorial has demonstrated the comprehensive capabilities of the PTIestimate structure:
#
# ### Forward Modeling Capabilities:
# 1. **Parameter-based interferogram synthesis** from phase, tilts, and background
# 2. **Ground truth generation** for algorithm development and testing
# 3. **Physical model validation** and sensitivity analysis
#
# ### Inverse Problem Capabilities:
# 4. **Automatic initialization** from measured interferogram data
# 5. **Unified data management** with the optional data field
# 6. **Iterative parameter estimation** with exact tilt refinement algorithms
# 7. **Clean API design** with convenience methods for common operations
#
# ### Key Advantages:
# - **Unified interface** for both forward and inverse problems
# - **Data-aware design** distinguishing between modeling and measurement scenarios
# - **Comprehensive state management** for iterative algorithms
# - **Scientific accuracy** with exact mathematical implementations
#
# The PTIestimate structure provides a robust foundation for Phase-Tilted Interferometry analysis, supporting both educational exploration and research applications.
#
# ## Next Steps
#
# - Apply the framework to real experimental interferogram data
# - Explore advanced tilt estimation algorithms and convergence analysis
# - Investigate performance with different noise levels and aperture geometries
# - Extend to more complex phase patterns and multi-wavelength applications
