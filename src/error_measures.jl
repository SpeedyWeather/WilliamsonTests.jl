
export ErrorMeasures

"""
$(TYPEDSIGNATURES)
Create a struct containing the error measures detailed in Williamson et al. 1992, J Computational Physics
`reference` is a function of time, longitude and latitude (both in degrees)
$(TYPEDFIELDS)"""
@kwdef mutable struct ErrorMeasures <: SpeedyWeather.AbstractCallback
    reference::Function
    output_file::String = "errormeasures.nc"

    # both of these are computed in the initialize! routine
    initial_mean::Float64 = 0.0
    initial_vorticity::Float64 = 0.0

    timestep_counter::Int = 0
    time::Vector{DateTime} = []
    l1error::Vector{Float64} = []
    l2error::Vector{Float64} = []
    l∞error::Vector{Float64} = []
    mean::Vector{Float64} = []
    variance::Vector{Float64} = []
    minimum::Vector{Float64} = []
    maximum::Vector{Float64} = []
end

"""
$(TYPEDSIGNATURES)
Error measures for the test cases of Williamson et al. 1992, J Computational Physics"""
function SpeedyWeather.initialize!(callback::ErrorMeasures,
                                   progn::PrognosticVariables,
                                   diagn::DiagnosticVariables,
                                   model::AbstractModel)

    n = progn.clock.n_timesteps + 1
    callback.time = zeros(DateTime, n)
    callback.l1error = zeros(n)
    callback.l2error = zeros(n)
    callback.l∞error = zeros(n)
    callback.mean = zeros(n)
    callback.variance = zeros(n)
    callback.minimum = zeros(n)
    callback.maximum = zeros(n)

    η0 = diagn.grid.pres_grid
    H = model.atmosphere.layer_thickness
    Hb = model.orography.orography
    h0 = @. η + H - Hb

    callback.initial_mean = ∬dA(h0)
    callback.initial_variance = ∬dA((h0 .- callback.initial_mean) .^ 2)

    callback.time[1] = progn.clock.time
    callback.timestep_counter = 1

    l1, l2, l∞, mean, var, min, max = error_measures(callback, diagn, model)

    callback.l1error[1] = l1
    callback.l2error[1] = l2
    callback.l∞error[1] = l∞
    callback.mean[1] = mean
    callback.variance[1] = var
    callback.minimum[1] = min
    callback.maximum[1] = max

end

"""
$(TYPEDSIGNATURES)
Error measures for the test cases of Williamson et al. 1992, J Computational Physics"""
function SpeedyWeather.callback!(callback::ErrorMeasures,
                                 progn::PrognosticVariables,
                                 diagn::DiagnosticVariables,
                                 model::AbstractModel)

    callback.timestep_counter += 1

    callback.time[i] = progn.clock.time

    l1, l2, l∞, mean, var, min, max = error_measures(callback, diagn, model)

    callback.l1error[i] = l1
    callback.l2error[i] = l2
    callback.l∞error[i] = l∞
    callback.mean[i] = mean
    callback.variance[i] = var
    callback.minimum[i] = min
    callback.maximum[i] = max
end

using NCDatasets

"""
$(TYPEDSIGNATURES)
Error measures for the test cases of Williamson et al. 1992, J Computational Physics"""
function SpeedyWeather.finalize!(callback::ErrorMeasures,
                                 progn::PrognosticVariables,
                                 diagn::DiagnosticVariables,
                                 model::AbstractModel)

    n_timesteps = callback.timestep_counter

    ds = NCDataset(joinpath(pwd(), callback.output_file), "c")

    defDim(ds, "time", n_timesteps)
    defVar(ds, "time",     callback.time,     ("time",))
    defVar(ds, "l1error",  callback.l1error,  ("time",))
    defVar(ds, "l2error",  callback.l2error,  ("time",))
    defVar(ds, "l∞error",  callback.l∞error,  ("time",))
    defVar(ds, "mean",     callback.mean,     ("time",))
    defVar(ds, "variance", callback.variance, ("time",))
    defVar(ds, "minimum",  callback.minimum,  ("time",))
    defVar(ds, "maximum",  callback.maximum,  ("time",))

    close(ds)

    return nothing
end

"""
$(TYPEDSIGNATURES)
Computes the surface integral of a function on the sphere
"""
function ∬dA(f, S::SpectralTransform)
    return real(transform(f, S)[1]) / S.norm_sphere
end

∬dA(f, model::AbstractModel) = ∬dA(f, model.spectral_transform)

"""
$(TYPEDSIGNATURES)
Error measures for the test cases of Williamson et al. 1992, J Computational Physics"""
function error_measures(callback::ErrorMeasures,
                        diagn::DiagnosticVariables,
                        model::AbstractModel)

    η = diagn.grid.pres_grid
    ηref = zeros(typeof(η), η.nlat_half)
    londs, latds = RingGrids.get_londlatds(η)
    rings = RingGrids.eachring(η)
    for (j, ring) in enumerate(rings)
        for ij in ring
            ηref[ij] = callback.reference(callback.time[callback.timestep_counter], londs[ij], latds[ij])
        end
    end

    R = model.spectral_grid.radius
    H = model.atmosphere.layer_thickness
    Hb = model.orography.orography

    h = @. η + H - Hb
    href = @. ηref + H - Hb
    l1error = ∬dA(abs.(h .- href), model) / ∬dA(abs.(href, model))
    l2error = sqrt(∬dA((h .- href) .^ 2)) / sqrt(∬dA(href .^ 2))
    l∞error = max.( abs.(h .- href) ) / max.( abs.(href) )
    hmean = ∬dA(h, model)
    hrefmean = ∬dA(href, model)
    mean = (hmean - hrefmean) / callback.initial_mean
    variance =  (∬dA((h .- hmean) .^ 2) - ∬dA((href .- hrefmean) .^ 2)) / callback.initial_variance
    minimum = (max.(h) - max.(href)) / callback.Δh
    maximum = (min.(h) - min.(href)) / callback.Δh

    return l1error, l2error, l∞error, mean, variance, minimum, maximum
end
