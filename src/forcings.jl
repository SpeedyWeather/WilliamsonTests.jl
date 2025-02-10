
export NoVorDivTendencies

"""
A forcing that removes vorticity and divergence tendences in each time step
so that only the advection is simulated."""
@kwdef struct NoVorDivTendencies{NF} <: SpeedyWeather.AbstractForcing
end

"""
$(TYPEDSIGNATURES)
A forcing that removes vorticity and divergence tendences in each time step."""
function SpeedyWeather.initialize!(forcing::NoVorDivTendencies{NF},
                                   model::AbstractModel) where NF
    return nothing
end

"""
$(TYPEDSIGNATURES)
A forcing that removes vorticity and divergence tendences in each time step."""
function SpeedyWeather.forcing!(diagn::DiagnosticVariables,
                                progn::PrognosticVariables,
                                forcing::NoVorDivTendencies{NF},
                                lf::Integer,
                                model::AbstractModel) where NF

    forcing!(diagn, forcing, model.spectral_transform)
end

"""
$(TYPEDSIGNATURES)
A forcing that removes vorticity and divergence tendences in each time step."""
function forcing!(diagn::DiagnosticVariables,
                  forcing::NoVorDivTendencies{NF},
                  spectral_transform::SpectralTransform) where NF

    diagn.tendencies.vor_tend .= 0
    diagn.tendencies.div_tend .= 0
    return nothing
end
