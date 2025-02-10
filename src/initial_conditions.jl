
export WilliamsonInicond1

"""
$(TYPEDSIGNATURES)
Initial conditions for the first test case of Williamson et al. 1992.
α defaults to 0.0 and should be set to the transport direction being tested.
$(TYPEDFIELDS)"""
@kwdef struct WilliamsonInicond1 <: SpeedyWeather.AbstractInitialConditions
   h0::Float64 = 1000 # m
   θc::Float64 = 0.0 # degrees
   λc::Float64 = 270.0 # degrees
   α::Float64 # degrees
end

"""
$(TYPEDSIGNATURES)
Initial conditions for the first test case of Williamson et al. 1992, J Computational Physics"""
function SpeedyWeather.initialize!(progn::PrognosticVariables,
                                   initial_conditions::WilliamsonInicond1,
                                   model::AbstractModel)

   R = model.spectral_grid.radius/3
   u0 = 2π*model.spectral_grid.radius/(12*24*60*60) # advecting wind velocity ≈ 40m/s

   u(θ, λ, σ) = u0*(cosd(θ)*cosd(α) + sind(θ)*cosd(λ)*sind(α))
   v(θ, λ, σ) = -u0*sind(λ)*sind(α)
   function η(θ, λ)
      rdist = R*acosd(sind(θc)*sind(θ) + cos(θc)*cosd(θ)*cos(λ-λc))
      if rdist < r
         (h0/2)*(1+cosd(180*rdist/r))
      else
          0.0
      end
   end
   set!(model, u=u, v=v, η=η)
end
