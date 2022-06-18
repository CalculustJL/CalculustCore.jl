#
struct SpectralGalerkin{Tsp::AbstractSpectralSpace}
    space::Tsp
end

@forward SpectralGalerkin.space (
                                 grid, domain, npoints,
                                 boundary_nodes
                                )
# implement vector calculus ops here

struct SpectralCollocation{Tsp::AbstractSpectralSpace}
    space::Tspace
end
                      
@forward SpectralCollocation.space (
                                    grid, domain, npoints,
                                    boundary_nodes
                                   )
#
