#
struct SpectralGalerkin{Tsp::AbstractSpectralSpace}
    space::Tsp
end

@forward SpectralGalerkin.space (
                                 get_grid, get_domain, numpoints,
                                 boundary_nodes
                                )
# implement vector calculus ops here

struct SpectralCollocation{Tsp::AbstractSpectralSpace}
    space::Tspace
end
                      
@forward SpectralCollocation.space (
                                    get_grid, get_domain, numpoints,
                                    boundary_nodes
                                   )


