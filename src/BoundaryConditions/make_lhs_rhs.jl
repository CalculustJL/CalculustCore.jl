#
###
# boundary masks
###

function dirichlet_mask(space::AbstractSpace, dom::AbstractDomain, indices, bc_dict)
    tags = boundary_tags(dom)
    x = points(space) |> first
    M = similar(x, Bool) * false .+ true |> vec

    for i in 1:num_boundaries(dom)
        tag = boundary_tag(dom, i)
        bc = bc_dict[tag]

        if bc isa DirichletBC
            idx = indices[i]
            set_val!(M, false, idx)
        end
    end

    DiagonalOperator(M)
end

function boundary_antimasks(space::AbstractSpace, dom::AbstractDomain, indices)
    glo_num = global_numbering(space)

    antimasks = []
    for i in 1:num_boundaries(dom)
        M = similar(glo_num, Bool) * false |> vec
        idx = indices[i]
        set_val!(M, true, idx)
        push!(antimasks, M)
    end

    DiagonalOperator.(antimasks)
end

###
# form LHS, RHS
###

function makeLHS(op::AbstractSciMLOperator, bc::AbstractBoundaryCondition)
    @unpack mask_dir, amask_dir = bc

    #TODO
    """
    how to empty dirichlet row? -
        https://github.com/vpuri3/PDEInterfaces.jl/issues/1

    is having the construct u := u_inhom + u_hom the only way?
    then would have to interpolate u_inhom into interior.
    what are the consequences?
    """

    mask_dir * op + amask_dir
end

function makeRHS(f, bc::AbstractBoundaryCondition)
    @unpack bc_dict, antimasks, mask_dir, space, discr = bc

    M = massOp(space, discr)
    b = M * f

    dirichlet = zero(b)
    neumann = zero(b)
    robin = zero(b)

    pts = points(space)
    dom = domain(space)

    for i in 1:num_boundaries(dom)
        tag = boundary_tag(dom, i)
        bc = bc_dict[tag]
        amask = antimasks[i]

        if bc isa DirichletBC
            dirichlet += amask * bc.f(pts...)
        elseif bc isa NeumannBC
            neumann += amask * bc.f(pts...)
        elseif bc isa RobinBC
            # TODO
        elseif bc isa PeriodicBC
            continue
        end
    end

    (mask_dir * b) + dirichlet - neumann + robin
end
#
#
