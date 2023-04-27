#=
# TODO Gordon Hall
function mesh_domain(dom1::BoxDomain{<:Number,D},
                     space::AbstractSpace{<:Number,D}
                    ) where{D}

    #function gordonHall(xrm,xrp,xsm,xsp,yrm,yrp,ysm,ysp,zr,zs)
    ze  = [-1,1]
    Jer = interpMat(zr,ze)
    Jes = interpMat(zs,ze)

    xv = [xrm[1]   xrp[1]
          xrm[end] xrp[end]]

    yv = [yrm[1]   yrp[1]
          yrm[end] yrp[end]]

    xv = ABu(Jes,Jer,xv)
    yv = ABu(Jes,Jer,yv)

    #display(mesh(xv,yv,0*xv,0,90))

    x = ABu([],Jer,vcat(xrm',xrp')) .+ ABu(Jes,[],hcat(xsm,xsp)) .- xv
    y = ABu([],Jer,vcat(yrm',yrp')) .+ ABu(Jes,[],hcat(ysm,ysp)) .- yv

    return
end
=#

###
# Conveniences
###

function default_tags(i::Integer)
    tag1 = Symbol("Lower$(i)")
    tag2 = Symbol("Upper$(i)")

    (tag1, tag2)
end

function GaussLobattoLegendreDomain(D; periodic_dirs = ())
    domain = BoxDomain()
    endpts = (-true, true)
    for i in 1:D
        interval = IntervalDomain(endpts...;
                                  periodic = i ∈ periodic_dirs,
                                  boundary_tags = default_tags(i))
        domain = domain ⊗ interval
    end
    domain
end

function ChebyshevDomain(D; periodic_dirs = ())
    domain = BoxDomain()
    endpts = (-true, true)
    for i in 1:D
        interval = IntervalDomain(endpts...;
                                  periodic = i ∈ periodic_dirs,
                                  boundary_tags = default_tags(i))
        domain = domain ⊗ interval
    end
    domain
end

function FourierDomain(D; periodic_dirs = 1:D)
    domain = BoxDomain()
    endpts = (-π, π)
    for i in 1:D
        interval = IntervalDomain(endpts...;
                                  periodic = i ∈ periodic_dirs,
                                  boundary_tags = default_tags(i))
        domain = domain ⊗ interval
    end
    domain
end

function polarMap(r, θ)
    x = @. r * cos(θ)
    y = @. r * sin(θ)
    x, y
end

function AnnulusDomain(rinner, router)
    intR = IntervalDomain(rinner, router, false, (:Inner, :Outer))
    intθ = IntervalDomain(-π, π, true, (:Periodic, :Periodic))

    dom = intR ⊗ intθ

    deform(dom, polar)
end
#
