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

function polarMap(r, θ)
    x = @. r * cos(θ)
    y = @. r * sin(θ)
    x, y
end

function reference_box(D; periodic_dirs=())
    domain = BoxDomain()
    for i=1:D
        periodic = i ∈ periodic_dirs
        tag1 = Symbol("Lower$(i)")
        tag2 = Symbol("Upper$(i)")
        tags = (tag1, tag2)
        interval = IntervalDomain(-true, true; periodic=periodic, boundary_tags=tags)
        domain = domain ⊗ interval
    end
    domain
end

function annulus_2D(r0, r1)
    intR = IntervalDomain(r0, r1, false, (:Inner, :Outer))
    intθ = IntervalDomain(-π,  π, true , (:Periodic, :Periodic))

    dom = intR ⊗ intθ

    deform(dom, polar)
end
#
