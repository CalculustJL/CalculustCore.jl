#
struct TensorProductSpace{T,D1+D2,
                          Tspace1<:AbstractSpace{<:Number, D1},
                          Tspace2<:AbstractSpace{<:Number, D2},
                         } <: AbstractTensorProductSpace{T,D1+D2}
    space1::Tspace1
    space2::Tspace2
end
