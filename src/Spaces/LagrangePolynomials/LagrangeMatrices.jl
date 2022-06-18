#
###
# Lagrange polynomial matrices
###

function lagrange_barycentric_weights(x)
    n = length(x)

    a = ones(1,n)
    for i=1:n
        for j=1:(i-1) a[i]=a[i]*(x[i]-x[j]) end
        for j=(i+1):n a[i]=a[i]*(x[i]-x[j]) end
    end

    a = 1 ./ a
end

"""
    Compute the Lagrange interpolation matrix from xi to xo.
    lagrange_poly_interp_mat(xₒ,xᵢ)
"""
function lagrange_interp_mat(xo,xi)

    no = length(xo)
    ni = length(xi)

    a = lagrange_barycentric_weights(xi)

    J = zeros(no,ni)
    s = ones(1,ni)
    t = ones(1,ni)

    for i=1:no
        x = xo[i]
        for j=2:ni
            s[j]      = s[j-1]    * (x-xi[j-1]   ) # left sum
            t[ni+1-j] = t[ni+2-j] * (x-xi[ni+2-j]) # right sum
        end
        J[i,:] = a .* s .* t
    end

    return J
end

"""
 Compute derivative matrix for lagrange
 interpolants on points [x]
"""
function lagrange_deriv_mat(x)
    
    n = length(x)
    a = lagrange_barycentric_weights(x)

    # diagonal elements
    D = x .- x'
    for i=1:n D[i,i] = 1. end
    D = 1 ./ D
    for i=1:n
        D[i,i] = 0.
        D[i,i] = sum(D[i,:])
    end

    # off-diagonal elements
    for j=1:n for i=1:n
        if(i!=j) D[i,j] = a[j] / (a[i]*(x[i]-x[j])) end
    end end

    return D
end
#
