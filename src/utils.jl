#
import SciMLOperators: AbstractSciMLOperator, _reshape, _vec
import LinearSolve: default_tol

function set_val!(M, val, idx)
    len = length(idx)

    if len < 1
        M
    elseif len == 1
        M[idx] = val
    else
        M[idx] = [val for i in 1:length(idx)]
    end

    M
end

