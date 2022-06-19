#
abstract type AbstractRegressionProblem end
abstract type AbstractRegressionCache end
abstract type RegressionAlgorithm end

struct RegressionProblem <: AbstractRegressionProblem
    """ target function """
    target_func
    """ function space """
    space
    """ interpolant """
    J
    """ model """
    model
    """ parameters """
    p
    """ kwargs """
    kwargs
    """ loss """
    loss
end

# dispatch to OptimizationProblem
function SciMLBase.init(prob::RegressionProblem)
end
