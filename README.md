## PDEInterfaces.jl

Tired of writing boilerplate code for PDE solvers? Just want to focus on discretizations? Want something that plays nice with `DiffEq` ecosystem? AbstractPDEs.jl contains separate abstract interaces for multidimensional domains, fields, operators, and function spaces. it is general enough that anybody can plug in their discretizations and start solving PDEs. in fact i aim to use it with my experimental NN discretizations. All you need to do is provide a gradient operator, and a mass operator (integration).

Once you plug in your discretizations, you can do a lot of cool things like apply any random deformations to the space. AbstractPDEs.jl translate all your vector calculus operators correctly. That means the same code could solve convection diffusion on a square as well as an annulus with no extra work and basically conserved accuracy.

After describing your problem, it should spit out the right BoundaryValueProblem  or ODEProblem  that you can solve using the correct `DiffEq` package.

Goals:
- [ ] Add `AbstractCalculus` object between function space, and `AbstractDiscretization`
- [ ] Abstract Domain interface
  - [ ] Is it possibe to just use `DomainSets.jl` and add some metadata info?
  - [X] Logically rectangular
  - [X] Deform domain
  - [X] Boundary tags
  - [ ] Interior tags
  - [ ] Gordon Hall interpolation (transfinite interpolation)
  - [ ] general meshes
  - [ ] signed distance geometries
- [ ] Abstract Field interface `<: AbstractVector`
  - [X] spectral
  - [ ] spectral element - overload inner product
- [ ] Operator interface `<: AbstractDiffEqOperator`
  - [X] linear algebra operations
  - [X] lazy composition
  - [X] can use array reductions
  - [ ] caching
  - [ ] Gather-Scatter operator using `NNlib`
  - [ ] General interpolation operator on element-meshes
- [ ] Spaces
  - [X] Deform space
  - [X] orthogonal polynomials
  - [ ] Spectral with transforms (Fourier, Cosine, Sin, Ultraspherical, Jacobi)
  - [ ] option to solve in transformed space
  - [ ] spectral elements
- [ ] Create a distinction between `Space`, and `Discretization`
  - [ ] Space is how to represent functions
  - [ ] Discretization is how you form operators
- [ ] Boundary Condition interface
  - [X] apply "this" boundary condition based on "that" domain tag (use `Dict`)
  - [ ] Dirichlet
  - [ ] Neumann
  - [ ] Robin
- [ ] Problems
  - [ ] Describe problem with `ModelingToolkit.jl`?
    - [ ] ref for input interface https://github.com/zoemcc/DFNExperiments.jl/blob/main/src/multi_dimensional_function.jl
  - [ ] Boundary Value Problems
    - [ ] move boundary information to RHS
    - [ ] dispatch to `LinearSolve.jl`, `NonlinearSolve.jl`
  - [ ] Method of Lines
    - [ ] play nice with `OrdinaryDiffEq.jl`

