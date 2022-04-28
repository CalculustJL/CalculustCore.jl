## AbstractPDEs.jl

Tired of writing boilerplate code for PDE solvers? Just want to focus on discretizations? Want something that plays nice with `DiffEq` ecosystem? AbstractPDEs.jl contains separate abstract interaces for multidimensional domains, fields, operators, and function spaces. it is general enough that anybody can plug in their discretizations and start solving PDEs. in fact i aim to use it with my experimental NN discretizations. All you need to do is provide a gradient operator, and a mass operator (integration).

Once you plug in your discretizations, you can do a lot of cool things like apply any random deformations to the space. AbstractPDEs.jl translate all your vector calculus operators correctly. That means the same code could solve convection diffusion on a square as well as an annulus with no extra work and basically conserved accuracy.

After describing your problem, it should spit out the right BoundaryValueProblem  or ODEProblem  that you can solve using the correct `DiffEq` package.


Goals:
- [ ] Domain interface
  - [X] Logically rectangular
  - [X] Deform domain
  - [ ] Gordon Hall interpolation
  - [ ] general meshes
  - [ ] signed distance geometries
- [ ] Field interface `<: AbstractVector`
  - [x] spectral
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
  - [ ] Spectral with transforms (Fourier, Cosine, Sin, Ultraspherical,)
  - [ ] option to solve in transformed space
  - [ ] spectral elements
- [ ] Create a distinction between `Space`, and `Discretization`
- [ ] Boundary Condition interface
  - [ ] apply "this" boundary condition based on "that" domain tag (use `Dict`)
- [ ] Problems
  - [ ] Describe problem with `ModelingToolkit.jl`?
  - [ ] Boundary Value Problems
    - [ ] move boundary information to RHS
    - [ ] dispatch to `LinearSolve.jl`, `NonlinearSolve.jl`
  - [ ] Method of Lines
    - [ ] play nice with `OrdinaryDiffEq.jl`


