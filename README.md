## AbstractPDEs.jl

Aim is to write a general API that any discretization can be plugged into

Goals:
- [ ] Domain interface
  - [X] Logically rectangular
  - [X] Deform domain
  - [ ] Gordon Hall interpolation
  - [ ] general meshes
- [ ] Field interface
  - [x] spectral
  - [ ] spectral element - overload inner product
- [ ] Operator interface
  - [X] `<: AbstractSciMLOperators`
  - [X] linear algebra operations
  - [ ] caching
  - [ ] Gather-Scatter operator using `NNlib`
  - [ ] General interpolation operator on element-meshes
- [ ] Spaces
  - [X] Deform space
  - [X] Spectral polynomials
  - [ ] Spectral with transforms (Fourier, Cosine, Sin, Ultraspherical,)
- [ ] Boundary Condition interface
  - [ ] apply "this" boundary condition based on "that" domain tag
- [ ] Problems
  - [ ] Boundary Value Problems
    - [ ] move boundary information to RHS
    - [ ] dispatch to `LinearSolve.jl`, `NonlinearSolve.jl`
  - [ ] Method of Lines
    - [ ] play nice with `OrdinaryDiffEq.jl`
