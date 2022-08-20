### `AbstractPDEInterfaces.jl`

A common interface for PDE solvers and discretizations.

PDEs are messy: weird boundary conditions, moving domains, singular operators, you name it! And solvers more so with their discretization schemes, stability charectoristics, and performance tradeoffs. Not to mention the layers of complications added by the software stack, or computer architecture, or pre/post-processing pipelines.

Numerical PDE problems drive high-performance computing. The biggest supercomputers on the planet run physics simulations on millions of MPI ranks over years of wall-time. An optimized workflow is the difference between having a solution, or going home emptyhanded, and every ounce of performance has to be squeezed. As such, highly tuned software packages that specialize on a set class of problems dominate the market. With specializations, however, generalizability and interpoeratibility take a hit.

`AbstractPDEInterfaces.jl` is written so that package authors won't need to write a method for the Laplacian (`Δ`) every time a new scheme for the gradient (`∇`) comes along. We provides abstractions over components of PDE solvers that reduce the amount of boilerplate code that any new solver would need, and improve interpoeratibility between solvers. Furthermore, a requisite for the developement of Machine Learning based PDE solvers is the ability to mix state of the art discretizations with large neural network models on a AD-compatible, accelerator-friendly test-bench.

Julia's multiple-dispatch based programming model allows for abstractly-typed composable code lets users just plug-and-play without sacrifising performace Package authors shouldn't waste their time reinventing the wheel, rather focus on what's new and interesting! 

We believe switching from 'Fourier-spectral collocation method' to `discontinuous Galerkin finite elements' should be a one-line change for the user. And users should not need to deal with inconsistent syntax between solvers for specifying boundary conditions, and forming vector-calculus operators.


workflow = physcical theory <--> PDE <--> geometry description <--> geometry discretization <--> 

current stacks:

D type of domain descriptions -> M types of meshes -> S types of spatial discr -> X types of schemes

# `AbstractSpace` interface

`AbstractPDEInterfacess.jl` contains separate abstract interaces for multidimensional domains, fields, operators, and function spaces. it is general enough that anybody can plug in their discretizations and start solving PDEs. In fact, I aim to use it with my experimental NN discretizations. All you need to do is provide a gradient operator, and a mass operator (integration).

Once you plug in your discretizations, you can do a lot of cool things like apply any random deformations to the space. AbstractPDEInterfacess.jl translate all your vector calculus operators correctly. That means the same code could solve convection diffusion on a square as well as an annulus with no extra work and basically conserved accuracy.

After describing your problem, it should spit out the right BoundaryValueProblem  or ODEProblem  that you can solve using the correct `DiffEq` package.

Goals:
- [ ] Abstract Domain interface
  - [X] Logically rectangular
  - [X] Deform domain
  - [X] Boundary tags
  - [ ] Interior tags (for multiphase flows, conjugate heat-transfer)
  - [ ] Gordon Hall interpolation (transfinite interpolation)
  - [ ] meshed domains
  - [ ] signed distance geometries
  - [ ] Time-varying domains
  - [ ] Is it possibe to just use `DomainSets.jl` and add some metadata info?
- [ ] Abstract Field interface `<: AbstractVector`
  - [X] Spectral polynomial
  - [ ] Transformation based spectral (fourier, some chebychev, etc)
  - [ ] Box/ full spectral elements -> create `SpectralElementSpaces.jl`
- [ ] Operator interface
  - [X] linear algebra operations
  - [X] lazy composition
  - [X] can use array reductions
  - [X] move as much as possible to `SciMLOperators`
  - [X] caching
  - [ ] Gather-Scatter operator using `NNlib`
  - [ ] General interpolation operator on element-meshes
- [ ] Spaces
  - [X] Deform space
  - [X] orthogonal polynomials
  - [X] option to solve in transformed space
  - [ ] Spectral with transforms (Fourier, Cosine, Sin, Ultraspherical, Jacobi)
  - [ ] spectral elements
- [X] Create a distinction between `Space`, and `Discretization`
  - [ ] Space is how to represent functions
  - [ ] Discretization is how you form operators
- [ ] Boundary Condition interface
  - [X] apply "this" boundary condition based on "that" domain boundary
  - [X] Dirichlet
  - [X] Neumann
  - [ ] Robin
- [ ] Problems
  - [ ] Describe problem with `ModelingToolkit.jl` frontend
  - [ ] Boundary Value Problems
    - [ ] move boundary information to RHS
    - [ ] dispatch to `LinearSolve.jl`, `NonlinearSolve.jl`
  - [ ] Method of Lines
    - [ ] play nice with `OrdinaryDiffEq.jl`

