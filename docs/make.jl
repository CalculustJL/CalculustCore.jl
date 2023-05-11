using CalculustCore
using Documenter

DocMeta.setdocmeta!(CalculustCore, :DocTestSetup, :(using CalculustCore);
                    recursive = true)

include("pages.jl")

makedocs(sitename = "CalculustCore.jl",
         authors = "Vedant Puri",
         modules = [
             CalculustCore,
             CalculustCore.Domains,
             CalculustCore.Spaces,
             CalculustCore.BoundaryConditions,
         ],
         clean = true,
         doctest = false, # TODO - set to true later
         linkcheck = true,
         strict = [
             :doctest,
             :linkcheck,
             :parse_error,
             :example_block,
             # Other available options are
             # :autodocs_block, :cross_references, :docs_block, :eval_block, :example_block, :footnote, :meta_block, :missing_docs, :setup_block
         ],
         # format = Documenter.HTML(
         #                          # assets = ["assets/favicon.ico"],
         #                          # canonical = "https://docs.calculust.dev/CalculustCore/stable/"
         #                         ),

         pages = pages)

deploydocs(;
           repo = "github.com/vpuri3/CalculustCore.jl",
           push_preview = true)
