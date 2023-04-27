using AbstractPDEInterfaces
using Documenter

DocMeta.setdocmeta!(AbstractPDEInterfaces, :DocTestSetup, :(using AbstractPDEInterfaces); recursive = true)

include("pages.jl")

makedocs(
         sitename = "AbstractPDEInterfaces.jl",
         authors = "Vedant Puri",
         modules = [
                    AbstractPDEInterfaces,
                    AbstractPDEInterfaces.Domains,
                    AbstractPDEInterfaces.Spaces,
                    AbstractPDEInterfaces.BoundaryConditions,
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
         #                          # canonical = "https://docs.sciml.ai/AbstractPDEInterfaces/stable/"
         #                         ),

         pages = pages
        )

deploydocs(;
           repo = "github.com/vpuri3/AbstractPDEInterfaces.jl",
           push_preview = true)
