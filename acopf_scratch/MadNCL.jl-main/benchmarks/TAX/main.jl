
using AmplNLReader
using LazyArtifacts

include(joinpath(@__DIR__, "..", "common.jl"))

tax_apml_path = artifact"tax"

tax_1d() = AmplModel(joinpath(tax_apml_path, "tax1D.nl"))
tax_2d() = AmplModel(joinpath(tax_apml_path, "tax2D.nl"))
tax_3d() = AmplModel(joinpath(tax_apml_path, "pTax3D.nl"))
tax_4d() = AmplModel(joinpath(tax_apml_path, "pTax4D.nl"))
tax_5d() = AmplModel(joinpath(tax_apml_path, "pTax5D.nl"))

TAX_INSTANCES = [
    (tax_1d, (), "tax_1d"),
    (tax_2d, (), "tax_2d"),
    (tax_3d, (), "tax_3d"),
    (tax_4d, (), "tax_4d"),
    (tax_5d, (), "tax_5d"),
]

function main(
    name::String="tax";
    tol::Float64=1e-6,
    max_iter::Int=1000,
    kkt_benchmark::Bool=false,
    madnlp::Bool=true,
    madncl::Bool=true,
    verbose::Bool=false
)

    results_dir = joinpath(@__DIR__, "..",  "results") |> normpath
    !isdir(results_dir) && mkpath(results_dir)

    prefix = kkt_benchmark ? "kkt_" : ""
    _madnlp = kkt_benchmark ? build_madnlp : solve_madnlp
    _madncl = kkt_benchmark ? build_madncl : solve_madncl
    _benchmark = kkt_benchmark ? run_kkt : run_benchmark

    all_solvers = []

    if madnlp
        append!(all_solvers, [
            (_madnlp, "madnlp", Ma27Solver, "ma27", MadNLP.SparseKKTSystem, "K2"),
            (_madnlp, "madnlp", Ma57Solver, "ma57", MadNLP.SparseKKTSystem, "K2"),
        ])
    end
    if madncl
        append!(all_solvers, [
            (_madncl, "madncl", Ma27Solver, "ma27", MadNLP.SparseKKTSystem, "K2"),
            (_madncl, "madncl", Ma27Solver, "ma27", MadNCL.K2rAuglagKKTSystem, "K2r"),
            (_madncl, "madncl", Ma27Solver, "ma27", MadNCL.K1sAuglagKKTSystem, "K1s"),
        ])
    end

    @info "Warm-start -- TAX"
    for (builder, algo, solver, solver_name, kkt, kkt_name) in all_solvers
        results = _benchmark("tax_1d", builder, TAX_INSTANCES; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=false)
    end

    for (builder, algo, solver, solver_name, kkt, kkt_name) in all_solvers
        @info "Benchmark $(algo)-$(solver_name)-$(kkt_name)"
        if name == "tax"
            results = _benchmark(builder, TAX_INSTANCES; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=verbose)
        else
            results = _benchmark(name, builder, TAX_INSTANCES; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=verbose)
        end
        path_results = joinpath(results_dir, "$(prefix)$(name)_$(algo)_$(kkt_name)_$(solver_name).txt")
        open(path_results, "w") do io
            writedlm(io, results)
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
