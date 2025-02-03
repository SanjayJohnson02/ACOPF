
using Pkg.Artifacts
using PowerModels
using ExaModels
using KernelAbstractions

using CUDA
using MadNLPGPU
CUDA.allowscalar(true)

pglib_path = joinpath(artifact"PGLib_opf", "pglib-opf-23.07")

include(joinpath(@__DIR__, "..", "common.jl"))
include(joinpath(@__DIR__, "acopf_model.jl"))

SELECTED_OPF_INSTANCES = [
    "pglib_opf_case1354_pegase",
    "pglib_opf_case1803_snem",
    "pglib_opf_case1888_rte",
    "pglib_opf_case1951_rte",
    "pglib_opf_case2000_goc",
    "pglib_opf_case2312_goc",
    "pglib_opf_case2383wp_k",
    "pglib_opf_case2736sp_k",
    "pglib_opf_case2737sop_k",
    "pglib_opf_case2742_goc",
    "pglib_opf_case2746wop_k",
    "pglib_opf_case2746wp_k",
    "pglib_opf_case2848_rte",
    "pglib_opf_case2853_sdet",
    "pglib_opf_case2868_rte",
    "pglib_opf_case2869_pegase",
    "pglib_opf_case3012wp_k",
    "pglib_opf_case3022_goc",
    "pglib_opf_case3120sp_k",
    "pglib_opf_case3375wp_k",
    "pglib_opf_case3970_goc",
    "pglib_opf_case4020_goc",
    "pglib_opf_case4601_goc",
    "pglib_opf_case4619_goc",
    "pglib_opf_case4661_sdet",
    "pglib_opf_case4837_goc",
    "pglib_opf_case4917_goc",
    "pglib_opf_case5658_epigrids",
    "pglib_opf_case6468_rte",
    "pglib_opf_case6470_rte",
    "pglib_opf_case6495_rte",
    "pglib_opf_case6515_rte",
    "pglib_opf_case7336_epigrids",
    "pglib_opf_case8387_pegase",
    "pglib_opf_case9241_pegase",
    "pglib_opf_case9591_goc",
    "pglib_opf_case10000_goc",
    "pglib_opf_case10192_epigrids",
    "pglib_opf_case10480_goc",
    "pglib_opf_case13659_pegase",
    "pglib_opf_case19402_goc",
    "pglib_opf_case20758_epigrids",
    "pglib_opf_case24464_goc",
    "pglib_opf_case30000_goc",
    "pglib_opf_case78484_epigrids",
]

ACOPF_INSTANCES = Dict(
    "cpu" => [],
    "cuda" => [],
)

# Build instances
for opf in SELECTED_OPF_INSTANCES
    instantiate = backend -> ac_power_model(joinpath(pglib_path, "$(opf).m"); prod=true, backend=backend)
    push!(ACOPF_INSTANCES["cpu"], (instantiate, (CPU(),), opf))
    if CUDA.has_cuda()
        push!(ACOPF_INSTANCES["cuda"], (instantiate, (CUDABackend(),), opf))
    end
end

function main(
    name::String="pglib";
    tol::Float64=1e-8,
    max_iter::Int=1000,
    kkt_benchmark::Bool=false,
    madnlp::Bool=true,
    madncl::Bool=true,
    hykkt::Bool=true,
    likkt::Bool=true,
    verbose::Bool=false,
)

    results_dir = joinpath(@__DIR__, "..",  "results") |> normpath
    !isdir(results_dir) && mkpath(results_dir)

    prefix = kkt_benchmark ? "kkt_" : ""
    _madnlp = kkt_benchmark ? build_madnlp : solve_madnlp
    _madncl = kkt_benchmark ? build_madncl : solve_madncl
    _hykkt = kkt_benchmark ? build_hykkt : solve_hykkt
    _likkt = kkt_benchmark ? build_likkt : solve_likkt
    _benchmark = kkt_benchmark ? run_kkt : run_benchmark

    all_solvers = []

    if madnlp
        append!(all_solvers, [
            (_madnlp, "madnlp", "cpu", Ma27Solver, "ma27", MadNLP.SparseKKTSystem, "K2"),
            (_madnlp, "madnlp", "cpu", Ma57Solver, "ma57", MadNLP.SparseKKTSystem, "K2"),
        ])
    end

    if madncl
        append!(all_solvers, [
            (_madncl, "madncl", "cpu", Ma27Solver, "ma27", MadNLP.SparseKKTSystem, "K2"),
            (_madncl, "madncl", "cpu", Ma27Solver, "ma27", MadNCL.K2rAuglagKKTSystem, "K2r"),
            (_madncl, "madncl", "cpu", Ma27Solver, "ma27", MadNCL.K1sAuglagKKTSystem, "K1s"),
            (_madncl, "madncl", "cuda", MadNLPGPU.CUDSSSolver, "cudss-ldl", MadNCL.K2rAuglagKKTSystem, "K2r"),
            (_madncl, "madncl", "cuda", MadNLPGPU.CUDSSSolver, "cudss-ldl", MadNCL.K1sAuglagKKTSystem, "K1s"),
        ])
    end

    if hykkt
        append!(all_solvers, [
            (_hykkt, "hykkt", "cuda", MadNLPGPU.CUDSSSolver, "cudss-ldl", HybridKKT.HybridCondensedKKTSystem, "K1"),
        ])
    end

    if likkt
        append!(all_solvers, [
            (_likkt, "likkt", "cuda", MadNLPGPU.CUDSSSolver, "cudss-ldl", MadNLP.SparseCondensedKKTSystem, "K1"),
        ])
    end

    @info "Warm-start -- OPF"
    for (builder, algo, backend, solver, solver_name, kkt, kkt_name) in all_solvers
        results = _benchmark("pglib_opf_case1354_pegase", builder, ACOPF_INSTANCES[backend]; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=false)
    end

    for (builder, algo, backend, solver, solver_name, kkt, kkt_name) in all_solvers
        @info "Benchmark $(algo)-$(backend)-$(solver_name)-$(kkt_name)"
        if name == "pglib"
            results = _benchmark(builder, ACOPF_INSTANCES[backend]; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=verbose)
        else
            results = _benchmark(name, builder, ACOPF_INSTANCES[backend]; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=verbose)
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
