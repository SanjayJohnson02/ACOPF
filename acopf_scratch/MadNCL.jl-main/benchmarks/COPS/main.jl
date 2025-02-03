
using COPSBenchmark
using NLPModelsJuMP
using ExaModels
using KernelAbstractions

using CUDA
using MadNLPGPU

# FIXME: we have to force allowscalar to true to compute dual residual at
# the optimum solution.
CUDA.allowscalar(true)

include(joinpath(@__DIR__, "..", "common.jl"))

SELECTED_COPS_INSTANCES = [
    (COPSBenchmark.bearing_model , (400, 400), "bearing" ),
    (COPSBenchmark.catmix_model  , (3200,)   , "catmix"  ),
    (COPSBenchmark.channel_model , (3200,)   , "channel" ),
    (COPSBenchmark.elec_model    , (400,)    , "elec"    ),
    (COPSBenchmark.gasoil_model  , (3200,)   , "gasoil"  ),
    (COPSBenchmark.marine_model  , (3200,)   , "marine"  ),
    # (COPSBenchmark.methanol_model, (3200,)   , "methanol"),
    (COPSBenchmark.minsurf_model , (400, 400), "minsurf" ),
    (COPSBenchmark.pinene_model  , (3200,)   , "pinene"  ),
    (COPSBenchmark.polygon_model , (100,)    , "polygon" ),
    (COPSBenchmark.robot_model   , (3200,)   , "robot"   ),
    (COPSBenchmark.steering_model, (3200,)   , "steering"),
    (COPSBenchmark.torsion_model , (200, 200), "torsion" ),
]

COPS_INSTANCES = Dict(
    "cpu" => [],
    "cuda" => [],
)

# Build instances
for (instance, params, name) in SELECTED_COPS_INSTANCES
    instantiate = backend -> ExaModels.ExaModel(instance(params...); prod=true, backend=backend)
    push!(COPS_INSTANCES["cpu"], (instantiate, (CPU(),), name))
    if CUDA.has_cuda()
        push!(COPS_INSTANCES["cuda"], (instantiate, (CUDABackend(),), name))
    end
end

function main(
    name::String="cops";
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

    @info "Warm-start -- COPS"
    for (builder, algo, backend, solver, solver_name, kkt, kkt_name) in all_solvers
        results = _benchmark("polygon", builder, COPS_INSTANCES[backend]; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=false)
    end

    for (builder, algo, backend, solver, solver_name, kkt, kkt_name) in all_solvers
        @info "Benchmark $(algo)-$(backend)-$(solver_name)-$(kkt_name)"
        if name == "cops"
            results = _benchmark(builder, COPS_INSTANCES[backend]; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=verbose)
        else
            results = _benchmark(name, builder, COPS_INSTANCES[backend]; kkt_system=kkt, linear_solver=solver, tol=tol, max_iter=max_iter, verbose=verbose)
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
