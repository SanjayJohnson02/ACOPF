using DelimitedFiles
using PlotlyJS

results_dir = joinpath(@__DIR__, "results")
eps_dir = joinpath(@__DIR__, "eps")
svg_dir = joinpath(@__DIR__, "svg")

# Create the folders "eps" and "svg" if needed
isdir(eps_dir) || mkdir(eps_dir)
isdir(svg_dir) || mkdir(svg_dir)

dict_solver = Dict("ma27" => "MA27", "ma57" => "MA57", "cudss-ldl" => "cuDSS")
dict_ipm = Dict("madnlp" => "MadNLP", "madncl" => "MadNCL", "hykkt" => "HybridKKT", "likkt" => "LiftedKKT")

name_tax = ("tax_1d", "tax_2d", "tax_3d", "tax_4d", "tax_5d")

name_cops = ("bearing", "chain", "catmix", "channel", "elec", "gasoil", "marine", "methanol", "minsurf",
             "pinene", "polygon", "robot", "steering", "torsion")

name_pglib = ("pglib_opf_case240_pserc", "pglib_opf_case4619_goc", "pglib_opf_case24464_goc", "pglib_opf_case4661_sdet",
              "pglib_opf_case24_ieee_rts", "pglib_opf_case4837_goc", "pglib_opf_case2736sp_k", "pglib_opf_case4917_goc",
              "pglib_opf_case2737sop_k", "pglib_opf_case500_goc", "pglib_opf_case2742_goc", "pglib_opf_case5658_epigrids",
              "pglib_opf_case10000_goc", "pglib_opf_case2746wop_k", "pglib_opf_case57_ieee", "pglib_opf_case10192_epigrids",
              "pglib_opf_case2746wp_k", "pglib_opf_case588_sdet", "pglib_opf_case10480_goc", "pglib_opf_case2848_rte",
              "pglib_opf_case5_pjm", "pglib_opf_case118_ieee", "pglib_opf_case2853_sdet", "pglib_opf_case60_c",
              "pglib_opf_case1354_pegase", "pglib_opf_case2868_rte", "pglib_opf_case6468_rte", "pglib_opf_case13659_pegase",
              "pglib_opf_case2869_pegase", "pglib_opf_case6470_rte", "pglib_opf_case14_ieee", "pglib_opf_case30000_goc",
              "pglib_opf_case6495_rte", "pglib_opf_case162_ieee_dtc", "pglib_opf_case300_ieee", "pglib_opf_case6515_rte",
              "pglib_opf_case179_goc", "pglib_opf_case3012wp_k", "pglib_opf_case7336_epigrids", "pglib_opf_case1803_snem",
              "pglib_opf_case3022_goc", "pglib_opf_case73_ieee_rts", "pglib_opf_case1888_rte", "pglib_opf_case30_as",
              "pglib_opf_case78484_epigrids", "pglib_opf_case19402_goc", "pglib_opf_case30_ieee", "pglib_opf_case793_goc",
              "pglib_opf_case1951_rte", "pglib_opf_case3120sp_k", "pglib_opf_case8387_pegase", "pglib_opf_case197_snem",
              "pglib_opf_case3375wp_k", "pglib_opf_case89_pegase", "pglib_opf_case2000_goc", "pglib_opf_case3970_goc",
              "pglib_opf_case9241_pegase", "pglib_opf_case200_activ", "pglib_opf_case39_epri", "pglib_opf_case9591_goc",
              "pglib_opf_case20758_epigrids", "pglib_opf_case3_lmbd", "pglib_opf_case2312_goc", "pglib_opf_case4020_goc",
              "pglib_opf_case2383wp_k", "pglib_opf_case4601_goc")

# kkt_benchmark = true
for (name_benchmark, name_instances) in [("tax", name_tax), ("cops", name_cops), ("pglib", name_pglib), ("maxopf", name_pglib)]
  for instance in name_instances
    x_name = String[]
    y_build = Float64[]
    y_factorization = Float64[]
    y_backsolve = Float64[]

    for ipm in ("madnlp", "madncl", "hykkt", "likkt")
      for kkt in ("K2", "K2r", "K1s", "K1")
        for solver in ("ma27", "ma57", "cudss-ldl")
          benchmark_file = joinpath(results_dir, "kkt_$(name_benchmark)_$(ipm)_$(kkt)_$(solver).txt")
          if isfile(benchmark_file)
            data, header = readdlm(benchmark_file, '\t', header=true)
            m, n = size(data)
            @assert n == 4
            k = 1
            found = false
            while !found && k ≤ m
              instance2 = data[k,1]
              if instance == instance2
                found = true
                timer_build = data[k,2]
                timer_factorization = data[k,3]
                timer_backsolve = data[k,4]
                push!(x_name, "$(dict_ipm[ipm]) + $kkt + $(dict_solver[solver])")
                push!(y_build, timer_build)
                push!(y_factorization, timer_factorization)
                push!(y_backsolve, timer_backsolve)
              else
                k = k+1
              end
            end
          else
            benchmark_file = joinpath(results_dir, "kkt_$(instance)_$(ipm)_$(kkt)_$(solver).txt")
            if isfile(benchmark_file)
              data, header = readdlm(benchmark_file, '\t', header=true)
              m, n = size(data)
              @assert data[1,1] == instance
              @assert m == 1
              @assert n == 4
              timer_build = data[1,2]
              timer_factorization = data[1,3]
              timer_backsolve = data[1,4]
              push!(x_name, "$(dict_ipm[ipm]) + $kkt + $(dict_solver[solver])")
              push!(y_build, timer_build)
              push!(y_factorization, timer_factorization)
              push!(y_backsolve, timer_backsolve)
            end
          end
        end
      end
    end

    if !isempty(x_name)
      data = [bar(name="build", x=x_name, y=y_build),
              bar(name="factorization", x=x_name, y=y_factorization),
              bar(name="backsolve", x=x_name, y=y_backsolve)]
      layout = Layout(; title=instance)
      p = plot(data, layout)
      relayout!(p, barmode="group")

      path_eps = joinpath(eps_dir, "kkt_$(instance).eps")
      path_svg = joinpath(svg_dir, "kkt_$(instance).svg")
      savefig(p, path_eps, format="eps")
      savefig(p, path_svg, format="svg")
    end
  end
end

# kkt_benchmark = false
for (name_benchmark, name_instances) in [("tax", name_tax), ("cops", name_cops), ("pglib", name_pglib), ("maxopf", name_pglib)]
  for instance in name_instances
    x_name = String[]
    y_timer = Float64[]

    for ipm in ("madnlp", "madncl", "hykkt", "likkt")
      for kkt in ("K2", "K2r", "K1s", "K1")
        for solver in ("ma27", "ma57", "cudss-ldl")
          benchmark_file = joinpath(results_dir, "$(name_benchmark)_$(ipm)_$(kkt)_$(solver).txt")
          if isfile(benchmark_file)
            data, header = readdlm(benchmark_file, '\t', header=true)
            m, n = size(data)
            @assert n == 13
            k = 1
            found = false
            while !found && k ≤ m
              instance2 = data[k,1]
              if instance == instance2
                found = true
                push!(x_name, "$(dict_ipm[ipm]) + $kkt + $(dict_solver[solver])")
                push!(y_timer, data[k,5])
              else
                k = k+1
              end
            end
          else
            benchmark_file = joinpath(results_dir, "$(instance)_$(ipm)_$(kkt)_$(solver).txt")
            if isfile(benchmark_file)
              data, header = readdlm(benchmark_file, '\t', header=true)
              m, n = size(data)
              @assert data[1,1] == instance
              @assert m == 1
              @assert n == 13
              push!(x_name, "$(dict_ipm[ipm]) + $kkt + $(dict_solver[solver])")
              push!(y_timer, data[1,5])
            end
          end
        end
      end
    end

    colors = [
      "#FF7F50",  # Corail
      "#9ACD32",  # Vert jaunâtre
      "#6495ED",  # Bleu clair
      "#FF69B4",  # Rose clair
      "#87CEEB",  # Bleu ciel
      "#FFD700",  # Or
      "#9370DB",  # Violet moyen
      "#48D1CC"   # Turquoise moyen
    ]

    if !isempty(x_name)
      data = [bar(name="Elapsed time", x=x_name, y=y_timer, marker_color=colors)]
      y_timer_int = Int.(round.(y_timer))

      annotations = [attr(text=string(y), x=x, y=y, xanchor="center", yanchor="bottom",
                          showarrow=false, font_size=12) for (x, y) in zip(x_name, y_timer_int)]

      layout = Layout(; title=instance, annotations=annotations, yaxis=attr(showticklabels=false))
      p = plot(data, layout)

      path_eps = joinpath(eps_dir, "benchmarks_$(instance).eps")
      path_svg = joinpath(svg_dir, "benchmarks_$(instance).svg")
      savefig(p, path_eps, format="eps")
      savefig(p, path_svg, format="svg")
    end
  end
end
