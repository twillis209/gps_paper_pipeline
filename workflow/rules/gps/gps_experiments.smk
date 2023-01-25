rule compute_per_chr_gps_for_trait_pair:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_per_chr/{chr}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}_per_chr/{chr}_{ecdf}_pert_{pert}_gps_value.tsv"
    params:
        epsilon_multiple = 2.0
    threads: lambda wildcards: 12 if wildcards.ecdf == 'naive' else 1
    resources:
        runtime = 20
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a P.A -b P.B -c {wildcards.effect_blocks_A} -d {wildcards.effect_blocks_B} -f {wildcards.ecdf} -n {threads} -p {wildcards.pert} -e {params.epsilon_multiple} -o {output}"

rule compute_per_chr_gps_with_intermediate_values:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_per_chr/{chr}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}_per_chr/{chr}_{ecdf}_pert_{pert}_gps_intermediates.tsv",
    threads: 1
    resources:
        runtime = 30
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/fitAndEvaluateEcdfsCLI -i {input} -a P.A -b P.B -f {wildcards.ecdf} -n {threads} -o {output}"

rule permute_per_chr_gps_for_trait_pair:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}_per_chr/{chr}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}_per_chr/{draws}_permutations/{chr}.tsv",
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}",
        no_of_pert_iterations = 1
    threads: 12
    resources:
        mem_mb = get_mem_mb,
        runtime = get_permute_time,
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {params.a_colname} -b {params.b_colname} -c {threads} -n {wildcards.draws} -p {params.no_of_pert_iterations}"

rule test_per_chr_gps:
    input:
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/seed_48801_tags_1-2_per_chr/chr22_naive_pert_0_gps_value.tsv"
#        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/seed_48801_tags_1-2_per_chr/chr22_naive_pert_0_gps_value.tsv"

rule compute_mean_stat_for_trait_pair:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}/mean_stat/{ecdf}_pert_{pert}_gps_value.tsv"
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}",
        epsilon_multiple = 2.0
    threads: lambda wildcards: 12 if wildcards.ecdf == 'naive' else 1
    resources:
        runtime = 5
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input} -a {params.a_colname} -b {params.b_colname} -c {wildcards.effect_blocks_A} -d {wildcards.effect_blocks_B} -f {wildcards.ecdf} -n {threads} -p {wildcards.pert} -e {params.epsilon_multiple} -o {output} -s mean"

rule permute_mean_stat_for_trait_pair:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}/mean_stat/{draws}_permutations/permutations.tsv",
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}",
        no_of_pert_iterations = 1
    threads: 12
    resources:
        mem_mb = get_mem_mb,
        runtime = get_permute_time,
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input} -o {output} -a {params.a_colname} -b {params.b_colname} -c {threads} -n {wildcards.draws} -p {params.no_of_pert_iterations} -s mean"

rule write_out_mean_stat_components_for_trait_pair:
    input:
        "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}-{tag_B}.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}/mean_stat/{ecdf}_pert_{pert}_components.tsv"
    params:
        a_colname = lambda wildcards: f"p.{wildcards.tag_A}",
        b_colname = lambda wildcards: f"p.{wildcards.tag_B}",
        epsilon_multiple = 2.0
    threads: lambda wildcards: 12 if wildcards.ecdf == 'naive' else 1
    resources:
        runtime = 5
    group: "gps"
    shell:
        "workflow/scripts/gps_cpp/build/apps/writeMeanStatComponentsCLI -i {input} -a {params.a_colname} -b {params.b_colname} -f {wildcards.ecdf} -n {threads} -p {wildcards.pert} -e {params.epsilon_multiple} -o {output}"

rule compute_pvalue_for_mean_stat_for_trait_pair:
    input:
        permutations = "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}/mean_stat/{draws}_permutations/permutations.tsv",
        mean_stat = "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}/mean_stat/{ecdf}_pert_{pert}_gps_value.tsv"
    output:
        "results/gps/simgwas/{no_reps}_reps/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_tags_{tag_A,\d+}-{tag_B,\d+}/mean_stat/{draws}_permutations/{ecdf}_pert_{pert}_pvalue.tsv",
    group: "sn_pvalue"
    script:
        "../scripts/compute_pvalue_for_mean_stat.R"

rule compute_mean_stat_for_ukbb_trait_pair:
    input:
      "resources/ukbb_sum_stats/{trait_A}.done",
      "resources/ukbb_sum_stats/{trait_B}.done",
      sum_stats_file = "resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv",
    output:
        "results/{join}/{snp_set}/window_{window}_step_{step}/mean_stat/{trait_A}-{trait_B}/{ecdf}_pert_{pert}_mean_value.tsv"
    threads: 1
    group: "permute"
    resources:
        runtime = 10
    shell:
      "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {wildcards.trait_A} -d {wildcards.trait_B} -f {wildcards.ecdf} -n {threads} -p {wildcards.pert} -o {output} -s mean"

rule permute_ukbb_trait_pair_for_mean_stat:
    input:
      "resources/ukbb_sum_stats/{trait_A}.done",
      "resources/ukbb_sum_stats/{trait_B}.done",
      sum_stats_file = "resources/pruned_sum_stats/{join}/{snp_set}/window_{window}_step_{step}/pruned_merged_sum_stats.tsv",
    output:
        "results/{join}/{snp_set,all_pruned_snps|sans_mhc}/window_{window}_step_{step}/mean_stat/{trait_A}-{trait_B}/{draws}_permutations/permutations.tsv"
    params:
        no_of_perturbations = 1
    threads: 12
    resources:
        mem_mb = get_mem_mb,
        runtime = get_permute_time,
    group: "permute"
    shell:
      "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {wildcards.trait_A} -b {wildcards.trait_B} -c {threads} -n {wildcards.draws} -p {params.no_of_perturbations} -s mean"

rule compute_pvalue_for_mean_stat_for_ukbb_trait_pair:
    input:
        mean_stat = "results/{join}/{snp_set}/window_{window}_step_{step}/mean_stat/{trait_A}-{trait_B}/{ecdf}_pert_{pert}_mean_value.tsv",
        permutations = "results/{join}/{snp_set,all_pruned_snps|sans_mhc}/window_{window}_step_{step}/mean_stat/{trait_A}-{trait_B}/{draws}_permutations/permutations.tsv"
    output:
        "results/{join}/{snp_set,all_pruned_snps|sans_mhc}/window_{window}_step_{step}/mean_stat/{trait_A}-{trait_B}/{draws}_permutations/{ecdf}_pert_{pert}_pvalue.tsv"
    resources:
        runtime = 5
    group: "permute"
    shell:
        "../scripts/compute_pvalue_for_mean_stat.R"

rule collate_mean_stat_pvalues_for_ukbb:
    input:
        pvalue_files = ["results/{join}/{snp_set}/window_{window}_step_{step}/mean_stat/%s/{draws}_permutations/{ecdf}_pert_{pert}_pvalue.tsv" % x for x in ukbb_trait_pairs]
    output:
        combined_pvalue_file = "results/gps/combined/{join}/{snp_set}/window_{window}_step_{step}/mean_stat/{draws}_permutations/pvalues_{ecdf}_pert_{pert}.tsv"
    resources:
        runtime = 20
    run:
        with open(output.combined_pvalue_file, 'w') as outfile:
            outfile.write(("\t".join(["Trait_A", "Trait_B", "xi", "omega", "alpha", "GPS", "pval"]))+"\n")
            for i,x in enumerate(input.pvalue_files):
                with open(x, 'r') as infile:
                    line = infile.readline()
                    line = infile.readline()

                outfile.write(line)

rule write_out_components_for_one_s400_simulation:
    input:
        "results/gps/simgwas/400_reps/randomised/10000_10000_10000_10000/s400_s400_s0/window_1000kb_step_50/seed_48801_tags_1-2/mean_stat/pp_pert_1_components.tsv"
