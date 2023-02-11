rule fit_gev_to_increasing_n:
    input:
        "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/10000_permutations/{trait_A}-{trait_B}.tsv"
    output:
        "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_1000-10000_permutations_estimates.tsv"
    shell:
      "Rscript workflow/scripts/ukbb/fit_gev_to_increasing_n.R -a {wildcards.trait_A} -b {wildcards.trait_B} -p {input} -n 500 1000 2000 3000 4000 5000 6000 7000 8000 9000 10000 -o {output}"

rule fit_gev_to_increasing_n_for_selected_ukbb_trait_pairs:
    input:
        ["results/ukbb_with_mhc/snps_only/window_1000kb_step_50_r2_0_2/%s_1000-10000_permutations_estimates.tsv" % x for x in trait_pairs_for_increasing_perm_fits]
    output:
        "results/ukbb_with_mhc/snps_only/window_1000kb_step_50_r2_0_2/1000-10000_combined_permutations_estimates.tsv"
    shell:
        """
        echo -e 'trait_A\ttrait_B\tn\tloc\tloc.sd\tscale\tscale.sd\tshape\tshape.sd' >> {output}
        for x in {input}; do
        cat <(tail -n +2 $x) >> {output}
        done
        """

rule plot_gev_estimates_for_increasing_n:
    input:
        "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_1000-10000_permutations_estimates.tsv"
    output:
        "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_1000-10000_permutations_estimates.png"
    shell:
        "Rscript workflow/scripts/ukbb/plot_gev_estimates_for_increasing_n.R -f {input} -o {output}"

rule plot_gev_estimates_for_increasing_n_for_selected_ukbb_trait_pairs:
    input:
        ["results/{snp_set}/{variant_set}/window_1000kb_step_50_r2_0_2/%s_1000-10000_permutations_estimates.png" % x for x in trait_pairs_for_increasing_perm_fits]

rule plot_gev_estimates_for_increasing_no_snps:
    input:
        [f"results/{{snp_set}}/{{variant_set}}/window_{{window}}_step_{{step}}_r2_{{r2}}/{x}_snps/{{trait_A}}-{{trait_B}}_{{draws}}_permutations_estimates.tsv" for x in [10000, 50000, 100000, 200000, 300000, 400000]]
    output:
        "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_variable_no_snps_estimates.png"
    params:
        fit_file_string = " ".join([f"results/{{snp_set}}/{{variant_set}}/window_{{window}}_step_{{step}}_r2_{{r2}}/{x}_snps/{{trait_A}}-{{trait_B}}_{{draws}}_permutations_estimates.tsv" for x in [10000, 50000, 100000, 200000, 300000, 400000]])
    shell:
      "Rscript workflow/scripts/ukbb/plot_gev_estimates_for_increasing_no_snps.R -i {params.fit_file_string} -n 10000 50000 100000 200000 300000 400000 -o {output}"

rule plot_gev_estimates_for_increasing_no_snps_for_selected_ukbb_traits:
    input:
        ["results/ukbb_with_mhc/snps_only/window_1000kb_step_50_r2_0_2/%s_3000_permutations_variable_no_snps_estimates.png" % x for x in trait_pairs_for_increasing_perm_fits]

rule plot_null_dists_to_compare:
    input:
        perm_file = "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/{trait_A}-{trait_B}.tsv",
        fit_file = "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_gps_pvalue.tsv"
    output:
        exp1_null = "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_exp1_null_dist.png",
        gev_null = "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_gev_null_dist.png",
        exp1_gev_combined = "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_null_dists.png"
    shell:
      "Rscript workflow/scripts/ukbb/plot_null_dists_to_compare.R -f {input.fit_file} -p {input.perm_file} --exp1_null {output.exp1_null} --gev_null {output.gev_null} --exp1_gev_combined {output.exp1_gev_combined}"

rule plot_all_null_dists_to_compare:
    input:
        ["results/ukbb_with_mhc/snps_only/window_1000kb_step_50_r2_0_2/%s_3000_permutations_null_dists.png" % x for x in ukbb_trait_pairs]

rule plot_gof_plots:
    input:
        perm_file = "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{draws}_permutations/{trait_A}-{trait_B}.tsv",
        fit_file = "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_gps_pvalue.tsv"
    output:
        "results/{snp_set}/{variant_set}/window_{window}_step_{step}_r2_{r2}/{trait_A}-{trait_B}_{draws}_permutations_gof_plots.png"
    shell:
      "Rscript workflow/scripts/ukbb/plot_gev_gof_plots.R -f {input.fit_file} -p {input.perm_file} -l {wildcards.trait_A}-{wildcards.trait_B} -o {output}"

rule plot_all_gof_plots:
    input:
        ["results/{snp_set}/{variant_set}/window_1000kb_step_50_r2_0_2/%s_3000_permutations_gof_plots.png" % x for x in ukbb_trait_pairs]

rule compile_naive_gps_values:
    input:
        [f"results/{{snp_set}}/{{variant_set}}/window_1000kb_step_50_r2_0_2/{trait_pair}_naive_gps_value.tsv" for trait_pair in ukbb_trait_pairs]
    output:
        "results/{snp_set}/{variant_set}/window_1000kb_step_50_r2_0_2/compiled_naive_gps_values.tsv"
    shell:
        """
        echo -e "Trait_A\tTrait_B\tGPS" >{output}

        for x in {input}; do
            tail -n 1 $x >>{output}
        done
        """
