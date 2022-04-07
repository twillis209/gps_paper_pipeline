import os

rule compute_gps_for_sim_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}{tag_B}.tsv.gz"
    output:
        temp("results/gps/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}{tag_B}_gps_value.tsv")
    params:
        a_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_A)+1),
        b_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_B)+1)
    shell:
        "workflow/scripts/gps_cpp/build/apps/computeGpsCLI -i {input.sum_stats_file} -a {params.a_colname} -b {params.b_colname} -c {wildcards.effect_blocks_A} -d {wildcards.effect_blocks_B} -o {output}"

rule permute_sim_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}{tag_B}.tsv.gz"
    output:
        "results/gps/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{draws,\d+}_permutations/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}{tag_B}.tsv"
    params:
        a_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_A)+1),
        b_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_B)+1)
    threads: 8
    resources:
        mem_mb = get_mem_mb,
        time = get_permute_time,
    shell:
        "workflow/scripts/gps_cpp/build/apps/permuteTraitsCLI -i {input.sum_stats_file} -o {output} -a {params.a_colname} -b {params.b_colname} -c {threads} -n {wildcards.draws}"

rule fit_gev_and_compute_gps_pvalue_for_sim_pair:
    input:
        gps_file = "results/gps/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}{tag_B}_gps_value.tsv",
        perm_file = "results/gps/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{draws,\d+}_permutations/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}{tag_B}.tsv"
    output:
        "results/gps/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}{tag_B}_gps_pvalue.tsv"
    shell:
      "Rscript workflow/scripts/fit_gev_and_compute_gps_pvalue.R -g {input.gps_file} -p {input.perm_file} -a {wildcards.effect_blocks_A} -b {wildcards.effect_blocks_B} -o {output}"

rule compute_hoeffdings_for_sim_pair:
    input:
        sum_stats_file = "results/simgwas/simulated_sum_stats/whole_genome_sum_stats/randomised/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}/window_{window}_step_{step}/seed_{seed}_pruned_sum_stats_tags_{tag_A}{tag_B}.tsv.gz"
    output:
        "results/hoeffdings/simgwas/randomised/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{shared_effect_blocks}_seed_{seed}_tags_{tag_A}{tag_B}_hoeffdings.tsv"
    params:
        a_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_A)+1),
        b_colname = lambda wildcards: "p.%d" % (tags.index(wildcards.tag_B)+1)
    shell:
    """
    Rscript workflow/scripts/compute_hoeffdings.R -i {input.sum_stats_file} -a {params.a_colname} -b {params.b_colname} -o {output} -nt 1
    sed -i 's/{params.a_colname}/{wildcards.effect_blocks_A}_{wildcards.shared_effect_blocks}/' {output}
    sed -i 's/{params.b_colname}/{wildcards.effect_blocks_B}_{wildcards.shared_effect_blocks}/' {output}
    """

#    output:
#        "results/gps/simgwas/window_1000kb_step_50/combined_gps_results.tsv"
#    run:
#        with open(output[0], 'w') as outfile:
#
#            for i,x in enumerate(input):
#                with open(x, 'r') as infile:
#                    lines = infile.readlines()
#
#                    if i == 0:
#                        outfile.write("blocks.A\tblocks.B\teffect_size\tncases.A\tncontrols.A\tncases.B\tncontrols.B\ttag_pair\t%s" % lines[0])
#
#                    head, tail = os.path.split(x)
#
#                    ncases_A, ncontrols_A, ncases_B, ncontrols_B = re.search("\d+_\d+_\d+_\d+", head).group().split('_')
#
#                    effect_size = re.search("[smlvhi]", tail).group()
#
#                    effect_blocks_A, effect_blocks_B = tail.split('_')[:2]
#
#                    tag_pair = tail.split('_')[4]
#
#                    outfile.write(f"{effect_blocks_A}\t{effect_blocks_B}\t{effect_size}\t{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{tag_pair}\t{lines[1]}")
#
#rule compute_hoeffdings_for_sim_pair:
#    input:
#        sum_stats_file = "results/simgwas/simulated_sum_stats/pruned/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_sum_stats.tsv"
#    output:
#        "results/hoeffdings/simgwas/window_{window}_step_{step}/{ncases_A}_{ncontrols_A}_{ncases_B}_{ncontrols_B}/{effect_blocks_A}_{effect_blocks_B}_{tag_pair}_hoeffdings.tsv"
#    params:
#        a_colname = lambda wildcards: "p.%d.A" % (tags.index(wildcards.tag_pair[0])+1),
#        b_colname = lambda wildcards: "p.%d.B" % (tags.index(wildcards.tag_pair[1])+1)
#    shell:
#        """
#        Rscript workflow/scripts/compute_hoeffdings.R -i {input.sum_stats_file} -a {params.a_colname} -b {params.b_colname} -o {output} -nt 1
#        sed -i 's/p\.1\.A/{wildcards.effect_blocks_A}/' {output}
#        sed -i 's/p\.2\.B/{wildcards.effect_blocks_B}/' {output}
#        """
#
#rule run_hoeffdings_sim_data:
#    input:
#        ["results/hoeffdings/simgwas/window_1000kb_step_50/%d_%d_%d_%d/%s_%s_hoeffdings.tsv" % (i, i, i, i, j, k) for i in sample_sizes for k in tag_pairs for j in effect_block_pairs]
#    output:
#        "results/hoeffdings/simgwas/window_1000kb_step_50/combined_hoeffdings_results.tsv"
#    run:
#        with open(output[0], 'w') as outfile:
#            outfile.write("ncases.A\tncontrols.A\tncases.B\tncontrols.B\todds_ratio.A\todds_ratio.B\tblocks.A\tblocks.B\tno_shared_blocks\ttag_pair\tn\tpval\n")
#            for x in input:
#                head, tail = os.path.split(x)
#                head_res = re.match("results/hoeffdings/simgwas/window_1000kb_step_50/(\d+)_(\d+)_(\d+)_(\d+)", head)
#                ncases_A = int(head_res.group(1))
#                ncontrols_A = int(head_res.group(2))
#                ncases_B = int(head_res.group(3))
#                ncontrols_B = int(head_res.group(4))
#
#                effect_blocks_wc_A, effect_blocks_wc_B, tag_pair = tail.split('_')[:3]
#
#                odds_ratios_A = []
#                odds_ratios_B = []
#
#                effect_blocks_A = []
#                effect_blocks_B = []
#
#                for y in effect_blocks_wc_A.split('+'):
#                    if effect_blocks_wc_A != 'null':
#                        block_match = re.match('^(\d+)-(.+)', y)
#
#                        chrom = int(block_match.group(1))
#
#                        range_match = re.match('([smlvhi])(\d+):(\d+)', block_match.group(2))
#
#                        odds_ratios_A.append(odds_ratio_dict[range_match.group(1)])
#
#                        effect_blocks_A += ["%d-%s%d" % (chrom, range_match.group(1), z) for z in range(int(range_match.group(2)), int(range_match.group(3))+1) if z in block_dict[chrom]]
#                    else:
#                        odds_ratios_A.append(odds_ratio_dict['n'])
#
#                for y in effect_blocks_wc_B.split('+'):
#                    if effect_blocks_wc_B != 'null':
#                        block_match = re.match('^(\d+)-(.+)', y)
#
#                        chrom = int(block_match.group(1))
#
#                        range_match = re.match('([smlvhi])(\d+):(\d+)', block_match.group(2))
#
#                        odds_ratios_B.append(odds_ratio_dict[range_match.group(1)])
#
#                        effect_blocks_B += ["%d-%s%d" % (chrom, range_match.group(1), z) for z in range(int(range_match.group(2)), int(range_match.group(3))+1) if z in block_dict[chrom]]
#                    else:
#                        odds_ratios_B.append(odds_ratio_dict['n'])
#
#                odds_ratios_A = ','.join(set([str(x) for x in odds_ratios_A]))
#                odds_ratios_B = ','.join(set([str(x) for x in odds_ratios_B]))
#
#                no_shared_blocks = len([z for z in effect_blocks_A if z in effect_blocks_B])
#
#                with open(x, 'r') as infile:
#                    line = infile.readline()
#                    line = infile.readline()
#
#                    n = line.split()[2]
#                    pval = line.split()[5]
#                    outfile.write(f"{ncases_A}\t{ncontrols_A}\t{ncases_B}\t{ncontrols_B}\t{odds_ratios_A}\t{odds_ratios_B}\t{effect_blocks_wc_A}\t{effect_blocks_wc_B}\t{no_shared_blocks}\t{tag_pair}\t{n}\t{pval}\n")
