from scipy.stats import chi2
import pandas as pd

rule process_ukbb_sum_stats:
    input:
        "resources/ukbb_sum_stats/{snp_set}/merged_ukbb_sum_stats.tsv.gz"
    output:
        "resources/ukbb_sum_stats/{snp_set}/{trait_code}.assoc"
    params:
        z_colname = lambda wildcards: '.'.join(['tstat', wildcards.trait_code]),
        sample_size_colname = lambda wildcards: '.'.join(['n_complete_samples', wildcards.trait_code])
    threads: 12
    resources:
        runtime = 8
    group: "sumher"
    script:
        "../../scripts/process_ukbb_sum_stats.R"

rule estimate_rg_with_ldak_thin_for_ukbb:
    input:
        wg_tagging_file = "results/ldak/ldak-thin/{snp_set}/whole_genome.tagging",
        sum_stats_file_A = "resources/ukbb_sum_stats/{snp_set}/{trait_A}.assoc",
        sum_stats_file_B = "resources/ukbb_sum_stats/{snp_set}/{trait_B}.assoc"
    output:
        cors_full_file = "results/ldak/ldak-thin/{snp_set}/rg/{trait_A}-{trait_B}.cors.full"
    log:
        log_file = "results/ldak/ldak-thin/{snp_set}/rg/{trait_A}-{trait_B}.log"
    params:
        output_stem = "results/ldak/ldak-thin/{snp_set}/rg/{trait_A}-{trait_B}"
    resources:
        runtime = 5
    group: "sumher"
    shell:
        """
        $ldakRoot/ldak --sum-cors {params.output_stem} --tagfile {input.wg_tagging_file} --summary {input.sum_stats_file_A} --summary2 {input.sum_stats_file_B} --allow-ambiguous YES --check-sums NO --cutoff 0.01 > {log.log_file}
        """

rule ukbb_sumher:
    input:
       sumher_files = [f"results/ldak/ldak-thin/{{snp_set}}/rg/{trait_pair}.cors.full" for trait_pair in ukbb_trait_pairs]
    output:
        "results/ldak/ldak-thin/{snp_set}/rg/compiled_ukbb_sumher_results.tsv"
    run:
        d = []

        for x in input.sumher_files:
            head, tail = os.path.split(x)

            trait_A, trait_B = re.match("(\w+)-(\w+).cors.full", tail).groups()

            with open(x, 'r') as infile:
                line = infile.readline()
                line = infile.readline()

            # Category Trait1_Her SD Trait2_Her SD Both_Coher SD Correlation SD
            _, h2_A, h2_A_se, h2_B, h2_B_se, cov, cov_se, rg, rg_se = line.split()

            rg_z = float(rg)/float(rg_se)

            rg_p = chi2.sf(rg_z**2, df = 1, loc = 0, scale = 1)

            d.append(
                {
                    'trait.A' : trait_A,
                    'trait.B' : trait_B,
                    'snp.set' : wildcards.snp_set,
                    'h2.A.obs.sr' : float(h2_A),
                    'h2.A.obs.se.sr' : float(h2_A_se),
                    'h2.B.obs.sr' : float(h2_B),
                    'h2.B.obs.se.sr' : float(h2_B_se),
                    'gcov.obs.sr' : float(cov),
                    'gcov.obs.se.sr' : float(cov_se),
                    'rg.sr' : float(rg),
                    'rg.se.sr' : float(rg_se),
                    'rg.z.sr' : rg_z,
                    'rg.p.sr' : rg_p
                }
            )

            print(x)

        pd.DataFrame(d).to_csv(output[0], sep = '\t', index = False)
