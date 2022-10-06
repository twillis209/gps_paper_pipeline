rule preprocess_ukbb_sum_stats_trait:
    input:
        flag_file = "resources/ukbb_sum_stats/{trait}.done",
        sum_stats_file = "resources/ukbb_sum_stats/{snp_set}/merged_ukbb_sum_stats.tsv.gz",
        snplist_file = "resources/ldsc/ldsc_snps.tsv.gz"
    params:
        pval_col = lambda wildcards: f"pval.{wildcards.trait}",
        tstat_col = lambda wildcards: f"tstat.{wildcards.trait}",
        n_col = lambda wildcards: f"n_complete_samples.{wildcards.trait}",
    output:
        temp("results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait}_preprocessed_sum_stats.tsv.gz")
    threads: 4
    resources:
        tmpdir = 'tmp'
    script: "../../scripts/ldsc/preprocess_ukbb_sum_stats_for_ldsc.R"

rule munge_ukbb_sum_stats:
    input:
        "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait}_preprocessed_sum_stats.tsv.gz"
    output:
        "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait}.tsv.sumstats.gz"
    params:
        output_filename = "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait}.tsv",
        # NB: The '0' below gives the null value for beta
        signed_sumstats_col = lambda wildcards: f"tstat.{wildcards.trait},0",
        pvalue_col = lambda wildcards: f"pval.{wildcards.trait}",
        n_col = lambda wildcards: f"n_complete_samples.{wildcards.trait}",
        id_col = "SNP"
    log:
        log = "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait}.log"
    threads: 1
    resources:
        runtime = 10
    conda:
        "envs/ldsc.yaml"
    shell:
        """
        python $ldsc/munge_sumstats.py --sumstats {input} --N-col {params.n_col} --snp {params.id_col} --out {params.output_filename} --signed-sumstats {params.signed_sumstats_col} --p {params.pvalue_col} --a1 a1 --a2 a2 --frq EUR;
        """

rule estimate_rg_for_ukbb_sum_stats:
    input:
        ["resources/ldsc/eur_w_ld_chr/%d.l2.ldscore.gz" % i for i in range(1,23)],
        ["resources/ldsc/eur_w_ld_chr/%d.l2.M_5_50" % i for i in range(1,23)],
        sum_stats_A = "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait_A}.tsv.sumstats.gz",
        sum_stats_B = "results/ldsc/munged_sum_stats/ukbb/{snp_set}/{trait_B}.tsv.sumstats.gz"
    output:
        "results/ldsc/rg/ukbb/{snp_set}/{h2_intercept,fixed|free}_h2_{rg_intercept,fixed|free}_rg_intercept/{trait_A}-{trait_B}.log"
    log:
        log = "results/ldsc/rg/ukbb/{snp_set}/{h2_intercept,fixed|free}_h2_{rg_intercept,fixed|free}_rg_intercept/{trait_A}-{trait_B}.actual_log"
    params:
        log_file_par = "results/ldsc/rg/ukbb/{snp_set}/{h2_intercept,fixed|free}_h2_{rg_intercept,fixed|free}_rg_intercept/{trait_A}-{trait_B}",
        # NB: Trailing '/' is needed in ld_score_root
        ld_score_root = "resources/ldsc/eur_w_ld_chr/",
        h2_intercept = lambda wildcards: "--intercept-h2 1,1" if wildcards.h2_intercept == "fixed" else "",
        rg_intercept = lambda wildcards: "--intercept-gencov 0,0" if wildcards.rg_intercept == "fixed" else ""
    resources:
        runtime = 2
    group: "ldsc_hoeffding_sumher_gps_sans_permutation"
    conda:
        "envs/ldsc.yaml"
    shell:
        # Hacky fix to retain 'log' file (the output is regrettably so-called), which we need whether or not the estimation process failed
        "python $ldsc/ldsc.py --rg {input.sum_stats_A},{input.sum_stats_B} --ref-ld-chr {params.ld_score_root} --w-ld-chr {params.ld_score_root} --out {params.log_file_par} {params.h2_intercept} {params.rg_intercept} >{log.log} || true"

rule ukbb_ldsc:
    input:
        [f"results/ldsc/rg/ukbb/{{snp_set}}/fixed_h2_free_rg_intercept/{trait_pair}.log" for trait_pair in ukbb_trait_pairs]
    output:
        "results/ldsc/rg/ukbb/{snp_set}/fixed_h2_free_rg_intercept/compiled_results.tsv"
    run:
        d = []

        h2_regex = r"Total Observed scale h2: (.+)\s+\((.+)\)"
        int_regex = r"Intercept: (.+)\s+\((.+)\)"

        gcov_regex = r"Total Observed scale gencov: (.+)\s+\((.+)\)"
        gcov_zprod_regex = r"Mean z1\*z2: (.+)"

        for x in input:
            m = re.match(r"results/ldsc/rg/ukbb/(?P<snp_set>\w+)/fixed_h2_free_rg_intercept/(?P<trait_A>\w+)-(?P<trait_B>\w+)\.log", x)

            with open(x, 'r') as infile:
                line = infile.readline()

                # TODO fix these for the null case
                while re.match(h2_regex, line) is None and re.match('ERROR', line) is None:
                    line = infile.readline()

                if re.match('ERROR', line):
                    d.append(
                        {
                            'trait.A' : m.group('trait_A'),
                            'trait.B' : m.group('trait_B'),
                            'snp.set' : wildcards.snp_set,
                            'h2.A.obs.ldsc' : nan,
                            'h2.A.obs.se.ldsc' : nan,
                            'h2.B.obs.ldsc' : nan,
                            'h2.B.obs.se.ldsc' : nan,
                            'gcov.obs.ldsc' : nan,
                            'gcov.obs.se.ldsc' : nan,
                            'rg.ldsc' : nan,
                            'rg.se.ldsc' : nan,
                            'rg.z.ldsc' : nan,
                            'rg.p.ldsc' : nan
                        }
                    )

                else:
                    h2_match_A = re.match(h2_regex, line)
                    h2_A = float(h2_match_A.group(1))
                    h2_A_se = float(h2_match_A.group(2))

                    line = infile.readline()
                    line = infile.readline()
                    line = infile.readline()

                    h2_int_A_match = re.match(int_regex, line)

                    if h2_int_A_match:
                        h2_int_A = float(h2_int_A_match.group(1))
                        h2_int_A_se = float(h2_int_A_match.group(2))
                    elif 'constrained to 1.' in line:
                        h2_int_A = 1.0
                        h2_int_A_se = nan
                    else:
                        raise Exception("No match for h2_B int_regex")

                    while re.match(h2_regex, line) is None:
                        line = infile.readline()

                    h2_match_B = re.match(h2_regex, line)
                    h2_B = float(h2_match_B.group(1))
                    h2_B_se = float(h2_match_B.group(2))

                    line = infile.readline()
                    line = infile.readline()
                    line = infile.readline()

                    h2_int_B_match = re.match(int_regex, line)

                    if h2_int_B_match:
                            h2_int_B = float(h2_int_B_match.group(1))
                            h2_int_B_se = float(h2_int_B_match.group(2))
                    elif 'constrained to 1.' in line:
                            h2_int_B = 1.0
                            h2_int_B_se = nan
                    else:
                            raise Exception("No match for h2_A int_regex")

                    while re.match(gcov_regex, line) is None:
                        line = infile.readline()

                    gcov_match = re.match(gcov_regex, line)
                    gcov = float(gcov_match.group(1))
                    gcov_se = float(gcov_match.group(2))

                    line = infile.readline()

                    gcov_zprod_match = re.match(gcov_zprod_regex, line)
                    gcov_zprod = float(gcov_zprod_match.group(1))

                    line = infile.readline()

                    gcov_int_match = re.match(int_regex, line)

                    if gcov_int_match:
                        gcov_int = float(gcov_int_match.group(1))
                        gcov_int_se = float(gcov_int_match.group(2))
                    elif 'constrained to 0.' in line:
                        gcov_int = 0.0
                        gcov_int_se = nan
                    else:
                        raise Exception("No match for gcov_int_regex")

                    line = infile.readline()

                    while re.match("^p1\s", line) is None:
                        line = infile.readline()

                    line = infile.readline()

                    rg, rg_se, rg_z, rg_p = [float(z) if z != 'NA' else nan for z in line.split()[2:6]]
                    d.append(
                        {
                            'trait.A' : m.group('trait_A'),
                            'trait.B' : m.group('trait_B'),
                            'snp.set' : wildcards.snp_set,
                            'h2.A.obs.ldsc' : h2_A,
                            'h2.A.obs.se.ldsc' : h2_A_se,
                            'h2.B.obs.ldsc' : h2_B,
                            'h2.B.obs.se.ldsc' : h2_B_se,
                            'gcov.obs.ldsc' : gcov,
                            'gcov.obs.se.ldsc' : gcov_se,
                            'rg.ldsc' : rg,
                            'rg.se.ldsc' : rg_se,
                            'rg.z.ldsc' : rg_z,
                            'rg.p.ldsc' : rg_p
                        }
                    )

        pd.DataFrame(d).to_csv(output[0], sep = '\t', index = False)
