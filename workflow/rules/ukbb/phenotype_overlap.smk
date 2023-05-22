import pandas as pd
import re
from itertools import chain

def get_trait_token(trait_wildcard):
    if '20002_' in trait_wildcard:
        return trait_wildcard.replace('20002_', '')
    elif trait_wildcard == 'E4_DM1':
        return 'E10'
    else:
        return trait_wildcard

rule match_data_fields_to_trait_codes:
    input:
        manifest = "resources/ukbb_sum_stats/phenotype_manifest_201807.tsv",
        trait_metadata = "resources/ukbb_sum_stats/trait_metadata.tsv"
    output:
        protected("results/metadata/ukbb/meta_on_manifest.tsv")
    localrule: True
    run:
        manifest = pd.read_csv(input.manifest, sep = '\t', header = 0)
        manifest = manifest.query('Sex == \'both_sexes\'')
        manifest = manifest.assign(field = manifest['UK Biobank Data Showcase Link'].str.split("field\.cgi\?id=", expand = True)[[1]])

        trait_metadata = pd.read_csv(input.trait_metadata, sep = '\t', header = 0)

        merged = manifest.merge(trait_metadata, left_on = 'Phenotype Code', right_on = 'code')

        merged.to_csv(output[0], sep = '\t', index = False)

# TODO check contents of file to make sure we are grepping for the right sort of thing, check numbers

rule slice_ukb_tab_file:
    input:
        "resources/ukbb/ukbb_tab_file.tsv"
    output:
        "results/metadata/ukbb/trait_columns.tsv"
    localrule: True
    resources:
        runtime = 20
    shell: """
        cut -f6019-6023,6762-6795,10124,15054-15133 -d$'\t' {input} >{output}
    """

rule recode_ukb_tab_file:
    input:
        "results/metadata/ukbb/trait_columns.tsv"
    output:
        "results/metadata/ukbb/trait_columns_recoded.tsv"
    localrule: True
    threads: 8
    script: "../../scripts/recode_ukbb_fields.R"

rule count_field_overlapping_trait_instances:
    input:
        "results/metadata/ukbb/trait_columns_recoded.tsv"
    output:
        "results/metadata/ukbb/overlaps/{trait_A}_and_{trait_B}.tsv"
    params:
        trait_A_token = lambda wildcards: get_trait_token(wildcards.trait_A),
        trait_B_token = lambda wildcards: get_trait_token(wildcards.trait_B)
    localrule: True
    run:
        a_cases = int(shell("tail -n +2 {input} | grep -c {params.trait_A_token}", read = True))
        b_cases = int(shell("tail -n +2 {input} | grep -c {params.trait_B_token}", read = True))
        ab_cases = int(shell("tail -n +2 {input} | grep {params.trait_A_token} | grep -c {params.trait_B_token}", read = True))
        a_controls = int(shell("tail -n +2 {input} | grep -v -c {params.trait_A_token}", read = True))
        b_controls = int(shell("tail -n +2 {input} | grep -v -c {params.trait_B_token}", read = True))
        ab_controls = int(shell("tail -n +2 {input} | grep -v {params.trait_A_token} | grep -v -c {params.trait_B_token}", read = True))

        with open(output[0], 'w') as out:
            out.write(f"{wildcards.trait_A}\t{wildcards.trait_B}\t{a_cases}\t{b_cases}\t{ab_cases}\t{a_controls}\t{b_controls}\t{ab_controls}\n")

traits_to_keep = ["20002_1464","20002_1381","20002_1111","22126","20002_1462","K51","20002_1465","20002_1473","K80","20002_1452","20002_1154","D25","6148_2","6148_5","20002_1226","E4_DM1","I42","N80"]

overlapping_instance_pairs = list(chain(*[[f"{traits_to_keep[i]}_and_{traits_to_keep[j]}" for j in range(i+1,len(traits_to_keep))] for i in range(len(traits_to_keep))]))

rule compile_overlapping_instances_files:
    input:
        [f"results/metadata/ukbb/overlaps/{x}.tsv" for x in overlapping_instance_pairs]
    output:
        temp("results/metadata/ukbb/overlaps/compiled_pairs.tsv")
    localrule: True
    shell:"""
        echo -e "trait_A\ttrait_B\tA_cases\tB_cases\tAB_cases\tA_controls\tB_controls\tAB_controls" >>{output}
        cat {input} >>{output}
    """

rule annotate_overlap_results:
    input:
        overlap = "results/metadata/ukbb/overlaps/compiled_pairs.tsv",
        metadata = "resources/ukbb_sum_stats/trait_metadata.tsv"
    output:
        "results/metadata/ukbb/overlaps/annotated_compiled_pairs.tsv"
    params:
        traits_to_keep = traits_to_keep
    localrule: True
    script: "../../scripts/ukbb/annotate_overlap_results.R"
