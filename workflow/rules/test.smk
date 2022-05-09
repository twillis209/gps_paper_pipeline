# Try groups next
rule test_one:
    output:
        "test_one_{number}.txt"
    group: "test"
    resources:
        time = 1
    shell:
        "touch {output}"

rule test_two:
    input:
        "test_one_{number}.txt"
    output:
        "test_two_{number}.txt"
    group: "test"
    resources:
        time = 1
    shell:
        "touch {output}"
