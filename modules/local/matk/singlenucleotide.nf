process MATK_SINGLENUCLEOTIDE {
    tag "$meta.id"

    input:
    tuple val(meta), path(ipbam), path(controlbam), path(bed)
    path  bit
    path  matk_jar

    output:
    tuple val(meta), path("*.m6A_sites.bed")         , emit: sites
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ip = ipbam instanceof Path ? ipbam.toString() : (ipbam as List).join(';')
    def control = controlbam instanceof Path ? controlbam.toString() : (controlbam as List).join(';')
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    java -jar $matk_jar \\
        -singleNucleotide \\
        -ip "${ip}" \\
        -input "${control}" \\
        -bed $bed \\
        -2bit ${bit} \\
        -out ${prefix}.m6A_sites.bed \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MATK: matk_jar
    END_VERSIONS
    """
}
