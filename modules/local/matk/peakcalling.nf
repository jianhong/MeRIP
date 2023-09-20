process MATK_PEAKCALLING {
    tag "$meta.id"

    input:
    tuple val(meta), path(ipbam), path(controlbam)
    path  gtf
    path  matk_jar

    output:
    tuple val(meta), path("*.bed")                   , emit: peak
    path  "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def gtf_opt  = gtf ? "-gtf $gtf" : ''
    def ip = ipbam instanceof Path ? ipbam.toString() : (ipbam as List).join(';')
    def control = controlbam instanceof Path ? controlbam.toString() : (controlbam as List).join(';')
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    java -jar $matk_jar \\
        -peakCalling \\
        -ip "${ip}" \\
        -input "${control}" \\
        -out ${prefix}.bed \\
        ${gtf_opt} \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MATK: matk_jar
    END_VERSIONS
    """
}
