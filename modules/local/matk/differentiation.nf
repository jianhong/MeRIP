process MATK_DIFFERENTIATION {
    tag "$meta.id"

    input:
    tuple val(meta), path(ipbam), path(controlbam), path(bed), path(ipbam2), path(controlbam2), path(bed2)
    path  gtf
    path  matk_jar

    output:
    tuple val(meta), path("*.m6A_differentiation.bed"), emit: diff
    path  "versions.yml"                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def ip = ipbam instanceof Path ? ipbam.toString() : (ipbam as List).join(';')
    def control = controlbam instanceof Path ? controlbam.toString() : (controlbam as List).join(';')
    def ip2 = ipbam2 instanceof Path ? ipbam2.toString() : (ipbam2 as List).join(';')
    def control2 = controlbam2 instanceof Path ? controlbam2.toString() : (controlbam2 as List).join(';')
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    java -jar $matk_jar \\
        -diff \\
        -control_ip "${ip2}" \\
        -control_input "${control2}" \\
        -treated_ip "${ip}" \\
        -treated_input "${control}" \\
        -control_bed ${bed2} \\
        -treated_bed ${bed} \\
        -gtf {$gtf} \\
        -out ${prefix}.m6A_differentiation.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        MATK: matk_jar
    END_VERSIONS
    """
}
