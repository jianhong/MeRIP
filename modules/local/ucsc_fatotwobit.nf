process UCSC_FATOTWOBIT {
    tag "${fa}"
    label 'process_low'

    conda "bioconda::ucsc-fatotwobit=447"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ucsc-fatotwobit:447--h954228d_0':
        'biocontainers/ucsc-fatotwobit:447--h954228d_0' }"

    input:
    path fa

    output:
    path "*.2bit"                      , emit: bit
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '447' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    faToTwoBit \\
        $fa \\
        ${fa.getSimpleName()}.2bit

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ucsc: $VERSION
    END_VERSIONS
    """
}
