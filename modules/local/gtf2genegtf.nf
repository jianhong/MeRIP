process GTF2GENEGTF {
    tag "$gtf"
    label 'process_low'

    conda "bioconda::gffutils=0.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gffutils:0.12--pyh7cba7a3_0' :
        'biocontainers/gffutils:0.12--pyh7cba7a3_0' }"

    input:
    path gtf

    output:
    path '*.withgene.gtf'       , emit: gtf
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python

    #######################################################################
    #######################################################################
    ## Created on April. 29, 2021 prepare the gtf with gene annotation for MATK
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    import gffutils
    db = gffutils.create_db('${gtf}', '${gtf.getSimpleName()}.db', keep_order=True, sort_attribute_values=True)
    with open('${gtf.getSimpleName()}.withgene.gtf', 'w') as fout:
     for gene in db.features_of_type('gene', order_by='start'):
         fout.write(str(gene))
         for i in db.children(gene, order_by='start'):
            fout.write(str(i) + '\n')
    """
}
