//
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF              } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_SPIKEIN_FASTA    } from '../../modules/nf-core/gunzip/main'

include { UNTAR as UNTAR_STAR_INDEX         } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_BOWTIE2_INDEX       } from '../../modules/nf-core/untar/main'

include { CUSTOM_GETCHROMSIZES              } from '../../modules/nf-core/custom/getchromsizes/main'
include { GFFREAD                           } from '../../modules/nf-core/gffread/main'
include { STAR_GENOMEGENERATE               } from '../../modules/nf-core/star/genomegenerate/main'
include { HISAT2_EXTRACTSPLICESITES         } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD                      } from '../../modules/nf-core/hisat2/build/main'
include { BOWTIE2_BUILD                     } from '../../modules/nf-core/bowtie2/build/main'
include { BWA_INDEX                         } from '../../modules/nf-core/bwa/index/main'
include { KHMER_UNIQUEKMERS                 } from '../../modules/nf-core/khmer/uniquekmers/main'

include { GTF2BED                           } from '../../modules/local/gtf2bed'
include { CAT_ADDITIONAL_FASTA              } from '../../modules/local/cat_additional_fasta'
include { STAR_GENOME_CHECK                 } from '../../modules/local/star_genome_check'
include { STAR_GENOMEGENERATE_IGENOMES      } from '../../modules/local/star_genomegenerate_igenomes'

workflow PREPARE_GENOME {
    main:
    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = Channel.value(file(params.fasta))
    }


    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = Channel.value(file(params.gtf))
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [:], params.gff ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = Channel.value(file(params.gff))
        }
        ch_gtf      = GFFREAD ( ch_gff ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

    //
    // Uncompress spikein fasta file and concatenate with reference fasta and gtf files
    //
    if (params.spikein_fasta) {
        if (params.spikein_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_SPIKEIN_FASTA ( [ [:], params.spikein_fasta ] ).gunzip.map { it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
        } else {
            ch_add_fasta = Channel.value(file(params.spikein_fasta))
        }
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_add_fasta, params.gencode ? "gene_type" : "gene_biotype" )
        ch_fasta    = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf      = CAT_ADDITIONAL_FASTA.out.gtf
        ch_versions = ch_versions.mix(CAT_ADDITIONAL_FASTA.out.versions)
    }


    //
    // Uncompress gene BED annotation file or create from GTF if required
    //
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], params.gene_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = Channel.value(file(params.gene_bed))
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta.map { [ [:], it ] } )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai.map { it[1] }
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes.map { it[1] }
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)


    //
    // Uncompress STAR index or generate from scratch if required
    //
    // Check if an AWS iGenome has been provided to use the appropriate version of STAR
    //
    is_aws_igenome = false
    if (params.fasta && params.gtf) {
        if ((file(params.fasta).getName() - '.gz' == 'genome.fa') && (file(params.gtf).getName() - '.gz' == 'genes.gtf')) {
            is_aws_igenome = true
        }
    }
    ch_star_index = Channel.empty()
    if ('star' == params.aligner) {
        // check index version
        is_compatible_star_index = false
        if (params.star_index) {
            compatible_star_index = STAR_GENOME_CHECK(Channel.value(file(params.star_index))).compatable
            if ("$compatible_star_index" == 'yes') {
                is_compatible_star_index = true
            }
        }
        if (is_compatible_star_index) {
            if (params.star_index.endsWith('.tar.gz')) {
                ch_star_index = UNTAR_STAR_INDEX ( [ [:], params.star_index ] ).untar.map { it[1] }
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.value(file(params.star_index))
            }
        } else {
            if (is_aws_igenome) {
                ch_star_index = STAR_GENOMEGENERATE_IGENOMES ( ch_fasta, ch_gtf ).index
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE_IGENOMES.out.versions)
            } else {
                ch_star_index = STAR_GENOMEGENERATE ( ch_fasta.map{[ [:], it ]}, ch_gtf.map{[ [:], it ]} ).index.map{ it[1] }
                ch_versions   = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
            }
        }
    }


    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()
    if ('hisat2' == params.aligner) {
        if (!params.splicesites) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf.map { [ [:], it ] } ).txt.map { it[1] }
            ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        } else {
            ch_splicesites = Channel.value(file(params.splicesites))
        }
        if (params.hisat2_index) {
            if (params.hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [:], params.hisat2_index ] ).untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
            } else {
                ch_hisat2_index = Channel.value(file(params.hisat2_index))
            }
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_fasta.map { [ [:], it ] }, ch_gtf.map { [ [:], it ] }, ch_splicesites.map { [ [:], it ] } ).index.map { it[1] }
            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }


    //
    // Uncompress bowtie2 index or generate from scratch if required
    //
    ch_bowtie2_index = Channel.empty()
    if (params.bowtie2_index) {
        if (params.bowtie2_index.endsWith('.tar.gz')) {
            ch_bowtie2_index = UNTAR_BOWTIE2_INDEX ( [ [:], params.bowtie2_index ] ).untar.map { it[1] }
            ch_versions     = ch_versions.mix(UNTAR_BOWTIE2_INDEX.out.versions)
        } else {
            ch_bowtie2_index = Channel.value(file(params.bowtie2_index))
        }
    } else {
        if ('bowtie2' == params.aligner) {
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_fasta ).index.map{ it[1] }
            ch_versions     = ch_versions.mix(BOWTIE2_INDEX.out.versions)
        }
    }


    //
    // Uncompress bwa index or generate from scratch if required
    //
    ch_bwa_index = Channel.empty()
    if (params.bwa_index) {
        ch_bwa_index = Channel.value(file(params.bwa_index))
    } else {
        if ('bwa' == params.aligner) {
            ch_bwa_index = BWA_INDEX ( ch_fasta.map{ [[:], it] } ).index.map{ it[1] }
            ch_versions     = ch_versions.mix(BWA_INDEX.out.versions)
        }
    }


    //
    // MODULE: Calculute genome size with khmer
    //
    if (!params.macs_gsize) {
        KHMER_UNIQUEKMERS (
            ch_fasta,
            params.read_length
        )
        ch_macs_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it.text.trim() }
    } else {
        ch_macs_gsize = params.macs_gsize
    }

    emit:
    fasta            = ch_fasta                  // channel: path(genome.fasta)
    gtf              = ch_gtf                    // channel: path(genome.gtf)
    fai              = ch_fai                    // channel: path(genome.fai)
    gene_bed         = ch_gene_bed               // channel: path(gene.bed)
    chrom_sizes      = ch_chrom_sizes            // channel: path(genome.sizes)
    splicesites      = ch_splicesites            // channel: path(genome.splicesites.txt)
    star_index       = ch_star_index             // channel: path(star/index/)
    hisat2_index     = ch_hisat2_index           // channel: path(hisat2/index/)
    bowtie2_index    = ch_bowtie2_index          // channel: path(bowtie2/index/)
    bwa_index        = ch_bwa_index              // channel: path(bwa/index/)
    macs_gsize       = ch_macs_gsize             // channel: val(gsize)

    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
}
