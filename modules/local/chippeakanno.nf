process BIOC_CHIPPEAKANNO {
    tag "$meta.id"
    label 'process_medium'
    label 'error_ignore'

    conda "bioconda::bioconductor-chippeakanno=3.34.1"
    container "${ workflow.containerEngine == 'singularity' &&
                    !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bioconductor-chippeakanno:3.34.1--r43hdfd78af_0' :
        'biocontainers/bioconductor-chippeakanno:3.34.1--r43hdfd78af_0' }"

    input:
    tuple val(meta), path(peak)
    path gtf

    output:
    tuple val(meta), path("${prefix}/*")          , emit: anno
    tuple val(meta), path("${prefix}/**.anno.csv"), emit: csv
    path "${prefix}/*.png", optional:true         , emit: png
    path "versions.yml"                           , emit: versions

    script:
    prefix   = task.ext.prefix ?: "$meta.id"
    """
    #!/usr/bin/env Rscript

    #######################################################################
    #######################################################################
    ## Created on April. 29, 2021 call ChIPpeakAnno
    ## Copyright (c) 2021 Jianhong Ou (jianhong.ou@gmail.com)
    ## This source code is licensed under the MIT license
    #######################################################################
    #######################################################################
    pkgs <- c("ChIPpeakAnno", "rtracklayer", "GenomicFeatures", "ggplot2")
    versions <- c("${task.process}:")
    for(pkg in pkgs){
        # load library
        library(pkg, character.only=TRUE)
        # parepare for versions.yml
        versions <- c(versions,
            paste0("    ", pkg, ": ", as.character(packageVersion(pkg))))
    }
    writeLines(versions, "versions.yml") # write versions.yml

    gtf <- "${gtf}"
    pf <- file.path("${prefix}")
    extensionPattern <- ".(csv|bed|peaks|txt|narrowPeak|broadPeak)\$"
    detbl <- dir(".", extensionPattern,
                recursive = TRUE, full.names = TRUE)
    detbl <- detbl[!grepl("anno.csv", detbl)] ## in case of re-run
    txdb <- makeTxDbFromGFF(gtf) ## create annotation data from gtf file
    gtf <- import(gtf)
    id2symbol <- function(gtf){ ## convert entriz id to gene symbol
        if(is.null(gtf\$gene_name)) return(NULL)
        x <- data.frame(id=gtf\$gene_id, symbol=gtf\$gene_name)
        x <- unique(x)
        x <- x[!duplicated(x\$id), ]
        x <- x[!is.na(x\$id), , drop=FALSE]
        if(nrow(x)==0) return(NULL)
        y <- x\$symbol
        names(y) <- x\$id
        y
    }
    id2symbol <- id2symbol(gtf)
    anno <- toGRanges(txdb)
    promoters <- promoters(anno, upstream=2000, downstream=500)
    resList <- list() # save annotation results to a list

    dir.create(pf, showWarnings = FALSE, recursive = TRUE)
    for(det in detbl){
        extension <- sub("^.*\\\\.(.*?)\$", "\\\\1", det)
        switch(extension,
            "csv"={
                DB <- read.csv(det)
                rownames(DB) <- paste0("p", seq.int(nrow(DB)))
                DB.gr <- with(DB, GRanges(chr, IRanges(start, end, name=rownames(DB))))
            },
            "peaks"={
                DB <- read.table(det, header=TRUE)
                rownames(DB) <- paste0("p", seq.int(nrow(DB)))
                DB.gr <- with(DB, GRanges(chr, IRanges(start, end, name=rownames(DB))))
            },
            "txt"={
                header <- tryCatch(
                    read.table(det, header=FALSE, nrow=1),
                    error = function(.e){
                        message(.e)
                    },
                    finally = NULL
                )
                if(length(header)>0){
                    hasHeader <- all(c("chr", "start", "end") %in%
                                        header[1, , drop=TRUE])
                    DB <- read.table(det, header = hasHeader,
                                        stringsAsFactors = FALSE)
                    if(!hasHeader){
                        colnames(DB)[1:3] <- c("chr", "start", "end")
                    }
                }
                rownames(DB) <- paste0("p", seq.int(nrow(DB)))
                DB.gr <- with(DB, GRanges(chr, IRanges(start, end, name=rownames(DB))))
            },
            {
                DB.gr <- import(det)
                names(DB.gr) <- paste0("p", seq_along(DB.gr))
                DB <- as.data.frame(DB.gr)
            })
        if(nrow(DB)<1) next

        groupName <- gsub(extensionPattern, "", basename(det))
        # Annotation
        DB.anno <- annotatePeakInBatch(DB.gr, AnnotationData = anno,
                                        output = "both",
                                        PeakLocForDistance = "middle",
                                        FeatureLocForDistance = "TSS",
                                        ignore.strand = TRUE)
        if(length(id2symbol)>0) DB.anno\$symbol[!is.na(DB.anno\$feature)] <- id2symbol[DB.anno\$feature[!is.na(DB.anno\$feature)]]
        resList[[groupName]] <- DB.anno
        ## unique annotation plots
        ## only take the nearest annotation
        DB.anno.srt <- DB.anno[order(abs(DB.anno\$distancetoFeature))]
        DB.anno.srt <- unique(DB.anno.srt)

        pff <- file.path(pf, sub(extensionPattern, "", det))
        dir.create(dirname(pff), recursive = TRUE, showWarnings = FALSE)
        DB.anno <- mcols(DB.anno)
        DB <- cbind(DB[DB.anno\$peak, ], DB.anno)
        pff <- file.path(pf, sub(extensionPattern, ".anno.csv", det))
        dir.create(dirname(pff), recursive = TRUE, showWarnings = FALSE)
        write.csv(DB, pff, row.names = FALSE)
    }


    if(packageVersion("ChIPpeakAnno")>="3.23.12"){
        if(length(resList)>0){
            if(is.list(resList)){
                resList <- GRangesList(resList[lengths(resList)>0])
            }
            if(length(resList)>0){
                out <- genomicElementDistribution(resList,
                                                TxDb = txdb,
                                                promoterRegion=c(upstream=2000, downstream=500),
                                                geneDownstream=c(upstream=0, downstream=2000),
                                                promoterLevel=list(
                                                # from 5' -> 3', fixed precedence 3' -> 5'
                                                    breaks = c(-2000, -1000, -500, 0, 500),
                                                    labels = c("upstream 1-2Kb", "upstream 0.5-1Kb",
                                                            "upstream <500b", "TSS - 500b"),
                                                    colors = c("#FFE5CC", "#FFCA99",
                                                            "#FFAD65", "#FF8E32")),
                                                plot = FALSE)
                saveRDS(out\$peaks, file.path(pf, "genomicElementDistribuiton.RDS"))
                ggsave(file.path(pf, "genomicElementDistribuiton.pdf"), plot=out\$plot, width=9, height=9)
                ggsave(file.path(pf, "genomicElementDistribuiton.png"), plot=out\$plot)
                out <- metagenePlot(resList, txdb)
                ggsave(file.path(pf, "metagenePlotToTSS.pdf"), plot=out, width=9, height=9)
                ggsave(file.path(pf, "metagenePlotToTSS.png"), plot=out)
            }
        }
    }
    """
}
