## Copyright (C) 2021  Roel Janssen <roel@gnu.org>

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

#' Run AneuFinder for all cells of a donor.
#'
#' @param output_directory         The directory to write output to.
#' @param donor                    The name of the donor to find cells for.
#' @param samplesheet              A data frame of the samplesheet.
#' @param numCPU                   The number of CPUs to use.
#' @param copyNumberCallingBinSize The bin size to use for copy number calling.
#' @param stepsize                 The step size to use for copy number calling.
#' @param blacklist.file           The blacklist file to use.
#' @param sequenceability.file     The sequenceability factors file to use.
#' @param correction.method        The correction method to apply.
#' @param plotting                 Whether plots should be generated.
#' @param reference.genome         The BSgenome to use as reference genome.
#'
#' @importFrom AneuFinder Aneufinder
#'
#' @export

runAneufinderForDonor <- function (output_directory,
                                   donor,
                                   samplesheet,
                                   numCPU,
                                   copyNumberCallingBinSize,
                                   stepsize,
                                   blacklist.file=NULL,
                                   sequenceability.file=NULL,
                                   correction.method=c("GC"),
                                   plotting=FALSE,
                                   reference.genome,
                                   autosomes,
                                   allosomes)
{
    outputTempFolder       <- paste0 (output_directory, "/.tmp_", donor)
    aneufinderOutputFolder <- paste0 (output_directory, "/", donor)

    dir.create (outputTempFolder, showWarnings = FALSE, recursive = TRUE)
    dir.create (aneufinderOutputFolder, showWarnings = FALSE, recursive = TRUE)

    for (filepath in samplesheet$filename)
    {
        name        <- basename (filepath)
        destination <- paste0(outputTempFolder, "/", name)

        file.symlink (filepath, destination)
        file.symlink (paste0 (filepath, ".bai"), paste0 (destination, ".bai"))
    }

    Aneufinder (inputfolder            = outputTempFolder,
                outputfolder           = aneufinderOutputFolder,
                numCPU                 = numCPU,
                binsizes               = copyNumberCallingBinSize,
                stepsizes              = stepsize,
                correction.method      = correction.method,
                chromosomes            = c(allosomes, autosomes),
                remove.duplicate.reads = TRUE,
                reads.store            = FALSE,
                blacklist              = blacklist.file,
                strandseq              = FALSE,
                GC.BSgenome            = reference.genome,
                states                 = c("zero-inflation", paste0(0:10, "-somy")),
                method                 = c("edivisive"),
                min.mapq               = 10,
                sequenceability.file   = sequenceability.file,
                plotting               = plotting)

    unlink (outputTempFolder, recursive = TRUE)
}

#' Create blacklist
#'
#' @param outputDirectory    The directory to write output to.
#' @param samplesheet        A data frame of the samplesheet.
#' @param blacklistBinSize   The bin size to use for backlisting regions.
#' @param genome             The BSgenome to use as reference.
#' @param allosomes          The allosomal chromosomes to include.
#' @param autosomes          The autosomal chromosomes to include.
#' @param numCPU             The number of CPUs to use.
#'
#' @return The file path to the blacklist regions file.
#'
#' @importFrom AneuFinder exportGRanges
#' @importFrom parallel   mclapply
#'
#' @export

createBlacklistFromSamplesheet <- function (outputDirectory,
                                            samplesheet,
                                            blacklistBinSize,
                                            genome,
                                            allosomes,
                                            autosomes,
                                            numCPU=16)
{
    blacklist.file <- paste0(outputDirectory, "/blacklist_", blacklistBinSize, ".bed.gz")
    if (! file.exists (blacklist.file))
    {
        chromosomeFilter  <- c(autosomes, allosomes)
        coverages         <- mclapply (samplesheet$filename,
                                       function (filename) {
                                           return (coveragePerBin (filename,
                                                                   genome,
                                                                   chromosomeFilter,
                                                                   blacklistBinSize))
                                       }, mc.cores = numCPU)

        totalBinCounts    <- mergeBinCounts (coverages)
        outlierBins       <- determineOutlierRegions (NULL,
                                                      genome,
                                                      blacklistBinSize,
                                                      autosomes,
                                                      allosomes,
                                                      totalBinCounts,
                                                      0.05, 0.95,
                                                      0.05, 0.95)

        ## ".bed.gz" will be appended by the exportGRanges function.
        blacklist.file <- paste0(outputDirectory, "/blacklist_", blacklistBinSize)

        AneuFinder::exportGRanges (outlierBins,
                                   filename = blacklist.file,
                                   header = FALSE,
                                   chromosome.format = "NCBI")

        ## Let the file path match the actually created path.
        blacklist.file <- paste0(blacklist.file, ".bed.gz")
    }

    return(blacklist.file)
}

#' Create sequenceability factors
#'
#' @param outputDirectory          The directory to write output to.
#' @param samplesheet              A data frame of the samplesheet.
#' @param copyNumberCallingBinSize The bin size to use for copy number calling.
#' @param reference.genome         The BSgenome to use as reference.
#' @param numCPU                   The number of CPUs to use.
#'
#' @return The file path to the sequenceability factors file.
#'
#' @importFrom Rsamtools BamFile
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom AneuFinder fixedWidthBins
#'
#' @export

createSequenceabilityFactorsFromSamplesheet <- function (outputDirectory,
                                                         samplesheet,
                                                         copyNumberCallingBinSize,
                                                         stepsize,
                                                         reference.genome,
                                                         allosomes,
                                                         autosomes,
                                                         numCPU=16)
{
    sequenceability.file <- paste0(outputDirectory,
                                   "/sequenceability.factors.",
                                   copyNumberCallingBinSize,
                                   ".gc.RData")

    if (! file.exists (sequenceability.file))
    {
        bamfile       <- samplesheet$filename[1]
        chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(bamfile))
        bins          <- fixedWidthBins (chrom.lengths = chrom.lengths, chromosomes=c(allosomes, autosomes),
                                         binsizes      = copyNumberCallingBinSize,
                                         stepsizes     = copyNumberCallingBinSize)

        sequenceability.folder <- paste0 (outputDirectory, "/.tmp_sequenceability_factors")
        runAneufinderForDonor (sequenceability.folder,
                               "Aneufinder",
                               samplesheet,
                               numCPU,
                               copyNumberCallingBinSize,
                               stepsize,
                               NULL,
                               NULL,
                               c("GC"),
                               plotting=FALSE,
                               reference.genome,
                               autosomes,
                               allosomes)

        sequenceability.folder  <- paste0 (sequenceability.folder, "/Aneufinder/binned-GC")
        sequenceability.factors <- determineSequenceabilityFactors (sequenceability.folder, bins)
        save(sequenceability.factors, file=sequenceability.file)
    }

    return(sequenceability.file)
}

#' Run AneuFinder for all samples in the specified samplesheet.
#'
#' @param outputDirectory              The directory to write output to.
#' @param samplesheet                  A data frame of the samplesheet.
#' @param blacklistBinSize             The bin size to use for backlisting regions.
#' @param copyNumberCallingBinSize     The bin size to use for copy number calling.
#' @param stepsize                     The step size to use for copy number calling.
#' @param genome                       The BSgenome to use as reference.
#' @param allosomes                    The allosomal chromosomes to include.
#' @param autosomes                    The autosomal chromosomes to include.
#' @param applySequenceabilityFactors  Whether to apply sequenceability factors.
#' @param numCPU                       The number of CPUs to use.
#'
#' @export

runAneufinderForSamplesheet <- function (outputDirectory,
                                         samplesheet,
                                         blacklistBinSize,
                                         copyNumberCallingBinSize,
                                         stepsize,
                                         genome,
                                         autosomes,
                                         allosomes,
                                         applySequenceabilityFactors = FALSE,
                                         numCPU = 16,
                                         plotting = FALSE)
{
		chromosomeFilter  <- c(autosomes, allosomes)
    sf_samplesheet    <- samplesheet[which (samplesheet$include_in_sf == 1),]

    ## -----------------------------------------------------------------------
    ## CREATE BLACKLIST
    ## -----------------------------------------------------------------------
    blacklist.file <- createBlacklistFromSamplesheet (outputDirectory,
                                                      samplesheet,
                                                      blacklistBinSize,
                                                      genome,
                                                      allosomes,
                                                      autosomes,
                                                      numCPU)

		## -----------------------------------------------------------------------
    ## CREATE SEQUENCEABILITY FACTORS
    ## -----------------------------------------------------------------------

		if (! applySequenceabilityFactors) {
			 sequenceability.file <- NULL
		} else {
        sequenceability.file <- createSequenceabilityFactorsFromSamplesheet (
            outputDirectory,
            sf_samplesheet,
            copyNumberCallingBinSize,
            stepsize,
            genome,
            allosomes,
            autosomes,
            numCPU)
    }

		## -----------------------------------------------------------------------
    ## RUN ANEUFINDER
    ## -----------------------------------------------------------------------

    donors <- unique(samplesheet$donor)
    mclapply (donors, function (donor) {
        donor_samplesheet <- samplesheet[which (samplesheet$donor == donor),]
        correction_method <- "GCSC"
        if (!applySequenceabilityFactors) { correction_method <- "GC" }
        runAneufinderForDonor (outputDirectory,
                               donor,
                               donor_samplesheet,
                               numCPU,
                               copyNumberCallingBinSize,
                               stepsize,
                               blacklist.file,
                               sequenceability.file,
                               correction.method=c(correction_method),
                               plotting=plotting,
                               reference.genome=genome,
                               autosomes,
                               allosomes)
    }, mc.cores=numCPU)
}

#' Gather quality metrics from AneuFinder output.
#'
#' @param base_directory  The directory in which AneuFinder output was written.
#' @param samplesheet     The samplesheet used to run the pipeline.
#'
#' @return A data frame containing various quality metrics for all cells
#'         associated with the specified donor.
#'
#' @export

gatherQualityInfoForSamplesheet <- function (base_directory, samplesheet)
{
    donor_samples         <- samplesheet[["sample_name"]]
    number_of_samples     <- length(donor_samples)

    ## Pre-allocate numeric arrays
    name          <- character(number_of_samples)
    num.segments  <- numeric(number_of_samples)
    bhattacharyya <- numeric(number_of_samples)
    spikiness     <- numeric(number_of_samples)
    entropy       <- numeric(number_of_samples)
    read.count    <- numeric(number_of_samples)

    ## Look up scores in the data files
    for (sample_index in 1:number_of_samples)
    {
        sample_name  <- donor_samples[sample_index]
        donor_name   <- samplesheet[which(samplesheet[["sample_name"]] == sample_name),][["donor"]]
        file_name    <- Sys.glob(paste0(base_directory, "/",
                                        donor_name ,"/MODELS",
                                        "/method-edivisive/", sample_name,
                                        "_dedup.bam_*.RData"))

        if (identical(file_name, character(0))) {
            name[sample_index]          <- sample_name
            num.segments[sample_index]  <- NA
            bhattacharyya[sample_index] <- NA
            entropy[sample_index]       <- NA
            spikiness[sample_index]     <- NA
            read.count[sample_index]    <- NA
        }
        else {
            sample       <- get(load(file_name))

            name[sample_index]          <- sample_name
            num.segments[sample_index]  <- sample$qualityInfo$num.segments
            bhattacharyya[sample_index] <- sample$qualityInfo$bhattacharyya
            entropy[sample_index]       <- sample$qualityInfo$entropy
            spikiness[sample_index]     <- sample$qualityInfo$spikiness
            read.count[sample_index]    <- sample$qualityInfo$total.read.count
        }
    }

    output <- data.frame (name, num.segments, bhattacharyya, entropy, spikiness, read.count)

    ## Exclude bin scores that have no read support.
    output$entropy[which(!is.finite(output$entropy))]             <- NA
    output$num.segments[which(!is.finite(output$num.segments))]   <- NA
    output$spikiness[which(!is.finite(output$spikiness))]         <- NA
    output$bhattacharyya[which(!is.finite(output$bhattacharyya))] <- NA

    return (output)
}


#' Create a single matrix with all copy number events in a samplesheet
#'
#' @param base_directory  The directory in which AneuFinder output was written.
#' @param samplesheet     The samplesheet used to run the pipeline.
#' @param maximum.size    The maximum event size to look for.
#' @param numCPU          The number of tasks to run in parallel.
#'
#' @importFrom BiocGenerics width
#'
#' @return A matrix of the size number-of-bins * number-of-cells.
#'
#'

copyNumberMatrixForSamplesheet <- function (base_directory, samplesheet, numCPU=16)
{
    ## ------------------------------------------------------------------------
    ## Build matrices with copy number states.
    ## ------------------------------------------------------------------------

    donor_samples     <- samplesheet[["sample_name"]]
    number_of_samples <- length(donor_samples)

    ## Determine the number of bins.
    sample_name       <- donor_samples[1]
    donor_name        <- samplesheet[which(samplesheet[["sample_name"]] == sample_name),][["donor"]]
    file_name         <- Sys.glob(paste0(base_directory, "/",
                                         donor_name ,"/MODELS",
                                         "/method-edivisive/", sample_name,
                                         "_dedup.bam_*.RData"))
    sample            <- get(load(file_name))
    bins              <- granges(sample[["bincounts"]][[1]])
    bin_size          <- width(bins[1])
    number_of_bins    <- length(bins)
    rm(sample)

    ## Not ideal, but reduce the memory footprint before forking.
    gc(full=TRUE)

    ## Create matrix with copy-number states and relative state changes.
    results <- mclapply (1:number_of_samples, function(sample_index)
    {
        sample_name   <- donor_samples[sample_index]
        donor_name    <- samplesheet[which(samplesheet[["sample_name"]] == sample_name),][["donor"]]
        file_name     <- Sys.glob(paste0(base_directory, "/",
                                         donor_name ,"/MODELS",
                                         "/method-edivisive/", sample_name,
                                         "_dedup.bam_*.RData"))
        sample        <- get(load(file_name))

        ## Determine the sample's base copy number state.
        base_cn_state <- median(sample$segments$copy.number)

        row <- sapply (seq_along(bins), function (range_index) {
            range       <- subsetByOverlaps(sample[["segments"]], bins[range_index])
            events_df   <- range$copy.number
            relative_df <- base_cn_state / range$copy.number

            return (c(events_df, relative_df))
        })

        rm(sample)
        return(list(base_cn_state, row))
    }, mc.cores=numCPU)

    ##events.df         <- matrix(ncol=number_of_bins, nrow=number_of_samples)
    relative.df       <- matrix(ncol=number_of_bins, nrow=number_of_samples)
    ##base.cn.state     <- numeric(number_of_samples)

    ## Assign the values gathered using mclapply in the pre-allocated matrix
    for (index in 1:number_of_samples) {
        ##base.cn.state[index] <- results[[index]][[1]]
        relative.df[index,]  <- results[[index]][[2]][2,]
        ##events.df[index,]    <- results[[index]][[2]][1,]
    }

    rownames(relative.df) <- samplesheet[,"sample_name"]
    colnames(relative.df) <- paste0(seqnames(bins), ":", ranges(bins))

    return(relative.df)
}


#' Summarize events per donor.
#'
#' @param base_directory  The directory in which AneuFinder output was written.
#' @param samplesheet     The samplesheet used to run the pipeline.
#' @param donor           The name of the donor to process.
#' @param numCPU          The number of tasks to run in parallel. (default=16)
#'
#' @importFrom BiocGenerics width
#'
#' @return A GRanges object containing the regions of recurrent events.
#'
#' @export

eventsForDonor <- function (base_directory, samplesheet, donor, numCPU=16)
{
    donor_samplesheet <- samplesheet[which (samplesheet[["donor"]] == donor),]
    cn.matrix         <- copyNumberMatrixForSamplesheet (base_directory, donor_samplesheet, numCPU)
    number_of_bins    <- ncol(cn.matrix)
    bin_state         <- integer (length=number_of_bins)
    bin_logical       <- logical (length=number_of_bins)
    threshold         <- round(nrow(donor_samplesheet) / 30)

    ## Decide for each bin whether it's a "population-wide" event or not.
    for (column_index in 1:number_of_bins)
        for (row_index in 1:nrow(cn.matrix))
            bin_state[column_index] = bin_state[column_index] +
                (cn.matrix[row_index,column_index] != 1.0)

    ## Make the scoring logical.
    for (column_index in 1:number_of_bins)
        for (row_index in 1:nrow(cn.matrix))
            bin_logical[column_index] = (bin_state[column_index] > 6)

    return (bin_state)
}


#' Find recurring copy-number events.
#'
#' @param base_directory  The directory in which AneuFinder output was written.
#' @param samplesheet     The samplesheet used to run the pipeline.
#' @param maximum.size    The maximum event size to look for.
#' @param numCPU          The number of tasks to run in parallel.
#'
#' @importFrom BiocGenerics width
#'
#' @return A GRanges object containing the regions of recurrent events.
#'
#' @export

findRecurringEvents <- function (base_directory, samplesheet, maximum.size=10000, numCPU=16)
{
    relative.df <- copyNumberMatrixForSamplesheet (base_directory,
                                                   samplesheet,
                                                   numCPU=16)

    ## ------------------------------------------------------------------------
    ## Determine recurring events.
    ## ------------------------------------------------------------------------

    bins_needed_for_exclusion <- ceiling (10000000 / bin_size, 1)

    is.recurrent      <- numeric(number_of_bins)
    for (column_index in 1:ncol(events.df)) {
        bin             <- relative.df[,column_index]
        occurrences     <- as.data.frame(table(bin))

        ## XXX: The criteria should be: Events that occur in more than 3 donors.
        ## If more than 90% of the samples deviate from their base CN...
        if ((occurrences$Freq[which(occurrences$bin==1)] / number_of_samples) < 0.1) {
            is.recurrent[column_index] <- 1
        } else {
            is.recurrent[column_index] <- 0
        }
    }

    ## Separate the bins per chromosome.
    bins.df           <- as.data.frame(bins)
    chromosomes       <- unique(bins.df$seqnames)
    chromosome.bins   <- sapply (chromosomes,
                                 function (chromosome) {
                                     ## We can use SUM here because TRUE is 1,
                                     ## and FALSE is 0.
                                     sum(bins.df$seqnames == chromosome)
                                 })

    ## Check adjacency of the bins in the same chromosome.
    chromosome.bins[[chromosomes[1]]]

    return (relative.df)
}

#' Gather quality metrics from AneuFinder output.
#'
#' @param base_directory  The directory in which AneuFinder output was written.
#' @param samplesheet     The samplesheet used to run the pipeline.
#' @param donor           The name of the donor to extract cells for.
#'
#' @return A data frame containing various quality metrics for all cells
#'         associated with the specified donor.
#'
#' @export

gatherQualityInfoForDonor <- function (base_directory, samplesheet, donor)
{
    donor_samples         <- samplesheet[(samplesheet$donor == donor),][["sample_name"]]
    number_of_samples     <- length(donor_samples)

    ## Pre-allocate numeric arrays
    name          <- character(number_of_samples)
    num.segments  <- numeric(number_of_samples)
    bhattacharyya <- numeric(number_of_samples)
    spikiness     <- numeric(number_of_samples)
    entropy       <- numeric(number_of_samples)
    read.count    <- numeric(number_of_samples)

    ## Look up scores in the data files
    for (sample_index in 1:number_of_samples)
    {
        sample_name  <- donor_samples[sample_index]
        file_name    <- Sys.glob(paste0(base_directory, "/",
                                        donor, "/MODELS/method-edivisive/",
                                        sample_name, "_dedup.bam_*.RData"))

        if (identical(file_name, character(0))) {
            name[sample_index]          <- sample_name
            num.segments[sample_index]  <- NA
            bhattacharyya[sample_index] <- NA
            entropy[sample_index]       <- NA
            spikiness[sample_index]     <- NA
            read.count[sample_index]    <- NA
        }
        else {
            sample       <- get(load(file_name))

            name[sample_index]          <- sample_name
            num.segments[sample_index]  <- sample$qualityInfo$num.segments
            bhattacharyya[sample_index] <- sample$qualityInfo$bhattacharyya
            entropy[sample_index]       <- sample$qualityInfo$entropy
            spikiness[sample_index]     <- sample$qualityInfo$spikiness
            read.count[sample_index]    <- sample$qualityInfo$total.read.count
        }
    }

    output <- data.frame (name, num.segments, bhattacharyya, entropy, spikiness, read.count)

    ## Exclude bin scores that have no read support.
    output$entropy[which(!is.finite(output$entropy))]             <- NA
    output$num.segments[which(!is.finite(output$num.segments))]   <- NA
    output$spikiness[which(!is.finite(output$spikiness))]         <- NA
    output$bhattacharyya[which(!is.finite(output$bhattacharyya))] <- NA

    return (output)
}

#' Exclude cells from the analysis after performing quality control.
#'
#' @param base_directory           The directory in which AneuFinder output was written.
#' @param samplesheet              The samplesheet used to run the pipeline.
#' @param donor                    The name of the donor to extract cells for.
#' @param bhattacharyya_threshold  Threshold for the Bhattacharyya.
#' @param spikiness_threshold      Threshold for the spikiness.
#' @param plotOverlap              Whether to make a Venn diagram to show the overlap
#'                                 between filter criteria.
#' @param minimum_number_of_reads  Exclude cells with less reads than this parameter.
#' @param threshold_bhattacharyya  Exclude the worst scoring N percentage of cells.
#' @param threshold_spikiness      Exclude the worst scoring N percentage of cells.
#'
#' @importFrom ggplot2     ggsave
#' @importFrom VennDiagram venn.diagram
#' @importFrom utils       head tail
#' @importFrom scales      alpha
#'
#' @export

excludedCellsForRun <- function (base_directory,
                                 samplesheet,
                                 donor,
                                 bhattacharyya_threshold=NULL,
                                 spikiness_threshold=NULL,
                                 plotOverlap=FALSE,
                                 minimum_number_of_reads=200000,
                                 threshold_bhattacharyya=10,
                                 threshold_spikiness=10)
{
    run_samplesheet             <- samplesheet[which (samplesheet$donor == donor),]
    scores.df                   <- gatherQualityInfoForDonor (base_directory, samplesheet, donor)

    ncells                      <- nrow(scores.df)
    all_cells                   <- scores.df[["name"]]
    nreads_after_filter         <- scores.df[which(scores.df$read.count > minimum_numer_of_reads),"name"]

    if (is.null(bhattacharyya_threshold)) {
        bhattacharyya_threshold <- tail(head(sort(scores.df[["bhattacharyya"]]), round(ncells / threshold_bhattacharyya)), 1)
    }
    if (is.null(spikiness_threshold)) {
        spikiness_threshold     <- head(tail(sort(scores.df[["spikiness"]]), round(ncells / threshold_spikiness)), 1)
    }

    bhattacharyya_after_filter  <- scores.df[which(scores.df$bhattacharyya > bhattacharyya_threshold),"name"]
    spikiness_after_filter      <- scores.df[which(scores.df$spikiness < spikiness_threshold),"name"]

    included_cells              <- Reduce(intersect, list(nreads_after_filter,
                                                          bhattacharyya_after_filter,
                                                          spikiness_after_filter))
    excluded_cells              <- setdiff(all_cells, included_cells)

    if (plotOverlap) {
        temp   <- venn.diagram(
            x               = list (nreads_after_filter,
                                    bhattacharyya_after_filter,
                                    spikiness_after_filter),
            category.names  = c ("Reads", "Bhattacharyya", "Spikiness"),
            filename        = NULL,
            width           = 5,
            height          = 5,
            lwd             = 3,
            lty             = 'solid',
            col             = c("#56B4E9", "#E69F00", "#009E73"),
            fill            = c(alpha("#56B4E9",0.4),
                                alpha('#E69F00',0.4),
                                alpha('#009E73',0.4)),
            cex             = .9,
            fontface        = "bold",
            fontfamily      = "sans",
            cat.cex         = .9,
            cat.default.pos = "text",
            cat.pos         = c(0, 0, 0),
            cat.fontfamily  = "sans");

        ggsave(paste0(donor, "_filter_overlap.svg"), temp, width=5, height=5,units="cm", dpi=300)
    }

    return(excluded_cells)
}

#' Exclude cells from the analysis after performing quality control.
#'
#' @param base_directory           The directory in which AneuFinder output was written.
#' @param samplesheet              The samplesheet used to run the pipeline.
#' @param bhattacharyya_threshold  Threshold for the Bhattacharyya.
#' @param spikiness_threshold      Threshold for the spikiness.
#' @param plotOverlap              Whether to make a Venn diagram to show the overlap
#'                                 between filter criteria.
#' @param minimum_number_of_reads  Exclude cells with less reads than this parameter.
#' @param threshold_bhattacharyya  Exclude the worst scoring N percentage of cells.
#' @param threshold_spikiness      Exclude the worst scoring N percentage of cells.
#'
#' @importFrom ggplot2     ggsave
#' @importFrom VennDiagram venn.diagram
#' @importFrom utils       head tail
#' @importFrom scales      alpha
#'
#' @export

excludedCells <- function (base_directory,
                           samplesheet,
                           bhattacharyya_threshold=NULL,
                           spikiness_threshold=NULL,
                           plotOverlap=FALSE,
                           minimum_number_of_reads=200000,
                           threshold_bhattacharyya=10,
                           threshold_spikiness=10)
{
    scores.df                   <- gatherQualityInfoForSamplesheet (base_directory, samplesheet)
    ncells                      <- nrow(scores.df)
    all_cells                   <- scores.df[["name"]]
    nreads_after_filter         <- scores.df[which(scores.df$read.count > minimum_number_of_reads),"name"]

    if (is.null(bhattacharyya_threshold)) {
        bhattacharyya_threshold <- tail(head(sort(scores.df[["bhattacharyya"]]), round(ncells / threshold_bhattacharyya)), 1)
    }
    if (is.null(spikiness_threshold)) {
        spikiness_threshold     <- head(tail(sort(scores.df[["spikiness"]]), round(ncells / threshold_spikiness)), 1)
    }

    bhattacharyya_after_filter  <- scores.df[which(scores.df$bhattacharyya > bhattacharyya_threshold),"name"]
    spikiness_after_filter      <- scores.df[which(scores.df$spikiness < spikiness_threshold),"name"]

    included_cells              <- Reduce(intersect, list(nreads_after_filter,
                                                          bhattacharyya_after_filter,
                                                          spikiness_after_filter))
    excluded_cells              <- setdiff(all_cells, included_cells)

    if (plotOverlap) {
        temp   <- venn.diagram(
            x               = list (nreads_after_filter,
                                    bhattacharyya_after_filter,
                                    spikiness_after_filter),
            category.names  = c ("Reads", "Bhattacharyya", "Spikiness"),
            filename        = NULL,
            width           = 5,
            height          = 5,
            lwd             = 3,
            lty             = 'solid',
            col             = c("#56B4E9", "#E69F00", "#009E73"),
            fill            = c(alpha("#56B4E9",0.4),
                                alpha('#E69F00',0.4),
                                alpha('#009E73',0.4)),
            cex             = .9,
            fontface        = "bold",
            fontfamily      = "sans",
            cat.cex         = .9,
            cat.default.pos = "text",
            cat.pos         = c(0, 0, 0),
            cat.fontfamily  = "sans");

        ggsave("filter_overlap.svg", temp, width=5, height=5,units="cm", dpi=300)
    }

    return(excluded_cells)
}

#' Remove AneuFinder output of samples.
#'
#' @param base_directory  The directory in which AneuFinder output was written.
#' @param samples         A vector of sample names to exclude.
#'
#' @return The samplesheet without the specified samples.
#'
#' @export

removeExcludedCellsFromOutput <- function (base_directory, samples)
{
    number_of_samples <- length(samples)
    for (sample_index in 1:number_of_samples)
    {
        sample_name   <- samples[sample_index]
        file_name     <- Sys.glob(paste0(base_directory, "/*/MODELS/method-edivisive/",
                                         sample_name, "_dedup.bam_*.RData"))
        if (! identical(file_name, character(0))) {
            file.remove(file_name, showWarnings=FALSE)
        }
    }

    return(TRUE)
}

#' Remove samples from samplesheets.
#'
#' @param samplesheet     The samplesheet used to run the pipeline.
#' @param samples         A vector of sample names to exclude.
#'
#' @return The samplesheet without the specified samples.
#'
#' @export

removeCellsFromSamplesheet <- function (samplesheet, samples)
{
    return(samplesheet[which(! samplesheet$sample_name %in% samples),])
}
