## Copyright (C) 2021  Roel Janssen <roel@gnu.org>

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

#' Return a data frame containing the chromosome name and its length in bp.
#'
#' @param genome             A BSgenome object.
#' @param chromosomeFilter   A vector containing the chromosome names to keep.
#'
#' @return A data frame containing the chromosome name and its length in bp.
#'
#' @importFrom GenomeInfoDb  seqlengths
#'
#' @export

chromosomeLengths <- function (genome, chromosomeFilter)
{
    chr_lengths        <- seqlengths(genome)
    names(chr_lengths) <- gsub(pattern = "chr", replacement = "", x = names(chr_lengths))
    chr_lengths        <- chr_lengths[which(names(chr_lengths) %in% chromosomeFilter)]
    chr_lengths        <- data.frame(seqnames= names(chr_lengths), chr_size = as.numeric(chr_lengths))
    return (chr_lengths)
}

#' Merge overlapping CNVs
#'
#' This function can be used to merge regions (CNVs) based on a percentage of
#' overlap between the regions (x and y, Genomicranges object).  This
#' percentage is adjusted for the size of the region. For example, it
#' increases for larger regions.  The function can also be used to filter
#' away overlapping regions, by setting ' mode = "Filter".
#'
#' This function was written by Sjors Middelkamp.
#'
#' @param  x                A GRanges object.
#' @param  y                A GRanges object.
#' @param  minimal_overlap  The minimum percentage of overlap between regions
#'                          to be considered "overlapping".
#' @param  resolution       The bin size used in AneuFinder.
#' @param  mode             Either "Merge" or "Filter".
#'                          With "Merge", overlapping regions are merged into
#'                          one.  With "Filter", ranges in 'x' that overlap
#'                          with 'y' are removed.
#'
#' @return A GRanges object with the overlapping regions between 'x' and 'y'.
#'
#' @importFrom GenomicRanges granges findOverlaps pintersect
#' @importFrom IRanges width
#' @importFrom S4Vectors queryHits subjectHits
#'
#' @export

mergeOverlappingCNVs <- function(x, y,
                                 minimal_overlap = 0.3,
                                 resolution = 10e6,
                                 mode = "Merge")
{
    x       <- granges(x)
    y       <- granges(y)
    hits    <- findOverlaps(x, y)
    xhits   <- x[queryHits(hits)]
    yhits   <- y[subjectHits(hits)]

    # Determine how much overlap there is between the two ranges:
    frac    <- width(pintersect(xhits, yhits)) / pmax(width(xhits), width(yhits))

    # Determine the minimal overlap that is necessary to merge the two ranges.
    # (This depends on the size of the largest fragment)
    minfrac <- minimal_overlap + floor(pmax(width(xhits), width(yhits)) / resolution) * 0.1

    # Why do we have a maximum_overlap?
    maximum_overlap <- 0.9
    minfrac[which(minfrac > maximum_overlap)] <- maximum_overlap
    merge   <- frac >= minfrac

    if (mode == "Merge")
    {
        # The CNVs that overlap with more than the minfrac will be merged together.
        merged <- c(reduce(c(xhits[merge], yhits[merge])),
                    xhits[!merge], yhits[!merge],
                    x[-queryHits(hits)], y[-subjectHits(hits)])

        # The merging generates duplicates, remove the duplicates:
        merged <- merged[!duplicated(merged)]
    }
    else if (mode == "Filter")
    {
        if (length(hits) > 0)
        {
            # Remove the entries in x that overlap with y
            merged <-  c(x[-queryHits(hits)], xhits[!merge])
        }
        else
        {
            merged <- x
        }
  }

  return (merged)
}

#' Create populations
#'
#' This function uses ‘mergeOverlappingCNVs’ to merge a list of GRanges.
#'
#' @param  population_set   A list of GRanges.
#' @return A GRanges object with the merge of the 'population_set'.
#'
#' @export

createPopulation <- function (population_set)
{
    population <- NULL
    for (cell in population_set)
    {
        population <-  if (is.null (population)) { cell$segments }
                       else { mergeOverlappingCNVs (population, cell$segments) }
    }

    return (population)
}

#' Obtain the median copy number state per cell
#'
#' @param cells.list  A list of GRanges objects to determine the copy number
#'                    state of.
#'
#' @return A data frame with two columns: the ID and the median copy number
#'         state of the cell.
#'
#' @importFrom stats median
#' @export

getMedianCopyNumberStates <- function (cells.list)
{
    # Pre-allocate the vectors.
    number_of_items <- length(cells.list)
    cn.state.median <- numeric(number_of_items)
    ID              <- character(number_of_items)

    for (index in 1:number_of_items)
    {
        cell                   <- cells.list[[index]]
        segments               <- cell$segments
        cn.state.median[index] <- median(segments$copy.number)
        ID[index]              <- cell$ID
    }

    # Convert the vectors to a single data frame.
    return (data.frame (ID, cn.state.median, stringsAsFactors=FALSE))
}

#' Coverage per bin
#'
#' @param compositeBam     The BAM file to determine unreliable regions from.
#' @param genome           A BSgenome object matching the reference genome
#'                         used in the 'compositeBam'.
#' @param chromosomeFilter A vector containing the chromosome names to keep.
#' @param binSize          The size of a single bin.
#' @param minimumMappingQuality  The minimum mapping quality to include.
#'                               This should be a value between 0 and 60,
#'                               and defaults to 0 (include all reads).
#'
#' @return A GRanges object containing binned regions and their coverage.
#'
#' @importFrom GenomicRanges tileGenome GRanges
#' @importFrom IRanges IRanges
#' @importFrom stats setNames
#'
#' @export

coveragePerBin <- function (compositeBam, genome, chromosomeFilter, binSize, minimumMappingQuality = 0)
{
    if (! file.exists (compositeBam))
    {
        cat (paste0 ("File '", compositeBam, "' does not exist.\n"))
        return (FALSE)
    }

    ## Generate a BAM index so that idxstats are generated efficiently.
    index_filename <- paste0 (compositeBam, ".bai")
    if (! file.exists (index_filename))
    {
        createBamIndex (compositeBam);
    }

    ## Determine chromosome lengths based on the BSgenome.
    chromosome_lengths <- chromosomeLengths (genome, chromosomeFilter)
    chromosomes        <- GRanges (seqnames = chromosome_lengths[,1],
                                   IRanges (start = 1,
                                            end = chromosome_lengths[,2]))

    ## Split the genome in bins.
    bins <- tileGenome (seqlengths             = setNames(chromosome_lengths$chr_size,
                                                          chromosome_lengths$seqnames),
                        tilewidth              = as.numeric (binSize),
                        cut.last.tile.in.chrom = TRUE)

    df_bins <- data.frame (bins, stringsAsFactors = F)[,1:3]
    number_of_bins <- length(bins)
    regions_vector <- character(number_of_bins)

    for (index in 1:number_of_bins)
    {
        regions_vector[index] <- paste0(df_bins$seqnames[index], ":",
                                        df_bins$start[index], "-",
                                        df_bins$end[index])
    }

    ## Count the number of reads per bin.
    read.count       <- readsInRegions (compositeBam, regions_vector, minimumMappingQuality)
    read.count.total <- data.frame (df_bins, read.count, stringsAsFactors=FALSE)

    return (read.count.total)
}


#' Determine list of unreliable regions
#'
#' This function divides the genome in bins and returns a list of the
#' 2% bins with the least number of reads, and of the 3% bins with the
#' most number of reads for autosomal chromosomes, and the least 3% and
#' most 5% for allosomal chromosomes.
#'
#' @param compositeBam  The BAM file to determine unreliable regions from.
#' @param genome        A BSgenome object matching the reference genome used
#'                      in the 'compositeBam'.
#' @param autosomes     A vector containing the autosomal chromosome names.
#' @param allosomes     A vector containing the sex chromosome names.
#' @param binSize       The size of a single bin.
#' @param coverage      A data frame like the output of ‘coveragePerBin’ or
#'                      ‘mergeBinCounts’.
#' @param autosomalCutoffLow (default=0.02)
#' @param autosomalCutoffHigh (default=0.97)
#' @param allosomalCutoffLow (default=0.03)
#' @param allosomalCutoffHigh (default=0.95)
#'
#' @return A GRanges object containing the regions to exclude from further
#'         analysis.
#'
#' @importFrom GenomicRanges makeGRangesFromDataFrame reduce
#' @importFrom stats quantile
#'
#' @export

determineOutlierRegions <- function (compositeBam,
                                     genome,
                                     binSize,
                                     autosomes,
                                     allosomes,
                                     coverage = NULL,
                                     autosomalCutoffLow=0.02,
                                     autosomalCutoffHigh=0.97,
                                     allosomalCutoffLow=0.03,
                                     allosomalCutoffHigh=0.95)
{
    chromosomeFilter <- c(autosomes, allosomes)

    if (is.null(coverage))
        coverage <- coveragePerBin (compositeBam, genome, binSize, chromosomeFilter)

    ## Determine the outlier bins.
    coverage_autosomal <- coverage[coverage[,1] %in% autosomes,]
    low_cutoff  <- quantile (coverage_autosomal[,4], autosomalCutoffLow)
    high_cutoff <- quantile (coverage_autosomal[,4], autosomalCutoffHigh)

    autosomal_excluded_bins <- coverage_autosomal[coverage_autosomal[,4] <= low_cutoff |
                                                  coverage_autosomal[,4] >= high_cutoff,]

    coverage_allosomal <- coverage[coverage[,1] %in% allosomes,]
    low_cutoff  <- quantile (coverage_allosomal[,4], allosomalCutoffLow)
    high_cutoff <- quantile (coverage_allosomal[,4], allosomalCutoffHigh)

    allosomal_excluded_bins <- coverage_allosomal[coverage_allosomal[,4] <= low_cutoff |
                                                  coverage_allosomal[,4] >= high_cutoff,]

    combined_excluded_bins <- rbind (autosomal_excluded_bins, allosomal_excluded_bins)

    names (combined_excluded_bins) <- c("seqnames", "start", "end", "number.of.reads")
    output <- makeGRangesFromDataFrame (combined_excluded_bins, keep.extra.columns=TRUE)

    ## TODO: Figure out whether min.gapwidth is equivalent to the commented
    ## out piece below.
    ## TODO: The call to reduce removes the 'number.of.reads' metadata column.
    output <- reduce (output, min.gapwidth = (binSize / 2))

    return (output)

    ## ## To merge bins together, look half a binSize to the left.
    ## shifted_start <- combined_excluded_bins[,2] - (binSize / 2)

    ## ## For the first bin, we would end up with a negative start position.
    ## ## Reset it to 1.
    ## shifted_start[which (shifted_start < 0)] <- 1

    ## shifted_end <- combined_excluded_bins[,3] - (binSize / 2)

    ## # Merge bins together when they are less one bin apart.
    ## combined <- GRanges (seqnames = combined_excluded_bins[,1],
    ##                      IRanges(start = shifted_start,
    ##                              end   = shifted_end))
    ## output   <- reduce (combined)

    ## start_output <- start (output) + (binSize / 2)
    ## start (output) <- start_output
    ## start (output) <- start (output) + (binSize / 2)
    ## end   (output) <- end   (output) + (binSize / 2)

    ## return (output)
}


#' Merge the results of multiple outputs of ‘coveragePerBin’, summing
#' the total number of reads per bin.
#'
#' @param cells.list  A list of data frames returned by 'coveragePerBin'.
#'
#' @return A data frame containing the seqname, start, end, and the sum of
#'         reads in a bin.
#'
#' @export

mergeBinCounts <- function (cells.list)
{
    numberOfCells <- length(cells.list)
    numberOfBins  <- nrow(cells.list[[1]])

    totalPerBin   <- numeric(numberOfBins)
    for (binIndex in 1:numberOfBins)
    {
        binCounts <- numeric(numberOfCells)
        for (cellIndex in 1:numberOfCells)
        {
            binCounts[cellIndex] <- cells.list[[cellIndex]]$read.count[binIndex]
        }
        totalPerBin[binIndex] <- sum(binCounts)
    }

    output <- data.frame(seqname    = cells.list[[1]]$seqname,
                         start      = cells.list[[1]]$start,
                         end        = cells.list[[1]]$end,
                         read.count = totalPerBin)

    return (output)
}


#' Calculate a correction factor for the sequencing depth of each bin.
#'
#' @param cells.list  A list of data frames returned by 'coveragePerBin'.
#'
#' @return A list with the correction factor per bin.
#'
#' @export

determineCorrectionFactorPerBin <- function (cells.list)
{
    numberOfCells <- length(cells.list)
    numberOfBins  <- nrow(cells.list[[1]])

    medianPerBin  <- numeric(numberOfBins)
    for (binIndex in 1:numberOfBins)
    {
        binCounts <- numeric(numberOfCells)
        for (cellIndex in 1:numberOfCells)
        {
            binCounts[cellIndex] <- cells.list[[cellIndex]]$read.count[binIndex]
        }
        medianPerBin[binIndex] <- median(binCounts)
    }

    averageBinCount <- mean(medianPerBin)
    output <- averageBinCount / medianPerBin

    return (output)
}

#' Plot the correction factors calculated by ‘determineCorrectionFactorPerBin’.
#'
#' @param cells.list         The same parameter passed to
#'                           ‘determineCorrectionFactorPerBin’.
#' @param correction.factors The output of ‘determineCorrectionFactorPerBin’.
#'
#' @importFrom ggplot2 ggplot ylim geom_point geom_line theme aes xlab ylab element_text rel
#'
#' @return A list of ggplot2 objects (one per chromosome).
#'
#' @export

plotCorrectionFactorPerBin <- function (cells.list, correction.factors)
{
    # Please the static-analysis tool.
    correction.factor <- NULL
    bin <- NULL

    df <- data.frame(bin = paste0(cells.list[[1]]$seqname, ":",
                                  cells.list[[1]]$start, "-",
                                  cells.list[[1]]$end),
                     correction.factor = correction.factors)

    chromosomeNames <- unique(cells.list[[1]]$seqname)

    plots <- lapply (chromosomeNames, function (chromosome)
    {
        chr_df <- df[grep (paste0("^", chromosome, ":"), df$bin, perl = TRUE),]
        plot <- ggplot (chr_df, aes(x=bin, y=correction.factor, group=1)) +
            xlab("Bins") +
            ylab("Correction factor") +
            ylim(0, 2) +
            geom_point() +
            geom_line() +
            theme(axis.text.x = element_text(size  = rel(0.5),
                                             angle = 90,
                                             vjust = 0.5,
                                             hjust = 1))
        return (plot)
    })

    return (plots)
}

#' Create a summary table with events per cell.
#'
#' @param base_directory  The directory in which AneuFinder output was written.
#' @param samplesheet     The samplesheet used to run the pipeline.
#'
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom IRanges subsetByOverlaps
#' @importFrom S4Vectors elementMetadata
#'
#' @return A data frame containing the summary counts on success,
#'         or NULL otherwise.
#' @export

summaryCountsTable <- function (base_directory, samplesheet)
{
    samples                <- samplesheet[["sample_name"]]
    number_of_samples      <- length(samples)

    ## Pre-allocate numeric arrays
    name                   <- character(number_of_samples)
    cell                   <- numeric(number_of_samples)
    donor                  <- numeric(number_of_samples)
    num.segments           <- numeric(number_of_samples)
    num.chromosomal.gains  <- numeric(number_of_samples)
    num.chromosomal.losses <- numeric(number_of_samples)
    num.gains              <- numeric(number_of_samples)
    num.losses             <- numeric(number_of_samples)
    num.reciprocal         <- numeric(number_of_samples)
    num.nonReciprocal      <- numeric(number_of_samples)

    ## This will be set when the first sample is loaded.
    chromosomes <- NULL

    ## Look up scores in the data files
    for (sample_index in 1:number_of_samples)
    {
        sample_name <- samples[sample_index]
        row         <- samplesheet[which(samplesheet$sample_name == sample_name),]

        file_name   <- Sys.glob(paste0(base_directory,
                                       "/*/MODELS/method-edivisive/",
                                       sample_name, "_dedup.bam_*.RData"))

        name[sample_index]  <- sample_name
        cell[sample_index]  <- row[["cell"]]
        donor[sample_index] <- row[["embryo"]]

        if (! identical(file_name, character(0))) {
            sample <- get(load(file_name))

            chromosome_ranges <- sample[["segments"]]

            ## The first sample to load is used to determine the chromosomes
            ## and their lengths.
            if (is.null(chromosomes)) {
                chromosomes               <- as.data.frame(seqinfo(chromosome_ranges))
                chromosomes[["seqnames"]] <- rownames(chromosomes)
            }

            losses <- chromosome_ranges[(elementMetadata(chromosome_ranges)[,"copy.number"] < 2)]
            gains  <- chromosome_ranges[(elementMetadata(chromosome_ranges)[,"copy.number"] > 2)]

            num.segments[sample_index]  <- sample$qualityInfo$num.segments
            num.losses[sample_index]    <- length(losses)
            num.gains[sample_index]     <- length(gains)

            ## Determine whether whole-chromosome events occurred.
            num.chromosomal.losses[sample_index] <- 0
            num.chromosomal.gains[sample_index]  <- 0
            for (chromosome in chromosomes[["seqnames"]]) {
                total_length <- chromosomes[chromosome,"seqlengths"]
                query        <- GRanges(seqnames = chromosome,
                                        ranges   = IRanges(start = 1,
                                                           end   = total_length))

                chromosome_losses <- subsetByOverlaps(losses, query)
                loss_area         <- sum(width(chromosome_losses))

                chromosome_gains  <- subsetByOverlaps(gains, query)
                gain_area         <- sum(width(chromosome_gains))

                ## Losing or gaining 90% of the total chromosome length is
                ## considered a whole-chromosome loss or gain.
                if (loss_area / total_length > 0.9) {
                    num.chromosomal.losses[sample_index] = num.chromosomal.losses[sample_index] + 1
                }
                else if (gain_area / total_length > 0.9) {
                    num.chromosomal.gains[sample_index] = num.chromosomal.gains[sample_index] + 1
                }
            }
        } else {
            num.segments[sample_index]           <- NA
            num.chromosomal.gains[sample_index]  <- NA
            num.chromosomal.losses[sample_index] <- NA
            num.gains[sample_index]              <- NA
            num.losses[sample_index]             <- NA
        }
    }

    output <- data.frame (name,
                          cell,
                          donor,
                          num.segments,
                          num.chromosomal.gains,
                          num.chromosomal.losses,
                          num.gains,
                          num.losses,
                          num.reciprocal,
                          num.nonReciprocal)

    return (output)
}


#' Get the number of reads in a region.
#'
#' @param bamFilename            The BAM file to look for reads.
#' @param region                 A region specification in the form
#'                               'chr:start-end'.
#' @param minimumMappingQuality  The minimum mapping quality to include.
#'                               This should be a value between 0 and 60,
#'                               and defaults to 0 (include all reads).
#'
#' @return The number of reads found in that region.
#'
#' @useDynLib USCDtools, .registration = TRUE
#'
#' @export

readsInRegion <- function (bamFilename, region, minimumMappingQuality = 0)
{
    .Call ("count_reads_for_range", bamFilename, region, as.integer(minimumMappingQuality))
}

#' Get the number of reads for a list of regions.
#'
#' @param bamFilename            The BAM file to look for reads.
#' @param regions                A vector of regions in the form
#'                               'chr:start-end'.
#' @param minimumMappingQuality  The minimum mapping quality to include.
#'                               This should be a value between 0 and 60,
#'                               and defaults to 0 (include all reads).
#'
#' @return A list with the number of reads in each region.
#'
#' @useDynLib USCDtools, .registration = TRUE
#'
#' @export

readsInRegions <- function (bamFilename, regions, minimumMappingQuality = 0)
{
    .Call ("count_reads_for_ranges", bamFilename, regions, as.integer(minimumMappingQuality))
}

#' Create a BAM index
#'
#' @param bamFilename  The BAM file to create an index for.
#'
#' @return TRUE on sucess, FALSE on failure.
#'
#' @useDynLib USCDtools, .registration = TRUE
#'
#' @export

createBamIndex <- function (bamFilename)
{
    .Call ("create_bam_index", bamFilename)
}

#' Create a symbolic link.
#'
#' @param target       The existing target path.
#' @param destination  The non-existing destination path.
#'
#' @return TRUE on success, FALSE on failure.
#'
#' @export

createSymlink <- function (target, destination)
{
    .Call ("create_symbolic_link", target, destination)
}
