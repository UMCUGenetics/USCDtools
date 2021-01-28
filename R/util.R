## Copyright (C) 2021  Roel Janssen <roel@gnu.org>

## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

#' SPARQL function for SPARQLing-genomics
#'
#' @param endpoint      The endpoint URI.
#' @param projectId     The SPARQLing-genomics project ID.
#' @param token         An authentication token.
#' @param query         A SPARQL query to execute.
#'
#' @return A data frame containing the query results.
#'
#' @importFrom RCurl    basicTextGatherer curlPerform
#' @importFrom jsonlite fromJSON
#'
#' @export

SGSPARQL <- function (endpoint, projectId, token, query)
{
    accumulator <- basicTextGatherer()
    accumulator$reset()

    curlPerform(url           = paste0(endpoint, "?project-id=", projectId),
                httpheader    = c("Accept"       = "application/json",
                                  "Cookie"       = paste0("SGSession=", token),
                                  "Content-Type" = "application/sparql-update"),
                customrequest = "POST",
                postfields    = query,
                writefunction = accumulator$update)

    jsonData    <- accumulator$value()
    data        <- fromJSON(jsonData)

    accumulator$reset()
    return (data)
}


#' Function to merge and sort BAM files.
#'
#' @param filelist     A list of input BAMs
#' @param output_path  The path to the output BAM filename.
#' @param threads      Number of threads to use to merge and sort.
#'
#' @return TRUE when all went well, or FALSE otherwise.
#'
#' @export

createCompositeBam <- function (filelist, output_path, threads=32)
{
    if (file.exists(output_path))
    {
        cat(paste0("Not overwriting ", output_path, ".\n"))
        return (FALSE)
    }
    else
    {
        mergeCommand <- paste0 ("samtools merge -@ 32 -O bam ",
                                output_path, " ",
                                paste(filelist, collapse=" "))

        cat(paste0("Merging ", length(filelist),
                   " files to ", output_path ,".\n"))
        tryCatch(system (mergeCommand), error = function (e) {
            cat(paste0("Merging failed.\n"))
            return (FALSE)
        })
    }

    ### The merge process already sorts it.
    ## sortCommand <- paste0 ("samtools sort -T ", dirname (output_path),
    ##                        " -@ ", threads ," -O bam -o ",
    ##                        str_replace (output_path, ".bam", ".sorted.bam"),
    ##                        " ", output_path)
    ##
    ## cat(paste0("Sorting ", output_path ,".\n"))
    ## tryCatch(system (sortCommand), error = function (e) {
    ##     cat(paste0("Sorting failed.\n"))
    ##     return (FALSE)
    ## })

    return (TRUE)
}

#' Return a data frame containing the chromosome name and its length in bp.
#'
#' @param genome             A BSgenome object.
#' @param chromosomeFilter   A vector containing the chromosome names to keep.
#'
#' @return A data frame containing the chromosome name and its length in bp.
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

### IDEA: Also include the bins that fall into telomere and centromere regions.
#' Determine list of unreliable regions
#'
#' This function divides the genome in bins and returns a list of the
#' 2% bins with the least number of reads, and of the 3% bins with the
#' most number of reads for autosomal chromosomes, and the least 3% and
#' most 5% for allosomal chromosomes.
#'
#' @param compositeBam     The BAM file to determine unreliable regions from.
#' @param genome           A BSgenome object matching the reference genome used
#'                         in the 'compositeBam'.
#' @param autosomes        A vector containing the autosomal chromosome names to keep.
#' @param allosomes        A vector containing the sex chromosome names to keep.
#' @param binSize          The size of a single bin.
#' @param threads          Number of threads to use for this process.
#'
#' @return A GRanges object containing the regions to exclude from further
#'         analysis.
#'
#' @importFrom GenomicRanges tileGenome GRanges
#' @importFrom tools         file_path_sans_ext
#'
#' @export

prepareRegionalBlacklist <- function (compositeBam,
                                      genome,
                                      binSize,
                                      autosomes,
                                      allosomes,
                                      threads=32)
{
    compositeBam_sans_ext <- tools::file_path_sans_ext(compositeBam)
    chromosomeFilter <- c(autosomes, allosomes)

    index_filename <- paste0 (compositeBam, ".bai")
    if (! file.exists (index_filename))
    {
        system (paste0 ("samtools index -@ ", threads, " ", compositeBam))
    }

    idxstats_filename <- paste0 (compositeBam, ".idxstats")
    idxstats_errors   <- paste0 (compositeBam, ".idxstats.errors")
    if (! file.exists (idxstats_filename))
    {
        system2 ('samtools', args = c('idxstats', compositeBam),
                 stdout = idxstats_filename,
                 stderr = idxstats_errors)
    }

    chromosome_lengths <- chromosomeLengths (genome, chromosomeFilter)
    chromosome_lengths[,1]
    chromosome_lengths[,2]
    chromosomes <- GRanges (seqnames = chromosome_lengths[,1],
                            IRanges (start = 1, end = chromosome_lengths[,2]))

    seqlengths (chromosomes) <- chromosome_lengths[,2]

    bins <- tileGenome (seqlengths             = chr_lengths,
                        tilewidth              = as.numeric (binSize),
                        cut.last.tile.in.chrom = TRUE)

    bins_output <- data.frame(bins, stringsAsFactors = F)[,1:3]

    # The coverage on the X chromosome is expected to be lower in males.
    # Therefore, we treat it seperately from the autosomes.

    bins_file <- paste0 (compositeBam_sans_ext, "_",
                         format(binSize, scientific=FALSE), "bp_bins.bed")

    write.table (bins_output,
                 file      = bins_file,
                 sep       = "\t",
                 quote     = FALSE,
                 col.names = FALSE,
                 row.names = FALSE)

    intersect_output_file   <- paste0(compositeBam_sans_ext, "_intersect.bed")
    intersect_output_errors <- paste0(compositeBam_sans_ext, "_intersect.errors")

    exit_code <- system2 ('bedtools',
                          args = c('intersect',
                                   '-a', bins_file,
                                   '-b', compositeBam,
                                   '-sorted', '-c', '-g', idxstats_filename),
                          stdout = intersect_output_file,
                          stderr = intersect_output_errors)

    if (exit_code != 0)
    {
        cat (paste0 ("Failed to run 'bedtools'.\n"))
        return (FALSE)
    }

    coverage_per_bin_total <- read.delim (intersect_output_file,
                                          stringsAsFactors = FALSE,
                                          header           = FALSE)

    coverage_per_bin_autosomal <- coverage_per_bin_total[coverage_per_bin_total[,1] %in% autosomes]
    low_cutoff  <- quantile (coverage_per_bin_autosomal[,4], 0.02)
    high_cutoff <- quantile (coverage_per_bin_autosomal[,4], 0.97)

    autosomal_excluded_bins <- coverage_per_bin_autosomal[coverage_per_bin_autosomal[,4] <= low_cutoff |
                                                          coverage_per_bin_autosomal[,4] >= high_cutoff,]


    coverage_per_bin_allosomal <- coverage_per_bin_total[coverage_per_bin_total[,1] %in% allosomes]
    low_cutoff  <- quantile (coverage_per_bin_allosomal[,4], 0.03)
    high_cutoff <- quantile (coverage_per_bin_allosomal[,4], 0.95)

    allosomal_excluded_bins <- coverage_per_bin_allosomal[coverage_per_bin_allosomal[,4] <= low_cutoff |
                                                          coverage_per_bin_allosomal[,4] >= high_cutoff,]

    combined_excluded_bins <- rbind (autosomal_excluded_bins, allosomal_excluded_bins)

    # Merge bins together when they are less one bin apart.
    combined <- GRanges (seqnames = combined_excluded_bins[,1],
                         IRanges(start = combined_excluded_bins[,2] - (binSize / 2),
                                 end   = combined_excluded_bins[,3] - (binSize / 2)))
    output   <- reduce (combined)

    start (output) <- start (output) + (binSize / 2)
    end   (output) <- end   (output) + (binSize / 2)

    return (output)
}
