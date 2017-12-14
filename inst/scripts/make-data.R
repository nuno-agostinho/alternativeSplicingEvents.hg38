require(plyr)
require(GenomicRanges)
require(S4Vectors)
source(system.file("scripts/events.R",
                   package="alternativeSplicingEvents.hg38"))

# Obtain alternative splicing annotation -------------------------------------

## The latest release of the annotation files for the Human genome (hg19
## assembly) from MISO and VAST-TOOLS were retrieved from the following links:
##
## - MISO: https://miso.readthedocs.io/en/fastmiso/annotation.html
## - VAST-TOOLS: http://vastdb.crg.eu/libs/
##
## The downloaded folders contain GFF3 files (MISO) or tab-delimeted tables
## (VAST-TOOLS) with the alternative splicing events for each event type.

## SUPPA and rMATS identify alternative splicing events and generate annotation
## files from a transcript annotation file in GTF format. The annotation file
## was retrieved from UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables)
## by selecting the GRCh37/hg19 assembly, "Genes and Gene Predictions" group,
## "UCSC Genes" track, "knownGene" table for all genome in the GTF format.
## Misleadingly, the "transcript_id" column contains gene identifiers. As such,
## the proper transcript identifiers were retrieved from UCSC Table Browser in
## TXT format and and the following steps were taken:

annotationGTF <- "ensembl_hg19.gtf"
annotationTXT <- "ensGene.txt"

# Create a match table between Gene ID and Transcript ID using the Ensembl
# annotation from UCSC (TXT file)
require(data.table) # faster to load data frames
txt <- fread(annotationTXT, data.table = FALSE)
idTable <- txt$V13          # Save gene ID
names(idTable) <- txt$V2    # Save transcript ID

# Retrieve transcript IDs from the GTF file
gtf <- fread(annotationGTF, data.table = FALSE)
transcriptid <- gsub(".*\\\"(.*)\\\".*", "\\1", gtf$V9)

# Get gene IDs for the respective transcript ID
geneid <- idTable[transcriptid]

# Correctly place the gene ID in the GTF annotation
# WARNING: this erases other GTF attributes
gtf$V9 <- sprintf("gene_id \"%s\"; transcript_id \"%s\";", geneid, transcriptid)

# Save as a tab-delimited file with no header and no quotes
outputFilename <- "ensembl_hg19_fixed.gtf"
write.table(gtf, file = outputFilename, quote = FALSE, col.names = FALSE,
            sep = "\t", row.names = FALSE)

## The resulting transcripts annotation file is used to generate the alternative
## splicing events with SUPPA (by running the command generateEvents for all
## event types) and rMATS (by running the main rMATS script with any two BAM or
## FASTQ samples; events generated are independent of the sample files used).
## The resulting alternative splicing events are split by files based on their
## event type (IOE files for SUPPA and TXT for rMATS -- including novel events
## in the case of rMATS).

# Parse alternative splicing annotation --------------------------------------

## After obtaining the resulting annotation files, all files are parsed to
## obtain the identifier, chromosome, strand and coordinates of each splicing
## event per event type. Some rMATS coordinates are incremented by one to be
## comparable to events from other annotations.
miso  <- getEventsFromMisoAnnotation(folder = "miso_annotation")
mats  <- getEventsFromMatsAnnotation(folder = "mats_annotation")
suppa <- getEventsFromSuppaAnnotation(folder = "suppa_annotation")

## Note the folder that matters in case of VAST-TOOLS is the TEMPLATES folder
vast  <- getEventsFromVastToolsAnnotation(folder = "Hsa/TEMPLATES")
events <- list(miso, mats, suppa, vast)

## The "chr" prefix from the chromosome field is removed
for (each in seq_along(events)) {
    events[[each]]$Chromosome <- gsub("chr", "", events[[each]]$Chromosome)
}

## Organise splicing events by event type and then by program in a list of list
## of dataframes
events <- rbind.fill(events)
events <- dlply(events, .(Event.type))
events <- lapply(events, dlply, .(Program))

## Sort coordinates for some event types (some programs differ the sorting of
## specific event types; to make them comparable, we can order all coordinates
## by increasing or decreasing order, depending on strand)
types <- names(events)
for (type in types) {
    coord <- sortingCoordinates(type)
    events[[type]] <- colsAsNumbers(type, events)
    if (!is.null(coord)) {
        for (program in names(events[[type]])) {
            print(paste(type, program))
            table <- events[[type]][[program]]
            plus <- table[["Strand"]] == "+"
            plusOrd <- apply(table[plus, coord], 1, sort)
            minusOrd <- apply(table[!plus, coord], 1, sort, decreasing = TRUE)
            events[[type]][[program]][plus, coord] <- t(plusOrd)
            events[[type]][[program]][!plus, coord] <- t(minusOrd)
        }
    }
}

## A full outer join is performed per event type. This allows to cross-reference
## alternative splicing events between the different programs.
join <- joinEventsPerType(events)

# Add 1st constitutive exon's end and 2nd constituve exon's start from SUPPA or
# rMATS to AFE and ALE events, respectively (other programs do not state these
# coordinates)
AFE.C2.start        <- join$AFE$SUPPA.C2.start
nulls               <- sapply(AFE.C2.start, is.null)
AFE.C2.start[nulls] <- join$AFE$MATS.C2.start[nulls]
nulls               <- sapply(AFE.C2.start, is.null)
AFE.C2.start[nulls] <- NA
AFE.C2.start        <- as.numeric(AFE.C2.start)

ALE.C1.end        <- join$ALE$SUPPA.C1.end
nulls             <- sapply(ALE.C1.end, is.null)
ALE.C1.end[nulls] <- join$ALE$MATS.C1.end[nulls]
nulls             <- sapply(ALE.C1.end, is.null)
ALE.C1.end[nulls] <- NA
ALE.C1.end        <- as.numeric(ALE.C1.end)

## Organise columns
annot <- lapply(names(join), function(i) {
    join[[i]][, c("Chromosome", "Strand", getSplicingEventCoordinates(i),
                  grep("Event.ID", names(join[[i]]), value = TRUE))]
})
names(annot) <- names(join)
annot$AFE[["C2.start"]] <- AFE.C2.start
annot$ALE[["C1.end"]]   <- ALE.C1.end
events <- annot

# Prepare combined annotation -----------------------------------------------

## Clean the annotation (make it ready for external use) and join everything in
## one organised table
types <- c(SE="Skipped exon", MXE="Mutually exclusive exon",
           A3SS="Alternative 3' splice site", A5SS="Alternative 5' splice site",
           AFE="Alternative first exon", ALE="Alternative last exon",
           RI="Retained intron", TandemUTR="Tandem UTR")

for (type in names(types)) {
    events[[type]] <- cbind("Event type"=types[[type]], events[[type]])
}
events <- rbind.fill(events)

coords <- c("C1.start"="Constitutive exon 1 start",
            "C1.end"="Constitutive exon 1 end",
            "A1.start"="Alternative exon 1 start",
            "A1.end"="Alternative exon 1 end",
            "A2.start"="Alternative exon 2 start",
            "A2.end"="Alternative exon 2 end",
            "C2.start"="Constitutive exon 2 start",
            "C2.end"="Constitutive exon 2 end")
names(events)[match(names(coords), names(events))] <- coords
events$Gene <- NA
eventId <- grep("Event.ID", names(events), value = TRUE)

# Order rows by event type, chromosome and the first exons coordinate and order
# columns
ord <- order(events$`Event type`, events$Chromosome,
             events$`Constitutive exon 1 start`,
             events$`Constitutive exon 1 end`,
             events$`Alternative exon 1 start`,
             events$`Alternative exon 1 end`)
events <- events[ord, c("Event type", "Chromosome", "Strand", "Gene",
                        coords, eventId)]

# Identify genes of alternative splicing events ------------------------------

## Identify the gene associated with each splicing event based on a file with
## coordinates of transcripts and their respective gene symbol

getGRangesFromCoordinates <- function(filename, byExon=TRUE) {
    gene_coordinates <- read.delim(filename, header=FALSE, comment.char="#")
    names(gene_coordinates) <- c("chr", "strand", "txStart", "txEnd",
                                 "exonStarts", "exonEnds", "gene")
    gene_coordinates <- unique(gene_coordinates)
    gene_coordinates$chr <- gsub("^chr", "", gene_coordinates$chr)

    if (byExon) {
        start <- gene_coordinates[["exonStarts"]]
        end   <- gene_coordinates[["exonEnds"]]
    } else {
        start <- gene_coordinates[["txStart"]]
        end   <- gene_coordinates[["txEnd"]]
    }

    start_split <- strsplit(as.character(start), ",")
    end_split   <- strsplit(as.character(end), ",")

    reps     <- sapply(start_split, length)
    chr_2    <- rep(gene_coordinates$chr, reps)
    strand_2 <- rep(gene_coordinates$strand, reps)
    gene_2   <- rep(gene_coordinates$gene, reps)

    start_split <- as.numeric(unlist(start_split))
    end_split   <- as.numeric(unlist(end_split))
    coordinates <- data.frame("chr"=chr_2, "strand"=strand_2,
                              "start"=start_split, "end"=end_split,
                              "gene"=gene_2)

    makeGRangesFromDataFrame(coordinates, keep.extra.columns=T)
}

assignGenesFromGRange <- function(events, pos) {
    # Get event type-specifc coordinates of alternative exons for gene
    # attribution
    df <- cbind(events[2:3], start=NA, end=NA)
    eventTypes <- events$`Event type`
    for (type in unique(eventTypes)) {
        print(type)
        rows <- eventTypes == type

        if (type %in% c("Skipped exon", "Mutually exclusive exon",
                        "Alternative first exon", "Alternative last exon")) {
            start <- events[rows, "Alternative exon 1 start"]
            end   <- events[rows, "Alternative exon 1 end"]
        } else if (type %in% "Alternative 3' splice site") {
            start <- events[rows, "Alternative exon 1 start"]
            end   <- events[rows, "Alternative exon 2 start"]
        } else if (type %in% "Alternative 5' splice site") {
            start <- events[rows, "Alternative exon 1 end"]
            end   <- events[rows, "Alternative exon 2 end"]
        } else if (type %in% "Retained intron") {
            start <- events[rows, "Constitutive exon 1 start"]
            end   <- events[rows, "Constitutive exon 1 end"]
        } else if (type %in% "Tandem UTR") {
            start <- events[rows, "Alternative exon 2 start"]
            end   <- events[rows, "Alternative exon 2 end"]
        } else {
            stop("Undefined event type:", type)
        }

        df[rows, "start"] <- ifelse(start < end, start, end)
        df[rows, "end"]   <- ifelse(start < end, end, start)
    }

    names(df) <- c("chr", "strand", "start", "end")
    eventsGR <- makeGRangesFromDataFrame(df)

    # Overlap coordinates to find gene names
    print("Looking for overlapping coordinates...")
    overlaps <- findOverlaps(eventsGR, pos)
    query <- as.character(pos$gene[subjectHits(overlaps)])
    names(query) <- queryHits(overlaps)

    # Remove duplicated gene names per alternative splicing event
    geneSymbol <- split(query, names(query))
    geneSymbol <- lapply(geneSymbol, unique)

    # Label events with no genes as "Hypothetical"
    seq <- as.character(1:length(eventsGR))
    noGene <- !seq %in% names(geneSymbol)
    geneSymbol[seq[noGene]] <- "Hypothetical"

    # Sort by numeric value and add to events
    ind <- as.character(sort(as.integer(names(geneSymbol))))
    events$Gene <- geneSymbol[ind]

    # Remove duplicate events
    print("Removing duplicate events...")
    return(unique(events))
}

## Retrieve coordinates from the UCSC table browser by selecting Human (hg19
## assembly), "Gene and Gene Predictions" group, "UCSC genes" track,
## "knownGene" table for all genome. In the output format, choose "selected
## fields from primary and related tables" and name the output file
## "gene_coordinates.txt". After clicking on "get output", choose fields
## "chrom", "strand", "txStart", "txEnd", "exonStarts" and "exonEnds". Also,
## select the field "geneSymbol" from the kgXref table.

# Get genes based on exon coordinates
exon_pos <- getGRangesFromCoordinates("gene_coordinates.txt", byExon=TRUE)
events   <- assignGenesFromGRange(events, exon_pos)

# Get remaining genes based on transcript coordinates
gene_pos <- getGRangesFromCoordinates("gene_coordinates.txt", byExon=FALSE)
events[events$Gene == "Hypothetical", ] <- assignGenesFromGRange(
    events[events$Gene == "Hypothetical", ], gene_pos)

# Convert from hg19 to hg38 coordinates -------------------------------------

## Convert hg19 to hg38 coordinates based on the UCSC chain file titled
## "hg19ToHg38.over.chain.gz"
require("rtracklayer")
chain <- "hg19ToHg38.over.chain.gz"
chain <- R.utils::gunzip(chain, remove=FALSE)
chain <- import.chain(chain)

convertToGRanges <- function(events, cols) {
    noNAs <- !is.na(events[[cols[1]]]) | !is.na(events[[cols[2]]])
    start <- events[[cols[1]]][noNAs]
    end   <- events[[cols[2]]][noNAs]

    # Avoid NAs
    startNA <- is.na(start)
    start   <- ifelse(startNA, end, start)
    endNA   <- is.na(end)
    end     <- ifelse(endNA, start, end)

    # Start must be a value less than or equal to end
    startLTend <- start < end
    start <- ifelse(startLTend, start, end)
    end   <- ifelse(startLTend, end, start)

    # Convert hg19 coordinates to hg38
    df <- data.frame(chr=paste0("chr", events$Chromosome[noNAs]),
                     start=start, end=end, strand=events$Strand[noNAs])
    lifted <- liftOver(makeGRangesFromDataFrame(df), chain)

    # Return coordinates as they were before
    start                      <- as.numeric(rep(NA, nrow(df)))
    startLifted                <- start(lifted)
    startLiftedLen             <- sapply(startLifted, length)
    start[startLiftedLen == 1] <- as.numeric(startLifted[startLiftedLen == 1])
    start[startLiftedLen != 1] <- NA

    end                    <- as.numeric(rep(NA, nrow(df)))
    endLifted              <- end(lifted)
    endLiftedLen           <- sapply(endLifted, length)
    end[endLiftedLen == 1] <- as.numeric(endLifted[endLiftedLen == 1])
    end[endLiftedLen != 1] <- NA

    start <- ifelse(startLTend, start, end)
    end   <- ifelse(startLTend, end, start)
    start[startNA] <- NA
    end[endNA]     <- NA

    # Replace hg19 coordinates with their respective hg38
    events[[cols[1]]][noNAs] <- start
    events[[cols[2]]][noNAs] <- end

    # Fix chromosomes if needed
    chr <- runValue(seqnames(lifted[startLiftedLen == 1]))
    chr <- as.character(unlist(chr))
    events$Chromosome[noNAs][startLiftedLen == 1] <- gsub("chr", "", chr,
                                                          fixed=TRUE)
    return(events)
}

constitutiveExon1 <- paste("Constitutive exon 1", c("start", "end"))
events <- convertToGRanges(events, constitutiveExon1)
constitutiveExon2 <- paste("Constitutive exon 2", c("start", "end"))
events <- convertToGRanges(events, constitutiveExon2)
alternativeExon1  <- paste("Alternative exon 1", c("start", "end"))
events <- convertToGRanges(events, alternativeExon1)
alternativeExon2  <- paste("Alternative exon 2", c("start", "end"))
events <- convertToGRanges(events, alternativeExon2)

# Remove events not properly converted to hg38 coordinates
for (type in unique(events$`Event type`)) {
    if (type == "Skipped exon") {
        coords <- c("Constitutive exon 1 end",
                    "Alternative exon 1 start", "Alternative exon 1 end",
                    "Constitutive exon 2 start")
    } else if (type == "Mutually exclusive exon") {
        coords <- c("Constitutive exon 1 end",
                    "Alternative exon 1 start", "Alternative exon 1 end",
                    "Alternative exon 2 start", "Alternative exon 2 end",
                    "Constitutive exon 2 start")
    } else if (type %in% c("Alternative 5' splice site",
                           "Alternative first exon")) {
        coords <- c("Constitutive exon 1 end", "Alternative exon 1 end",
                    "Constitutive exon 2 start")
    } else if (type %in% c("Alternative 3' splice site",
                           "Alternative last exon")) {
        coords <- c("Constitutive exon 1 end",
                    "Alternative exon 1 start", "Constitutive exon 2 start")
    } else if (type == "Retained intron") {
        coords <- c("Constitutive exon 1 start", "Constitutive exon 1 end",
                    "Constitutive exon 2 start", "Constitutive exon 2 end")
    } else if (type == "TandemUTR") {
        coords <- c("Constitutive exon 1 start", "Constitutive exon 1 end",
                    "Alternative exon 1 end")
    }

    print(type)
    # Check if any coordinates are missing values
    coordsOkay <- Reduce("&", lapply(coords, function(i) is.na(events[[i]])))

    eventType <- events$`Event type` == type
    events <- events[!(eventType & coordsOkay), ]
}

# Final clean up ------------------------------------------------------------

# Organise by event type and remove columns with NAs only
events <- split(events[-1], events$`Event type`)
for (type in names(events)) {
    naCols <- apply(events[[type]], 2, function(col) all(is.na(col)))
    events[[type]] <- events[[type]][!naCols]
}

# Save variable in rda
alternativeSplicingEvents.hg38 <- events
save(alternativeSplicingEvents.hg38, file="alternativeSplicingEvents.hg38.rda")
