source(system.file("scripts/events.R",
                   package="alternativeSplicingEvents.hg19"))

## The latest release of the annotation files for the Human genome (hg19
## assembly) from MISO and VAST-TOOLS were retrieved from the following links:
##
## - MISO: https://miso.readthedocs.io/en/fastmiso/annotation.html
## - VAST-TOOLS: http://vastdb.crg.eu/libs/
##
## The downloaded folders contain GFF3 files (MISO) or tab-delimeted tables
## (VAST-TOOLS) with the alternative splicing events for each event type.

## SUPPA and rMATS identify alternative splicing events and generate annotation
## files from a transcripts annotation file in GTF format. The annotation file
## was retrieved from UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables)
## by selecting the GRCh37/hg19 assembly, "Genes and Gene Predictions" group,
## "UCSC Genes" track, "knownGene" table for all genome in the GTF format.
## Unfortunately, this file has a bug where the "transcript_id" is the same as
## the "gene_id". To fix the transcripts' identifier from the GTF file, the
## exact same information was retrieved from UCSC Table Browser in TXT format
## and the following steps were taken:

annotationGTF <- "ensembl_hg19.gtf"
annotationTXT <- "ensGene.txt"

# Create a match table between Gene ID and Transcript ID using the Ensembl
# annotation from UCSC (TXT file)
library(data.table) # faster to load data frames
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
## FASTQ samples). Just as with the other programs, the resulting alternative
## splicing events are divided in files by event type (IOE files for SUPPA and
## TXT for rMATS -- including novel events in the case of rMATS).

## After obtaining the resulting annotation files, all files are parsed to
## obtain the identifier, chromosome, strand and coordinates of each splicing
## event per event type. Some rMATS coordinates increment one unit to be
## compared to the annotation of other programs.
events <- list(
    miso  = getEventsFromMisoAnnotation(folder = "miso_annotation"),
    mats  = getEventsFromMatsAnnotation(folder = "mats_output"),
    suppa = getEventsFromSuppaAnnotation(folder = "suppa_output"),
    vast  = getEventsFromVastToolsAnnotation(folder = "vast_annotation")
)

## The "chr" prefix from the chromosome field is removed
for (each in seq_along(events)) {
    events[[each]]$Chromosome <- gsub("chr", "", events[[each]]$Chromosome)
}

## Organise splicing events by event type and then by program in a list of list
## of dataframes
events <- plyr::rbind.fill(events)
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
suppaAFE <- join$AFE$SUPPA.C2.start
matsAFE  <- join$AFE$MATS.C2.start
AFE.C2.start <- as.numeric( ifelse(
    sapply(suppaAFE, is.null),
    ifelse(sapply(matsAFE, is.null), NA, unlist(matsAFE)),
    unlist(suppaAFE)))

suppaALE <- join$ALE$SUPPA.C1.end
matsALE  <- join$ALE$MATS.C1.end
ALE.C1.end <- as.numeric( ifelse(
    sapply(suppaALE, is.null),
    ifelse(sapply(matsALE, is.null), NA, unlist(matsALE)),
    unlist(suppaALE)))

## Organise columns
annot <- lapply(names(join), function(i) {
    join[[i]][, c("Chromosome", "Strand", getSplicingEventCoordinates(i),
                  grep("Event.ID", names(join[[i]]), value = TRUE))]
})
names(annot) <- names(join)
annot$AFE[["C2.start"]] <- AFE.C2.start
annot$ALE[["C1.end"]]   <- ALE.C1.end
events <- annot

## Clean the annotation (make it ready for external use) and join everything in
## one organised table
types <- c(SE="Skipped exon", MXE="Mutually exclusive exon",
           A3SS="Alternative 3' splice site", A5SS="Alternative 5' splice site",
           AFE="Alternative first exon", ALE="Alternative last exon",
           RI="Retained intron", TandemUTR="Tandem UTR")

for (type in names(types)) {
    events[[type]] <- cbind("Event type"=types[[type]], events[[type]])
}
events <- plyr::rbind.fill(events)

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

# Order rows by event type, chromosome and the first exons coordinates and order
# columns
ord <- order(events$`Event type`, events$Chromosome,
             events$`Constitutive exon 1 start`,
             events$`Constitutive exon 1 end`,
             events$`Alternative exon 1 start`,
             events$`Alternative exon 1 end`)
events <- events[ord, c("Event type", "Chromosome", "Strand", "Gene",
                        coords, eventId)]

## Identify the gene associated with each splicing event by loading a file with
## coordinates of transcripts and the respective gene symbol, as retrieved from
## the UCSC table browser by selecting Human (hg19 assembly), "Gene and Gene
## Predictions" group, "UCSC genes" track, "known genes" table for all genome.
## In the output format, choose "selected fields from primary and related
## tables" and name the output file "gene_coordinates.txt". After clicking on
## "get output", choose fields "chrom", "strand", "txStart" and "txEnd". Also,
## select the field "geneSymbol" from the kgXref table.

coords <- function(FUN, ...) unname(apply(events[ , 5:12], 1, FUN, ...))
mini <- coords(min, na.rm=TRUE)
maxi <- coords(max, na.rm=TRUE)
df <- cbind(events[2:3], mini, maxi)
names(df) <- c("chr", "strand", "start", "end")
eventsGR <- GenomicRanges::makeGRangesFromDataFrame(df)

gene_coordinates <- read.delim("gene_coordinates.txt", header=FALSE,
                               comment.char="#")
names(gene_coordinates) <- c("chr", "strand", "start", "end", "gene")
gene_coordinates <- unique(gene_coordinates)
gene_coordinates$chr <- gsub("^chr", "", gene_coordinates$chr)
gene_pos <- GenomicRanges::makeGRangesFromDataFrame(gene_coordinates,
                                                    keep.extra.columns = T)

# Overlap coordinates to find gene names
overlaps <- GenomicRanges::findOverlaps(eventsGR, gene_pos)
query <- as.character(gene_pos$gene[S4Vectors::subjectHits(overlaps)])
names(query) <- S4Vectors::queryHits(overlaps)

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

# Filter events with more than 4 genes (probably annotation made by error)
fewGenes <- vapply(events$Gene, length, numeric(1)) <= 4
events <- events[fewGenes, ]

# Organise by event type and remove columns with NAs only
events <- split(events[-1], events$`Event type`)
for (type in names(events)) {
    naCols <- apply(events[[type]], 2, function(col) all(is.na(col)))
    events[[type]] <- events[[type]][!naCols]
}

# Save variable in rda
alternativeSplicingEvents.hg19 <- events
save(alternativeSplicingEvents.hg19, file="alternativeSplicingEvents.hg19.rda")
