insideFile <- function(...)
    source(system.file(..., package="alternativeSplicingEvents.hg38"))

insideFile("scripts/events-miso.R")
insideFile("scripts/events-mats.R")
insideFile("scripts/events-suppa.R")
insideFile("scripts/events-vast-tools.R")

#' Creates a template of alternative splicing junctions
#'
#' @param nrow Integer: Number of rows
#' @param program Character: Program used to get the junctions
#' @param event.type Character: Event type of the respective events
#' @param chromosome Character: Chromosome of the junctions
#' @param strand Character: positive ("+") or negative ("-") strand of the event
#' @param id Character: events' ID
#'
#' @return A data frame with the junctions coordinate names pre-filled with NAs
#'
#' @examples
#' psichomics:::createJunctionsTemplate(nrow = 8)
createJunctionsTemplate <- function(nrow, program = character(0),
                                    event.type = character(0),
                                    chromosome = character(0),
                                    strand = character(0),
                                    id = character(0)) {
    ## TODO(NunoA): only accept a "+" or a "-" strand
    parsed <- as.data.frame(matrix(NA, nrow = nrow, ncol = 8),
                            stringsAsFactors = FALSE)
    names(parsed) <- c("C1.start", "C1.end",
                       "A1.start", "A1.end",
                       "A2.start", "A2.end",
                       "C2.start", "C2.end")

    if (length(program) > 0)    parsed[["Program"]] <- "MISO"
    if (length(event.type) > 0) parsed[["Event.type"]] <- event.type
    if (length(chromosome) > 0) parsed[["Chromosome"]] <- chromosome
    if (length(strand) > 0)     parsed[["Strand"]] <- strand
    if (length(id) > 0)         parsed[["Event.ID"]] <- id
    return(parsed)
}

#' Get events from alternative splicing annotation
#'
#' @param folder Character: path to folder
#' @param types Character: type of events to retrieve (depends on the program of
#' origin)
#' @param genome Character: genome of interest (for instance, "hg19"; depends on
#' the program of origin)
#'
#' @importFrom utils read.delim
#' @importFrom plyr rbind.fill
#'
#' @return Retrieve data frame with events based on a given alternative splicing
#' annotation
getEventsFromMisoAnnotation <- function(
    folder,
    types=c("SE", "AFE", "ALE", "MXE", "A5SS", "A3SS", "RI", "TandemUTR"),
    genome="hg19") {

    cat("Retrieving MISO annotation...", fill=TRUE)
    typesFile <- file.path(folder, paste0(types, ".", genome, ".gff3"))
    annot <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                    comment.char="#", header=FALSE)

    ## TODO: ALE events are baldy formatted, they have two consecutive gene
    ## lines... remove them for now
    annot[[3]] <- annot[[3]][-c(49507, 49508), ]

    cat("Parsing MISO annotation...", fill=TRUE)
    events <- lapply(annot, parseMisoEvent)
    events <- plyr::rbind.fill(events)
    return(events)
}

#' @rdname getEventsFromMisoAnnotation
getEventsFromSuppaAnnotation <- function(
    folder,
    types=c("SE", "AF", "AL", "MX", "A5", "A3", "RI"),
    genome="hg19") {

    cat("Retrieving SUPPA annotation...", fill=TRUE)
    typesFile <- file.path(folder, paste0(genome, "_", types, ".ioe"))
    annot <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                    comment.char="#", header=TRUE)

    cat("Parsing SUPPA annotation...", fill=TRUE)
    eventsID <- lapply(annot, "[[", "event_id")
    events <- lapply(eventsID, parseSuppaEvent)
    events <- plyr::rbind.fill(events)
    return(events)
}

#' @rdname getEventsFromMisoAnnotation
getEventsFromMatsAnnotation <- function(
    folder,
    types=c("SE", "AFE", "ALE", "MXE", "A5SS", "A3SS", "RI"),
    genome="fromGTF") {

    cat("Retrieving rMATS annotation...", fill=TRUE)
    typesFile <- file.path(folder, paste(
        genome, c(types, paste0("novelEvents.", types)), "txt", sep = "."))
    names(typesFile) <- rep(types, 2)
    annot <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                    comment.char="#", header=TRUE)

    cat("Parsing rMATS annotation...", fill=TRUE)
    types <- names(annot)
    events <- lapply(seq_along(annot), function(i)
        if (nrow(annot[[i]]) > 0)
            return(parseMatsEvent(annot[[i]], types[[i]])))
    events <- plyr::rbind.fill(events)

    # Sum 1 position to the start/end of MATS events (depending on the strand)
    matsNames <- names(events)
    plus <- events$Strand == "+"
    # Plus
    start <- matsNames[grep(".start", matsNames)]
    events[plus, start] <- events[plus, start] + 1
    # Minus
    end <- matsNames[grep(".end", matsNames)]
    events[!plus, end] <- events[!plus, end] + 1
    return(events)
}

#' @rdname getEventsFromMisoAnnotation
getEventsFromVastToolsAnnotation <- function(
    folder,
    types=c("ALT3", "ALT5", "COMBI", "IR", "MERGE3m", "MIC", "EXSK", "MULTI"),
    genome="Hsa") {

    cat("Retrieving VAST-TOOLS annotation...", fill=TRUE)
    typesFile <- file.path(folder,
                           sprintf("%s.%s.Template%s.txt", genome, types,
                                   c(rep("", 6), rep(".2", 2))#, rep(".2", 2))
    ))
    names(typesFile) <- types

    annot <- lapply(typesFile, read.delim, stringsAsFactors = FALSE,
                    comment.char="#", header=TRUE)

    cat("Parsing VAST-TOOLS annotation...", fill=TRUE)
    types <- names(annot)
    events <- lapply(seq_along(annot),
                     function(i) {
                         cat(types[i], fill=TRUE)
                         a <- annot[[i]]
                         if (nrow(a) > 0)
                             return(parseVastToolsEvent(a))
                     })
    events <- plyr::rbind.fill(events)
    events <- unique(events)
    return(events)
}

#' Returns the coordinates of interest for a given event type
#' @param type Character: alternative splicing event type
#' @return Coordinates of interest according to the alternative splicing event
#' type
getSplicingEventCoordinates <- function(type) {
    switch(type,
           "SE"   = c("C1.end", "A1.start", "A1.end", "C2.start"),
           "A3SS" = c("C1.end", "A2.start", "A1.start"),
           "A5SS" = c("A2.end", "C2.start", "A1.end"),
           "AFE"  = c("A2.start", "A2.end", "A1.start", "A1.end"),
           "ALE"  = c("A1.start", "A1.end", "A2.start", "A2.end"),
           "RI"   = c("C1.start", "C1.end", "C2.start", "C2.end"),
           "MXE"  = c("C1.end", "A1.start", "A1.end",
                      "A2.start", "A2.end", "C2.start"),
           "TandemUTR" = c("A2.start", "A2.end", "A1.end"))
}

#' Convert a column to numeric if possible and ignore given columns composed
#' of lists
#'
#' @param table Data matrix: table
#' @param by Character: column names of interest
#' @param toNumeric Boolean: which columns to convert to numeric (FALSE by
#' default)
#'
#' @return Processed data matrix
#' @examples
#' event <- read.table(text = "ABC123 + 250 300 350
#'                             DEF456 - 900 800 700")
#' names(event) <- c("Event ID", "Strand", "C1.end", "A1.end", "A1.start")
#'
#' # Let's change one column to character
#' event[ , "C1.end"] <- as.character(event[ , "C1.end"])
#' is.character(event[ , "C1.end"])
#'
#' event <- psichomics:::getNumerics(event, by = c("Strand", "C1.end", "A1.end",
#'                                   "A1.start"),
#'                                   toNumeric = c(FALSE, TRUE, TRUE, TRUE))
#' # Let's check if the same column is now integer
#' is.numeric(event[ , "C1.end"])
getNumerics <- function(table, by = NULL, toNumeric = FALSE) {
    # Check which elements are lists of specified length
    bool <- TRUE
    for (each in by)
        bool <- bool & vapply(table[[each]], length, integer(1)) == 1

    table <- table[bool, ]
    # Convert elements to numeric
    conv <- by[toNumeric]
    table[conv] <- as.numeric(as.character(unlist(table[conv])))
    return(table)
}

#' Full outer join all given events based on select columns
#' @param events Data frame or matrix: alternative splicing events
#' @param types Character: alternative splicing types
#' @return List of events joined by alternative splicing event type
joinEventsPerType <- function(events, types) {
    if (missing(types)) types <- names(events)
    joint <- lapply(types, function(type, events) {
        cat(type, fill=TRUE)
        # Create vector with comparable columns
        id <- c("Strand", "Chromosome", "Event.type")
        by <- c(id, getSplicingEventCoordinates(type))
        toNumeric <- !by %in% id

        # Convert given columns to numeric if possible
        tables <- lapply(events[[type]], getNumerics, by, toNumeric)

        # Make the names of non-comparable columns distinct
        cols <- lapply(names(tables), function(k) {
            ns <- names(tables[[k]])
            inBy <- ns %in% by
            ifelse(inBy, ns, paste(k, ns, sep="."))
        })

        # Full join all the tables
        res <- Reduce(function(x, y) dplyr::full_join(x, y, by), tables)
        names(res) <- unique(unlist(cols))

        # Remove equal rows
        return(unique(res))
    }, events)
    names(joint) <- types
    return(joint)
}

#' Write events to a file
#'
#' @param jointEvents List of lists of data frame
#' @param eventType Character: type of event
#' @param filename Character: path to the file
#' @param showID Boolean: show the events' ID? FALSE by default
#' @param rds Boolean: write to a RDS file? TRUE by default; otherwise, write to
#' TXT
#'
#' @importFrom utils write.table
#'
#' @return Invisible TRUE if successful
writeAnnotation <- function(jointEvents, eventType,
                            filename = paste0("data/annotation_",
                                              eventType, ".txt"),
                            showID = FALSE, rds = TRUE) {
    res <- jointEvents[[eventType]]
    # Show the columns Chromosome, Strand and coordinates of interest
    by <- c("Chromosome", "Strand", getSplicingEventCoordinates(eventType))
    ord <- 0

    # Show the events' ID if desired
    if (showID) {
        cols <- grep("Event.ID", names(res), value = TRUE)
        by <- c(cols, by)
        ord <- length(cols)
    }
    res <- subset(res, select = by)

    ## TODO(NunoA): clean this mess
    # Order by chromosome and coordinates
    orderBy <- lapply(c(1 + ord, (3 + ord):ncol(res)),
                      function(x) return(res[[x]]))
    res <- res[do.call(order, orderBy), ]

    res <- unique(res)

    if (rds)
        saveRDS(res, file = filename)
    else
        write.table(res, file = filename, quote = FALSE, row.names = FALSE,
                    sep = "\t")
    return(invisible(TRUE))
}

#' Returns the coordinates of interest to sort for a given event type
#' @export
sortingCoordinates <- function(type) {
    switch(type,
           "A3SS" = c("A2.start", "A1.start"),
           "A5SS" = c("A2.end", "A1.end"),
           "AFE"  = c("A1.start", "A1.end", "A2.start", "A2.end"),
           "ALE"  = c("A1.start", "A1.end", "A2.start", "A2.end"),
           "MXE"  = c("A1.start", "A1.end", "A2.start", "A2.end"),
           "TandemUTR" = c("A1.end", "A2.end"))
}

colsAsNumbers <- function(type, annotation) {
    # Create vector with comparable columns
    id <- c("Strand", "Chromosome", "Event.type")
    by <- c(id, getSplicingEventCoordinates(type))
    toNumeric <- !by %in% id

    # Convert given columns to numeric if possible
    tables <- lapply(annotation[[type]], getNumerics, by, toNumeric)
    return(tables)
}
