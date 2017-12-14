\name{alternativeSplicingEvents.hg38.rda}
\alias{alternativeSplicingEvents.hg38.rda}
\docType{data}
\title{
Alternative splicing event annotation for Human (hg38)
}
\description{
Data frame containing alternative splicing events per event type. Each splicing
event is characterised by its chromosome, strand, splice junction coordinates,
associated gene and original event identifier from their source.

The following event types are available:
\itemize{
    \item{Skipped Exon (SE)}
    \item{Mutually Exclusive Exon (MXE)}
    \item{Alternative First Exon (AFE)}
    \item{Alternative Last Exon (ALE)}
    \item{Alternative 3' Splice Site (A3SS)}
    \item{Alternative 5' Splice Site (A5SS)}
    \item{Retained Intron (RI)}
    \item{Tandem UTR}
}

Alternative splicing events were compiled from the input annotation files used
by the splicing quantification tools MISO, VAST-TOOLS, SUPPA and rMATS,
followed by a lift-over of the Human annotation hg19 coordinates to hg38 (see
details).
}
%\usage{data("data/alternativeSplicingEvents.hg38.RDS")}
\format{
%  The format is: chr "data/alternativeSplicingEvents.hg38.RDS"
}
\details{
The compiled annotation of human alternative splicing events was based on
annotation files used as input to MISO, VAST-TOOLS, rMATS, and SUPPA.

Annotation files from MISO and VAST-TOOLS are provided in their respective
websites (see "Sources"), whereas SUPPA and rMATS identify alternative splicing
events and generate such annotation files based on a given isoform-centred
transcript annotation in the GTF format.

As such, the transcript annotation file was retrieved from UCSC Table Browser by
selecting the GRCh37/hg19 assembly, "Genes and Gene Predictions" group, "UCSC
Genes" track and "knownGene"" table for all genome in GTF and TXT formats, so
gene identifiers in the GTF file (misleadingly identical to transcript
identifiers) could be replaced with proper ones from the TXT version.

After obtaining the resulting annotation files, all files were parsed to obtain
the identifier, chromosome, strand and coordinates of each splicing event per
event type. Then, a full outer join was performed on the annotation of the
different programs per event type using \code{dplyr::full_join}.

Note that before the outer join, some rMATS coordinates must be incremented by
one to be comparable to other annotation files and the "chr" prefix in the
chromosome field of some annotation files also must be removed.

Finally, Ensembl gene symbols (retrieved from UCSC Table Browser) are assigned
to alternative splicing events based on overlapping exon and gene coordinates.
Splicing events not assigned to any genes are marked as "Hypothetical".

For more details, check the \code{make-data.R} file.

The coordinates from the combined annotation were lifted-over from Human
annotation hg19 to hg38 using \code{rtracklayer::liftOver} based on an UCSC
chain file.
}
\source{
    \itemize{
        \item{\strong{MISO:}
            \url{https://miso.readthedocs.io/en/fastmiso/annotation.html}}
        \item{\strong{rMATS:}
            \url{http://rnaseq-mats.sourceforge.net/user_guide.htm}}
        \item{\strong{SUPPA:}
            \url{https://bitbucket.org/regulatorygenomicsupf/suppa}}
        \item{\strong{VAST-TOOLS:}
            \url{http://vastdb.crg.eu/libs/}}
    }

    Transcript annotation files used to generate events from SUPPA and rMATS
    and annotation with gene symbols were retrieved from UCSC Table Browser
    (\url{https://genome.ucsc.edu/cgi-bin/hgTables}).
}
\references{
    \itemize{
        \item{\strong{MISO:}
            Katz Y, Wang ET, Airoldi EM, Burge CB. Analysis and design of RNA
            sequencing experiments for identifying isoform regulation. Nature
            methods. 2010;7(12):1009-1015. doi:10.1038/nmeth.1528.}
        \item{\strong{rMATS:}
            Shen S, Park JW, Lu Z, et al. rMATS: Robust and flexible detection
            of differential alternative splicing from replicate RNA-Seq data.
            Proceedings of the National Academy of Sciences of the United States
            of America. 2014;111(51):E5593-E5601. doi:10.1073/pnas.1419161111.}
        \item{\strong{SUPPA:}
            Alamancos GP, Pagès A, Trincado JL, Bellora N, Eyras E. Leveraging
            transcript quantification for fast computation of alternative
            splicing profiles. RNA. 2015;21(9):1521-1531.
            doi:10.1261/rna.051557.115.}
        \item{\strong{VAST-TOOLS:}
            Irimia M, Weatheritt RJ, Ellis J, et al. A highly conserved program
            of neuronal microexons is misregulated in autistic brains. Cell.
            2014;159(7):1511-1523. doi:10.1016/j.cell.2014.11.035.}
    }
}
\examples{
library(AnnotationHub)
hub <- AnnotationHub()

## Load the alternative splicing events for Human (hg38)
events <- query(hub, "alternativeSplicingEvents.hg38")[[1]]
}
\keyword{datasets}
