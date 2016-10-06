metadata <- data.frame(
    Title="Alternative Splicing Annotation for Homo sapiens (Human)",
    Description=paste("List of data frames containing alternative splicing",
                      "events per event type. Each splicing event is",
                      "characterised by its chromosome, strand, splice",
                      "junctions coordinates relevant to the type of event and",
                      "associated gene."),
    Species="Homo sapiens",
    TaxonomyId=9606,
    Genome="hg19",
    Maintainer="Nuno Agostinho <nunodanielagostinho@gmail.com>",
    RDataClass="list",
    DispatchClass="list",
    SourceUrl=paste("https://miso.readthedocs.io/en/fastmiso/annotation.html",
                    "http://rnaseq-mats.sourceforge.net/user_guide.htm",
                    "https://bitbucket.org/regulatorygenomicsupf/suppa",
                    "http://vastdb.crg.eu/libs/", sep=", "),
    SourceType="tab",
    SourceVersion=NA_character_,
    DataProvider=NA_character_,
    BiocVersion=package_version("3.4"),
    Coordinate_1_based=TRUE,
    RDataDateAdded=as.POSIXct(Sys.time()),
    Recipe=NA_character_,
    ResourceName="alternativeSplicingEvents.hg19.rda",
    Tags=I(list(c("Human", "Alternative", "Splicing", "Events", "Annotation", "hg19")))
)

write.csv(metadata, file="inst/extdata/metadata.csv", row.names=FALSE)
