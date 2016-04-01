library(Rsamtools)
library(XML)
library(biomaRt)
#' @export getAHGAXML


getXML <- function(db = "allpops_auto_maf0.005geno0.1", refresh = FALSE) {
  indiv_labelXMLroot <- attr(getXML, "indiv_labelXMLroot")
  pop_geoXMLroot <- attr(getXML, "pop_geoXMLroot")
  fileURL <- paste("http://genome-data.cri.uchicago.edu/ahga/", db, sep = "")
  if (refresh || is.null(indiv_labelXMLroot)) {
    indiv_labelXMLfile <- XML::xmlTreeParse(paste(fileURL, ".indiv_label.xml", sep = ""))
    pop_geoXMLfile <- XML::xmlTreeParse(paste(fileURL, ".pop_geo.xml", sep = ""))
    indiv_labelXMLroot <- XML::xmlRoot(indiv_labelXMLfile)
    pop_geoXMLroot <- XML::xmlRoot(pop_geoXMLfile)
  }
  attr(getXML, "indiv_labelXMLroot") <<- indiv_labelXMLroot
  attr(getXML, "pop_geoXMLroot") <<- pop_geoXMLroot
  return(list(indiv_labelXMLroot, pop_geoXMLroot))
}

getSNP <- function(the.snp = "rs12913832", chr = NULL, pos = NULL, build = "B37") {
  message(paste0("looking up SNP ", the.snp, " in build ", build))
  snp.db <- biomaRt::useEnsembl(biomart = "snp", dataset = "hsapiens_snp", GRCh = 37)
  nt.biomart <- biomaRt::getBM(c("refsnp_id", "allele", "chr_name", "chrom_start",
                                 "chrom_strand", "allele_1"),
                               filters="snp_filter", values = the.snp, mart = snp.db)
  return(nt.biomart)
}

parseSNPID <- function(snpid, build) {
  snp <- getSNP(snpid)
  print(snp)
  pos <- as.numeric(snp["chrom_start"])
  chr <- as.character(snp["chr_name"])
  alleles <- unlist(strsplit(as.character(snp["allele"]), "/"))
  #anc <- as.character(snp["allele_1"])
  
  anc <- ancestralState(chr, pos, build)
  #print("----------------")
  #print(anc)
  anc.der <- alleles
  
  #print(anc.der)
  if (anc %in% alleles)
    anc.der <- alleles[which(alleles == anc):which(alleles != anc)]
  names(anc.der) <- c("anc", "der")
  message(paste0("snp ID: ", snpid))
  message(paste0("pos B37: ", chr, ":", pos))
  message(paste0("ancestral: ", anc.der["anc"], " derived: ", anc.der["der"]))
  snpinfo <- list("snpid" = snpid, "chr" = chr, "pos" = pos, "anc.der" = anc.der)
  #print(snpinfo)
  return(snpinfo)
}

ancestralState <- function(chr, pos, build = "B37") {
  strchr <- paste0("chr",as.character(chr))
  site <- GenomicRanges::GRanges(seqnames = strchr, IRanges::IRanges(pos, pos))
  server <- "http://genome-data.cri.uchicago.edu/ahga/"
  dir <- "human_ancestor_GR37_e59/"
  if (build == 'B36')
    dir <- "Ancestral_hg18/"
  prefix <- "human_ancestor_"
  suffix <- ".fa"
  fafile=paste0(server, dir,prefix,as.character(chr),suffix)
  fa <- Rsamtools::FaFile(fafile)
  a <- Rsamtools::getSeq(fa,param=site)
  return (as.character(a))
}

getPopLevels <- function(refresh = FALSE) {
  plevels <- attr(getPopLevels, "plevelAttr")
  if (is.null(plevels) || refresh) {
    res <- getXML()
    plevels <- length(getNodeSet(getNodeSet(res[[1]], "//Sample")[[1]], "//Population"))
  }
  attr(getPopLevels, "plevelAttr") <<- plevels
  return (plevels)
}

getSamples <- function(refresh = FALSE) {
  samp_pops <- attr(getSamples, "sampAttr")
  if (is.null(samp_pops) || refresh) {
    res <- getXML()
    plevels <- getPopLevels(refresh)
    samp_ids <- sapply(getNodeSet(res[[1]], "//Label"), function(x){xmlValue(x)})
    samp_pops <- NULL
    for (i in 1:plevels) {
      path <- paste("//Sample/Population[@level='",i,"']", sep = "")
      p <- sapply(getNodeSet(res[[1]], path), function(x){xmlValue(x)})
      samp_pops <- cbind(samp_pops, p)
    }
    colnames(samp_pops) <- paste("pop_level", 1:plevels, sep = "") 
    rownames(samp_pops) <- samp_ids
    attr(getSamples, "sampAttr") <<- samp_pops
  }
  return(samp_pops)
}

getCoords <- function() {
  res <- getXML()
  c <- res[[2]]$children
  pop_names <- sapply(c, function(x){xmlValue(x[["PopLabel"]])})
  pop_lats <- as.numeric(sapply(c, function(x){xmlValue(x[["Latitude"]])}))
  pop_longs <- as.numeric(sapply(c, function(x){xmlValue(x[["Longitude"]])}))
  coords <- rbind(pop_lats, pop_longs)
  rownames(coords) <- c("lat", "long")
  colnames(coords) <- pop_names
  return(coords)
}



test <- function(snpid, chr, pos) {
  snp = parseSNPID(snpid, build="B37")
  chrom = paste("chr", snp$chr, sep = "")
  tbx <- TabixFile("http://genome-data.cri.uchicago.edu/ahga/archdata/allpops_auto_maf0.005geno0.1.vcf.gz")
  param <- GRanges(chrom, IRanges(snp$pos, width=1))
  res=scanTabix(tbx, param=param)
  r2=gsub("1/1","2",res)
  r2=gsub("1/0","1",r2)
  r2=gsub("0/1","1",r2)
  r2=gsub("0/0","0",r2)
  r2=gsub("./.","NA",r2)

  vcf=read.csv(textConnection(r2[[1]]),sep="\t",header=F)
  
  colnames(vcf) = unlist(parseVCFHeader(tbx))
  attr(getVCF, "vcfAttr") <<- vcf
  s=getSamples()
  genos=vcf[1,rownames(s)]

  k=data.frame(s, t(genos), stringsAsFactors = TRUE)
  cols = colnames(k)
  cols[length(cols)]= "Genotype"
  colnames(k)=cols
  geojson <- buildGeojson(k, getCoords())
  clusterpie(geojson)
  return(k)
}

getVCF <- function() {
  vcf <- attr(getVCF, "vcfAttr")
  return(vcf)
}

parseAlleles <- function(vcf = getVCF()) {
  alleles <- vcf[c("REF", "ALT")]
  print(alleles)
}

parseVCFHeader <- function(tbx) {
  header <- attr(parseVCFHeader, "headerAttr")
  tbxAttr <- attr(parseVCFHeader, "tbxAttr")
  if (is.null(header) || is.null(tbxAttr) || !identical(tbxAttr, tbx)) {
    h <- headerTabix(tbx)[["header"]]
    header <- read.csv(textConnection(h[length(h)]), sep = "\t", header = FALSE)
    attr(parseVCFHeader, "headerAttr") <<- header
    attr(parseVCFHeader, "tbxAttr") <<- tbx
  }
  return(header)
}



k=test("rs12913832")

buildGeojson <- function(samp_data, coords) {
  ancder = "DUMMY"
  rs = "DUMMY"
  pname = "1kgenome_imputed_HGDP"
  geom <- NULL
  alleles <- parseAlleles()
  for (i in 1:nrow(samp_data)) {
    al <- alleles
    pop <- as.character(samp_data[i, "pop_level1"])
    lat <- coords["lat", pop]
    lon <- coords["long", pop]
    c <- paste0("[", lon, ", ", lat, "]")
    gen <- samp_data[i, "Genotype"]
    if (is.na(gen)) next
    if (gen == 2) 
      al["REF"] <- al["ALT"]
    else if (gen == 0)
      al["ALT"] <- al["REF"]
    g <- list(geometry=list(type="Point",coordinates = c), type="Feature", properties=c(popname=pop,allele=as.character(al[1, "REF"])))
    geom[[length(geom)+1]] <- g
    g <- list(geometry=list(type="Point",coordinates = c), type="Feature", properties=c(popname=pop,allele=as.character(al[1, "ALT"])))
    geom[[length(geom)+1]] <- g
  }
    geojson <- list(type="FeatureCollection", features=geom,
                    properties=list(fields=list(allele=list(lookup=ancder, name="Ancestral State"),
                                                popname=pname, attribution="XXXXX", description=rs)))
    return(RJSONIO::toJSON(geojson))
}

freqTable <- function(samp_data, level = 1) {
  if (level > getPopLevels()) {
    warning(paste0("The level exceeds the number of population levels (", getPopLevels(),")."))
    return(-1)
  }
  lev <- paste("pop_level", level, sep = "") 
  num <- tapply(samp_data$Genotype, samp_data[[lev]], sum, na.rm = TRUE)
  denom <- tapply(samp_data$Genotype, samp_data[[lev]], function(x){return(2*sum(!is.na(x)))})
  return (num/denom)
}

getFreqTable <- function(vcf, level="Population") {
  samp_pops <- getSamples()
  coords <- getCoords()
  pop_names <- colnames(coords)
  alleles <- vcf[c("REF", "ALT")]
  freq <- matrix(nrow=length(pop_names), ncol = 2)
  for (i in 1:length(pop_names)) {
    n <- names(which(samp_pops[level,] == pop_names[i]))
    hets <- sum(vcf[n] == "1/0" || vcf[n] == "0/1")
    freq[i, 1] <- 2*sum(vcf[n] == "0/0") + hets
    freq[i, 2] <- 2*sum(vcf[n] == "1/1") + hets
  }
  rownames(freq) <- pop_names
  colnames(freq) <- unlist(alleles)
  return(freq)
}

getRelFreqTable <- function(vcf) {
  freq <- getFreqTable(vcf)
  sums <- apply(freq, 1, sum)
  freq <- freq/sums
  return(freq)
}

clusterpie <- function(geojson, server = NULL, width = NULL, height = NULL) {
  if (is.null(server)) {
    data <- geojson
  }
  
  # forward options using x
  x = list(
    geojsonfile = geojson,
    data = data,
    server = server
  )
  
  # create widget
  htmlwidgets::createWidget(
    name = 'clusterpie4widget',
    x,
    width = width,
    height = height,
    package = 'RGGV'
  )
}

# Widget output function for use in Shiny
#
#
clusterpie4widgetOutput <- function(outputId, width = '100%', height = '400px'){
  shinyWidgetOutput(outputId, 'clusterpie4widget', width, height, package = 'RGGV')
}

# Widget render function for use in Shiny
#
#
renderClusterpie4widget <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  shinyRenderWidget(expr, clusterpie4widgetOutput, env, quoted = TRUE)
}
