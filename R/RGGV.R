

#' Map variant in GGV
#'
#' @description
#' \code{ggv} takes either a string \code{rs} or a numeric physical position specified by
#' arguments \code{chr} and \code{pos}. If neither identifier is specified, a random variant
#' is mapped. Variants from either the HGDP or 1000 genomes can be retrieved by setting
#' \code{db = "hgdp"} (default) or \code{db = "1kgenomes"}.
#'
#' @details
#' Creates a cluster pie map from variant data in the GGV. Variants can be
#' retrieved by rs ID or physical position. If neither rs ID or physical
#' position is specified, a random SNP will be retrieved. The variant data
#' can be accessed from either the HGDP or 1kgenome databases. To produce a table
#' instead of a map, set \code{table = TRUE}.
#'
#' @param rs A SNP ID (string)
#' @param chr Chromosome (numeric)
#' @param pos Physical position (numeric)
#' @param db Database "hgdp" (default) or "1kgenomes"
#' @param output "map" or "table" (string)
#' @examples
#' ggv("rs1834640")
#' ggv()
#' ggv(, 6, 130099903, db = "1kgenomes")
#' @export ggv
ggv <- function(rs = NULL, chr = NULL, pos = NULL, db = "HGDP", output = "map") {
  table = FALSE
  if (output == "table")
    table = TRUE
  if (db == "1kgenomes" || db == "1000genomes")
    db = "1000genomes"
  api = 'http://popgen.uchicago.edu/ggv_api/freq_table?data="'
  if (!is.null(rs)) {
    api <- paste0(api, db, '_table\"&rsID=', rs)
  }
  else if (!is.null(chr) && !is.null(pos)) {
    api = paste0(api, db, '_table\"&chr=', chr, "&pos=", pos)
    rs <- paste0("chr",as.character(chr),":",as.character(pos))
  }
  else {
    message(paste0("searching for random SNP in ", db))
    api = paste0(api, db, '_table"&random_snp=True')
  }

  ggvjson <- try(RJSONIO::fromJSON(api), silent=TRUE)
  
    if (inherits(ggvjson, "try-error")) {
    p = testggv()
    if (p == 0)
      message(paste0("Variant ", rs, " not found in ", toupper(db)))
  }
  else if (table) {
    tab <- json2table(ggvjson)
  
    if (is.null(rs)) 
      rs = paste0("chr", ggvjson[[1]]$chrom_pos)
    message(paste0("found SNP ", rs))
    return(tab)
  }
  else {
    geojson <- json2geojson(ggvjson, rs, db)
    if (!is.null(geojson))
      #print(geojson)
      clusterpie(geojson)
  }
}

#' @export testSNPs
testSNPs <- function() {
  s <- select.list(c("rs1426654","rs1834640"))
  if (s == "rs1426654")
    ggv("rs1426654", db="1kgenomes")
  else if (s == "rs1834640")
    ggv("rs1834640")
}


json2table <- function(json) {
   freqtable <- NULL
   a1 <- json[[1]]$alleles[1]
  for (i in 1:length(json)) {
    freq <- json[[i]]$rawfreq
    pop <- json[[i]]$pop
    coord <- unlist(json[[i]]$pos)
    nobs <- json[[i]]$nobs
    freqtable <- rbind(freqtable, c(pop, coord[2], coord[1], nobs, round(freq[1],4)))
  }
  freqtable <- data.frame(freqtable, stringsAsFactors=FALSE)

  colnames(freqtable) <- c("Pop", "lat", "long" , "nobs", paste0("freq_",a1))
  freqtable[,c("lat","long","nobs",paste0("freq_",a1))] <- sapply( freqtable[,c("lat","long","nobs",paste0("freq_",a1))], as.numeric)
                                                                           
  return(freqtable)
}



json2geojson <- function(json, rs, db) {
  geom <- NULL
  pos <- json[[1]]$chrom_pos
  alleles <- json[[1]]$alleles
  #  alleles = c("A", "G")
  splpos <- unlist(strsplit(pos,":"))
  ancstate <- ancestralState(splpos[1], as.numeric(splpos[2]), db)
  message(paste0("SNP pos: chr", splpos[1], ":", splpos[2]))
  message(paste0("Ancestral state: ", ancstate))
  for (i in 1:length(json)) {
    nobs <- as.numeric(json[[i]]$nobs)
    xobs <- as.numeric(json[[i]]$xobs)
    pop <- json[[i]]$pop
    coord <- as.numeric(json[[i]]$pos)

    #print(coord)
    for (j in 1:nobs) {
      al <- alleles[2]
      if (xobs > 0){
        al <- alleles[1]
        xobs <- xobs - 1
      }
      g <- list(geometry=list(type="Point",coordinates=coord), type="Feature", properties=c(popname=pop,allele=al))
      geom[[length(geom)+1]] <- g
    }
  }

  ancder <- c("Ancestral", "Derived")
  if (ancstate == alleles[2])
    ancder <- rev(ancder)
  else if (ancstate != alleles[1] && ancstate != alleles[2])
    ancder <- c("State1", "State2")
  names(ancder) <- alleles
  pname <- c(name=toupper(db))
  if (is.null(rs)) rs <- pos
  geojson <- list(type="FeatureCollection", features=geom,
                  properties=list(fields=list(allele=list(lookup=ancder, name="Ancestral State"),
                                              popname=pname, attribution="XXXXX", description=rs)))
  #writeLines(toJSON(geojson), "jnk.geojson")
  return(RJSONIO::toJSON(geojson))

}

parseSNPID <- function(snpid) {
  snp <- getSNP(snpid)
  pos <- as.numeric(snp["chrom_start"])
  chr <- as.character(snp["chr_name"])
  alleles <- unlist(strsplit(as.character(snp["allele"]), "/"))
  #anc <- as.character(snp["allele_1"])
  anc <- ancestralState(chr, pos, db = "")
  anc.der <- alleles[which(alleles == anc):which(alleles != anc)]
  names(anc.der) <- c("anc", "der")
  message(paste0("snp ID: ", snpid))
  message(paste0("pos B37: ", chr, ":", pos))
  message(paste0("ancestral: ", anc.der["anc"], " derived: ", anc.der["der"]))
  snpinfo <- list("snpid" = snpid, "chr" = chr, "pos" = pos, "anc.der" = anc.der)
  return(snpinfo)
}

getSNP <- function(the.snp = "rs12913832", chr=NULL, pos = NULL, build = "B37") {
  message(paste0("looking up SNP ", the.snp, " in build ", build))
  ####CHANGE BACK TO GRCh=37
  snp.db <- biomaRt::useEnsembl(biomart = "snp", dataset = "hsapiens_snp", GRCh=37)
  nt.biomart <- biomaRt::getBM(c("refsnp_id", "allele", "chr_name", "chrom_start",
                        "chrom_strand", "allele_1"),
                      filters="snp_filter", values = the.snp, mart = snp.db)
  return(nt.biomart)
}


flipAllele <- function(allele) {
  flip <- "X"
  if (allele == "A" || allele == "a")
    flip <- "T"
  else if (allele == "T" || allele == "t")
    flip <- "A"
  else if (allele == "G" || allele == "g")
    flip <- "C"
  else if (allele == "C" || allele == "c")
    flip <- "G"
  return(flip)
}

testggv <- function() {
  server = "popgen.uchicago.edu"
  p=suppressWarnings(system2("ping",server,stderr=F,stdout=F))
  if (p == 0)
    message("GGV server is alive.")
  else {
    message("GGV server is not responding.")
    message(paste0("Please contact: ", maintainer("RGGV")))
  }
  return(p)
}

ancestralState <- function(chr, pos, db = "") {
  strchr <- paste0("chr",as.character(chr))
  site <- GenomicRanges::GRanges(seqnames = strchr, IRanges::IRanges(pos, pos))
  server <- "http://genome-data.cri.uchicago.edu/ahga/"
  dir <- "human_ancestor_GR37_e59/"
  if (db == 'hgdp')
    dir <- "Ancestral_hg18/"
  prefix <- "human_ancestor_"
  suffix <- ".fa"
  fafile=paste0(server, dir,prefix,as.character(chr),suffix)
  fa <- Rsamtools::FaFile(fafile)
  a <- Rsamtools::getSeq(fa,param=site)
  return (as.character(a))
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




####################################################

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

getPopLevels <- function(refresh = FALSE) {
  plevels <- attr(getPopLevels, "plevelAttr")
  if (is.null(plevels) || refresh) {
    res <- getXML()
    plevels <- length(XML::getNodeSet(XML::getNodeSet(res[[1]], "//Sample")[[1]], "//Population"))
  }
  attr(getPopLevels, "plevelAttr") <<- plevels
  return (plevels)
}

getSamples <- function(refresh = FALSE) {
  samp_pops <- attr(getSamples, "sampAttr")
  if (is.null(samp_pops) || refresh) {
    res <- getXML()
    plevels <- getPopLevels(refresh)
    samp_ids <- sapply(XML::getNodeSet(res[[1]], "//Label"), function(x){XML::xmlValue(x)})
    samp_pops <- NULL
    for (i in 1:plevels) {
      path <- paste("//Sample/Population[@level='",i,"']", sep = "")
      p <- sapply(XML::getNodeSet(res[[1]], path), function(x){XML::xmlValue(x)})
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
  pop_names <- sapply(c, function(x){XML::xmlValue(x[["PopLabel"]])})
  pop_lats <- as.numeric(sapply(c, function(x){XML::xmlValue(x[["Latitude"]])}))
  pop_longs <- as.numeric(sapply(c, function(x){XML::xmlValue(x[["Longitude"]])}))
  coords <- rbind(pop_lats, pop_longs)
  rownames(coords) <- c("lat", "long")
  colnames(coords) <- pop_names
  return(coords)
}

#' @export testVCF
testVCF <- function(snpid = NULL, chr, pos) {
  snp <- NULL
  if (!is.null(snpid))
    snp = parseSNPID(snpid)
  else {
      snpid = paste0(chr, ":", pos)
      snp$chr <- chr
      snp$pos <- pos
  }
  chrom = paste("chr", snp$chr, sep = "")
  tbx <- Rsamtools::TabixFile("http://genome-data.cri.uchicago.edu/ahga/archdata/allpops_auto_maf0.005geno0.1.vcf.gz")
  param <- GenomicRanges::GRanges(chrom, IRanges::IRanges(snp$pos, width=1))
  res <- Rsamtools::scanTabix(tbx, param=param)
  if (length(res[[1]]) == 0) {
    stop(paste0("SNP ", snpid, " not in database."), call. = FALSE)
    }
  #print(res)
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
  geojson <- buildGeojson(k, getCoords(), snpid)
  clusterpie(geojson)
  #print(geojson)
  #return(k)
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
    h <- Rsamtools::headerTabix(tbx)[["header"]]
    header <- read.csv(textConnection(h[length(h)]), sep = "\t", header = FALSE)
    attr(parseVCFHeader, "headerAttr") <<- header
    attr(parseVCFHeader, "tbxAttr") <<- tbx
  }
  return(header)
}





buildGeojson <- function(samp_data, coords, snpid) {
  ancder = c("Ancestral", "Derived")
  names(ancder) = c("A", "G")
  rs = snpid
  pname = "1kgenome_imputed_HGDP"
  geom <- NULL
  alleles <- parseAlleles()
  poplevlist <- NULL
  for (i in 1:getPopLevels())
    poplevlist[[length(poplevlist)+1]] <- paste0("pop_level", as.character(i))
  poplev <- select.list(poplevlist)
  for (i in 1:nrow(samp_data)) {
    al <- alleles
    id <- rownames(samp_data[i,])
    pop <- as.character(samp_data[i, poplev])
    lat <- coords["lat", pop]
    lon <- coords["long", pop]
    #c <- paste("[", lon, ",", lat, "]")
    c <- c(lon, lat)
    gen <- samp_data[i, "Genotype"]
    if (is.na(gen)) next
    if (gen == 2)
      al["REF"] <- al["ALT"]
    else if (gen == 0)
      al["ALT"] <- al["REF"]
    g <- list(geometry=list(type="Point",coordinates = c), type="Feature", properties=c(id=id,popname=pop,allele=as.character(al[1, 1])))
    geom[[length(geom)+1]] <- g
    g <- list(geometry=list(type="Point",coordinates = c), type="Feature", properties=c(id=id,popname=pop,allele=as.character(al[1, 2])))
    geom[[length(geom)+1]] <- g
  }
  geojson <- list(type="FeatureCollection", features=geom,
                  properties=list(fields=list(allele=list(lookup=ancder, name="Ancestral State"),
                                              popname=pname, attribution="XXXXX", description=rs)))
  return(RJSONIO::toJSON(geojson))
}

freqTable <- function(samp_data, level = 1) {
  if (level > getPopLevels()) {
    stop(paste0("The level exceeds the number of population levels (", getPopLevels(),")."), call. = FALSE)
  }
  lev <- paste("pop_level", level, sep = "")
  num <- tapply(samp_data$Genotype, samp_data[[lev]], sum, na.rm = TRUE)
  denom <- tapply(samp_data$Genotype, samp_data[[lev]], function(x){return(2*sum(!is.na(x)))})
  return (num/denom)
}

seasonAve <- function(var='air.sig995', season="winter", ext='min', hemi="south") {
  
  vardata=NCEP.gather(variable=var,level='surface',years.minmax=c(1960,1960),lat.southnorth=c(-70,77),lon.westeast=c(0,360),reanalysis2 = F,status.bar=F,months.minmax=c(1,12))
  lat=as.numeric(dimnames(vardata)[[1]])
  latS <- which(lat < 0)
  latN <- which(lat >= 0)
  var_S_D=NCEP.aggregate(vardata[latS,,],HOURS = FALSE, fxn=ext)
  var_S_MAve=NCEP.aggregate(var_S_D,DAYS=FALSE,HOURS=FALSE, fxn="mean")
  var_S_MAve=NCEP.aggregate(var_S_MAve,DAYS=FALSE,HOURS=FALSE, fxn="mean")
  var_S_MAve=NCEP.aggregate(var_S_MAve,YEARS=FALSE,DAYS=FALSE,HOURS=FALSE, fxn="mean")
  
  var_N_D=NCEP.aggregate(vardata[latN,,],HOURS = FALSE, fxn=ext)
  var_N_MAve=NCEP.aggregate(var_N_D,DAYS=FALSE,HOURS=FALSE, fxn="mean")
  var_N_MAve=NCEP.aggregate(var_N_MAve,DAYS=FALSE,HOURS=FALSE, fxn="mean")
  var_N_MAve=NCEP.aggregate(var_N_MAve,YEARS=FALSE,DAYS=FALSE,HOURS=FALSE, fxn="mean")
  mnthS=c(6,7,8)
  mnthN=c(1,2,12)
  if (season=="summer") {
    mnthS=c(1,2,12)
    mnthN=c(6,7,8)
  }
  var_S_MAve=as.data.frame.table(NCEP.aggregate(var_S_MAve[,,mnthS],YEARS=FALSE,MONTHS=F,DAYS=FALSE,HOURS = FALSE, fxn="mean"), stringsAsFactors=FALSE)[,c(1,2,4)]
  var_N_MAve=as.data.frame.table(NCEP.aggregate(var_N_MAve[,,mnthN],YEARS=FALSE,MONTHS=F,DAYS=FALSE,HOURS = FALSE, fxn="mean"), stringsAsFactors=FALSE)[,c(1,2,4)]
  mat=rbind(var_S_MAve,var_N_MAve)
  mat[,1]=as.numeric(mat[,1])
  mat[,2]=as.numeric(mat[,2])
  colnames(mat)=c("lat","lon","var")
    envmap(mat)
}

envmap <- function(mat) {
  
  nogrid <- theme(
    axis.text = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.title = element_blank()
  )

    world <-map_data('world') %>% data.table()
  world <- world[region!='Antarctica',]
  g <- ggplot(world, aes(long,lat)) + 
    geom_polygon(aes(group=group),fill="white",colour="black",size=0.1) +
    coord_equal() + 
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) + theme_bw() + nogrid 
  
  mat2 <- mat
  ml <- which(mat2[,"lon"]>180)
  mat2[ml,"lon"] <- mat2[ml,"lon"] - 360
    varname <- colnames(mat2)[3]
   v <- g + geom_tile(data=mat2,aes(x=lon, y=lat, fill=var),alpha=.75) + scale_fill_gradientn(colours=rev(rainbow(8)))
  return(v)
}
