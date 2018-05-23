

#' Map variant in GGV
#'
#' @description
#' Retrieves GGV allele frequency data and creates interactive cluster map or table.
#'
#' @details
#' Downloads frequency data from one of the Geography of Genetic Variants (GGV) 
#' databases and creates a cluster pie map by default. Variants can be
#' retrieved by rs ID or physical position. If neither rs ID or physical
#' position is specified, a random SNP will be retrieved. Available variant 
#' databases are the HGDP, 1000 Genome, POPRES_Euro, and ExAC. To produce a table
#' instead of a map, set \code{table = TRUE}.
#' 
#' Genotype data can be imported directly by setting \code{vcf = URL}, where the 
#' URL is a string giving the .
#'
#' @param rs a string representing the SNP ID (NULL)
#' @param chr the chromosome number (NULL)
#' @param pos the physical position (NULL)
#' @param db a character string for the database ("HGDP" | "1kgenomes" | "POPRES" | "ExAC") (NULL)
#' @param output a character string for the output type ("map" | "table") ("map")
#' @param vcf a character string for the VCF filename (NULL)
#' @examples
#' ggv("rs1834640")
#' ggv()
#' ggv(, 6, 130099903, db = "1kgenomes")
#' @export ggv
ggv <- function(rs = NULL, chr = NULL, pos = NULL, db = NULL, output = "map", vcf = NULL) {
  table = FALSE
  if (output == "table")
    table = TRUE
 
  if (!is.null(vcf)) {
    getVCF(snpid=rs, chr, pos, vcffile=vcf, output=output)
  }
  
  else {
    ggvdb = c("HGDP", "1000genomes", "POPRES", "ExAC")
    if (db == "1kgenomes" || db == "1000Genomes")
      db = "1000genomes"
    if (is.null(db) || !(db %in% ggvdb)) {
      db <- select.list(ggvdb, title="GGV Databases")
    }
    
    api = 'http://popgen.uchicago.edu/ggv_api/freq_table?data="'
    searchstr <- '_table"&random_snp=True'
    if (!is.null(rs)) {
      searchstr <- paste0('_table\"&rsID=', rs)
    }
    else if (!is.null(chr) && !is.null(pos)) {
      searchstr <- paste0('_table\"&chr=', chr, "&pos=", pos)
      rs <- paste0("chr",as.character(chr),":",as.character(pos))
    }
    else {
      message(paste0("searching for random SNP in ", db))
    }
    api <- paste0(api, db, searchstr)
    ggvjson <- try(RJSONIO::fromJSON(api), silent=TRUE)
    if (inherits(ggvjson, "try-error")) {
      p <- testggv()
      if (p == 0)
        message(paste0("Variant ", rs, " not found in ", toupper(db)))
    }
    else if (table) {
      tab <- json2table(ggvjson)
  
      if (is.null(rs)) 
        rs <- paste0("chr", ggvjson[[1]]$chrom_pos)
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
   a2 <- json[[1]]$alleles[2]
  for (i in 1:length(json)) {
    freq <- json[[i]]$rawfreq
    pop <- json[[i]]$pop
    coord <- unlist(json[[i]]$pos)
    nobs <- json[[i]]$nobs
    freqtable <- rbind(freqtable, c(pop, coord[2], coord[1], nobs, round(freq[1],4), round(1-freq[1],4)))
  }
  freqtable <- data.frame(freqtable, stringsAsFactors=FALSE)

  colnames(freqtable) <- c("Pop", "lat", "long" , "nobs", paste0("freq_",a1),paste0("freq_",a2))
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
    package = 'RGGV',
    sizingPolicy = htmlwidgets::sizingPolicy(
                viewer.suppress = TRUE,
                browser.fill = TRUE

                )
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

getVCFMetaData <- function(vcffile = "allpops_auto_maf0.005geno0.1.vcf.gz", refresh = FALSE) {
  vcfpre <- strsplit(vcffile, "\\.")[[1]]
  vcfpre <- cat(vcfpre[1:length(vcfpre)-2], sep=".")
  indiv_labelXMLroot <- attr(getVCFMetaData, "indiv_labelXMLroot")
  pop_geoXMLroot <- attr(getVCFMetaData, "pop_geoXMLroot")
  fileURL <- paste0("http://genome-data.cri.uchicago.edu/ahga/", vcfpre)
  if (refresh || is.null(indiv_labelXMLroot)) {
    indiv_labelXMLfile <- XML::xmlTreeParse0(paste0(fileURL, ".indiv_label.xml"))
    pop_geoXMLfile <- XML::xmlTreeParse(paste0(fileURL, ".pop_geo.xml"))
    indiv_labelXMLroot <- XML::xmlRoot(indiv_labelXMLfile)
    pop_geoXMLroot <- XML::xmlRoot(pop_geoXMLfile)
  }
  attr(getVCFMetaData, "indiv_labelXMLroot") <<- indiv_labelXMLroot
  attr(getVCFMetaData, "pop_geoXMLroot") <<- pop_geoXMLroot
  return(list(indiv_labelXMLroot, pop_geoXMLroot))
}

getPopLevels <- function(vcffile, refresh = FALSE) {
  plevels <- attr(getPopLevels, "plevelAttr")
  if (is.null(plevels) || refresh) {
    res <- getVCFMetaData(vcffile)
    plevels <- length(XML::getNodeSet(XML::getNodeSet(res[[1]], "//Sample")[[1]], "//Population"))
  }
  attr(getPopLevels, "plevelAttr") <<- plevels
  return (plevels)
}

getSamples <- function(vcffile, refresh = FALSE) {
  samp_pops <- attr(getSamples, "sampAttr")
  if (is.null(samp_pops) || refresh) {
    res <- getVCFMetaData(vcffile)
    plevels <- getPopLevels(vcffile, refresh)
    samp_ids <- sapply(XML::getNodeSet(res[[1]], "//Label"), function(x){XML::xmlValue(x)})
    samp_pops <- NULL
    for (i in 1:plevels) {
      path <- paste0("//Sample/Population[@level='",i,"']")
      p <- sapply(XML::getNodeSet(res[[1]], path), function(x){XML::xmlValue(x)})
      samp_pops <- cbind(samp_pops, p)
    }
    colnames(samp_pops) <- paste0("pop_level", 1:plevels)
    rownames(samp_pops) <- samp_ids
    attr(getSamples, "sampAttr") <<- samp_pops
  }
  return(samp_pops)
}

getCoords <- function(vcffile) {
  res <- getVCFMetaData(vcffile)
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
testVCF <- function(snpid = NULL, chr, pos, vcfile = "http://genome-data.cri.uchicago.edu/ahga/archdata/allpops_auto_maf0.005geno0.1.vcf.gz") {
  getVCF(snpid, chr, pos, vcffile)
}

getVCF <- function(snpid = NULL, chr, pos, vcffile = "", output="map") {
  snp <- NULL
  if (!is.null(snpid))
    snp = parseSNPID(snpid)
  else {
    snpid = paste0(chr, ":", pos)
    snp$chr <- chr
    snp$pos <- pos
  }
  chrom <- paste0("chr", snp$chr)
  tbx <- Rsamtools::TabixFile(vcffile)
  param <- GenomicRanges::GRanges(chrom, IRanges::IRanges(snp$pos, width=1))
  res <- Rsamtools::scanTabix(tbx, param=param)
  if (length(res[[1]]) == 0) {
    stop(paste0("SNP ", snpid, " not in database."), call. = FALSE)
  }
  #print(res)  
  res=gsub("1/1|1\\|1","2",res)
  res=gsub("1/0|0/1|1\\|0|0\\|1","1",res)
  res=gsub("0/0|0\\|0","0",res)
  res=gsub("\\./\\.|\\.\\|\\.","NA",res)
 # print(res)
  vcf <- read.csv(textConnection(res[[1]]),sep="\t",header=F)
  colnames(vcf) <- unlist(parseVCFHeader(tbx))
 # print(vcf)
  attr(getVCF, "vcfAttr") <<- vcf
  s <- getSamples(vcffile)
  genos <- vcf[1,rownames(s)]
  k <- data.frame(s, t(genos), stringsAsFactors = TRUE)
  cols <- colnames(k)
  cols[length(cols)] <- "Genotype"
  colnames(k) <- cols
  if (output == "table") 
    geno2FrqTable(k, getCoords(vcffile), snpid)
  else {
    geojson <- buildGeojson(k, getCoords(vcffile), snpid)
    clusterpie(geojson)
  }
}

geno2FrqTable <- function(samp_data, coords, snpid, poplev = NULL) {
  if (is.null(poplev)) {
    poplevlist <- NULL
    for (i in 1:getPopLevels())
      poplevlist[[length(poplevlist)+1]] <- paste0("pop_level", as.character(i))
    poplev <- menu(poplevlist)
  }
  frq <- aggregate(Genotype ~ samp_data[,poplev], data=samp_data, FUN=frq)
  pops <- as.character(frq[,1])
  nobs <- aggregate(Genotype ~ samp_data[,poplev], data=samp_data, FUN=nobs)
  alleles <- parseAlleles()
  frqtable <- data.frame(cbind(pops, coords["lat", pops], coords["long", pops], nobs[,2],round(1-frq[,2], 4), round(frq[,2], 4)))
  colnames(frqtable) <- c("Pop", "lat", "long", "nobs", paste0("freq_", alleles[1,"REF"]), paste0("freq_", alleles[1,"ALT"]))
  rownames(frqtable) <- NULL
  return(frqtable)
}

frq <- function(x) {
  sum(x, na.rm=TRUE)/nobs(x)
}

nobs <- function(x) {
  2*sum(!is.na(x))
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





buildGeojson <- function(samp_data, coords, snpid, poplev = NULL) {
  ancder = c("Ancestral", "Derived")
  names(ancder) = c("A", "G")
  rs = snpid
  pname = "1kgenome_imputed_HGDP"
  geom <- NULL
  alleles <- parseAlleles()
  if (is.null(poplev)) {
    poplevlist <- NULL
    for (i in 1:getPopLevels())
      poplevlist[[length(poplevlist)+1]] <- paste0("pop_level", as.character(i))
    poplev <- menu(poplevlist)
  }
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
  lev <- paste0("pop_level", level)
  num <- tapply(samp_data$Genotype, samp_data[[lev]], sum, na.rm = TRUE)
  denom <- tapply(samp_data$Genotype, samp_data[[lev]], function(x){return(2*sum(!is.na(x)))})
  return (num/denom)
}

testseasonAve <- function(var='air.sig995', level="surface", season="winter", ext='min', hemi="south") {
  
  vardata=NCEP.gather(variable=var,level=level,years.minmax=c(1960,1960),lat.southnorth=c(-70,77.5),lon.westeast=c(0,360),reanalysis2 = F,status.bar=F,months.minmax=c(1,12))
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
  return(NCEP.aggregate(var_S_MAve[,,mnthS],YEARS=FALSE,MONTHS=F,DAYS=FALSE,HOURS = FALSE, fxn="mean"))
}

#' Map NCEP-NCAR environmental data
#'
#' @description
#' Creates a choropleth map displaying seasonal 
#' environmental data 
#'
#' @details
#' Creates a choropleth map displaying seasonal environmental data from the NCEP-NCAR database.
#'
#' @param envvar Enviornmental variable name (string)
#' @param season Season [winter|summer] (string)
#' @param level Data level (string)
#' @param ext Function for aggregating data ['min'|'max'|'mean'] (string)
#' @examples
#' envmap()
#' envmap(season="summer")
#' envmap(envvar='dswrf.surface',season='winter',level='gaussian',ext='max')
#' @export envmap
envmap <- function(grid=NULL, interior=FALSE, ncol=18, minlat=-61, scale=TRUE, col=NULL, title=NULL) {
  library(maps)
  if (is.null(col))
    col=sp::bpy.colors(ncol)
  else 
    ncol <- length(col)
  if (!is.null(grid)) {
    grid2 <- grid
    p <- which(grid$y >= minlat)
    grid2$y <- grid$y[p]
    grid2$z <- grid$z[,p]
    sgdf <- sgdf(grid2)
    zmin <- min(grid2$z)
    zmax <- max(grid2$z)
    what <- "both"
    if (!scale)
      what <- "image"
    if (is.null(title))
      title <- grid2$var
    br <- seq(zmin, zmax, (zmax-zmin)/ncol)
    sp::plot(sgdf, main=title, breaks=br, at=br, col=col, axis.pos=1, scale.shrink=.6, what=what)
    ymin=max(sp::bbox(sgdf)[2,"min"], -63)
    ymax=min(sp::bbox(sgdf)[2,"max"], 88)
    maps::map("world", add=TRUE, interior=interior,ylim=c(ymin,ymax))
  }
  else {
    maps::map("world", add=FALSE, fill=T, interior=interior, col="gray", lty=0, ylim=c(minlat,88))
  }
}

#' Add allele frequency pie charts to choropleth map
#'
#' @description
#' \code{addpie} Adds allele frequency pie charts for all populations in the  
#'  frequency dataframe.
#'
#' @details
#' Adds allele frequency pie charts for all populations in the frequency dataframe.
#'
#' @param f an allele frequency table generated by ggv(output="table") 
#' @param radius a number for the radius of the pie charts. (4)
#' @param ggvscaling a logical indicating whether to use pie chart scaling (TRUE)
#' @param cols a character vector giving the color of the pie chart slices (c("black", "orange"))
#' @alpha alpha a number indicating alpha transparency of pie chart slices (1)
#' @examples
#' freq_table <- ggv()
#' grid <- seasonAve()
#' envmap(grid)
#' addpie(freq_table)
#' @export addpie
addpie <- function(f, radius=4, ggvscaling = TRUE, cols=c("black", "orange"), alpha=1) {
  i = 1
  f2 <- f[order(f$long),]
  prevlong = -9999
  prevlat = -9999
  altcol= adjustcolor(cols[2], alpha.f=alpha)
  bkcol = adjustcolor(cols[1], alpha.f=alpha)
  if (ggvscaling && (max(f2[,5]) > 0.01 && max(f2[,5]) <= 0.10 || min(f2[,5]) < 0.99 && min(f2[,5]) >= 0.9))
    altcol = "green"
  else if (ggvscaling && (max(f2[,5]) > 0.00 && max(f2[,5]) <= 0.010 || min(f2[,5]) < 1.00 && min(f2[,5]) >= 0.99))
    altcol = "red"
  while (i <= nrow(f)) {
    currlong = f2[i,"long"]
    currlat = f2[i,"lat"]
        if (sqrt((currlong - prevlong)^2 + (currlat - prevlat)^2) < 2*radius)
      currlong <- prevlong + 2*radius
    col1 = altcol
    col2 = bkcol
    freq <- f2[i,5]
    if (altcol == "green") { 
      if (f2[i,5] <= 0.10) {
        col1 = altcol
        col2 = bkcol
        freq = 10 * f2[i,5]
      }
      else if (1-f2[i,5] <= 0.10) {
        col1 = bkcol
        col2 = altcol
        freq = 10 * (1 - f2[i,5])
      }
    }
    else if (altcol == "red") {
      if (f2[1,5] <= 0.01) {
        col1 = altcol
        col2 = bkcol
        freq <- 100 * f2[i,5]
      }
      else if (1-f2[i,5] <= 0.01) {
        col1 = bkcol
        col2 = altcol
        freq <- 100 * (1-f2[i,5])
      }
    }
  
    mapplots::add.pie(z=c(freq,1 - freq), x=currlong, y=f2[i,"lat"],radius=radius,col=c(col1,col2),labels="",init.angle=45,border="white")
    i <- i+1
    prevlong <- currlong
    prevlat <- currlat
  }
}

#' @export seasonAve
seasonAve <- function(var='air.sig995', season="winter", level='surface', ext='min', savedata = TRUE) {
  
  datafile=paste0(var,"_",level,".Rdata")
  if (file.exists(datafile)) {
    load(datafile)
  }
  else {
    vardata=RNCEP::NCEP.gather(variable=var,level=level,years.minmax=c(1960,1960),lat.southnorth=c(-63,88),lon.westeast=c(-179,180),reanalysis2 = F,status.bar=F,months.minmax=c(1,12))
     if (savedata) save(vardata, file=datafile) 
    }
    lat=as.numeric(dimnames(vardata)[[1]])
    latS <- which(lat < 0)
    latN <- which(lat >= 0)
    var_S_D <- RNCEP::NCEP.aggregate(vardata[latS,,],HOURS = FALSE, fxn=ext)
    var_S_MAve <- RNCEP::NCEP.aggregate(var_S_D,DAYS=FALSE,HOURS=FALSE, fxn="mean")
    var_S_MAve <- RNCEP::NCEP.aggregate(var_S_MAve,DAYS=FALSE,HOURS=FALSE, fxn="mean")
    var_S_MAve <- RNCEP::NCEP.aggregate(var_S_MAve,YEARS=FALSE,DAYS=FALSE,HOURS=FALSE, fxn="mean")
  
    var_N_D <- RNCEP::NCEP.aggregate(vardata[latN,,],HOURS = FALSE, fxn=ext)
    var_N_MAve <- RNCEP::NCEP.aggregate(var_N_D,DAYS=FALSE,HOURS=FALSE, fxn="mean")
    var_N_MAve <- RNCEP::NCEP.aggregate(var_N_MAve,DAYS=FALSE,HOURS=FALSE, fxn="mean")
    var_N_MAve <- RNCEP::NCEP.aggregate(var_N_MAve,YEARS=FALSE,DAYS=FALSE,HOURS=FALSE, fxn="mean")
      
  mnthS=c(6,7,8)
  mnthN=c(1,2,12)
  if (season=="summer") {
    mnthS=c(1,2,12)
    mnthN=c(6,7,8)
  }
  var_S_MAve <- RNCEP::NCEP.array2df(RNCEP::NCEP.aggregate(var_S_MAve[,,mnthS],YEARS=FALSE,MONTHS=F,DAYS=FALSE,HOURS = FALSE, fxn="mean"))[,2:4]
  var_N_MAve <- RNCEP::NCEP.array2df(RNCEP::NCEP.aggregate(var_N_MAve[,,mnthN],YEARS=FALSE,MONTHS=F,DAYS=FALSE,HOURS = FALSE, fxn="mean"))[,2:4]
  mat=rbind(var_S_MAve,var_N_MAve)

  #setorderv(mat,c("latitude","longitude"),order=c(-1,1))
  mat <- mat[order(-mat$latitude, mat$longitude),]
  mat <- unique(mat)

  longgridres <- 2*length(unique(mat$longitude))
  latgridres <- 2*length(unique(mat$latitude))
  t=tgp::interp.loess(mat[,"longitude"],mat[,"latitude"],mat[,"variable1"],gridlen=c(longgridres,latgridres),span=.001)
  t=c(t,"var"=paste0(var,"_",season,"_",ext))
  return(t)
}

sgdf <- function(t) {
  sg=sp::SpatialGrid(sp::GridTopology(c(min(t$x),min(t$y)), c(t$x[2]-t$x[1],t$y[2]-t$y[1]), c(length(t$x), length(t$y))),proj4string=sp::CRS("+proj=longlat"))
  m=t$z[,rev(seq_len(ncol(t$z)))]
  dim(m)<-c(length(m),1)
  sgdf=sp::SpatialGridDataFrame(grid=sg,data=as.data.frame(m))

 # sg=SpatialGrid(points2grid(mat),proj4string=CRS("+proj=longlat"))
  #  sgdf=SpatialGridDataFrame(grid=sg,data=as.data.frame(mat$variable1))
  return(sgdf)#' 
}


#' Interpolate environmental variables to population coordinates
#'
#' @description
#' \code{interpvar} Interpolates the environmental variable to the    
#' population coordinates
#'
#' @details
#' Interpolates the environmental variable to the population coordinates.
#'
#' @param t Environmental data grid. (dataframe)
#' @param f Frequency table. (dataframe)
#' @examples
#' 
#' 
#' 
#' @export interpvar
interpvar <- function(t,f) {
  longv <- sapply(f$long,nearestLong,t=t)
  latv <- sapply(f$lat, nearestLat,t=t)
  variable1 <- t$z[cbind(longv,latv)]
  f <- cbind(f, variable1)
  names(f)[length(f)] <- t$var
  return(f)
}

nearestLat <- function(t, lat) {
  return(which.min(abs(t$y-lat)))
}

nearestLong <- function(t, long) {
  return(which.min(abs(t$x-long)))
}



