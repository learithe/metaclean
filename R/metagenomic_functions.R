#' RemoveTaxa Function
#'
#' Function to remove unwanted taxa from a physeq object
#'
#' @param physeq phyloseq object to prune
#' @param badTaxa list of taxa (OTU IDs) to remove
#' @param prune remove samples that have no remaining taxa (default FALSE)
#' @keywords badtaxa, removetaxa
#' @export
#' @examples
#' ps_noblanks <- RemoveTaxa(ps, blank_taxa, prune=TRUE)
RemoveTaxa = function(physeq, badTaxa, prune=FALSE){
  allTaxa = phyloseq::taxa_names(physeq)
  keepTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
  physeq <- phyloseq::prune_taxa(keepTaxa, physeq)
  if(prune){
    physeq <- phyloseq::prune_samples(phyloseq::sample_sums(physeq)>0, physeq)
  }
  return(physeq)
}


#' SubsetSamples Function
#'
#' Function to subset a phyloseq object to a specific subset of samples,
#' and remove taxa not found in those samples from the phyloseq object.
#' Especially useful for separating out blank samples
#'
#' @param physeq phyloseq object to subset
#' @param samplelist character vector of samples to keep
#' @keywords subsettaxa, filter, blanks
#' @export
#' @examples
#' ps_blanks <- SubsetSamples(ps, blank_samples)
SubsetSamples = function(physeq, samplelist){
  
  #sanity check, because prune_samples error msg is too vague
  is_subset = samplelist %in% rownames(otu_table(physeq))
  if( !any(is_subset) ){
    warning("Fatal Error in SubsetSamples(): no ID's in samplelist were found in the physeq sample IDs")
    return(NULL)
  }

  ps <- phyloseq::prune_samples(samplelist, physeq)  #subset to just the chosen samples
  ps <- phyloseq::filter_taxa(ps, function(x) sum(x) > 0, TRUE) #subset to just taxa that contain reads in these samples

  return(ps)
}


#' CreateRelabund Function
#'
#' Function to generate relative-abundance table versions of a phyloseq object
#' Returns a list of phyloseq objects:
#' ret$orig = the original phyloseq object submitted
#' ret$r = the relative-abundance version
#' ret$fr = the filtered relative abundance version (if specified)
#'
#' @param physeq phyloseq object to convert
#' @param filt minimum mean value to filter to (default is FALSE (no filtering))
#' @keywords relative abundance, phyloseq
#' @export
#' @examples
#' ps <- CreateRelabund(physeq_object)
#' ps$orig; ps$r
#' ps <- CreateRelabund(physeq_object, filt=1e-5)
#' ps$orig; ps$r; ps$fr
CreateRelabund = function(physeq, filt=FALSE){

  #generate relative abundance version
  psr  = phyloseq::transform_sample_counts(physeq, function(x) x / sum(x) )

  #replace NaNs generated where sum(x) = 0 with 0's
  tmp <- phyloseq::otu_table(psr)  
  tmp[is.na(tmp)] <- 0
  phyloseq::otu_table(psr) <- tmp

  #prepare list of ps objects to return
  ps <- list(orig=physeq, r=psr)

  #add filtered relative abundance table if desired
  if(filt){
    psfr = phyloseq::filter_taxa(psr, function(x) mean(x) > filt, TRUE)
    ps$fr = psfr
  }
  return(ps)
}



#' ReplaceSeqIDs Function
#'
#' Replace dada2 sequence strings with seqence ID values
#' and generate a .fasta file to match these to their sequences.
#' Use to make working with dada2 output a bit easier
#'
#' @param physeq phyloseq object to convert
#' @param name of fasta file to write to
#' @keywords phyloseq, dada2, otu_table
#' @export
#' @examples
#' ret <- ReplaceSeqIDs(ps, "seqids.fa")
#' ps2 <- ret$physeq
#' otu_map <- ret$otu_map
#' system("head -4 seqids.fa", intern=TRUE) #peek at the fasta file
ReplaceSeqIDs = function(physeq, fastafile){

  otus <- phyloseq::otu_table(physeq)
  taxa <- phyloseq::tax_table(physeq)

  #create "OTU ID's" & generage a id -> seq table
  seq <- colnames(otus)
  otu_id <- as.character( lapply(1:length(seq), function(x) paste("seq.",x,sep="")) )
  otu_mapping <- cbind(otu_id,seq)
  otu_mapping <- as.data.frame(otu_mapping)
  colnames(otu_mapping) <- c("SeqID","Seq")

  lookup_id = function(seq,omap){
    r <- as.character( omap[which(omap$Seq == seq), 1] )
  }

  #create a "mapping" otu_id -> sequence list in fasta format
  fasta <- as.data.frame(otu_mapping) %>%
    mutate(fasta = paste(">",SeqID,"\n",Seq,"\n",sep="") )
  fasta <- fasta$fasta
  cat(sapply(fasta, toString), file=fastafile, sep="")

  #modify the phyloseq object
  colnames(otus) <- otu_id
  rownames(taxa) <- unlist( lapply(rownames(taxa), function(x) lookup_id(x, otu_mapping)) )
  ps <- phyloseq(otu_table(otus, taxa_are_rows=FALSE),
                 sample_data(physeq),
                 tax_table(taxa))

  return( list(physeq=ps, otu_map=otu_mapping) )

}



#' ExtendPlastidTax Function
#'
#' The dada2 silva_v128 taxonomy database only classifies mitochondria to
#' Family and chloroplasts to Class level. This is problematic when these taxa 
#' have (deliberately) not been removed from the dataset, are present in 
#' significant amounts, and one wants to look at distributions at lower
#' taxonomy levels.
#' 
#' This function takes a phyloseq object with this problem, and fills in 
#' the taxonomy table for these plastids down the the lowest level
#' (either Genus or Species)
#'
#' @param physeq phyloseq object to correct
#' @keywords phyloseq, dada2, taxonomy
#' @export
#' @examples
#' ps <- ExtendPlastids(ps)
#' #to check that it worked:
#' df <- as.data.frame(tax_table(ps))
#' filter( df[which(df$Family == "Mitochondria"), ] )[1:3, ]
ExtendPlastidTax = function(ps){
  
  tax <- as.data.frame(tax_table(ps))
  names <- rownames(tax)
  tax[] <- lapply(tax, as.character)
  
  tax <- tax %>%
    mutate(Genus=replace(Genus, Family=="Mitochondria", "Mitochondria")) %>%
    mutate(Genus=replace(Genus, Class=="Chloroplast", "Chloroplast")) %>%
    mutate(Family=replace(Family, Class=="Chloroplast", "Chloroplast")) %>%
    mutate(Order=replace(Order, Class=="Chloroplast", "Chloroplast")) %>%
    as.data.frame()
  
  if("Species" %in% colnames(tax)){
    tax <- tax %>%
      mutate(Species=replace(Species, Family=="Mitochondria", "Mitochondria")) %>%
      mutate(Species=replace(Species, Class=="Chloroplast", "Chloroplast")) %>%
      as.data.frame()     
  }  
  
  rownames(tax) <- names
  tax <- as.matrix(tax)
  
  ps2 <- phyloseq( otu_table(ps), tax_table(tax), sample_data(ps))
  return(ps2)
}
