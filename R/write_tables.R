#' CreateTables Function
#'
#' Function to convert phyloseq objects into a set of
#' human-readable OTU tables analagous to QIIME output
#' will also append sequences to master table if requested
#'
#' @param physeq phyloseq object to convert
#' @param metadata metadata table to merge
#' @param cnames char vector of metadata column names to include in output
#' @param otu_map optional: df["Seq","SeqID"] pairing sequence ID to its sequence
#' @param key name of key variable (default = "SampleID")
#' @keywords qiime, phyloseq, dada2, otu_table
#' @export
#' @examples
#' cnames = c("SampleID","SampleType","TimePoint","Concentration","Batch","Description")
#' ret <- WriteTables(ps$orig, Metadata, cnames, key="SampleID")
#' write.csv(ret$master, "master_table.csv")
#' write.csv(ret$sum, "summary_table.csv")
#' write.csv(ret$otu, "otu_table.csv")
#' write.csv(ret$tax, "tax_table.csv")
#' cat(sapply(ret$fasta, toString), file="otu_mapping.fasta", sep="")
CreateTables = function (physeq, metadata, cnames, key = "SampleID", otu_map = NA)
{
  otus <- phyloseq::otu_table(physeq)
  metadata <- suppressWarnings(as.matrix(metadata))
  nsamples <- nrow(phyloseq::sample_data(physeq)[, key])
  if (nsamples != nrow(otus)) {
    otus <- t(otus)
  }
  sum_table <- as.data.frame(otus)[, 0]
  sum_table <- as.data.frame(cbind(sum_table, row.names(sum_table)))
  names(sum_table) <- c(key)
  df <- as.data.frame(t(as.data.frame(otus)))
  counts <- t(summarise_all(df, funs(sum)))
  counts <- as.data.frame(cbind(counts, row.names(counts)))
  names(counts) <- c("ClassifiedReads", key)
  sum_table <- merge(sum_table, counts, by = key)
  meta <- metadata[, cnames]
  sum_table <- merge(sum_table, meta, by = key)
  classic_otu_table <- as.data.frame(otus)
  classic_tax_table <- as.data.frame(phyloseq::tax_table(physeq))
  otu_tax_table <- cbind(t(classic_otu_table), classic_tax_table)
  if (!is.na(otu_map)) {
    lookup_id = function(seqid, omap) {
      r <- as.character(omap[which(omap$SeqID == seqid),
                             2])
    }
    print(rownames(otu_tax_table))
    seqs <- unlist(lapply(rownames(otu_tax_table), function(x) lookup_id(x,
                                                                         otu_map)))
    otu_tax_table <- cbind(otu_tax_table, data.frame(Sequence = seqs))
  }
  ret <- list(sum = sum_table, otu = classic_otu_table, tax = classic_tax_table,
              master = otu_tax_table)
  return(ret)
}




#' CreateTablesSeqIDs Function
#'
#' Function to convert phyloseq objects of DADA2 output
#' into a set of human-readable OTU tables analagous to QIIME output,
#' including creating sequence ID's & a fasta mapping for the taxa found
#'
#' @param physeq phyloseq object to convert
#' @param metadata metadata table to merge
#' @param cnames char vector of metadata column names to include in output
#' @param key name of key variable (default = "SampleID")
#' @keywords qiime, phyloseq, dada2, otu_table
#' @export
#' @examples
#' cnames = c("SampleID","SampleType","TimePoint","Description")
#' ret <- WriteTablesSeqIDs(ps$orig, sample_data(ps$orig), cnames, key="SampleID")
#' write.csv(ret$master, "master_table.csv")
#' write.csv(ret$sum, "summary_table.csv")
#' write.csv(ret$otu, "otu_table.csv")
#' write.csv(ret$tax, "tax_table.csv")
#' cat(sapply(ret$fasta, toString), file="otu_mapping.fasta", sep="")
CreateTablesSeqIDs = function(physeq, metadata, cnames, key="SampleID"){

  otus <- phyloseq::otu_table(physeq)
  metadata <- suppressWarnings(as.matrix(metadata)) #remove any lingering sample_data() object stuff

  #check whether the matrix needs to be transposed
  nsamples <- nrow( phyloseq::sample_data(physeq)[ ,key] )
  if(nsamples != nrow(otus)){
    otus <- t(otus)
  }

  #initialise the summary table
  sum_table <- as.data.frame(otus)[,0]
  sum_table <- as.data.frame( cbind(sum_table,row.names(sum_table)) )
  names(sum_table) <- c(key)

  #add column of total classified read counts for each sample
  df <- as.data.frame( t( as.data.frame(otus) ) )
  counts <- t( summarise_all(df, funs(sum)) )
  counts <- as.data.frame( cbind(counts,row.names(counts)) )
  names(counts) <- c("ClassifiedReads",key)
  sum_table <- merge(sum_table, counts, by=key)

  #merge the two together & add the metadata
  meta <- metadata[ ,cnames]
  sum_table <- merge(sum_table, meta, by=key)

  #create "OTU ID's" & generage a id -> seq table
  seq <- colnames(otus)
  otu_id <- as.character( lapply(1:length(seq), function(x) paste("seq.",x,sep="")) )
  otu_mapping <- cbind(otu_id,seq)

  #create a "classic" otu table with the new seq ids
  #for some reason I can't get rownames_to_column() to work here
  classic_otu_table <- as.data.frame(otus)
  colnames(classic_otu_table) <- otu_id
  #print(sort(rownames(classic_otu_table)))
  #print(length(rownames(classic_otu_table)))
  #print(length(unique(rownames(classic_otu_table))))
  #print(nrow(classic_otu_table))
  #classic_otu_table <- as.data.frame(classic_otu_table)
  #classic_otu_table <- tibble::rownames_to_column(as_tibble(classic_otu_table), var=key)

  #create a "classic" taxonomy table with the new seq ids
  classic_tax_table <- as.data.frame(phyloseq::tax_table(physeq))
  rownames(classic_tax_table) <- otu_id
  classic_tax_table <- tibble::rownames_to_column(classic_tax_table, "SeqID")

  #create a "mapping" otu_id -> sequence list in fasta format
  fasta <- as.data.frame(otu_mapping) %>%
    mutate(fasta = paste(">",otu_id,"\n",seq,"\n",sep="") )
  fasta <- fasta$fasta

  #create a master OTU table with counts, taxonomy, and sequence all together
  otu_tax_table <- cbind(t(classic_otu_table),classic_tax_table)
  otu_tax_table <- cbind(otu_tax_table,otu_mapping[ ,"seq"])
  colnames(otu_tax_table)[length(colnames(otu_tax_table))] = "Sequence"

  ret <- list(sum=sum_table, otu=classic_otu_table, tax=classic_tax_table,
              fasta=fasta, master=otu_tax_table)

  return(ret)

}
