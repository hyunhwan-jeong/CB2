fa_path <- "inst/extdata/toydata/small_sample.fasta"
ids <- readLines(fa_path)[c(TRUE, FALSE)]
seqs <- readLines(fa_path)[c(FALSE, TRUE)]

message("Performing sanity checks.")

valid_ids <- startsWith(ids, ">")
valid_ids[3] <- FALSE
if(!all(valid_ids)) {
  stop("Something is not correct with the guide IDs. Please check the line ",
       which(!valid_ids)[1] * 2 -1)
  
}

valid_seqs <- sapply(seqs, function(x) nchar(x) == nchar(seqs[1]))

if(!all(valid_seqs)) {
  stop("Something is not correct with the guide RNAs. Please check the line ",
       which(!valid_seqs)[1] * 2)
}

ids <- stringr::str_remove(ids, ">")  

data.frame(
  id = ids,
  sgRNA = seqs
) %>% readr::write_csv("inst/extdata/toydata/small_sample.csv")

