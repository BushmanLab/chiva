library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

r1over <- function(seq) {
  paste0("TG",as.character(ShortRead::reverseComplement(Biostrings::DNAString(substr(seq,nchar(seq)-18,nchar(seq))))))
}

write.table(read.table(args[1], sep = ",", header = TRUE, stringsAsFactors = FALSE) %>%
              rowwise() %>%
              mutate(R1Over=r1over(LTR)) %>%
              ungroup() %>% as.data.frame,
            paste0(args[1],".new"), sep = ",", row.names = FALSE, quote = FALSE)
