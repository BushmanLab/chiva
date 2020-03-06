# This script is meant to concatenate sequence files within a directory together
# into a single file. 
# 
# Usage Rscript concat_seq_files.R {input.directory} {output.fasta}
# 
# Files within the input directory should be in fasta or fasta.gz formats.
# 

options(stringsAsFactors = FALSE, scipen = 99, width = 999)


code_dir <- dirname(sub(
  pattern = "--file=", 
  replacement = "", 
  x = grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
))

desc <- yaml::yaml.load_file(
  file.path(code_dir, "descriptions/trim.desc.yml")
)

#' Set up and gather command line arguments
parser <- argparse::ArgumentParser(
  description = desc$program_short_description,
  usage = "Rscript trim.R <seqFile> [-h/--help, -v/--version] [optional args]"
)


parser$add_argument(
  "-o", "--output", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--read2", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--genome", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--bowtie", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--build", nargs = 1, type = "character", help = desc$output
)

parser$add_argument(
  "--index", nargs = 1, type = "character", help = desc$output
)

args <- parser$parse_args(commandArgs(trailingOnly = TRUE))

# Perform initial alignment & variant calling
i <- 1
nextFasta <- paste0(args$output,".",i,".fasta")
nextIndex <- substr(nextFasta,1,nchar(nextFasta)-6)
tempBam <- paste0(args$output,".",i,".bam")
tempVcf <- paste0(args$output,".",i,".vcf.gz")

bowtieCommand <- paste(args$bowtie,"-x",args$index,
                       "-U",args$read2, "-k 2 --local -N 1",
                       "| /home/kevin/miniconda3/envs/chiva/bin/samtools sort",
                       "| /home/kevin/miniconda3/envs/chiva/bin/samtools view -Sb",
                       "> ", tempBam,
                       sep = " ")
freebayesCommand <- paste("freebayes","-f",args$genome,tempBam,
                          "| /home/kevin/miniconda3/envs/chiva/bin/bcftools view -i '%QUAL>=20'",
                          "| bgzip",
                          "> ",tempVcf, sep=" ")
tabixCommand <- paste("tabix",tempVcf, sep=" ")
bcftoolsCommand <- paste("cat",args$genome,
                         "| /home/kevin/miniconda3/envs/chiva/bin/bcftools consensus",
                         tempVcf,">",nextFasta,
                         sep=" ")

print(bowtieCommand)
system(bowtieCommand)
print(freebayesCommand)
system(freebayesCommand)
print(tabixCommand)
system(tabixCommand)
print(bcftoolsCommand)
system(bcftoolsCommand)

diffCommand <- paste("diff",args$output,args$genome)
diffTest <- system(diffCommand, intern=FALSE, ignore.stdout = TRUE)
lastFasta <- nextFasta
lastIndex <- nextIndex


while(diffTest>0 & i < 9) {
  i <- i + 1
  print("genomes different, continuing loop")
  
  nextFasta <- paste0(args$output,".",i,".fasta")
  nextIndex <- substr(nextFasta,1,nchar(nextFasta)-6)
  tempBam <- paste0(args$output,".",i,".bam")
  tempVcf <- paste0(args$output,".",i,".vcf.gz")
  
  
  bowtieIndexCommand <- paste(args$build,"-q",lastFasta,lastIndex, sep=" ")
  bowtieCommand <- paste(args$bowtie,"-x",lastIndex,
                         "-U",args$read2, "-k 2 --local -N 1",
                         "| /home/kevin/miniconda3/envs/chiva/bin/samtools sort",
                         "| /home/kevin/miniconda3/envs/chiva/bin/samtools view -Sb",
                         "> ", tempBam,
                         sep = " ")
  freebayesCommand <- paste("freebayes","-f",lastFasta,tempBam,
                            "| /home/kevin/miniconda3/envs/chiva/bin/bcftools view -i '%QUAL>=20'",
                            "| bgzip",
                            "> ",tempVcf, sep=" ")
  tabixCommand <- paste("tabix",tempVcf, sep=" ")
  bcftoolsCommand <- paste("cat",lastFasta,
                           "| /home/kevin/miniconda3/envs/chiva/bin/bcftools consensus",
                           tempVcf,">",nextFasta,
                           sep=" ")

  print(bowtieIndexCommand)
  system(bowtieIndexCommand)
  print(bowtieCommand)
  system(bowtieCommand)
  print(freebayesCommand)
  system(freebayesCommand)
  print(tabixCommand)
  system(tabixCommand)
  print(bcftoolsCommand)
  system(bcftoolsCommand)
  
  diffCommand <- paste("diff",lastFasta,nextFasta)
  print(diffCommand)
  diffTest <- system(diffCommand, intern=FALSE, ignore.stdout = TRUE)

  lastFasta <- nextFasta
  lastIndex <- nextIndex
  nextFasta <- paste0(args$output,".",i,".fasta")
  nextIndex <- gsub(".fasta","",nextFasta)
}
print(paste0("Genome converged after ",i," iterations"))
replaceFastaCommand <- paste("cp",nextFasta,args$output, sep=" ")
print(replaceFastaCommand)
system(replaceFastaCommand)
allTempFiles <- paste0(args$output,".*")
removeTempFilesCommand <- paste("rm",allTempFiles, sep=" ")
print(removeTempFilesCommand)
system(removeTempFilesCommand)
bowtieIndex <- substr(args$output,1,nchar(nextFasta)-6)
bowtieIndexCommand <- paste(args$build,"-q",args$output,bowtieIndex, sep=" ")
print(bowtieIndexCommand)
system(bowtieIndexCommand)
samtoolsIndexCommand <- paste("samtools","faidx",args$output, sep=" ")
print(samtoolsIndexCommand)
system(samtoolsIndexCommand)


