library(annotatr)

args = commandArgs(trailingOnly=TRUE)
bed_path <- args[1]
tmp <- read.csv(bed_path, header = FALSE, sep = "\t")

# # bed_path <- $1
df <-read_annotations(bed_path)
df$name <- tmp$V4
annots = c('hg38_cpgs')
annotations = build_annotations(genome = 'hg38', annotations = annots)
dm_annotated = annotate_regions(
    regions = df,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
df_dm_annotated = data.frame(dm_annotated)
df_dm_annotated <- df_dm_annotated[df_dm_annotated$annot.type == 'hg38_cpg_islands', ]
# # See the GRanges column of dm_annotaed expanded
write.csv(df_dm_annotated, args[2])