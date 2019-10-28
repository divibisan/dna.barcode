#' Make table of DNA sequence for graphing
#' @param seqs DNAStringSet containing sequences to align
#' @param trim logical, should leading or trailing gaps be trimmed?
#'   All sequences will be trimmed equally
#' @importFrom magrittr %>%
#' @return A data.frame in long-format, with 1 row per sequence and base
make_dna_vectors <- function(seqs, trim = FALSE) {
  align <- msa(seqs, method = 'Muscle')
  align_list = msaConvert(align, type = "bios2mds::align")

  if (trim) {
    trim_end <- max(sapply(align_list, function(s) sum(rev(cumsum(rev(!(s == '-')))) == 0)))
    trim_start <- max(sapply(align_list, function(s) sum(cumsum(!(s == '-')) == 0)))
    align_list <- lapply(align_list, function(s) s[trim_start:(length(s) - trim_end)])
  }

  lapply(align_list, function(s) {
    data.frame(base = s) %>%
      rowid_to_column()
  }) %>% `names<-`(names(align_list)) %>%
    bind_rows(.id = 'species')
}

#' Make DNA Barcode Plot
#' @param sequences DNAStringSet cnotaining sequences to be aligned and plotted.
#'  Should be loaded from a FASTA file with \code{sequences <- Biostrings::readDNAStringSet('seqs.fasta')}
#' @param species_sub character vector, the names of the sequences to align and plot.
#'  Must match `rownames(sequences)`
#' @param trim logical, should leading or trailing gaps be trimmed?
#'  All sequences will be trimmed equally
#' @return A ggplot containing a DNA barcode plot
#' @importFrom magrittr %>%
#' @export
make_barcode_plot <- function(sequences, species_sub, trim = FALSE) {
  dna_vectors <- make_dna_vectors(sequences[species_sub,], trim) %>%
    mutate(species = factor(species, levels = rev(species_sub)))

  ggplot(dna_vectors, aes(x = rowid, y = species, fill = base)) +
    geom_tile(height = 0.9) +
    scale_fill_manual(name = 'Base',
                      values = c('-' = 'grey',
                                 'A' = rgb(20,42,250, max = 255),
                                 'C' = rgb(104, 243, 18, max = 255),
                                 'T' = rgb(230, 46, 37, max = 255),
                                 'G' = rgb(0,0,0, max = 255))) +
    scale_x_continuous(name = 'Sequence Position',
                       expand = c(0,0),
                       breaks = seq(0,max(dna_vectors$rowid), by = 100)) +
    scale_y_discrete(name = NULL, expand = c(0,0),
                     labels = function(x) gsub('_', ' ', x)) +
    theme_minimal() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
}
