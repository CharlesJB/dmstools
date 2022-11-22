#' Produce list of expected sequences
#'
#' When doing a miniscreen, there is only a small subset of all possible
#' sequences that are expected. For instance, we do not expect more than 1
#' codon to be mutated in a typical DMS library.
#'
#' This function will list all those possible sequences compared to the
#' provided reference sequence using the 63 alternative codons.
#'
#' @param reference The WT sequence of the region of interest. Must be a
#' character string containing only A, C, G or T and must have a length that is
#' a multiple of 3.
#' @param codon_offset The offset of the region of interest compared to the
#' real start site. For example, if the region of interest starts at the 4th
#' codon of the gene, the \code{codon_offset} would be 3. Default: 0
#' @param no.init.codon Should the first codon considered not to be an
#' initialization codon? Default: \code{TRUE}
#'
#' @examples
#' expected_sequences <- produce_expected_sequences("ACGTGCAAC")
#'
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings translate
#' @importFrom tibble tibble
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom dplyr mutate
#' @importFrom dplyr if_else
#' @importFrom purrr map_dfr
#' @importFrom stringr str_sub
#'
#' @export
produce_expected_sequences <- function(reference, codon_offset = 0,
                                       no.init.codon = TRUE) {
    wt_df <- tibble::tibble(expected_sequence = reference,
                    codon = NA,
                    wt_nucl = NA,
                    variant_nucl = NA,
                    wt_aa = NA,
                    variant_aa = NA,
                    is_wt = TRUE,
                    is_synonymous = FALSE)
    expected_seqs <- purrr::map_dfr(1:(nchar(reference) / 3),
                                    produce_mutant_seq,
                                    reference,
                                    codon_offset,
                                    no.init.codon)
    rbind(expected_seqs, wt_df)
}

# Produce all possible mutant sequence 
#
# This function will produce all the possible mutant sequences from a wt
# sequence of reference at a specific codon position.
#
# @param pos The position where the mutant sequence should be produced.
# @param wt_seq The wt sequence to use as a reference.
# @param offset The offset to add to the relative codon position. For instance,
# the first codon in the sequence could be the tenth codon in the gene. In that
# case, the offset would be 1.
# @param no.init.codon Should the first codon considered not to be an
# initialization codon? Default: \code{TRUE}
#
# @return A \code{
produce_mutant_seq <- function(pos, wt_seq, offset, no.init.codon) {
    print(pos)
    nucl <- c("A", "C", "G", "T")
    x <- expand.grid(nucl, nucl, nucl)
    all_codons <- paste0(x$Var1, x$Var2, x$Var3)

    codons <- purrr::map(1:(nchar(wt_seq) / 3), ~ get_codon(wt_seq, .x))

    if (pos == 1) {
        before <- ""
    } else {
        before <- paste(codons[1:(pos-1)], collapse = "")
    }
    if (pos == length(codons)) {
        after <- ""
    } else {
        after <- paste(codons[(pos+1):length(codons)], collapse = "")
    }
    current_wt_seq <- codons[[pos]]

    non_wt_codons  <- all_codons[all_codons != codons[pos]]
    seqs <- paste0(before, non_wt_codons, after)
    variant_seqs <- non_wt_codons

    res <- tibble::tibble(expected_sequence = seqs,
                          codon = pos + offset,
                          wt_nucl = current_wt_seq,
                          variant_nucl = variant_seqs,
                          wt_aa = translate_sequence(current_wt_seq, no.init.codon),
                          variant_aa = translate_sequence(variant_seqs, no.init.codon)) %>%
        dplyr::mutate(is_wt = wt_nucl == variant_nucl,
               is_synonymous = dplyr::if_else(!is_wt, wt_aa == variant_aa, FALSE))
    res
}

# Extract the nucleotide corresponding to a specific codon
#
# This function will extract the nucleotide sequence corresponding to a condon
# of interest. For instance, the codon 1 corresponds to the nucleotide 1 to 3
# in the sequence.
#
# @param seq The nucleotide sequence from which the codon will be extracted
# @param pos The codon to extract
#
# @return The nucleotide sequence of the codon
get_codon <- function(seq, pos) {
    start <- ((pos - 1) * 3) + 1
    end <- start + 2
    stringr::str_sub(seq, start, end)
}

# Translate a sequence
#
# This function will take as an input a \code{character} and will output the
# translated sequence also in \code{character} format.
#
# @param seq The sequence to translate, in \code{character} format.
# @param no.init.codon Should the first codon considered not to be an
# initialization codon? Default: \code{TRUE}
#
# @return A \code{character} sequence corresponding to the sequence
# translation.
translate_sequence <- function(seq, no.init.codon = TRUE) {
    Biostrings::DNAStringSet(seq) %>%
        Biostrings::translate(no.init.codon = no.init.codon) %>%
        as.character
}
