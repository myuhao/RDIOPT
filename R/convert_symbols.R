#' Convert symbols to entrezID
#' @rdname covert_symbol
#'
#' @param symbols Gene Symbols
#'
#' @return A vector of mapped entrez ID
#'
#' @details This function will return a vector of the same length as
#' the input. If the input doesn't mtach, will return `NA` instead
#'
#' @keywords internal
#'
#' @seealso [AnnotationDbi::mapIds()] and [org.Hs.eg.db::org.Hs.eg.db]
convert_to_entrezID = function(
    symbols
    ) {

  if (all(is.na(symbols))) {
    return(c("NA" = NA_character_))
  }

  if (any(is.na(symbols))) {
    stop("There are NA in the query, fix them first\n")
  }

  ids = tryCatch(
    suppressMessages({
      AnnotationDbi::mapIds(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = symbols,
        keytype = "SYMBOL",
        column = "ENTREZID"
      )
    }),
    error = function(cond) {
      out = rep(NA_character_, length = length(symbols))
      names(out) = symbols
      out
    }
  )

  ids

}


#' Check if an input is enterzID
#'
#' @param input A vector of gene symbol, or name
#'
#' @importFrom stringr str_detect
#'
#' @return Is the input already enterzID
#'
#' @keywords internal
is_entrezID = function(input, silent = TRUE) {

  not_entrez = str_detect(input, "\\D") | is.na(input)

  if (any(not_entrez) && !silent) {
    warning(input[not_entrez], " might not be a entrezID")
  }

  return(!any(not_entrez))
}
