#' @title Find orthologs between species
#'
#' @description Convert entrez gene ID from one species to another using the
#' [DIOPT](https://www.flyrnai.org/diopt) service.
#' This function is a wrapper of the DIOPT
#' [API](https://www.flyrnai.org/tools/diopt/web/api).
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map imap pmap walk
#' @importFrom tibble enframe
#' @importFrom tidyr unnest
#' @importFrom dplyr select
#'
#' @param genes A chr vector of genes. See @details
#' @param from Source species. See @details
#' @param to Target species. See @details
#' @param filter Filter of the results. See @details.
#' @param version What version should be used. Currently v9 is still in beta.
#'
#' @return A tibble with requested information.
#'
#' @details
#' The `genes` argument can either be a enterzID, which will be used as is.
#' Or it can be a gene symbol, which will be converted to entrezID using
#' the [org.Hs.eg.db::org.Hs.eg.db] database.
#'
#' The `from` and `to` arguments are the NCBI Taxonomy ID.
#' DIOPT currently supports the following IDs:
#'
#'    - "all":  All Species (Avaiable for Target only)
#'    - "4896": Schizosaccharomyces pombe (Fission yeast)
#'    - "4932": Saccharomyces cerevisiae (Yeast)
#'    - "6239": Caenorhabditis elegans (Worm)
#'    - **"7227"**: Drosophila melanogaster (Fly)
#'    - "7955": Danio rerio (Zebrafish)
#'    - "8364": Xenopus tropicalis (Western clawed frog)
#'    - **"9606"**: Homo sapiens (Human)
#'    - "10090": Mus musculus (Mouse)
#'    - "10116": Rattus norvegicus (Rat)
#'    - "3702": Arabidopsis thaliana (Thale cress)
#'
#' If `from`/`to` is no the same length as the input `genes` or 1,
#' will warning and recycle to match the length of gene.
#'
#' The `filter` argument filter based on score. Options are:
#' `none`/`best_match`/`exclude_score_less_1`/`exclude_score_less_2`.
#' See [API Documentation](https://www.flyrnai.org/tools/diopt/web/api)
#'
#' The function is not vectorized, because DIOPT API is not vectorized.
#'
#'
#' @export
find_orthologs = function(genes, from = 9606, to = 7227, filter = "none", version = 8) {
  # Convert to Chr
  from = as.character(from)
  to = as.character(to)
  filter = as.character(filter)
  version = as.character(version)
  # Checkling parameters
  list(from, to) %>%
    walk(
      check_arg,
      allowed = c(
        "all", "4896", "4932", "6239",
        "7227", "7955", "8364",
        "9606", "10090", "10116", "3702"
        )
      )
  check_arg(filter, c("none", "best_match", "exclude_score_less_1", "exclude_score_less_2"))
  check_arg(version, c("8", "9"))

  # Conversion
  if (!is_entrezID(genes, silent = TRUE)) {
    rlang::check_installed("org.Hs.eg.db", reason = "Package \"org.Hs.eg.db\" is required for the gene symbol -> entrezID conversion")
    genes = convert_to_entrezID(genes)

    # After conversion clean up:
    if (length(genes) < 1) {
      stop("The input genes can't be converted to entrezIDs\n")
    }
  }

  # Fix parameters
  n_genes = length(genes)
  params = list(
    from = from,
    to = to,
    filter = filter,
    version = version
  ) %>%
    imap(
      function(arg, arg_name) {
        n_args = length(arg)
        if (n_args < n_genes) {
          # Recycle
          if (length(arg) != 1) {
            warning(arg_name, " is shorter than number of genes, will recycle...")
          }
          n_rep = n_genes %/% n_args + 1
          return(rep(arg, n_rep)[1:n_genes])
        } else if (n_args == n_genes) {
          # Perfect
          return(arg)
        } else {
          warning(arg_name, " is longer than number of genes. Use the first n_genes.")
          return(arg[1:n_genes])
        }
      }
    )

  # Run the query
  c(
    list(id = genes),
    params
  ) %>%
    pmap(submit_ids) %>%
    map(parse_response) %>%
    enframe() %>%
    unnest(value) %>%
    select(-name)

}


#' Check if all elements in the vecter is allowed.
#'
#' @importFrom stringr str_c
#' @importFrom dplyr setdiff
#'
#' @return The vector we are checking.
#'
#' @keywords internal
check_arg = function(
    vec,
    allowed
  ) {
  if (!all(vec %in% allowed)) {
    not_allowed = setdiff(vec, allowed)
    not_allowed = unique(not_allowed)
    if (length(not_allowed) > 3) {
      not_allowed = not_allowed[1:3]
    }
    not_allowed = str_c(not_allowed, collapse = ", ")
    stop("Parameters ", not_allowed, " are not supported. Please check function documentation")
  }
  invisible(vec)
}


#' Function to submit a query to DIOPT
#' @importFrom stringr str_c
#' @importFrom httr GET content
#' @importFrom magrittr %>%
#'
#' @param id A single entrez gene id
#' @param from Source species.
#' @param to Target species.
#' @param filter Filter of the results, see API [API](https://www.flyrnai.org/tools/diopt/web/api).
#' @param version What version should be used. Currently v9 is still in beta.
#'
#' @keywords internal
submit_ids = function(id, from = 9606, to = 7227, filter = "none", version = 8) {
  if (is.na(id)) {
    return(list(is_NA = TRUE))
  }
  response = str_c(
    "https://www.flyrnai.org/tools/diopt/web/diopt_api/",
    "v", version, "/",
    "get_orthologs_from_entrez/",
    from, "/",
    id, "/",
    to, "/",
    filter
  ) %>%
    httr::GET()

  if (response$status_code != 200) {
    warning("HTTP status warning: get status code: ", response$status_code)
  }
  response %>%
    httr::content()
}

#' Helper function to parse the query
#' @importFrom tibble as_tibble
#' @importFrom tidyr unnest_wider
#' @importFrom dplyr rename select
#' @importFrom magrittr %>%
#'
#' @param raw_res The content of a httr response
#'
#' @keywords internal
parse_search_detail = function(raw_res) {
  raw_res$search_details %>%
    as_tibble() %>%
    unnest_wider(gene_details) %>%
    rename(
      from_entrez_geneid = geneid
    ) %>%
    mutate(chromosome = as.character(chromosome)) %>%
    select(
      date,
      search_gene_entrez,
      filter,
      from_entrez_geneid,
      symbol,
      description,
      chromosome,
      gene_type,
      ensembl_geneid,
      diopt_version
    )
}

#' Helper function to parse the results that has ortholog
#' @importFrom tibble enframe
#' @importFrom purrr map_chr
#' @importFrom tidyr unnest_wider unnest_longer
#' @importFrom dplyr rename mutate relocate last_col
#' @importFrom magrittr %>%
#'
#' @param raw_res The content of a httr response
#'
#' @keywords internal
parse_results = function(raw_res) {
  raw_res$results %>%
    enframe(name = "from_id") %>%
    unnest_longer(value, indices_to = "to_id") %>%
    unnest_wider(value) %>%
    mutate(
      methods = map_chr(methods, str_c, collapse = "; ")
    ) %>%
    rename(
      to_symbol = symbol,
      to_entrez_geneid = geneid,
    ) %>%
    relocate(to_id, to_symbol, confidence, score, .after = from_id) %>%
    relocate(methods, .after = last_col())
}

#' Helper function to parse the results that has no ortholog
#'
#' @keywords internal
results_no_ortholog = function() {
  tibble::tibble(
    from_id = NA,
    to_id = NA_character_,
    to_symbol = "No Ortholog Data",
    confidence = NA,
    score = NA,
    best_score = NA,
    max_score = NA,
    best_score_rev = NA,
    best_score_count = NA,
    mist_ppi = NA,
    mist_genetic = NA,
    to_entrez_geneid = NA,
    species_id = NA,
    species_specific_geneid = NA,
    species_specific_geneid_type = NA,
    count = NA,
    methods = NA
    )
}

#' Helper function to parse the results when failed to convert to entrezID
#'
#' @importFrom magrittr %>%
#'
#' @keywords internal
results_NA_id = function() {
  tibble::tibble(
    date = NA,
    search_gene_entrez = NA,
    filter = NA,
    from_entrez_geneid = NA,
    symbol = NA,
    description = NA,
    chromosome = NA_character_,
    gene_type = NA,
    ensembl_geneid = NA,
    diopt_version = NA
  ) %>%
    tibble::add_column(results_no_ortholog()) %>%
    dplyr::mutate(to_symbol = NA)
}


#' Function to parse the response.
#' @importFrom stringr str_detect
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest
#' @importFrom magrittr %>%
#'
#' @param raw_res The content of a httr response
#'
#' @keywords internal
parse_response = function(raw_res) {
  if (isTRUE(raw_res$is_NA)) {
    return(results_NA_id())
  }

  query_info = parse_search_detail(raw_res)
  results_info = NULL
  if (str_detect(raw_res$message, "No Ortholog Data")) {
    results_info = results_no_ortholog()
  } else {
    results_info = parse_results(raw_res)
  }
  query_info %>%
    mutate(results = list(results_info)) %>%
    unnest(cols = results)
}


