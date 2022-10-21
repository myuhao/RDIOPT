#' @title Find orthologs between species
#'
#' @description Convert entrez gene ID from one species to another using the
#' [DIOPT](https://www.flyrnai.org/diopt) service.
#' This function is a wrapper of the DIOPT
#' [API](https://www.flyrnai.org/tools/diopt/web/api).
#'
#' @importFrom magrittr %>%
#' @importFrom purrr map imap pmap
#' @importFrom tibble enframe
#' @importFrom tidyr unnest
#' @importFrom dplyr select
#'
#' @param genes A single entrez gene id
#' @param from Source species. See @details
#' @param to Target species. See @details
#' @param filter Filter of the results. See @details.
#' @param version What version should be used. Currently v9 is still in beta.
#'
#' @return A tibble with requested information.
#'
#' @details The `from` and `to` arguments are the NCBI Taxonomy ID.
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
#'    - "3702": Arabidopsis thaliana (Thale cress)
#'
#' If `from`/`to` is no the same length as the input `genes` or 1,
#' will warning and recycle to match the length of gene.
#'
#' The `filter` argument filter based on score. Options are:
#' `none`/`best_match`/`exclude_score_less_1`/`exclude_score_less_2`.
#' See [API Documentation](https://www.flyrnai.org/tools/diopt/web/api)
#'
#' @export
convert_genes = function(genes, from = 9606, to = 7227, filter = "none", version = 8) {
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
submit_ids = function(id, from = 9606, to = 7227, filter = "none", version = 8) {
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
parse_search_detail = function(raw_res) {
  raw_res$search_details %>%
    as_tibble() %>%
    unnest_wider(gene_details) %>%
    rename(
      from_entrez_geneid = geneid
    ) %>%
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

#' Helper function to parse the results
#' @importFrom tibble enframe
#' @importFrom purrr map_chr
#' @importFrom tidyr unnest_wider unnest_longer
#' @importFrom dplyr rename mutate relocate last_col
#' @importFrom magrittr %>%
#'
#' @param raw_res The content of a httr response
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


#' Function to parse the response.
#' @importFrom stringr str_detect
#' @importFrom tibble tibble
#' @importFrom dplyr mutate
#' @importFrom tidyr unnest
#' @importFrom magrittr %>%
#'
#' @param raw_res The content of a httr response
parse_response = function(raw_res) {
  query_info = parse_search_detail(raw_res)
  results_info = NULL
  if (str_detect(raw_res$message, "No Ortholog Data")) {
    results_info = tibble(to_id = NA_character_, to_symbol = "No Ortholog Data")
  } else {
    results_info = parse_results(raw_res)
  }
  query_info %>%
    mutate(results = list(results_info)) %>%
    unnest(cols = results)
}
