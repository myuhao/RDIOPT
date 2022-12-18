
# RDIOPT

<!-- badges: start -->
<!-- badges: end -->

RDIOPT is a simple wrapper around the
[DIOPT](https://www.flyrnai.org/diopt) Ortholog Finder. The intention
for this package is to allow querying orthologs from R directly, without
copy and paste.

RDIOPT relies on the [DIOPT
API](https://www.flyrnai.org/tools/diopt/web/api) to query and search
results.

## Installation

You can install the development version of RDIOPT from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("myuhao/RDIOPT")
```

## Example

Call `find_orthologs()` with some of your gene of interests. The
function should give you a nice `tibble` with all the relevant
information.

``` r
# Entrez ID always work
gene_of_interests = c(1232, 7316)
# if you have AnnotationDbi installed, the function converts it for you
if (require(org.Hs.eg.db)) {
  gene_of_interests = c("CCR3", "UBC")
}
orthologs = RDIOPT::find_orthologs(gene_of_interests, from = 9606, to = 7227, filter = "none")
colnames(orthologs)
#>  [1] "date"                         "search_gene_entrez"          
#>  [3] "filter"                       "from_entrez_geneid"          
#>  [5] "symbol"                       "description"                 
#>  [7] "chromosome"                   "gene_type"                   
#>  [9] "ensembl_geneid"               "diopt_version"               
#> [11] "from_id"                      "to_id"                       
#> [13] "to_symbol"                    "confidence"                  
#> [15] "score"                        "best_score"                  
#> [17] "max_score"                    "best_score_rev"              
#> [19] "best_score_count"             "mist_ppi"                    
#> [21] "mist_genetic"                 "to_entrez_geneid"            
#> [23] "species_id"                   "species_specific_geneid"     
#> [25] "species_specific_geneid_type" "count"                       
#> [27] "methods"
```

The result table has many relevant information, but it can be summarized
in to three categories:

1.  Function parameters.
2.  Source gene information.
3.  Target gene information.

## References

1.  [DIOPT](https://www.flyrnai.org/diopt)
2.  [API](https://www.flyrnai.org/tools/diopt/web/api#)
