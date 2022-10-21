
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

<div>

> **Note**
>
> Note that currently DIOPT API only supports **Entrez ID**. You will
> need to convert your genes accordingly.

</div>

``` r
library(RDIOPT)

gene_of_interests = c(1232, 7316)
orthologs = find_orthologs(gene_of_interests, from = 9606, to = 7227, filter = "none")
colnames(orthologs)
#>  [1] "date"                         "search_gene_entrez"          
#>  [3] "filter"                       "from_entrez_geneid"          
#>  [5] "symbol"                       "description"                 
#>  [7] "chromosome"                   "gene_type"                   
#>  [9] "ensembl_geneid"               "diopt_version"               
#> [11] "to_id"                        "to_symbol"                   
#> [13] "from_id"                      "confidence"                  
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

#### Function parameters

DIOPT recorded the following parameter used for the query:

    #> # A tibble: 7 × 4
    #>   date              search_gene_entrez filter diopt_version
    #>   <chr>                          <int> <chr>          <int>
    #> 1 22-10-20 22:17:56               1232 none               8
    #> 2 22-10-20 22:17:56               7316 none               8
    #> 3 22-10-20 22:17:56               7316 none               8
    #> 4 22-10-20 22:17:56               7316 none               8
    #> 5 22-10-20 22:17:56               7316 none               8
    #> 6 22-10-20 22:17:56               7316 none               8
    #> 7 22-10-20 22:17:56               7316 none               8

#### Source gene information

DIOPT will try to find information of the gene you submitted.

    #> # A tibble: 7 × 6
    #>   from_entrez_geneid symbol description                  chrom…¹ gene_…² ensem…³
    #>                <int> <chr>  <chr>                          <int> <chr>   <chr>  
    #> 1               1232 CCR3   C-C motif chemokine recepto…       3 protei… ENSG00…
    #> 2               7316 UBC    ubiquitin C                       12 protei… ENSG00…
    #> 3               7316 UBC    ubiquitin C                       12 protei… ENSG00…
    #> 4               7316 UBC    ubiquitin C                       12 protei… ENSG00…
    #> 5               7316 UBC    ubiquitin C                       12 protei… ENSG00…
    #> 6               7316 UBC    ubiquitin C                       12 protei… ENSG00…
    #> 7               7316 UBC    ubiquitin C                       12 protei… ENSG00…
    #> # … with abbreviated variable names ¹​chromosome, ²​gene_type, ³​ensembl_geneid

#### Source gene information

Here is the actual ortholog information. Fileds are left `NA` in case no
orthologs can be found.

    #> # A tibble: 7 × 17
    #>   to_id  to_symbol from_id confi…¹ score best_…² max_s…³ best_…⁴ best_…⁵ mist_…⁶
    #>   <chr>  <chr>     <chr>   <chr>   <int> <chr>     <int> <chr>     <int> <chr>  
    #> 1 <NA>   No Ortho… <NA>    <NA>       NA <NA>         NA <NA>         NA <NA>   
    #> 2 38456  Ubi-p63E  7316    high       11 Yes          15 Yes           2 -      
    #> 3 326237 Ubi-p5E   7316    modera…     5 No           15 Yes           1 -      
    #> 4 31564  CG11700   7316    low         2 No           15 No            0 -      
    #> 5 33629  RpL40     7316    low         1 No           15 No            0 -      
    #> 6 34420  RpS27A    7316    low         1 No           15 No            0 -      
    #> 7 35151  Nedd8     7316    low         1 No           15 No            0 -      
    #> # … with 7 more variables: mist_genetic <chr>, to_entrez_geneid <int>,
    #> #   species_id <int>, species_specific_geneid <chr>,
    #> #   species_specific_geneid_type <chr>, count <int>, methods <chr>, and
    #> #   abbreviated variable names ¹​confidence, ²​best_score, ³​max_score,
    #> #   ⁴​best_score_rev, ⁵​best_score_count, ⁶​mist_ppi

## Reference
