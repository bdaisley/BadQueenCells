function_rollup_import <- function(input.list = input.list, rollup.list= rollup.list){
  collector.rollup <- list()
  for (x in rollup.list){
    d.b.collector <- data.frame(concat_column="")
    for(i in input.list){
      d.b.rollup <- tryCatch({
        df <- readr::read_tsv(paste0("data/metacerberus_results/final/rollup/", "rollup_prodigal_", i, "_trim_nohost-", x, ".tsv"))
        columns_to_concat <- colnames(df)[!colnames(df) %in% "Count"]
        df <- df %>% mutate(concat_column = apply(dplyr::select(., all_of(columns_to_concat)), 1, paste, collapse = "__")) %>%
          mutate(concat_column = gsub(" ", "_", concat_column)) %>% dplyr::select(concat_column, Count) %>% dplyr::rename(!!i := Count)
      }, error = function(e) {
        warning(paste("File not found or failed to read:", i, "_nohost-", x))
        return(NULL)
      })
      if (!is.null(d.b.rollup)) {
        d.b.collector <- merge(d.b.collector, d.b.rollup, by = "concat_column", sort = FALSE, all = TRUE) %>%
          mutate(across(-concat_column, ~ ifelse(is.na(.), 0, .))) %>% filter(concat_column != "")
      }
      if (is.null(d.b.rollup)) {
        d.b.rollup <- data.frame(concat_column="") %>% mutate(!!i := NA)
        d.b.collector <- merge(d.b.collector, d.b.rollup, by = "concat_column", sort = FALSE, all = TRUE) %>%
          mutate(across(-concat_column, ~ ifelse(is.na(.), 0, .))) %>% filter(concat_column != "")
      }
    }
    collector.rollup[[x]] <- d.b.collector
  }
  return(collector.rollup)
}