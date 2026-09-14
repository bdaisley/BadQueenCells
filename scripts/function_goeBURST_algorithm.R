library(dplyr)
library(readr)
library(stringr)
library(igraph)
library(ggraph)
library(tidyr)
library(ggraph)
library(scatterpie)
library(scales)
library(V8)
library(jsonlite)

goeBURST_algorithm <- function(
    profiles,
    id_col = 1,
    locus_cols = NULL,
    remove_duplicate_profiles = TRUE
) {
  if (!is.data.frame(profiles)) {
    stop("profiles must be a data.frame.")
  }
  
  # ------------------------------------------------------------
  # Identify ID and locus columns
  # ------------------------------------------------------------
  
  if (is.numeric(id_col)) {
    id_col <- names(profiles)[id_col]
  }
  
  if (!id_col %in% names(profiles)) {
    stop("id_col was not found in profiles.")
  }
  
  if (is.null(locus_cols)) {
    locus_cols <- setdiff(names(profiles), id_col)
  }
  
  if (length(locus_cols) == 0) {
    stop("No locus columns were selected.")
  }
  
  if (!all(locus_cols %in% names(profiles))) {
    stop("One or more locus_cols were not found.")
  }
  
  ids <- as.character(profiles[[id_col]])
  
  if (anyNA(ids) || any(trimws(ids) == "")) {
    stop("Profile identifiers cannot be missing or empty.")
  }
  
  if (anyDuplicated(ids)) {
    stop("Profile identifiers must be unique.")
  }
  
  allele_data <- profiles[, locus_cols, drop = FALSE]
  
  # Treat allele numbers as categorical identifiers
  allele_matrix <- do.call(
    cbind,
    lapply(allele_data, as.character)
  )
  
  allele_matrix <- as.matrix(allele_matrix)
  colnames(allele_matrix) <- locus_cols
  
  if (anyNA(allele_matrix) || any(trimws(allele_matrix) == "")) {
    stop(
      "Missing allele values are present. ",
      "They must be handled before running goeBURST."
    )
  }
  
  # ------------------------------------------------------------
  # Collapse identical profiles
  #
  # PHYLOViZ Online retains the first identifier associated with
  # an allelic profile and treats later identical profiles as
  # duplicates.
  # ------------------------------------------------------------
  
  profile_key <- apply(
    allele_matrix,
    1,
    paste,
    collapse = "\034"
  )
  
  first_occurrence <- !duplicated(profile_key)
  
  representative_index <- match(
    profile_key,
    profile_key[first_occurrence]
  )
  
  duplicate_map <- data.frame(
    ID = ids,
    representative_ID = ids[first_occurrence][representative_index],
    identical_profile = duplicated(profile_key),
    stringsAsFactors = FALSE
  )
  
  profile_frequency <- vapply(
    profile_key[first_occurrence],
    function(key) sum(profile_key == key),
    integer(1)
  )
  
  if (!remove_duplicate_profiles && anyDuplicated(profile_key)) {
    stop(
      "Exact PHYLOViZ goeBURST requires unique allelic profiles. ",
      "Set remove_duplicate_profiles = TRUE."
    )
  }
  
  ids <- ids[first_occurrence]
  allele_matrix <- allele_matrix[first_occurrence, , drop = FALSE]
  
  number_profiles <- nrow(allele_matrix)
  number_loci <- ncol(allele_matrix)
  
  nodes <- data.frame(
    name = ids,
    input_index = seq_len(number_profiles),
    input_index_zero_based = seq_len(number_profiles) - 1L,
    frequency = profile_frequency,
    stringsAsFactors = FALSE
  )
  
  if (number_profiles == 1) {
    return(
      list(
        nodes = nodes,
        edges = data.frame(
          source = character(),
          target = character(),
          distance = integer()
        ),
        distance_matrix = matrix(
          0L,
          nrow = 1,
          ncol = 1,
          dimnames = list(ids, ids)
        ),
        locus_variant_counts = matrix(
          0L,
          nrow = 1,
          ncol = number_loci,
          dimnames = list(
            ids,
            paste0(seq_len(number_loci), "LV")
          )
        ),
        duplicate_map = duplicate_map,
        graph = NULL
      )
    )
  }
  
  # ------------------------------------------------------------
  # Pairwise Hamming distances and nLV counts
  #
  # lvs[i, d] = number of profiles differing from profile i
  #             at exactly d loci.
  # ------------------------------------------------------------
  
  distance_matrix <- matrix(
    0L,
    nrow = number_profiles,
    ncol = number_profiles,
    dimnames = list(ids, ids)
  )
  
  locus_variant_counts <- matrix(
    0L,
    nrow = number_profiles,
    ncol = number_loci,
    dimnames = list(
      ids,
      paste0(seq_len(number_loci), "LV")
    )
  )
  
  for (i in seq_len(number_profiles - 1L)) {
    for (j in seq.int(i + 1L, number_profiles)) {
      distance <- sum(
        allele_matrix[i, ] != allele_matrix[j, ]
      )
      
      if (distance == 0L) {
        stop(
          "Identical profiles remained after duplicate removal."
        )
      }
      
      distance_matrix[i, j] <- distance
      distance_matrix[j, i] <- distance
      
      locus_variant_counts[i, distance] <-
        locus_variant_counts[i, distance] + 1L
      
      locus_variant_counts[j, distance] <-
        locus_variant_counts[j, distance] + 1L
    }
  }
  
  # ------------------------------------------------------------
  # Edge comparator
  #
  # Returns:
  #    -1 when edge e iss preferred
  #   0 when = equivalent
  #   1 when edge f is preferred
  # ------------------------------------------------------------
  
  compare_edges <- function(e, f) {
    e_distance <- distance_matrix[e[1], e[2]]
    f_distance <- distance_matrix[f[1], f[2]]
    
    # Rule 0: shortest allelic distance first
    if (e_distance < f_distance) return(-1L)
    if (e_distance > f_distance) return(1L)
    
    # Rules 1...n:
    # compare max endpoint nLV count, then min value
    for (level in seq_len(number_loci)) {
      e_max <- max(
        locus_variant_counts[e[1], level],
        locus_variant_counts[e[2], level]
      )
      
      f_max <- max(
        locus_variant_counts[f[1], level],
        locus_variant_counts[f[2], level]
      )
      
      if (e_max > f_max) return(-1L)
      if (e_max < f_max) return(1L)
      
      e_min <- min(
        locus_variant_counts[e[1], level],
        locus_variant_counts[e[2], level]
      )
      
      f_min <- min(
        locus_variant_counts[f[1], level],
        locus_variant_counts[f[2], level]
      )
      
      if (e_min > f_min) return(-1L)
      if (e_min < f_min) return(1L)
    }
    
    # Final PHYLOViZ Online tie-break
    ##--------------------------------
    # The javascript source uses zero-based input indices:
    # first lower max endpoint index, then higher min
    # endpoint index
    e_zero <- e - 1L
    f_zero <- f - 1L
    
    e_max_index <- max(e_zero)
    f_max_index <- max(f_zero)
    
    if (e_max_index < f_max_index) return(-1L)
    if (e_max_index > f_max_index) return(1L)
    
    e_min_index <- min(e_zero)
    f_min_index <- min(f_zero)
    
    if (e_min_index > f_min_index) return(-1L)
    if (e_min_index < f_min_index) return(1L)
    
    0L
  }
  
  # ------------------------------------------------------------
  # Prim implementation
  # ------------------------------------------------------------
  
  in_tree <- rep(FALSE, number_profiles)
  parent <- rep(NA_integer_, number_profiles)
  
  # PHYLOViZ starts from the first profile
  in_tree[1] <- TRUE
  parent[seq.int(2L, number_profiles)] <- 1L
  
  selected_edges <- vector(
    "list",
    number_profiles - 1L
  )
  
  edge_number <- 0L
  
  while (sum(in_tree) < number_profiles) {
    remaining <- which(!in_tree)
    
    best_vertex <- remaining[1]
    
    if (length(remaining) > 1L) {
      for (vertex in remaining[-1]) {
        candidate_edge <- c(parent[vertex], vertex)
        current_edge <- c(
          parent[best_vertex],
          best_vertex
        )
        
        if (
          compare_edges(
            candidate_edge,
            current_edge
          ) < 0L
        ) {
          best_vertex <- vertex
        }
      }
    }
    
    best_parent <- parent[best_vertex]
    
    edge_number <- edge_number + 1L
    
    selected_edges[[edge_number]] <- data.frame(
      source = ids[best_parent],
      target = ids[best_vertex],
      distance = distance_matrix[
        best_parent,
        best_vertex
      ],
      source_index = best_parent,
      target_index = best_vertex,
      stringsAsFactors = FALSE
    )
    
    in_tree[best_vertex] <- TRUE
    
    # Update best incoming edge for every remaining profile
    remaining <- which(!in_tree)
    
    for (vertex in remaining) {
      new_edge <- c(best_vertex, vertex)
      old_edge <- c(parent[vertex], vertex)
      
      if (compare_edges(new_edge, old_edge) < 0L) {
        parent[vertex] <- best_vertex
      }
    }
  }
  
  edges <- do.call(rbind, selected_edges)
  rownames(edges) <- NULL
  
  graph <- NULL
  
  if (requireNamespace("igraph", quietly = TRUE)) {
    graph <- igraph::graph_from_data_frame(
      edges[, c("source", "target", "distance")],
      directed = FALSE,
      vertices = nodes
    )
  }
  
  list(
    nodes = nodes,
    edges = edges,
    distance_matrix = distance_matrix,
    locus_variant_counts = locus_variant_counts,
    duplicate_map = duplicate_map,
    graph = graph
  )
}

equal_edge_tree_layout <- function(
    graph,
    edge_length = 1,
    root = NULL,
    fan_width = pi * 0.95
) {
  if (!igraph::is_tree(graph)) {
    stop("The graph must be a connected tree.")
  }
  
  number_nodes <- igraph::vcount(graph)
  number_edges <- igraph::ecount(graph)
  
  if (is.null(root)) {
    root <- which.max(igraph::degree(graph))
  }
  
  if (is.character(root)) {
    root <- match(root, igraph::V(graph)$name)
    
    if (is.na(root)) {
      stop("The requested root was not found.")
    }
  }
  
  if (length(edge_length) == 1) {
    edge_length <- rep(edge_length, number_edges)
  }
  
  if (length(edge_length) != number_edges) {
    stop(
      "edge_length must have length 1 or equal the number of edges."
    )
  }
  
  # Edge-length lookup matrix
  edge.ends <- igraph::ends(
    graph,
    igraph::E(graph),
    names = FALSE
  )
  
  edge.length.matrix <- matrix(
    NA_real_,
    nrow = number_nodes,
    ncol = number_nodes
  )
  
  edge.length.matrix[
    cbind(edge.ends[, 1], edge.ends[, 2])
  ] <- edge_length
  
  edge.length.matrix[
    cbind(edge.ends[, 2], edge.ends[, 1])
  ] <- edge_length
  
  # Root the tree using breadth-first traversal
  adjacency <- igraph::as_adj_list(
    graph,
    mode = "all"
  )
  
  parent <- rep(NA_integer_, number_nodes)
  parent[root] <- 0L
  
  queue <- root
  traversal.order <- integer()
  
  while (length(queue) > 0) {
    vertex <- queue[1]
    queue <- queue[-1]
    
    traversal.order <- c(
      traversal.order,
      vertex
    )
    
    neighbours <- as.integer(
      adjacency[[vertex]]
    )
    
    for (neighbour in neighbours) {
      if (is.na(parent[neighbour])) {
        parent[neighbour] <- vertex
        queue <- c(queue, neighbour)
      }
    }
  }
  
  children <- lapply(
    seq_len(number_nodes),
    function(vertex) {
      which(parent == vertex)
    }
  )
  
  # Number of terminal descendants in each subtree
  subtree.weight <- rep(1, number_nodes)
  
  for (vertex in rev(traversal.order)) {
    child.vertices <- children[[vertex]]
    
    if (length(child.vertices) > 0) {
      subtree.weight[vertex] <- sum(
        subtree.weight[child.vertices]
      )
    }
  }
  
  coordinates <- matrix(
    0,
    nrow = number_nodes,
    ncol = 2
  )
  
  direction <- rep(0, number_nodes)
  
  for (vertex in traversal.order) {
    child.vertices <- children[[vertex]]
    
    if (length(child.vertices) == 0) {
      next
    }
    
    if (vertex == root) {
      lower.angle <- 0
      upper.angle <- 2 * pi
    } else {
      lower.angle <- direction[vertex] - fan_width / 2
      upper.angle <- direction[vertex] + fan_width / 2
    }
    
    weights <- subtree.weight[child.vertices]
    cumulative <- c(
      0,
      cumsum(weights / sum(weights))
    )
    
    child.angles <- lower.angle +
      (
        head(cumulative, -1) +
          tail(cumulative, -1)
      ) / 2 *
      (upper.angle - lower.angle)
    
    for (index in seq_along(child.vertices)) {
      child <- child.vertices[index]
      
      current.edge.length <-
        edge.length.matrix[vertex, child]
      
      coordinates[child, ] <-
        coordinates[vertex, ] +
        current.edge.length *
        c(
          cos(child.angles[index]),
          sin(child.angles[index])
        )
      
      direction[child] <- child.angles[index]
    }
  }
  
  rownames(coordinates) <-
    igraph::V(graph)$name
  
  coordinates
}