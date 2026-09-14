# Direct port of TEToolkit/Counting.py::resolve_annotation_ambiguity(): for
# each read whose single locus overlaps multiple features, split its weight
# proportionally to those features' already-resolved counts (processed in
# order, as in the Python version -- later ambiguous reads see earlier ones'
# contributions). Not vectorized across reads (the dependency is sequential
# by construction), but it only runs over the ambiguous subset, which is
# normally tiny relative to total reads.
.resolve_ambiguity <- function(counts, leftover_ids_list, leftover_weight) {
  for (k in seq_along(leftover_ids_list)) {
    ids <- unique(leftover_ids_list[[k]])
    w <- leftover_weight[k]
    vals <- counts[ids]
    total <- sum(vals)
    if (total > 0) {
      counts[ids] <- counts[ids] + w * vals / total
    } else {
      counts[ids] <- counts[ids] + w / length(ids)
    }
  }
  counts
}
