## ---------------------------------------------------------------------------
## Workaround for an AnnotationHub 3.6.0 / dbplyr 2.5.0 incompatibility.
##
## AnnotationHub (and BiocFileCache) call `dplyr::collect(tbl, Inf)`, passing the
## row limit POSITIONALLY. In dbplyr 2.x the signature is
##     collect.tbl_sql(x, ..., n = Inf, warn_incomplete = TRUE)
## so `Inf` is captured by `...` instead of `n`. dbplyr's `db_collect()` then
## calls rlang::check_dots_used(), which aborts with
##     "Arguments in `...` must be used. Problematic argument: ..1 = Inf"
## surfacing to the user as "Failed to collect lazy table" and, from sesameData,
## the opaque "ExperimentHub Caching fails".
##
## Because `n = Inf` is also the default, simply ignoring the stray dot restores
## the intended behaviour. Neutralising rlang's dots checks for the lifetime of
## this R session is therefore both sufficient and semantically safe here.
##
## Trap: use `body(f) <- parse(text = txt)[[1]]`, never `eval(parse(...))`.
## ---------------------------------------------------------------------------
## Note: assignInNamespace("check_dots_used", ..., ns = "rlang") is NOT enough.
## dbplyr imports the symbol at install time, so its imports environment holds a
## separate (locked) binding to the original closure. The imports env of every
## consuming package must be rebound as well.
patch_ah_collect <- function() {
  noop <- function(...) invisible(NULL)
  nms  <- c("check_dots_used", "check_dots_empty", "check_dots_empty0")
  n <- 0L
  for (nm in nms) {
    ok <- try(assignInNamespace(nm, noop, ns = "rlang"), silent = TRUE)
    if (!inherits(ok, "try-error")) n <- n + 1L
  }
  for (pkg in c("dbplyr", "dplyr", "AnnotationHub", "BiocFileCache",
                "ExperimentHub", "AnnotationDbi")) {
    if (!nzchar(system.file(package = pkg))) next
    ns <- try(asNamespace(pkg), silent = TRUE)
    if (inherits(ns, "try-error")) next
    imp <- parent.env(ns)
    for (nm in nms) {
      if (!exists(nm, envir = imp, inherits = FALSE)) next
      if (bindingIsLocked(nm, imp)) unlockBinding(nm, imp)
      assign(nm, noop, envir = imp)
      lockBinding(nm, imp)
      n <- n + 1L
    }
  }
  message("patched rlang dots checks (no-op) at ", n, " binding(s)")
  invisible(n)
}
