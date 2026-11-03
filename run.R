#!/usr/bin/env Rscript
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
repo_root <- dirname(normalizePath(script, mustWork = TRUE))
scripts <- sort(list.files(file.path(repo_root, "analyses"), pattern = "\\.R$", full.names = TRUE))
ids <- sub("\\.R$", "", basename(scripts))
args <- commandArgs(TRUE)
if (!length(args) || identical(args, "--list")) {
  cat(paste(ids, collapse = "\n"), "\n")
  quit(status = 0)
}
if (!identical(args, "all")) {
  wanted <- unique(unlist(lapply(args, function(a) which(startsWith(ids, a)))))
  if (!length(wanted) || any(!vapply(args, function(a) any(startsWith(ids, a)), logical(1)))) stop("Unknown analysis. Use --list.")
  scripts <- scripts[wanted]; ids <- ids[wanted]
}
source(file.path(repo_root, "R", "common.R"))
receipt <- list()
for (i in seq_along(scripts)) {
  out_dir <- file.path(output_root, ids[i])
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  message("Running ", ids[i])
  started <- Sys.time()
  env <- new.env(parent = globalenv())
  status <- tryCatch({
    sys.source(scripts[i], envir = env)
    env$run()
    "PASS"
  }, error = function(e) {message(conditionMessage(e)); paste0("FAIL: ", conditionMessage(e))})
  receipt[[i]] <- data.frame(analysis = ids[i], status = status,
    elapsed_seconds = round(as.numeric(difftime(Sys.time(), started, units = "secs")), 2),
    pdf_files = length(list.files(out_dir, pattern = "\\.pdf$")))
  write.table(do.call(rbind, receipt), file.path(output_root, "run_receipt.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
capture.output(sessionInfo(), file = file.path(output_root, "sessionInfo.txt"))
paths <- as.list(used_inputs)
if (length(paths)) write.table(data.frame(key = names(paths), path = unlist(paths),
  md5 = vapply(paths, function(p) if (file.info(p)$isdir) NA_character_ else unname(tools::md5sum(p)), character(1))),
  file.path(output_root, "private_input_receipt.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (any(vapply(receipt, function(x) x$status != "PASS", logical(1)))) quit(status = 1)
