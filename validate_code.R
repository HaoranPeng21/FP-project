#!/usr/bin/env Rscript
# Read-only code checks. No private data or external R packages are needed.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
root <- dirname(normalizePath(script, mustWork = TRUE))
files <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE)
files <- files[!file.info(files)$isdir]
relative <- substring(files, nchar(root) + 2)
keep <- !grepl("(^|/)\\.git(/|$)", relative)
files <- files[keep]; relative <- relative[keep]
allowed <- grepl("\\.(R|md)$", relative) | relative == ".gitignore"
if (any(!allowed)) stop("Non-code files found: ", paste(relative[!allowed], collapse = ", "))
if (any(nzchar(Sys.readlink(files)))) stop("Symlinks are not allowed in the code-only release")
if (any(grepl("(^|/)config\\.local\\.R$", relative))) stop("Remove private local configuration before release")
rfiles <- files[grepl("\\.R$", files)]
invisible(lapply(rfiles, parse))
forbidden_roots <- paste0("/", c("scratch", "home[0-9]*"), "/", collapse = "|")
for (p in rfiles) {
  text <- readLines(p, warn = FALSE)
  if (any(grepl(forbidden_roots, text))) stop("Machine-specific path in ", basename(p))
  if (any(grepl("system2?\\s*\\(", text))) stop("External executable dependency in ", basename(p))
}
cat(length(rfiles), "R files parsed; code-only file allowlist, path and symlink checks passed.\n")
