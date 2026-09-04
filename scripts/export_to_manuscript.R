# Copies the generated tables and figures into the manuscript directory, so that
# yeast.tex can \input them instead of carrying hand-transcribed copies.
#
#   Rscript scripts/export_to_manuscript.R [manuscript_directory]
#
# Run make_tables.R, make_retail_tables.R and make_figures.R first. This script
# only copies; it never edits an artifact. If a caption or a label is wrong, fix
# the generator and regenerate, so that the repository stays the single source.

DEFAULT_MANUSCRIPT_DIR = "../arXiv-2406.16523v5"

args = commandArgs(trailingOnly = TRUE)
manuscript_dir = if (length(args) >= 1) args[1] else DEFAULT_MANUSCRIPT_DIR

if (!dir.exists(manuscript_dir)) {
  stop(sprintf("manuscript directory not found: %s", manuscript_dir))
}

TABLES = list.files("tables", pattern = "[.]tex$", full.names = TRUE)
FIGURES = "power-4.png"

missing = FIGURES[!file.exists(FIGURES)]
if (length(TABLES) == 0 || length(missing) > 0) {
  stop(sprintf("run the generators first; missing: %s",
               paste(c(if (length(TABLES) == 0) "tables/*.tex", missing),
                     collapse = ", ")))
}

target_tables = file.path(manuscript_dir, "tables")
dir.create(target_tables, showWarnings = FALSE, recursive = TRUE)

for (artifact in TABLES) {
  file.copy(artifact, file.path(target_tables, basename(artifact)),
            overwrite = TRUE)
  cat(sprintf("copied %s -> %s\n", artifact, target_tables))
}
for (artifact in FIGURES) {
  file.copy(artifact, file.path(manuscript_dir, artifact), overwrite = TRUE)
  cat(sprintf("copied %s -> %s\n", artifact, manuscript_dir))
}
