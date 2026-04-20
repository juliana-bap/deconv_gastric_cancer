# Setup — Environment installation

This folder contains everything needed to bootstrap the R and Python environments for the gastric cancer deconvolution pipeline.

---

## 1. R packages

Run the install script once after cloning the repository:

```bash
Rscript scripts/setup/install_r_deps.R
```

This installs all CRAN, Bioconductor, and GitHub packages required by scripts `00_*` to `09c_*`. The script is idempotent — re-running it skips already-installed packages.

### Required R version

- **R ≥ 4.4** recommended
- **Seurat ≥ 5.0** (the pipeline uses the v5 layer architecture)

### Optional: lock versions with `renv`

For full reproducibility, after the first successful pipeline run, initialize `renv` from inside R:

```r
install.packages("renv")
renv::init()        # scans the project and creates renv.lock
renv::snapshot()    # commit the lockfile
```

Other contributors can then restore the exact same environment:

```r
renv::restore()
```

The resulting `renv.lock` and `renv/activate.R` should be committed to the repository. The `renv/library/` folder is automatically gitignored.

---

## 2. Python environment (CellTypist)

Script `09b_celltypist.py` requires Python with `scanpy` and `celltypist`. Create the environment from the provided `environment.yml` (at the repository root):

```bash
conda env create -f environment.yml
conda activate deconv_gc_py
```

To update later:

```bash
conda env update -f environment.yml --prune
```

### Telling R where Python is

Scripts `09a_pre_annotation_singler.R` and `09b_run_celltypist.R` use **reticulate** to talk to Python. Point reticulate at the conda env by adding this line to `~/.Renviron` (and restarting R):

```
RETICULATE_PYTHON=/path/to/miniconda3/envs/deconv_gc_py/bin/python
```

You can find the path with:

```bash
conda activate deconv_gc_py
which python
```

---

## 3. Quarto (for notebooks)

Install from <https://quarto.org/docs/get-started/>. Tested with Quarto ≥ 1.4.

To render a notebook:

```bash
cd notebooks/<topic>
quarto render <notebook>.qmd
```

---

## 4. HPC server setup (PBS)

On the cluster, the pipeline uses a conda environment named `R_servidor`. PBS submission scripts under `scripts/sc_pre_proc/pbs/` already activate it. Make sure that env contains the same R packages listed above (you can run `install_r_deps.R` once on a login node to populate it).

---

## Quick check

After installation, verify everything is in place:

```bash
# R side
Rscript -e 'library(Seurat); library(SingleR); library(ROGUE); cat("R OK\n")'

# Python side
conda activate deconv_gc_py
python -c "import scanpy, celltypist; print('Python OK')"
```

If both lines print `OK`, the environment is ready and you can start running the pipeline (see the main `README.md` for the run order).
