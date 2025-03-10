# nf-core/mag protocol supplement

This is a complementary step by step guide with explicit commands for running nf-core/mag for ancient metagenomes.

## Requirements

This tutorial assumes you have:

<!-- TODO -->

- A Unix based machine (Linux, OSX)
- At least XYZ memory
- At least XYZ CPUs
- At least XYZ gb harddrive space
- An internet connection

## Installation and Setup

1. Set up directory structure for tutorial by cloning this repository and creating the required subdirectories

   ```bash
    git clone https://github.com/paleobiotechnology/nfcore-mag-adna-protocol.git

   ## Change into the cloned repo, and set this as root for the remainder of the tutorial
   cd nfcore-mag-adna-protocol
   TUTORIAL_DIR=$(pwd)
   ```

> [!NOTE]
> We assume you will want to remove the entire contents of this tutorial on completion.
> If you wish to retain certain files (e.g. downloaded databases), make sure to place them in a safe location outside the tutorial directory, and update file paths accordingly.

2. (if not already installed) Install conda through miniforge

   ```bash
   ## Download the installation script
   curl -L -o ${TUTORIAL_DIR}/bin/Miniforge3-$(uname)-$(uname -m).sh "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

   ## Run installer
   bash ${TUTORIAL_DIR}/bin/Miniforge3-$(uname)-$(uname -m).sh -b -s -p ${TUTORIAL_DIR}/bin/Miniforge3

   ## Configure miniconda - replace 'YOUR_SHELL_NAME' in the `eval` command in lower case - e.g. shell.bash
   eval "$($TUTORIAL_DIR/bin/Miniforge3/bin/conda shell.YOUR_SHELL_NAME hook)"
   conda init
   conda config --set auto_activate_base false
   ```

> [!NOTE]
> Miniforge is the preferred distribution of conda, as it does not come with any restrictive usage licenses as it does not include the Anaconda Inc.'s paid default channel.

3. Exit your terminal, and make a new window. Change back into the tutorial directory.

   ```bash
   cd $TUTORIAL_DIR
   ```

   > You might have to redefine the `TUTORIAL_DIR` environment variable

4. Create environment (`-y` is specified to automatically accept proposed dependencies, remove if you wish to check)

   ```bash
   conda create -y -n nextflow -c bioconda nextflow=24.10.4
   ```

5. Load environment, and set the NXF_HOME to allow efficient cleanup at end of tutorial

   ```bash
   conda activate nextflow

   ## Set home
   mkdir -p bin/nextflow/assets
   NXF_HOME=${TUTORIAL_DIR}/bin/nextflow/
   NXF_ASSETS=${TUTORIAL_DIR}/bin/nextflow/assets
   ```

6. Download nf-core/mag

   ```bash
   nextflow pull nf-core/mag -r 3.3.1
   ```

7. Deactivate conda environement

   ```bash
   conda deactivate
   ```

8. Create a configuration file

   Visit nf-core/launch, and select the nf-core/mag pipeline, with the version 3.3.1.
   Alternatively, a `nf-params.json` file is already available in [`analysis/nf-params.json`](analysis/nf-params.json)

## Database Downloading

1. Make conda environment for tools that require a specific tool to download the database

   - GUNC (`-y` is specified in the conda command to automatically accept proposed dependencies, remove if you wish to check)

   ```bash
   conda create -y -n gunc -c bioconda gunc=1.0.6
   conda activate gunc
   gunc download_db $TUTORIAL_DIR/cache/database/gunc_db
   conda deactivate
   ```

2. Download databases for each tool into cache directory

   - GTDB:

   ```bash
   wget -O $TUTORIAL_DIR/cache/database/gtdbtk_r220_data.tar.gz https://data.gtdb.ecogenomic.org/releases/release220/220.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r220_data.tar.gz
   tar -xzf $TUTORIAL_DIR/cache/database/gtdbtk_r220_data.tar.gz gtdbtk_r220_data.tar.gz -C $TUTORIAL_DIR/cache/database/gtdbtk_r220
   ```

   - CheckM:

   ```bash
   wget -O $TUTORIAL_DIR/cache/database/checkm_data_2015_01_16.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
   tar -xzf $TUTORIAL_DIR/cache/database/checkm_data_2015_01_16.tar.gz -C $TUTORIAL_DIR/cache/database/checkm_data_2015_01_16
   ```

3. Activate main Nextflow conda environment

```bash
conda activate nextflow
```

## Pipeline setup and run

There are multiple different ways to specify nf-core/mag parameters, we'll highlight the two main ones.

### The CLI way

Via the CLI paramaters directly. Note that a backslash character can be used to break up a long single command into multiple lines.

```bash
nextflow run nf-core/mag -r 3.3.1 \
-profile conda \
--input samplesheet.csv \
--outdir ./results \
--reads_minlength 30 \
--host_genome GRCh37 \
--centrifuge_db false \
--kraken2_db false \
--skip_krona \
--krona_db false \
--skip_metaspades \
--skip_metaeuk \
--ancient_dna \
--bbnorm \
--binning_map_mode own \
--min_contig_size 500  \
--save_assembly_mapped_reads \
--exclude_unbins_from_postbinning \
--binqc_tool CheckM \
--checkm_db $TUTORIAL_DIR/ancientdna-nfcoremag-tutorial/cache/database/checkm_data_2015_01_16
--run_gunc \
--gunc_db $TUTORIAL_DIR/ancientdna-nfcoremag-tutorial/cache/database/gunc_db \
--postbinning_input both \
--cat_db false \
--gtdb_db $TUTORIAL_DIR/ancientdna-nfcoremag-tutorial/cache/database/gtdbtk_r220
```

> [!CAUTION]
> TODO check bbnorm usage

### The JSON way

However, setting up each parameter directly from the CLI can quickly become cumbersome.
To remedy this, another way to specify parameters is to use a parameters JSON file, here saved as `nf-params.json`.
This JSON file can either be prepared manually in a text editor, or with the help of the online interactivate nf-core/launch tool ([https://nf-co.re/launch/](https://nf-co.re/launch/)). One parameters have been specified, nf-core tools provides a `nf-params.json` to download.
![nf-core-launch](https://hackmd.io/_uploads/H1PFLMsqkx.png)

This is the content of the `nf-params.json` file

```json
{
  "input": "samplesheet.csv",
  "outdir": "results",
  "reads_minlength": 30,
  "host_genome": "GRCh37",
  "bbnorm": true,
  "skip_krona": true,
  "gtdb_db": "ancientdna-nfcoremag-tutorial/cache/database/gtdbtk_r220",
  "skip_spades": true,
  "skip_spadeshybrid": true,
  "skip_metaeuk": true,
  "binning_map_mode": "own",
  "min_contig_size": 500,
  "save_assembly_mapped_reads": true,
  "binqc_tool": "checkm",
  "checkm_db": "ancientdna-nfcoremag-tutorial/bin ancientdna-nfcoremag-tutorial/cache/database/db/checkm_db",
  "refine_bins_dastool": true,
  "postbinning_input": "both",
  "run_gunc": true,
  "gunc_db": "ancientdna-nfcoremag-tutorial/cache/database/gunc_db",
  "ancient_dna": true
}
```

This `nf-params.json` is then used to specify the parameters on the command line like so:

```bash
nextflow run nf-core/mag -r 3.3.1 -params-file analysis/nf-params.json
```

## Clean up

To remove the entire tutorial directory

```bash
conda init --reverse $(basename SHELL)

## Triple check this command before executing that it definitely points to the tutorial directory, -rf can be dangerous!
cd ../
rm -rf ./ancientdna-nfcoremag-tutorial/
```
