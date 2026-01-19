# nf-core/mag protocol supplement

This is a complementary step by step guide with explicit commands for running nf-core/mag for ancient metagenomes.

## Requirements

This tutorial assumes you have:

- A UNIX based machine (Linux, OSX)
- At least 720.GB memory
- At least 64 CPUs
- At least 175 GB of free harddrive space
- An internet connection

## Installation and Setup

1. Set up directory structure for tutorial by cloning this repository and creating the required subdirectories

   ```bash
    git clone https://github.com/paleobiotechnology/nfcore-mag-adna-protocol.git

   ## Change into the cloned repo, and set this as root for the remainder of the tutorial
   cd nfcore-mag-adna-protocol
   export TUTORIAL_DIR=$(pwd)
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
   conda create -y -n nextflow -c bioconda nextflow=25.04.2
   ```

5. Load environment, and set the NXF_HOME to allow efficient cleanup at end of tutorial

   ```bash
   conda activate nextflow

   ## Set home
   mkdir -p bin/nextflow/assets
   export NXF_HOME=${TUTORIAL_DIR}/bin/nextflow/
   export NXF_ASSETS=${TUTORIAL_DIR}/bin/nextflow/assets
   ```

6. Download nf-core/mag

   ```bash
   nextflow pull nf-core/mag -r 5.3.0
   ```

7. Deactivate conda environement

   ```bash
   conda deactivate
   ```

## Database Downloading

1. Make conda environment for tools that require a specific tool to download the database
   - GUNC (`-y` is specified in the conda command to automatically accept proposed dependencies, remove if you wish to check)

   ```bash
   conda create -y -n gunc -c bioconda gunc=1.0.6
   conda activate gunc
   mkdir $TUTORIAL_DIR/cache/database/gunc_db
   gunc download_db $TUTORIAL_DIR/cache/database/gunc_db
   conda deactivate
   cd $TUTORIAL_DIR/
   ```

> [!TIP]
> If you get errors about `ModuleNotFoundError: No module named 'pkg_resources'`
> While _within_ the GUNC conda environment, run the following command:
>
> ```bash
> conda install -c conda-forge setuptools
> ```

2. Download databases for each tool into cache directory
   - GTDB:

   ```bash
   screen -R gdtb_download ## to allow disconnection from server while running, ust ctrl + a + d to detach, and screen -r gtdb_download to re-attach
   wget -O $TUTORIAL_DIR/cache/database/gtdbtk_r226_data.tar.gz https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
   tar -xzf $TUTORIAL_DIR/cache/database/gtdbtk_r226_data.tar.gz gtdbtk_r226_data.tar.gz -C $TUTORIAL_DIR/cache/database/gtdbtk_r226
   cd $TUTORIAL_DIR/
   ```

   > [!WARNING]
   > This is very large >110GB and take a long time to download.
   > We recommend re-using an already downloaded database if possible.
   > In this case symlink the `gtdbtk_r226/` directory to the cache directory:
   >
   > ```bash
   > ln -s /<your>/<path>/<to>/gtdbtk_r226/ $TUTORIAL_DIR/cache/database/gtdbtk_r226
   > ```
   - CheckM:

     ```bash
     wget -O $TUTORIAL_DIR/cache/database/checkm_data_2015_01_16.tar.gz https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
     mkdir $TUTORIAL_DIR/cache/database/checkm_data_2015_01_16/
     tar -xzf $TUTORIAL_DIR/cache/database/checkm_data_2015_01_16.tar.gz -C $TUTORIAL_DIR/cache/database/checkm_data_2015_01_16
     cd $TUTORIAL_DIR/
     ```

## Example data retrieval

For the purposes of this tutorial, we will use an Iberian ancient human dental calculus sample 'ECO014.B', originally published in [Fellows Yates et al (2021, PNAS)](https://doi.org/10.1073/pnas.2021655118), that comes from a Mesolithic site in Valencia.
This sample was analysed in a study to reconstruct the taxonomic profiles of ancient oral microbiomes from in ancient hominins.

It has three sequencing runs from two libraries: two shallow-sequenced (one non-UDG treated, one-full UDG treated to keep and full remove damage respectively), and one deep sequenced with half-UDG treated. All we sequenced paried-end on an Illumina NextSeq 500 platform.

We can retrieve the raw sequencing data and prepared input sheets from via the tool `amdirt` ([Borry et al. 2024, F1000Research](https://doi.org/10.12688/f1000research.134798.2)) from the European Nucleotide Archive (ENA).

We can install `amdirt` in a separate conda environment:

```bash
   conda create -y -n amdirt -c bioconda amdirt=1.7.0
   conda activate amdirt
```

> [!NOTE]
> You could also use the `amdirt` graphical user interface filter for the correct rows in the table, however we assume you're running on a remote server due to the resource requirements for the pipeline.

First we download the AncientMetagenomeDir host-associated libraries table:

```bash
mkdir -p data/amdirt/
amdirt download  --output $TUTORIAL_DIR/data/amdirt/ -t ancientmetagenome-hostassociated -y samples -r v25.12.0
amdirt download  --output $TUTORIAL_DIR/data/amdirt/ -t ancientmetagenome-hostassociated -y libraries -r v25.12.0
```

We can then use basic `bash` tools to filter to the required sequencing runs.

```bash
for i in samples libraries; do
   grep -e 'project_name' -e 'ECO004.B' $TUTORIAL_DIR/data/amdirt/ancientmetagenome-hostassociated_${i}_v25.12.0.tsv > $TUTORIAL_DIR/data/amdirt/ancientmetagenome-hostassociated_${i}_v25.12.0_filtered.tsv
done
```

We can then convert these samplesheets using `amdirt` to different files for preparing the input for nf-core/mag.

```bash
## Raw data
amdirt convert --libraries $TUTORIAL_DIR/data/amdirt/ancientmetagenome-hostassociated_libraries_v25.12.0_filtered.tsv -o data/raw_data --curl $TUTORIAL_DIR/data/amdirt/ancientmetagenome-hostassociated_samples_v25.12.0_filtered.tsv ancientmetagenome-hostassociated

## Samplesheet
amdirt convert --libraries $TUTORIAL_DIR/data/amdirt/ancientmetagenome-hostassociated_libraries_v25.12.0_filtered.tsv -o analysis/mag/ --mag --bibliography $TUTORIAL_DIR/data/amdirt/ancientmetagenome-hostassociated_samples_v25.12.0_filtered.tsv ancientmetagenome-hostassociated
conda deactivate
```

We can then run the `curl` shell script to download the raw data from ENA, and update the samplesheet with the correct paths to the raw data.

```bash
screen -R eco-curl-dl
cd $TUTORIAL_DIR/data/raw_data/
bash AncientMetagenomeDir_curl_download_script.sh
cd $TUTORIAL_DIR/
```

Once downloaded, we can add the full paths in the samplesheet.

```bash
## To verify the correct paths
sed "s#,ERR#,$TUTORIAL_DIR/data/raw_data/ERR#g" analysis/mag/AncientMetagenomeDir_nf_core_mag_input_paired_table.csv

## To update the samplesheet
sed -i "s#,ERR#,$TUTORIAL_DIR/data/raw_data/ERR#g" analysis/mag/AncientMetagenomeDir_nf_core_mag_input_paired_table.csv
```

We will also need to do a little but of repair work because `amdirt` currently supports an older version of nf-core/mag

```bash
sed -i "s/group,/group,short_reads_platform,/g" analysis/mag/AncientMetagenomeDir_nf_core_mag_input_paired_table.csv
sed -i "s/ERS3774460,/ERS3774460,ILLUMINA,/" analysis/mag/AncientMetagenomeDir_nf_core_mag_input_paired_table.csv
```

## Change column header from 'library_id' to 'sample_id'

## Pipeline setup and run

There are multiple different ways to specify nf-core/mag parameters, we'll highlight the two main ones.

For both methods, first activate main Nextflow conda environment

```bash
conda activate nextflow
```

For preparing our command we will skip a few steps that are not necessary for the purposes of this tutorial for reasons of speed, such as skipping metaeuk (can be slow, for eukaryotic contig detection), SPAdes (which is very slow and requires large amounts of computational resources), and Prodigal annotation (which will instead be performed at the bin level with Prokka).
We also exclude unbinned contigs from post-binning to reduce run time.

We will then tweak some relevant parameters that should be adjusted due to the nature of ancient DNA samples, such as setting a short minimum read length, host genome to GRCh37 for host DNA removal, as well as removing contigs less than 500 bp (which we expect many of because of the fragmented nature of aDNA).

We will also run co-assembly as we have three related libraries from the same sample.

Finally we also run `--ancient_dna` for automatic setting of suitable short-read-to-contig mapping parameters, damage authentication with pyDamage, and damage correction with freeBayes.

> [!WARNING]
> This command or `params.yml` method of configuration assumes you are running on a server with internet access, and specifically to the public AWS iGenome bucket for downloading pre-made host genome indices
> If you do not have internet access, you will need to download the host genome indices yourself, and specify the path to the local genome FASTA file with the `--host_fasta` parameter.
> Furthermore, if you have pre-built indices of the genome, you can additionally specify the `--host_fasta_bowtie2index` parameter to point to the a directory containing bowtie2 index files of the genome.

### The CLI way

We can run the pipeline via a full command line execution.

First make sure to create a `screen` session to allow disconnection from the server while running the pipeline, and everything you need to run the pipeline is loaded.

```bash
screen -R mag_run
export TUTORIAL_DIR=$(pwd)
cd analysis/mag
## Don't forget any other dependencies you need to run the pipeline on your specific infrastructure, such as `module load` commands!
```

Then you can run the pipeline with the following command.

Note that a backslash character can be used to break up a long single command into multiple lines.

> [!TIP]
> Alternatively, replace `-profile conda` with an institutional existing profile e.g. from nf-core/configs or a custom one with -c
> For example, in the below, we use a custom config file `custom.config` that is located in the same directory as the samplesheet to increase memory requirements for both the CHECKM_LINEAGEWF and GTDBTK_LINEAGEWF steps of the pipeline.

```bash
nextflow run nf-core/mag -r 5.3.0 \
-profile eva_grace,conda \
--input $TUTORIAL_DIR/analysis/mag/AncientMetagenomeDir_nf_core_mag_input_paired_table.csv \
--outdir  $TUTORIAL_DIR/analysis/mag/results \
--reads_minlength 30 \
--igenomes_base 's3://ngi-igenomes/igenomes/' \
--host_genome GRCh37 \
--skip_spades \
--skip_spadeshybrid \
--skip_prodigal \
--skip_metaeuk \
--binning_map_mode "own" \
--min_contig_size 500  \
--save_assembly_mapped_reads \
--exclude_unbins_from_postbinning \
--run_checkm \
--run_busco false \
--checkm_db $TUTORIAL_DIR/cache/database/checkm_data_2015_01_16 \
--refine_bins_dastool \
--postbinning_input "refined_bins_only" \
--run_gunc \
--gunc_db $TUTORIAL_DIR/cache/database/gunc_db/gunc_db_progenomes2.1.dmnd \
--postbinning_input "refined_bins_only" \
--gtdb_db $TUTORIAL_DIR/cache/database/release226 \
--ancient_dna \
-c custom.conf
```

> [!TIP]
> If you get an error message of `ERROR ~ /analysis`, check that the `$TUTORIAL_DIR environment variable is set correctly.

### The JSON way

However, setting up each parameter directly from the CLI can quickly become cumbersome.
To remedy this, another way to specify parameters is to use a parameters JSON file, here saved as `nf-params.json`.
This JSON file can either be prepared manually in a text editor, or with the help of the online interactivate nf-core/launch tool ([https://nf-co.re/launch/](https://nf-co.re/launch/)). One parameters have been specified, nf-core tools provides a `nf-params.json` to download.
![nf-core-launch](https://hackmd.io/_uploads/H1PFLMsqkx.png)

This is the content of the `nf-params.json` file

<!-- TO BE SYNCED WITH COMMAND ABOVE -->

```json
{
  "input": "samplesheet.csv",
  "outdir": "results",
  "reads_minlength": 30,
  "host_genome": "GRCh37",
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
  "gtdb_db": "ancientdna-nfcoremag-tutorial/cache/database/gtdbtk_r226",
  "ancient_dna": true
}
```

This `nf-params.json` is then used to specify the parameters on the command line like so:

```bash
nextflow run nf-core/mag -r 5.3.0 -params-file analysis/nf-params.json
```

## Evaluating the results

<!-- TODO -->

## Clean up

To remove the entire tutorial directory

```bash
conda init --reverse $(basename SHELL)

## Triple check this command before executing that it definitely points to the tutorial directory, -rf can be dangerous!
cd ../
rm -rf ./ancientdna-nfcoremag-tutorial/
```
