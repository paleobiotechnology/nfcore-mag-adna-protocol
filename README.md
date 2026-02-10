# nf-core/mag protocol supplement

This is a complementary step by step guide with explicit commands for running nf-core/mag for ancient metagenomes.

## Requirements

This tutorial assumes you have:

- A UNIX based machine (Linux, OSX)
- At least 768.GB memory
- At least 64 CPUs
- At least 300 GB of free harddrive space
- An internet connection

> [!NOTE]
> If you do not have sufficient resources, you can still follow the tutorial.
> Skipping to [Evaluating the results](#evaluating-the-results), where the necessary results files are available premade for you.

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

2. (⚠️ if not already installed) Install conda through miniforge

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
   nextflow pull nf-core/mag -r 5.4.0
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
     screen -R gdtb_download ## to allow disconnection from server while running, use ctrl + a + d to detach, and `screen -r gtdb_download` to re-attach
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

From this sample three libraries have been generated and sequenced across three sequencing runs.
All three were built into a double stranded library with the PfuTurbo Cx Hotstart DNA polymerase that is tolerant of ancient DNA damage and sequenced on the Illumina NextSeq 500 platform.

| sample_name | library_name       | archive | archive_project | archive_sample_accession | archive_data_accession | strand_type | library_polymerase       | library_treatment | instrument_model | library_layout | read_count |
| ----------- | ------------------ | ------- | --------------- | ------------------------ | ---------------------- | ----------- | ------------------------ | ----------------- | ---------------- | -------------- | ---------- |
| ECO004.B    | ECO004.B0101       | ENA     | PRJEB34569      | ERS3774460               | ERR3579731             | double      | PfuTurbo Cx Hotstart DNA | none              | NextSeq 500      | PAIRED         | 13,625,140 |
| ECO004.B    | ECO004.B0103.SG1.1 | ENA     | PRJEB34569      | ERS3774460               | ERR3579732             | double      | PfuTurbo Cx Hotstart DNA | full-udg          | NextSeq 500      | PAIRED         | 67,123,133 |
| ECO004.B    | ECO004.B0102       | ENA     | PRJEB55583      | ERS3774460               | ERR10114849            | double      | PfuTurbo Cx Hotstart DNA | half-udg          | NextSeq 500      | PAIRED         | 21,007,806 |

The primary difference between the three libraries is each have different 'UDG' treatment applied, that removes increasingly more damage from the DNA molecules, and the sequencing depth.
The two shallow-sequenced libraries had no or only partial damage removal, whereas the deeper sequenced library had full-UDG treatment and thus has had all damage removed.

We can retrieve the raw sequencing data and prepared input sheets from via the tool `amdirt` ([Borry et al. 2024, F1000Research](https://doi.org/10.12688/f1000research.134798.2)) from the European Nucleotide Archive (ENA).

We can install `amdirt` in a separate conda environment:

```bash
   conda create -y -n amdirt -c bioconda amdirt=1.7.1
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
nextflow run nf-core/mag \
-r 5.4.0 \
-profile conda \
--input $TUTORIAL_DIR/analysis/mag/AncientMetagenomeDir_nf_core_mag_input_paired_table.csv \
--outdir $TUTORIAL_DIR/analysis/mag/results \
--reads_minlength 30 \
--igenomes_base 's3://ngi-igenomes/igenomes/' \
--host_genome GRCh37 \
--skip_spades \
--skip_spadeshybrid \
--skip_prodigal \
--skip_metaeuk \
--binning_map_mode own \
--min_contig_size 500 \
--save_assembly_mapped_reads \
--exclude_unbins_from_postbinning \
--run_checkm \
--run_busco false \
--checkm_db $TUTORIAL_DIR/cache/database/checkm_data_2015_01_16 \
--refine_bins_dastool \
--refine_bins_dastool_threshold 0.3
--postbinning_input refined_bins_only \
--run_gunc \
--gunc_db $TUTORIAL_DIR/cache/database/gunc_db/gunc_db_progenomes2.1.dmnd \
--gtdb_db $TUTORIAL_DIR/cache/database/release226 \
--ancient_dna \
-c $TUTORIAL_DIR/analysis/mag/custom.conf
```

> [!TIP]
> If you get an error message of `ERROR ~ /analysis`, check that the `$TUTORIAL_DIR environment variable is set correctly.

### The JSON way

However, setting up each parameter directly from the CLI can quickly become cumbersome.

Alternatively, another way to specify parameters is to use a parameters JSON file, here saved as `nf-params.json`.
This JSON file can either be prepared manually in a text editor, or with the help of the online interactivate nf-core/launch tool ([https://nf-co.re/launch/](https://nf-co.re/launch/)).
One parameters have been specified, nf-core tools provides a `nf-params.json` to download.

![Screenshot of the nf-core-launch page for nf-core/mag 5.4.0, with a few input parameters shown (free text and radio button option ones)](assets/tutorial-figure-launchpage.png)

This is the content of the `nf-params.json` file

```json
{
  "reads_minlength": 30,
  "host_genome": "GRCh37",
  "skip_spades": true,
  "skip_spadeshybrid": true,
  "skip_prodigal": true,
  "skip_metaeuk": true,
  "binning_map_mode": "own",
  "min_contig_size": 500,
  "save_assembly_mapped_reads": true,
  "exclude_unbins_from_postbinning": true,
  "run_checkm": true,
  "run_busco": false,
  "refine_bins_dastool": true,
  "refine_bins_dastool_threshold": 0.3,
  "postbinning_input": "refined_bins_only",
  "run_gunc": true,
  "ancient_dna": true
}
```

This `nf-params.json` is then used to specify the parameters on the command line like so:

```bash
nextflow run nf-core/mag \
-r 5.4.0 \
-profile conda \
--input $TUTORIAL_DIR/analysis/mag/AncientMetagenomeDir_nf_core_mag_input_paired_table.csv \
--outdir $TUTORIAL_DIR/analysis/mag/results_params \
--igenomes_base 's3://ngi-igenomes/igenomes/' \
--checkm_db $TUTORIAL_DIR/cache/database/checkm_data_2015_01_16 \
--gunc_db $TUTORIAL_DIR/cache/database/gunc_db/gunc_db_progenomes2.1.dmnd \
--gtdb_db $TUTORIAL_DIR/cache/database/release226 \
-c $TUTORIAL_DIR/analysis/mag/custom.conf \
-params-file $TUTORIAL_DIR/analysis/mag/nf-params.json
```

> [!TIP]
> Using a `params.json` also helps promotes reproducible science, allowing yourself but also other researchers to reproduce your results without risk of typo-derived errors!
>
> We recommend using the JSON method for supplying parameters, and include this file as supplementary data in your publications (make sure not to include file paths, however - these should stay on the command line!).

For us on a shared(!) SLURM HPC cluster, the run (including the creation of fresh conda environments) took:

```log
Duration    : 14h 2m 23s
CPU hours   : 318.4
```

## Evaluating the results

Once our pipeline run has completed, before we dive into analysing each MAG separately, we want to get an overview on the quality of the data, assembly, and final bins.

For this, there are two primary files you should inspect:

- `<--outdir>/multiqc/multiqc_report.tsv`
- `<--outdir>/GenomeBinning/bin_summary.tsv`

### Raw data, preprocessing and assembly evaluation

To firstly evaluate the raw data (assuming not already performed prior running nf-core/mag), preprocessing. and initial assembly, we want to inspect the `multiqc_report.tsv` file.

> [!NOTE]
> If you have not executed the pipeline or still waiting for it to finish, you can find already made files in `$TUTORIAL_DIR/analysis/mag/data/premade_mag_results`.
>
> It does not matter which execution run you look at!

In the MultiQC report, the primary sections you will want to evaluate are:

- FastQC (before and after)
- fastp
- Bowtie2: PhiX removal
- Bowtie2: host removal
- QUAST: assembly

The remaining sections will be covered in [Binning evaluation](#binning-evaluation).

Following the summary of recommendations in Box 1 of Fellows Yates et al. 2026, we can use this report to evaluate the quality of data production, sequencing depth, and other preprocessing steps such as host removal.

#### FastQC (raw and after preprocessing)

As a general rule, evaluating the FastQC results follows typical NGS sequencing quality metrics.

In the context of running nf-core/mag with ancient DNA we should in particular check that in the first FastQC section 'Sequence Counts', we have the number of reads requested from the sequencing facility.
As a general rule, in metagenomic _de novo_ assembly the more the better - therefore you should double check you recieved what you expected (Recommendation 2 of Box 1 of Fellows Yates et al. 2026).

You can also check the 'Sequence Length Distribution' section for a distribution with a peak somewhere in the range of 20-70 bp, which is the typical length of ancient DNA reads.
If you have no peak around here, this _may_ (but not always!) indicate a badly preserved sample.

You should also compare the 'FastQC: preprocessing' results section with the first 'raw reads' section to make sure all metrics have improved as expected (e.g. adapter content has reduced - however in this tutorial as we are using public data, this has already been removed).

#### fastp

As a general rule, evaluating the fatp results follows typical NGS sequencing quality metrics.

You should see a general improvment in the reduction of N content and sequencing quality.

As with FastQC, the 'insert Size Distribution' sectuin can give you an indication of the presence of 'true' fragmented aDNA reads, with an expected peak - in this case of the ECO004 libraries around 30-40 bp.
If you see a loss of this short reads, it may imply your read filtering settings are incorrect (Recommendation 7 of Box 1 of Fellows Yates et al. 2026).

#### Bowtie2 (PhiX and host removal)

In the Bowtie2 sections, you expect to see a certain amount of PhiX (assuming it's used as a control by your sequencing facility) and host DNA (in the case of this tutorial, _Homo spaiens_) being removed.

The host DNA removed metrics should be particularly paid close attenion to, to ensure that you don't result in broken and/or chimeric contigs containing host DNA that can match contaminated reference genomes during taxonomic classification of contigs (Recommendation 4 of Box 1 of Fellows Yates et al. 2026).

The exact amount will depend on the sample type, and sequencing depth.
If you expected to see more, you should review your settings to ensure you have set alignment parameters suitable for aligning ultra-short DNA sequences (Recommendation 2 of Box 1 of Fellows Yates et al. 2026).

In the case of the tutorial, no PhiX spike-in was used during sequencing of the ECO004 libraries, and only a very small amount of Human DNA is expected (primarily from modern DNA contaminatiion during handling), as the sample type is from dental calculus, which is primarily a microbial biofilm that forms on the surface of teeth - rather than being a human tissue itself.

#### QUAST

We can finally use the MultiQC report to do an initial evaluation of the assemblies themselves, using typical _de novo_ assembly metrics (N50, Largest contig ).

We can already get a rough idea of the RATIO <!-- TODO --> (Recommendation 10 of Box 1 of Fellows Yates et al. 2026).

### Binning evaluation

Against summary of recommendations

1. MULTIQC: Did we have appropriate data production --> see input data
   - We can see the difference in the the retrieved bin likely derived from the same species
   - MEGAHIT-SemiBin2Refined-ERR10114849.1.fa; MEGAHIT-SemiBin2Refined-ERR3579731.1.fa; MEGAHIT-SemiBin2Refined-ERR3579732.110.fa
   - all hit against 'Clostridiales' with BUSCO and GTDBTL reports the latter two to be `d__Bacteria;p__Bacillota;c__Clostridia;o__Peptostreptococcales;f__Anaerovoracaceae;g__RGIG7111;s__RGIG7111 sp036839975`
2. MULTIQC: Sequence to a sufficient depth -> see data input
   - Confidence increases with more depth (in the full UDG One) for the `s__RGIG7111 sp036839975`
3. MULTIQC: State of the data -> see MultiQC
   1.
4. MULTIQC: Remove host sequences -> see MultiQC
5. Select appropriate tools and databases -> done by MAG
6. Co-assemble related samples -> did not do in this case to demonstrate difference
7. PARAMETERS Use aDNA tools and parameters -> done by MAG (`--ancient_dna`)
8. PARAMETERS Discard ultra-short contigs -> done by MAG (`--min_contig_size`)
9. Evaluate the prescne of chimeric MAGs -> ???
10. Check ratio of short to long contigs -> need to do manually
11. BINSUMMARY Validate low quality MAGs -> reduced bin score in DAS tool (` --refine_bins_dastool_threshold`)
12. BINSUMMARY Authetnicate contigs and bins -> check pyDamage results in `bin_summary.tsv`
13. BINSUMMARY Correct for aDNA damage -> done by MAG automatically

For pyDamage evaluation:

Alex Hübner The two important values that I look at generally are the predicted_accuracy and the qvalue.
The predicted_accuracy summarises how much data were available to compare them to the model predictions and usually correlates with the number of aligned reads. The rough guide is that you should have at least 50% predicted accuracy but Maxime wrote that there is a function that can evaluate the best cut off for the particular experiment. However, I have never used it. Second, I check the qvalue and simply ask whether it is < 0.05 or < 0.01 depending how accurate I would like to be.
The other summary statistics, such as null_model_p* or damage_model_p* , summarise the fit of each of the two models to the data. So a researcher could make use of these if they want to be very sure about their results. The column RSME summarises this to as it states how far the distance between the points and the curves inferred by the models are.
But to be far, I would use the latter only on edge cases when I am super interested in an individual contig and the predicted_accuracy and qvalue are not super indicative, whether there is ancient DNA damage. And this happens usually never.
[2:09 PM]James Fellows Yates Hahaha ok ok. This then tracks with what I thought I think.

So would you say we could maybe advise in the paper then for general rule of thumb 'evaluation':

You can be reasonable comfortable reporting an ancient bin if:

predicted_accuracy >= 0.5 (or 50%)q_value <= 0.05 (or whatever you feel is a suitable statistical significance threshold)(If you want 'classical' damage smiley plot information) Decreasing values from CtoT-0 to CtoT-12 columns according to the damage treatment type (if any)[2:09 PM]Would that be reasonable?
[2:11 PM]Alex Hübner Yes, it is. I would maybe add that you should take contigs with < 1,000 aligned reads with caution because of the stochasticity of the data. That's it.
[2:18 PM]James Fellows Yates OK awesome (even if our test data violates the last one, we aren't toooo bad when looking at the quality of the data)
[2:18 PM]Thank you! I'll add that to the text and to the tutorial

## Clean up

To remove the entire tutorial directory

```bash
conda init --reverse $(basename SHELL)

## Triple check this command before executing that it definitely points to the tutorial directory, -rf can be dangerous!
cd ../
rm -rf ./ancientdna-nfcoremag-tutorial/
```
