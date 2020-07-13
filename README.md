# dasatinib-TALL
## Overview

![Dasatinib project
random2](https://user-images.githubusercontent.com/19508721/87256975-c2a74980-c45c-11ea-9bb3-6194d09bee07.png)

### Part1 Data\_collection

-----

#### \-00.rsem\_RNAseq\_pipeline and rsemQCreport

We set up an RNAseq data analysis pipeline using **RSEM** package. First
**BAM** files from Discovery, Validation and TARGET cohorts were
converted into FASTQ files by **PICARD** package. The alignment report
were extracted and FPKM were used as the quantification result.  
\* See scripts: 00.rsem\_RNAseq\_pipeline.sh and 00.rsemQCreport.R

#### \-01.Discovery, Validation and TARGET data set

After run rsem\_RNAseq\_pipeline the raw data esets were created for the
Discovery, Validationa and TARGET cohorts.  
\* See scripts: 01.Discovery\_FPKM.eset.R; 02.Validation\_FPKM.eset.R;
03.TARGET\_FPKM.eset.R

#### \-02.mouse T cell development eset

T cell development expression profile was generated from **GSE15907**
data set. **23** T cell development stages were selected.  
\* See script:04.Mouse.eset.R

#### \-03.T-ALL PDX phosphoproteomic data eset

5 PDX samples (3 sensitive cases and 2 resistant cases) were treated
with or without Dasatinib (10nM, 1h) before multiplexed
tandem-mass-tab(TMT) labeling. log2 transformed peptide intensity were
used to create eset.  
\* See script:
05.PDX\_phosphoproteomic.eset.R

### Part2 T-ALL network buidling

-----

#### \-01.generate TARGET T-ALL cohort eset with published data from HTseq

We used the HTseq raw count data from [Liu’s
paper](https://pubmed.ncbi.nlm.nih.gov/28671688/). Each sample raw
counts were scaled up to 50M before log2 transfomation (with psudocount
of one). None informative features were filted out before building T-ALL
specific network. \* See scripts: 01.SJARACNeInputPreparation.R

#### \-02.Build T-ALL network

We used [SJARACNe](https://github.com/jyyulab/SJARACNe) to build T-ALL
specific transcriptional and signaling networks seperately.  
\* See scripts: 02.SJARACNEnetwork.sh

#### \-03.Create T-ALL network object

[NetBID](https://github.com/jyyulab/NetBID) were used to extract
transcriptional and signaling networks from SJARACNe results and two
networks were merged into one.  
\* See scripts: 03.CreateTALLnetwork.R

### Part3 Analysis

-----

> Steps:  
> 1.1 Discovery and TARGET datasets QC, feature selection and infer
> activity  
> 1.2 Combine Discovery, Validation and TARGET cohort dataset, QC ,
> remove batch effect and infer activity  
> 1.3 mouse dataset infer activity  
> 1.4 Calculate phosphorylate dataset kinase activity  
> 2.1 Discovery cohort master regulator discovery  
> 2.2 Kinase differential activity analysis  
> 2.3 Get biomarker scores in Discovery, TARGET, Validation and mouse
> development datasets  
> 2.4 TARGET cohort survival analysis  
> 2.5 mouse dataset development stage analysis

-----

