# hlavbseq

HLAVBSeq is software to estimate the most likely HLA types from high-throughput sequencing data.

## Overview

## Dependencies

* [hla-vbseq 1](http://nagasakilab.csml.org/hla/HLAVBSeq.jar)
* [bwa 0.7.17](https://github.com/lh3/bwa/archive/0.7.17.tar.gz)


## Usage

### Cromwell
```
java -jar cromwell.jar run hlavbseq.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputBam`|File|BWA aligned reads.
`inputBai`|File|index for BWA aligned reads.
`bwaMem.runBwaMem_bwaRef`|String|The reference genome to align the sample with by BWA
`bwaMem.runBwaMem_modules`|String|Required environment modules
`bwaMem.readGroups`|String|Complete read group header line


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`outputFileNamePrefix`|String|""|Output prefix, customizable. Default is the first file's basename.
`add_params`|String|"-M -P -L 10000 -a"|HLAVBSeq uses additional parameters for BWA which need to be specified.
`filterTag`|Int|256|To address the issue with some data, we use 256 (secondary alignment) as a default
`bwa_threads`|Int|8|Threads to use with BWA aligner.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`extractReads.jobMemory`|Int|8|Memory allocated to the task.
`extractReads.timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`extractReads.hlaIntervals`|String|"chr6:29907037-29915661 chr6:31319649-31326989 chr6:31234526-31241863 chr6:32914391-32922899 chr6:32900406-32910847 chr6:32969960-32979389 chr6:32778540-32786825 chr6:33030346-33050555 chr6:33041703-33059473 chr6:32603183-32613429 chr6:32707163-32716664 chr6:32625241-32636466 chr6:32721875-32733330 chr6:32405619-32414826 chr6:32544547-32559613 chr6:32518778-32554154 chr6:32483154-32559613 chr6:30455183-30463982 chr6:29689117-29699106 chr6:29792756-29800899 chr6:29793613-29978954 chr6:29855105-29979733 chr6:29892236-29899009 chr6:30225339-30236728 chr6:31369356-31385092 chr6:31460658-31480901 chr6:29766192-29772202 chr6:32810986-32823755 chr6:32779544-32808599 chr6:29756731-29767588"|Intervals where HLA-alleles are sitting
`extractReads.modules`|String|"samtools/1.9"|Names and versions of required modules.
`indexBam.jobMemory`|Int|18|Memory allocated to the task.
`indexBam.overhead`|Int|6|Ovrerhead for calculating heap memory, difference between total and Java-allocated memory
`indexBam.timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`indexBam.bamNameIndexJar`|String|"$HLA_VBSEQ_ROOT/bin/bamNameIndex.jar"|Jar file for bamNameInder
`indexBam.modules`|String|"hla-vbseq/1 hlavbseq-bwa-index/2.0"|Names and versions of required modules.
`makeFastq.jobMemory`|Int|24|Memory allocated to the task.
`makeFastq.overhead`|Int|6|Ovrerhead for calculating heap memory, difference between total and Java-allocated memory
`makeFastq.timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`makeFastq.bamNameIndexJar`|String|"$HLA_VBSEQ_ROOT/bin/bamNameIndex.jar"|Jar file for bamNameInder
`makeFastq.picardParams`|String|"VALIDATION_STRINGENCY=LENIENT"|Additional parameters for picard SamToFastq, Default is VALIDATION_STRINGENCY=LENIENT
`makeFastq.modules`|String|"samtools/1.9 picard/2.21.2"|Names and versions of required modules.
`bwaMem.adapterTrimmingLog_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimmingLog_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.indexBam_timeout`|Int|48|Hours before task timeout
`bwaMem.indexBam_modules`|String|"samtools/1.9"|Modules for running indexing job
`bwaMem.indexBam_jobMemory`|Int|12|Memory allocated indexing job
`bwaMem.bamMerge_timeout`|Int|72|Hours before task timeout
`bwaMem.bamMerge_modules`|String|"samtools/1.9"|Required environment modules
`bwaMem.bamMerge_jobMemory`|Int|32|Memory allocated indexing job
`bwaMem.runBwaMem_timeout`|Int|96|Hours before task timeout
`bwaMem.runBwaMem_jobMemory`|Int|32|Memory allocated for this job
`bwaMem.adapterTrimming_timeout`|Int|48|Hours before task timeout
`bwaMem.adapterTrimming_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.adapterTrimming_addParam`|String?|None|Additional cutadapt parameters
`bwaMem.adapterTrimming_modules`|String|"cutadapt/1.8.3"|Required environment modules
`bwaMem.slicerR2_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR2_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR2_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.slicerR1_timeout`|Int|48|Hours before task timeout
`bwaMem.slicerR1_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.slicerR1_modules`|String|"slicer/0.3.0"|Required environment modules
`bwaMem.countChunkSize_timeout`|Int|48|Hours before task timeout
`bwaMem.countChunkSize_jobMemory`|Int|16|Memory allocated for this job
`bwaMem.numChunk`|Int|1|number of chunks to split fastq file [1, no splitting]
`bwaMem.doTrim`|Boolean|true|if true, adapters will be trimmed before alignment
`bwaMem.trimMinLength`|Int|1|minimum length of reads to keep [1]
`bwaMem.trimMinQuality`|Int|0|minimum quality of read ends to keep [0]
`bwaMem.adapter1`|String|"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"|adapter sequence to trim from read 1 [AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC]
`bwaMem.adapter2`|String|"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"|adapter sequence to trim from read 2 [AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT]
`filterBam.jobMemory`|Int|12|Memory allocated to the task.
`filterBam.timeout`|Int|12|Timeout in hours, needed to override imposed limits.
`filterBam.modules`|String|"samtools/1.9"|Names and versions of required modules.
`runHlaVbSeq.jobMemory`|Int|32|Memory allocated to the task.
`runHlaVbSeq.timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`runHlaVbSeq.overhead`|Int|6|Difference between memory allocated and memory used as Java heap
`runHlaVbSeq.alphaZero`|Float|0.01|Hyperparameter, as described in the paper it is not recommended to change the default (0.01)
`runHlaVbSeq.hlavbseqJar`|String|"$HLA_VBSEQ_ROOT/bin/HLAVBSeq.jar"|path to HLAVBSeq jar file
`runHlaVbSeq.hlaFasta`|String|"$HLAVBSEQ_BWA_INDEX_ROOT/hla_all_v2.fasta"|fasta with HLA sequences
`runHlaVbSeq.modules`|String|"hla-vbseq/1 hlavbseq-bwa-index/2.0"|Names and versions of required modules.
`parseResults.jobMemory`|Int|8|Memory allocated to the task.
`parseResults.timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`parseResults.parsingScript`|String|"$HLA_VBSEQ_ROOT/bin/parse_result.pl"|Path to the parsing script
`parseResults.alleleFile`|String|"$HLAVBSEQ_BWA_INDEX_ROOT/Allelelist.txt"|File with allele information
`parseResults.modules`|String|"hla-vbseq/1 hlavbseq-bwa-index/2.0"|Names and versions of required modules.
`callHlaDigits.jobMemory`|Int|8|Memory allocated to the task.
`callHlaDigits.timeout`|Int|20|Timeout in hours, needed to override imposed limits.
`callHlaDigits.meanReadLength`|Int|90|Mean read length, default is 90
`callHlaDigits.resolution`|Int|4|Resolution for HLA allele may be 4, 6 or 8
`callHlaDigits.callingScript`|String|"$HLA_VBSEQ_ROOT/bin/call_hla_digits.py"|Path to the HLA allele calling script
`callHlaDigits.alleleFile`|String|"$HLAVBSEQ_BWA_INDEX_ROOT/Allelelist.txt"|File with allele information
`callHlaDigits.modules`|String|"hla-vbseq/1 hlavbseq-bwa-index/2.0"|Names and versions of required modules.


### Outputs

Output | Type | Description | Labels
---|---|---|---
`hlaVbSeqresult`|File|Results as they come from the tool, all HLA alleles are listed|vidarr_label: hlaVbSeqresult
`parsedResults`|File|Parsed Results, HLA alleles with non-zero signal|vidarr_label: parsedResults
`resultsFile`|File|File with the final HLA Allele calls|vidarr_label: resultsFile


## Commands
This section lists command(s) run by HLA-VBseq workflow
 
* Running HLA-VBseq
 
At this point, workflow uses a filtering step for alignments (removing secondary hits) which is not a part of the standard way to run this tool. It may be resolved in a future, but since the authors did not publish the source code it may require direct communication with them to clarify some issues. The workflow produces reports as expected, though filtering of reads may distort the assessment of HLA alleles.
 
### Extracting reads which overlap HLA alleles:
 
```
  samtools view INPUT_BAM HLA_INTERVALS | awk '{print $1}' | sort | uniq > [OUTPUT_PREFIX]_HLA-VBSeq_reads.txt
 
```
 
### Indexing bam file with HLA-VBseq and extracting HLA reads into a .sam file
 
```
  set -euo pipefail
  unset _JAVA_OPTIONS
  mkdir data
  ln -s INPUT_BAM -t data
  java -Xmx18G -jar BAM_NAME_INDEX_JAR index data/FILE_NAME  --indexFile data/FILE_NAME.idx
  java -Xmx18G -jar BAM_NAME_INDEX_JAR search data/FILE_NAME --name HLA_READS --output [FILE_NAME]_partial.sam
```
 
### Extract reads into fastq files, both for HLA-specific and unmapped reads. Merge
 
```
  set -euo pipefail
  unset _JAVA_OPTIONS
  java -Xmx24G -jar picard.jar SamToFastq I=PARTIAL_SAM F=PARTIAL_1.fastq F2=PARTIAL_2.fastq
  samtools view -bh -f 12 INPUT_BAM > [FILE_NAME].sorted_unmapped.bam
  java -Xmx24G -jar picard.jar SamToFastq I=[FILE_NAME].sorted_unmapped.bam F=UNMAPPED_1.fastq F2=UNMAPPED_2.fastq
  cat PARTIAL_1.fastq UNMAPPED_1.fastq | gzip -c > [FILE_NAME]_part_1.fastq.gz
  cat PARTIAL_2.fastq UNMAPPED_2.fastq | gzip -c > [FILE_NAME]_part_2.fastq.gz
```
 
### Filter alignments (the default is to remove secondary alignments)
 
```
  samtools view -h -F FILTER_TAG INPUT_BAM -b > FILTERED_BAM
```
 
### Run the analysis with HLA-VBseq, call HLA alleles
 
```
  set -euo pipefail
  unset _JAVA_OPTIONS
  java -jar -Xmx12g -Xms12g HLAVBSEQ_JAR \
                            HLA_FASTA \
                            INPUT_BAM \
                            [FILE_NAME]_HLA-VBSeq_results.txt  \
                            --alpha_zero ALPHA_ZERO \
                            --is_paired 
```

### Parsing raw calls from HLA-VBseq
 
```
  ~{parsingScript} ALLELE_FILE \
                   RESULTS_FILE > \
                   [FILE_NAME]_HLA-VBSeq_prediction.txt
```
 
### Post-processing of the results
 
```
  python3 CALLING_SCRIPT \
      -v RESULTS_FILE \
      -a ALLELE_FILE \
      -r MEAN_READ_LENGTH \
      -d RESOLUTION \
      --ispaired > \
      [FILE_NAME]_HLA-VBSeq_call.txt
```

## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
