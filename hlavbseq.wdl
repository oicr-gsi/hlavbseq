version 1.0

import "imports/pull_bwaMem.wdl" as bwaMem

# ================================================================================
# Workflow accepts two fastq files for paired-end sequencing, with R1 and R2 reads
# ================================================================================
workflow hlavbseq {
input {
  File inputBam 
  File inputBai
  String outputFileNamePrefix = ""
  String add_params = "-M -P -L 10000 -a"
  Int filterTag = 256
  Int bwa_threads = 8
}

String sampleID = if outputFileNamePrefix=="" then basename(inputBam, ".bam") else outputFileNamePrefix

call extractReads {
  input:
    inputBam = inputBam,
    inputBai = inputBai,
    outputFileNamePrefix = sampleID
}

call indexBam {
  input:
    inputBam = inputBam,
    hlaReads = extractReads.extractedReads,
    outputFileNamePrefix = sampleID
}

call makeFastq {
  input:
    inputBam = inputBam,
    partialSam = indexBam.partialAlignments,
    hlaReads = extractReads.extractedReads,
    outputFileNamePrefix = sampleID
}

call bwaMem.bwaMem {
  input:
    fastqR1 = makeFastq.partialF1,
    fastqR2 = makeFastq.partialF2,
    runBwaMem_threads = bwa_threads,
    runBwaMem_addParam = add_params,
    outputFileNamePrefix = sampleID
}

if (filterTag > 0) {
  call filterBam {
    input:
      inputBam = bwaMem.bwaMemBam,
      filterTag = filterTag
  }
}

call runHlaVbSeq {
  input:
    inputBam = select_first([filterBam.filteredBam, bwaMem.bwaMemBam]),
    outputFileNamePrefix = outputFileNamePrefix
}

call parseResults {
  input:
    resultsFile = runHlaVbSeq.resultsFile,
    outputFileNamePrefix = outputFileNamePrefix
}

call callHlaDigits {
  input:
    resultsFile = runHlaVbSeq.resultsFile,
    outputFileNamePrefix = outputFileNamePrefix
}

parameter_meta {
  inputBam: "BWA aligned reads."
  filterTag: "To address the issue with some data, we use 256 (secondary alignment) as a default"
  runBwaMem_threads: "Threads to use with BWA aligner."
  runBwaMem_addParam: "HLAVBSeq uses additional parameters for BWA which need to be specified."
  outputFileNamePrefix: "Output prefix, customizable. Default is the first file's basename."
}

meta {
    author: "Peter Ruzanov"
    email: "peter.ruzanov@oicr.on.ca"
    description: "HLAVBSeq is software to estimate the most likely HLA types from high-throughput sequencing data."
    dependencies: [
      {
        name: "hla-vbseq/1",
        url: "http://nagasakilab.csml.org/hla/HLAVBSeq.jar"
      },
      {
        name: "bwa/0.7.17",
        url: "https://github.com/lh3/bwa/archive/0.7.17.tar.gz"
      }
    ]
    output_meta: {
      hlaVbSeqresult: "Results as they come from the tool, all HLA alleles are listed",
      parsedResults: "Parsed Results, HLA alleles with non-zero signal",
      resultsFile: "File with the final HLA Allele calls"
    }
}

output {
  File hlaVbSeqresult = runHlaVbSeq.resultsFile
  File parsedResults = parseResults.parsedResults
  File resultsFile = callHlaDigits.finalResults
}
}


# ===============================================
#  1 of 6: Extract reads specific to HLA alleles
# ===============================================
task extractReads {
input {
  Int  jobMemory = 8
  Int  timeout   = 20
  File inputBam
  File inputBai
  String hlaIntervals = "chr6:29907037-29915661 chr6:31319649-31326989 chr6:31234526-31241863 chr6:32914391-32922899 chr6:32900406-32910847 chr6:32969960-32979389 chr6:32778540-32786825 chr6:33030346-33050555 chr6:33041703-33059473 chr6:32603183-32613429 chr6:32707163-32716664 chr6:32625241-32636466 chr6:32721875-32733330 chr6:32405619-32414826 chr6:32544547-32559613 chr6:32518778-32554154 chr6:32483154-32559613 chr6:30455183-30463982 chr6:29689117-29699106 chr6:29792756-29800899 chr6:29793613-29978954 chr6:29855105-29979733 chr6:29892236-29899009 chr6:30225339-30236728 chr6:31369356-31385092 chr6:31460658-31480901 chr6:29766192-29772202 chr6:32810986-32823755 chr6:32779544-32808599 chr6:29756731-29767588"
  String outputFileNamePrefix
  String modules = "samtools/1.9"
}

command <<<
 samtools view ~{inputBam} ~{hlaIntervals} | awk '{print $1}' | sort | uniq > ~{outputFileNamePrefix}_HLA-VBSeq_reads.txt
>>>

parameter_meta {
 inputBam: "Input bam file, BWA-aligned reads"
 inputBai: "index of the Input bam file"
 hlaIntervals: "Intervals where HLA-alleles are sitting"
 outputFileNamePrefix: "Output prefix for the result file"
 jobMemory: "Memory allocated to the task."
 modules: "Names and versions of required modules."
 timeout: "Timeout in hours, needed to override imposed limits."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File extractedReads = "~{outputFileNamePrefix}_HLA-VBSeq_reads.txt"
}
}

# ==========================================================================
#  2 of 6: Index bam file using HLA-VBSeq tools and subset reads
# ==========================================================================
task indexBam {
input {
  Int jobMemory = 18
  Int overhead = 6
  Int timeout = 20
  File inputBam
  File hlaReads
  String bamNameIndexJar = "$HLA_VBSEQ_ROOT/bin/bamNameIndex.jar"
  String outputFileNamePrefix
  String modules = "hla-vbseq/1 hlavbseq-bwa-index/2.0"
}

Int javaMemory = jobMemory - overhead
String fileName = basename(inputBam)

command <<<
 set -euo pipefail
 unset _JAVA_OPTIONS
 mkdir data
 ln -s ~{inputBam} -t data
 java -Xmx~{javaMemory}G -jar ~{bamNameIndexJar} index data/~{fileName}  --indexFile data/~{fileName}.idx
 java -Xmx~{javaMemory}G -jar ~{bamNameIndexJar} search data/~{fileName} --name ~{hlaReads} --output ~{outputFileNamePrefix}_partial.sam
>>>

parameter_meta {
 inputBam: "Input bam file, BWA-aligned reads"
 hlaReads: "Reads to use (HLA-allele specific), output from the previous step"
 bamNameIndexJar: "Jar file for bamNameInder"
 outputFileNamePrefix: "Output prefix for the result file"
 jobMemory: "Memory allocated to the task."
 overhead: "Ovrerhead for calculating heap memory, difference between total and Java-allocated memory"
 modules: "Names and versions of required modules."
 timeout: "Timeout in hours, needed to override imposed limits."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File partialAlignments = "~{outputFileNamePrefix}_partial.sam"
}
}


# ====================================
#   3 of 6: Make fastq files
# ====================================
task makeFastq {
input {
  Int jobMemory = 24
  Int overhead = 6
  Int timeout = 20
  File inputBam
  File partialSam
  File hlaReads
  String bamNameIndexJar = "$HLA_VBSEQ_ROOT/bin/bamNameIndex.jar"
  String outputFileNamePrefix
  String modules = "samtools/1.9 picard/2.21.2"
}

Int javaMemory = jobMemory - overhead

command <<<
 set -euo pipefail
 unset _JAVA_OPTIONS
 java -Xmx~{javaMemory}G -jar $PICARD_ROOT/picard.jar SamToFastq I=~{partialSam} F=PARTIAL_1.fastq F2=PARTIAL_2.fastq
 samtools view -bh -f 12 ~{inputBam} > ~{outputFileNamePrefix}.sorted_unmapped.bam
 java -Xmx~{javaMemory}G -jar $PICARD_ROOT/picard.jar SamToFastq I=~{outputFileNamePrefix}.sorted_unmapped.bam F=UNMAPPED_1.fastq F2=UNMAPPED_2.fastq
 cat PARTIAL_1.fastq UNMAPPED_1.fastq | gzip -c > ~{outputFileNamePrefix}_part_1.fastq.gz
 cat PARTIAL_2.fastq UNMAPPED_2.fastq | gzip -c > ~{outputFileNamePrefix}_part_2.fastq.gz

>>>

parameter_meta {
 inputBam: "Input bam file, BWA-aligned reads"
 hlaReads: "Reads to use (HLA-allele specific), output from the previous step"
 bamNameIndexJar: "Jar file for bamNameInder"
 outputFileNamePrefix: "Output prefix for the result file"
 jobMemory: "Memory allocated to the task."
 overhead: "Ovrerhead for calculating heap memory, difference between total and Java-allocated memory"
 modules: "Names and versions of required modules."
 timeout: "Timeout in hours, needed to override imposed limits."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File partialF1 = "~{outputFileNamePrefix}_part_1.fastq.gz"
  File partialF2 = "~{outputFileNamePrefix}_part_2.fastq.gz"
}
}

# ==================================================
#  3a of 6: FILTERING (OPTIONAL, ENABLED BY DEFAULT)
# ==================================================
task filterBam {
input {
  Int jobMemory = 12
  Int timeout   = 12
  Int filterTag = 256
  File inputBam
  String modules = "samtools/1.9"
}

String filteredBam = basename(inputBam, ".bam") + "filtered.bam"

command <<<
 samtools view -h -F ~{filterTag} ~{inputBam} -b > ~{filteredBam}
>>>

parameter_meta {
 jobMemory: "Memory allocated to the task."
 modules: "Names and versions of required modules."
 timeout: "Timeout in hours, needed to override imposed limits."
 inputBam: "Input bam file from the alignment with BWA"
 filterTag: "Binary Tag to use for filtering"
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File filteredBam = basename(inputBam, ".bam") + ".filtered.bam"
}
}


# =========================
#  4 of 6: RUNNING THE TOOL
# =========================
task runHlaVbSeq {
input {
  Int jobMemory = 32
  Int timeout   = 20
  Int overhead = 6
  File inputBam
  Float alphaZero = 0.01
  String hlavbseqJar = "$HLA_VBSEQ_ROOT/bin/HLAVBSeq.jar"
  String hlaFasta = "$HLAVBSEQ_BWA_INDEX_ROOT/hla_all_v2.fasta"
  String outputFileNamePrefix
  String modules = "hla-vbseq/1 hlavbseq-bwa-index/2.0"
}

Int javaMemory = jobMemory - overhead

command <<<
 set -euo pipefail
 unset _JAVA_OPTIONS
 java -jar -Xmx~{javaMemory}g -Xms~{javaMemory}g ~{hlavbseqJar} \
                           ~{hlaFasta} \
                           ~{inputBam} \
                           ~{outputFileNamePrefix}_HLA-VBSeq_results.txt  \
                           --alpha_zero ~{alphaZero} \
                           --is_paired 
>>>

parameter_meta {
 jobMemory: "Memory allocated to the task."
 modules: "Names and versions of required modules."
 timeout: "Timeout in hours, needed to override imposed limits."
 overhead: "Difference between memory allocated and memory used as Java heap"
 inputBam: "Input bam file from the alignment with BWA"
 alphaZero: "Hyperparameter, as described in the paper it is not recommended to change the default (0.01)"
 hlavbseqJar: "path to HLAVBSeq jar file"
 hlaFasta: "fasta with HLA sequences"
 outputFileNamePrefix: "Output prefix for the result file"
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File resultsFile = "~{outputFileNamePrefix}_HLA-VBSeq_results.txt"
}
}


# =========================
#  5 of 6: PARSING RESULTS
# =========================
task parseResults {
input {
  Int  jobMemory = 8
  Int  timeout   = 20
  File resultsFile
  String parsingScript = "$HLA_VBSEQ_ROOT/bin/parse_result.pl"
  String alleleFile = "$HLAVBSEQ_BWA_INDEX_ROOT/Allelelist.txt"
  String outputFileNamePrefix
  String modules = "hla-vbseq/1 hlavbseq-bwa-index/2.0"
}

command <<<
 ~{parsingScript} ~{alleleFile} \
                  ~{resultsFile} > \
                  ~{outputFileNamePrefix}_HLA-VBSeq_prediction.txt
>>>

parameter_meta {
 resultsFile: "Input is the results file from the previous step"
 parsingScript: "Path to the parsing script"
 alleleFile: "File with allele information" 
 outputFileNamePrefix: "Output prefix for the result file"
 jobMemory: "Memory allocated to the task."
 modules: "Names and versions of required modules."
 timeout: "Timeout in hours, needed to override imposed limits."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File parsedResults = "~{outputFileNamePrefix}_HLA-VBSeq_prediction.txt"
}
}


# ===========================================
#  6 of 6: Quantify HLA Allele Representation
# ===========================================
task callHlaDigits {
input {
  Int  jobMemory = 8
  Int  timeout   = 20
  Int meanReadLength = 90
  Int resolution = 4
  File resultsFile
  String callingScript = "$HLA_VBSEQ_ROOT/bin/call_hla_digits.py"
  String alleleFile = "$HLAVBSEQ_BWA_INDEX_ROOT/Allelelist.txt"
  String outputFileNamePrefix
  String modules = "hlavbseq/1 hlaminer-bwa-index/0.7.17"
}

command <<<
 python3 ~{callingScript} \
     -v ~{resultsFile} \
     -a ~{alleleFile} \
     -r ~{meanReadLength} \
     -d ~{resolution} \
     --ispaired > \
     ~{outputFileNamePrefix}_HLA-VBSeq_call.txt
>>>

parameter_meta {
 resultsFile: "Input is the results file from the previous step"
 callingScript: "Path to the HLA allele calling script"
 alleleFile: "File with allele information"
 outputFileNamePrefix: "Output prefix for the result file"
 meanReadLength: "Mean read length, default is 90"
 resolution: "Resolution for HLA allele may be 4, 6 or 8"
 jobMemory: "Memory allocated to the task."
 modules: "Names and versions of required modules."
 timeout: "Timeout in hours, needed to override imposed limits."
}

runtime {
  memory:  "~{jobMemory} GB"
  modules: "~{modules}"
  timeout: "~{timeout}"
}

output {
  File finalResults = "~{outputFileNamePrefix}_HLA-VBSeq_call.txt"
}
}
