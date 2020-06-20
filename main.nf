#! /usr/bin/env nextflow

//vim: syntax=groovy -*- mode: groovy;-*-

// Copyright (C) 2020 IRB Barcelona

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

log.info ""
log.info "-------------------------------------------------------------------------"
log.info "  SNPs-ML: machine learning models for SNPs data                         "
log.info "-------------------------------------------------------------------------"
log.info "Copyright (C) IRB Barcelona"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "-------------------------------------------------------------------------"
log.info ""

params.help = null

if (params.help) {
    log.info ''
    log.info '--------------------------------------------------'
    log.info '  USAGE              '
    log.info '--------------------------------------------------'
    log.info ''
    log.info 'Usage: '
    log.info 'nextflow run main.nf --train_vcf train.vcf.bgz --apply_vcf apply.vcf.bgz'
    log.info ''
    log.info 'Mandatory arguments:'
    log.info '    --mc3_file             FILE      Input MC3 file from GDC.'
    log.info '    AND --most_mutated_genes  FILE  Tab delimited file containing for each cancer type the.'
    log.info ''
    log.info '    OR --mc3_folder        FOLDER    Input folder containing one MC3 file per cancer type.'
    log.info '                                     most mutated genes (2 columns).'
    log.info '    --SNPs                 FILE      Input txt file containing SNPs in column, samples in rows.'
    log.info '    --samples_table        FILE      Input txt file "ID" and "sample_id" columns.'

    log.info 'Optional arguments:'
    log.info '    --output_folder        FILE      Output folder (default: SNPs_ML_output).'
    log.info '    --mem                  INT       Memory used (default=8, in Gb).'
    log.info '    --cpu                  INT       Number of CPU used in SVM (default=2).'
    log.info 'Flags:'
    log.info '    --no_cancer_type                 Run the pipeline on all cancer types merged (more power)'
    log.info '    --help                           Display this message'
    log.info ''
    exit 0
}

params.SNPs = null
params.mc3_folder = null
params.most_mutated_genes = null
params.samples_table = null
params.output_folder = "oncogenes_PRS_output"
params.mem = 8
params.cpu = "2"

if(params.SNPs == null & params.mc3_folder == null 1 params.samples_table){
  exit 1, "Please specify each of the following parameters: --SNPs, --samples_table and --mc3_folder"
}

SNPs = file(params.SNPs)
samples_table = file(params.samples_table)
mc3_files = Channel.fromPath(params.mc3_folder + "/*.txt")

process translate {
  cpus params.cpu
  memory params.mem+'G'

  publishDir params.output_folder+"/translate", mode: 'copy'

  input:
  file table from mc3_files
  file SNPs
  file samples_table

  output:
  file "*translate.txt" into snps_genes_table

  shell:
  '''
  Rscript !{baseDir}/bin/translate_oncogene.r --SNPs=!{SNPs} --mut_table=!{table} --samples_table=!{samples_table} --nb_cpu=!{params.cpu}
  '''
}

process preselection_svm {
  cpus params.cpu
  memory params.mem+'G'

  tag {tag}

  publishDir params.output_folder+"/PRS", mode: 'copy'

  input:
  file snps_gene from snps_genes_table

  output:
  set val(tag), file("PRS_*.Rdata") into prs_out

  shell:
  tag=snps_gene.baseName.replace("_translate.txt","").replace("MC3_","")
  '''
  Rscript !{baseDir}/bin/PRS.r --translate_SNPs=!{snps_gene} --nb_cpu=!{params.cpu} --tag=!{tag}
  '''
}


