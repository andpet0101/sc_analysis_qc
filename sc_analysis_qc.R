#!/usr/bin/env Rscript 

suppressPackageStartupMessages(library("optparse"))

# get script dir and working dir
arguments = commandArgs(trailingOnly = F)
script_index = grep("--file",arguments)
script_dir = normalizePath(dirname(sub("--file=","",arguments[script_index])))
script = basename(sub("--file=","",arguments[script_index]))
def_working_dir = normalizePath(getwd())

# get rid of leading arguments
arguments_length = length(arguments)
arguments_index = grep("--args",arguments)[1]
if(is.na(arguments_index)){
  arguments=c()
} else{
  arguments_index=arguments_index+1
  arguments = arguments[arguments_index:arguments_length]
}

# parse options
parser = OptionParser(usage = "usage: %prog [options] dataset1 [dataset2 ... ]",prog=script,
                      description="Assesses the quality of single cell data. One or more datasets can analysed where a dataset can be a featureCounts file, a kallisto directory or a 10x directory. When multiple datasets are analysed, dataset arguments can be named with a prefix 'name:'. Otherwise, the datasets will be prefixed with 'd1', 'd2', ... .")
parser = add_option(parser, c("-o", "--output"), action="store",type="character",dest="output_prefix",help="Output prefix path for the file.")
parser = add_option(parser, c("-b", "--bfx"), action="store",type="character",dest="bfx",help="BFX id of the project")
parser = add_option(parser, c("-a", "--annotation"), action="store",type="character",dest="annotation",help="Either Ensembl annotation version or a file with custom annotations which must contain at the least the columns feature_id and feature_symbol")
parser = add_option(parser, c("-s", "--species"), action="store",type="character",dest="species",help="Species name. When using Ensembl annotation, the name must be Ensembl-compatible.")
parser = add_option(parser, c("-w", "--workingdir"), action="store",type="character",dest="workingdir",default=def_working_dir,help="Working directory [%default]")
parser = add_option(parser, c("-r", "--source_dir"), action="store",type="character",dest="sourcedir",default=script_dir,help="Directory with the R and Rmd code [%default]")

parsed_args = parse_args(parser, args = arguments,positional_arguments = TRUE)

if(!"output_prefix" %in% names(parsed_args$options)) stop(paste0(script,": Please specify a output prefix!"))
if(!"bfx" %in% names(parsed_args$options)) stop(paste0(script,": Please specify a bfx project!"))
if(!"annotation" %in% names(parsed_args$options)) stop(paste0(script,": Please specify an annotation!"))
if(!"species" %in% names(parsed_args$options)) stop(paste0(script,": Please specify a species!"))

# check validity
datasets_arg = parsed_args$args
if(length(datasets_arg)==0) stop(paste0(script,": Please specify at least one dataset!"))

source_dir = normalizePath(file.path(normalizePath(parsed_args$options$sourcedir)),mustWork = FALSE)
if(!dir.exists(source_dir)) stop(paste0(script,": The source directory ",source_dir," with the R and Rmd files does not exist!"))

main_rmd_file = file.path(source_dir,'Rmd','ha.Rmd')
if(!file.exists(main_rmd_file)) stop(paste0(script,": The main Rmd file ",main_rmd_file," for rendering the report does not exist!"))

working_dir = normalizePath(file.path(parsed_args$options$workingdir),mustWork = FALSE)
if(!dir.exists(working_dir)) stop(paste0(script,": The working directory ",working_dir," does not exist!"))

output_prefix = normalizePath(file.path(parsed_args$options$output_prefix),mustWork = FALSE)
if(!dir.exists(dirname(output_prefix))) stop(paste0(script,": Could not find the directory for ",output_prefix,"!"))

list_of_datasets = c()
i = 1
if(length(datasets_arg)>1){
  for(a in strsplit(datasets_arg,":")){
    if(length(a)>1){
      list_of_datasets[a[1]] = a[2]
    }else{
      list_of_datasets[paste0("d",i)] = a[1]
      i = i + 1
    }
  }
}

print(list_of_datasets)

#rmarkdown::render(main_rmd_file,
#                  output_format = "pdf",
#                  output_file = basename(output_prefix),
#                  output_dir = dirname(output_prefix),
#                  intermediates_dir = working_dir,
#                  knit_root_dir = working_dir,
#                  clean=TRUE,
#                  run_pandoc=TRUE,
#                  quiet=TRUE,
#                  params = list(
#                    bfx = parsed_args$options$bfx,
#                    annotation = parsed_args$options$annotation,
#                    species = parsed_args$options$species,
#                    workingdir = working_dir,
#                    sourcedir = 
#                    )
#))