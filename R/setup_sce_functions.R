#' Creates a sce object from one or more countsdata objects.
#' 
#' @param countsdata  A counts data object. Can be featureCounts file, 10x matrix directory or kallisto directory. A list or named list can be given in order to load several countsdata  objects. If a named vector is given, the cell names will be prefixed with the names.
#' @param drop_zero_genes Drop gene rows with zero counts.
#' @param use_sparse_matrix Use a sparse matrix for the counts.
#' @return A SingleCellExperiment object.
#' @examples
#' sce = create_sce_from_countsdata("featurecounts_file.txt")
#' sce = create_sce_from_countsdata(list(Data1 ="featurecounts_file.txt",Data2="10x_matrix_dir",Data3="kallisto_dir"))
create_sce_from_countsdata = function(countsdata,drop_zero_genes=T,use_sparse_matrix=F){
  nsets = length(countsdata)
  full_data = vector("list", nsets)
  countobject_type = c()
  countobject_name = c()
  
  # go through countsdata objects
  for(i in seq_len(nsets)){
    file_info = file.info(countsdata[[i]])
    if(is.na(file_info$size)) flog.error("The countsdata object %s does not exist!",countsdata[[i]])
    
    # load countsdata object
    if(file_info$isdir){
      countdata_mtx = file.path(matrix_dir[[i]],'matrix.mtx')

      if(file.exists(countdata_mtx)){
        # 10X
        flog.info("Countsdata object %s is a directory with a 'matrix.mtx' file and therefore interpreted as 10x data.",countsdata[[i]])
        full_data[[i]] = create_matrix_from_10x(countsdata[[i]],use_sparse_matrix=use_sparse_matrix)
        data_types[[i]] = "10X"
        
      }else{
        # kallisto
        flog.info("Countsdata object %s is a directory without a 'matrix.mtx' file and therefore interpreted as kallisto data.",countsdata[[i]])
        full_data[[i]] = create_matrix_from_kallisto(countsdata[[i]],use_sparse_matrix=use_sparse_matrix)
        data_types[[i]] = "kallisto"
      }
      
    }else{
      # featureCounts
      flog.info("Countsdata object %s is a file and therefore interpreted as featureCounts data.",countsdata[[i]])
      full_data[[i]] = create_matrix_from_featurecounts(countsdata[[i]],use_sparse_matrix=use_sparse_matrix)
      data_types[[i]] = "featureCounts"
    }
    
    # add sample prefix (by countsdata object)
    if(!is.null(names(countsdata))){
      mod_sample_names = paste(names(countsdata)[[i]],colnames(full_data[[i]]),sep="-")
    }else{
      
    }
      
      colnames(full_data[[i]]) = paste(names(countsdata)[[i]],colnames(full_data[[i]]),sep="-")
    
    
  }
  
  # drop zero rows
  nonzero_rows = unique(unlist(lapply(full_data,function(x){rownames(x)[rowSums(x)>0]})))
  full_data = lapply(my_l,function(x){x[nonzero_rows,]})
  
  # add missing data as zero counts and reorder rows 
  all_rownames = sort(unique(unlist(lapply(full_data,row.names))))
  
  expand_full_matrix = function(matrix,all_row_names){
    missing_row_names = all_row_names[!all_row_names %in% rownames(matrix)]
    full_matrix = rbind(matrix,Matrix(0,nrow=length(missing_row_names),ncol=ncol(matrix),dimnames = list(missing_row_names,colnames(matrix)),sparse=class(matrix)=="dgCMatrix"))
    full_matrix[all_row_names,]
  }
  
  full_data = lapply(full_data,expand_full_matrix,all_rownames)
  
  # now combine in large matrix
  full_data = do.call(cbind, full_data)
  if(any(duplicated(colnames(countstable)))) flog.fatal("Multiple cells have the same name. Please rename the cells or - when using multiple datasets with same cell names - use a prefix for each dataset.")
  
  # create sce object
  sce = SingleCellExperiment(list(counts = full_data))
  
  # 
  
  # clear countstable and full_data and call gc to free memory
  rm(full_data)
  gc(verbose=F)

  sce  
}

#' Creates a matrix object from featureCounts table.
#' 
#' @param featurecounts_file A featureCounts file.
#' @param drop_zero_genes Drop gene rows with zero counts.
#' @param use_sparse_matrix Use a sparse matrix for counts.
#' @return A matrix.
#' #' @seealso create_matrix_from_10x,create_matrix_from_kallisto
#' @examples
#' counts_mtx = create_matrix_from_featurecounts(featurecounts_file)
create_matrix_from_featurecounts = function(featurecounts_file,drop_zero_genes=T,use_sparse_matrix=F){
  require(readr)
  require(Matrix)
  require(futile.logger)
  
  if(!file.exists(featurecounts_file)) flog.error("The featureCounts file %s does not exist!",featurecounts_file)
  
  countdata = read_table2(featurecounts_file,col_names=T,comment = "#")
  if(colnames(countdata)[1]!="Geneid") flog.error("First column in featureCounts file %s must be named 'Geneid'!",featurecounts_file)
  rownames(countdata) = countdata$Geneid
  
  countdata$Geneid = NULL
  countdata$Chr = NULL
  countdata$Start = NULL
  countdata$End = NULL
  countdata$Strand = NULL
  countdata$Length = NULL
  
  countdata = as.matrix(countdata)
  if(use_sparse_matrix) countdata = as(countdata, "dgCMatrix")
  
  if(drop_zero_genes){
    keep = rowSums(countdata > 0) > 0
    countdata = countdata[keep, ]
  }
  
  gc(verbose=F)

  return(sce)  
}

#' Creates a matrix object from 10X data.
#' 
#' @param matrix_dir A directory containing the matrix.mtx, genes.tsv, and barcodes.tsv files provided by 10x.
#' @param drop_zero_genes Drop gene rows with zero counts.
#' @param use_sparse_matrix Use a sparse matrix for counts.
#' @return A matrix.
#' #' @seealso create_matrix_from_featurecounts,create_matrix_from_kallisto
#' @examples
#' sce = create_matrix_from_10x(matrix_dir_from_cellranger)
create_matrix_from_10x = function(matrix_dir,drop_zero_genes=T){
  require(Matrix)
  require(futile.logger)

  
  if(!dir.exists(matrix_dir)) flog.error("Specified matrix directory",matrix_dir,"does not exist!")
  countdata_mtx = file.path(matrix_dir,'matrix.mtx')
  countdata_cols = file.path(matrix_dir,'barcodes.tsv')
  countdata_rows = file.path(matrix_dir,'genes.tsv')
  
  if(!file.exists(countdata_mtx) | !file.exists(countdata_cols) | file.exists(countdata_rows)) flog.error("Either matrix.mtx or barcodes.tsv or genes.tsv is missing for matrix directory ",matrix_dir,"!")
  
  countdata = readMM(countdata_mtx)
  crows = readLines(countdata_rows)[1]
  rownames(countdata) = crows
  ccols = readLines(countdata_cols)[1]
  colnames(countdata) = ccols
  
  countdata = as.matrix(countdata)
  if(!use_sparse_matrix) countdata = as.matrix(countdata)
  
  if(drop_zero_genes){
    keep = rowSums(countdata > 0) > 0
    countdata = countdata[keep, ]
  }
  
  gc(verbose=F)
}
  
  for(d in matrix_dir){
    if(!file.exists(d)) stop(paste("Specified matrix directory",d,"does not exist!"))
    countdata_mtx = file.path(d,'matrix.mtx')
    countdata_cols = file.path(d,'barcodes.tsv')
    countdata_rows = file.path(d,'genes.tsv')
    
  }
  
  nsets = length(matrix_dir)
  full_data = vector("list", nsets)
  for(i in seq_len(nsets)){
    countdata_mtx = file.path(matrix_dir[[i]],'matrix.mtx')
    countdata_cols = file.path(matrix_dir[[i]],'barcodes.tsv')
    countdata_rows = file.path(matrix_dir[[i]],'genes.tsv')
    
    if(!file.exists(countdata_mtx) | !file.exists(countdata_cols) | file.exists(countdata_rows)){
      stop(paste("Either matrix.mtx or barcodes.tsv or genes.tsv is missing for",matrix_dir[[i]],"!"))
    }
    
    full_data[[i]] = readMM(countdata_mtx)
    rownames(full_data[[i]]) = read.table(countdata_rows,sep="\t",header=F,stringsAsFactors = F)$V1
    colnames(full_data[[i]]) = read.table(countdata_cols,sep="\t",header=F,stringsAsFactors = F)$V1
    
    if(use_sparse_matrix) full_data[[i]] = as(full_data[[i]], "dgCMatrix")
    if(!is.null(names(matrix_dir))) colnames(full_data[[i]]) = paste(names(matrix_dir)[[i]],colnames(full_data[[i]]),sep="-")
  
  }
  
  # combine all datasets in one table
  if(any(sapply(full_data,function(x){any(rownames(full_data[[1]])!=rownames(x))}))) stop("Gene information is not the same for all 10x matrix dirs!")
  
  countstable = do.call(cbind, full_data)
  if(any(duplicated(colnames(countstable)))) stop("Multiple cells have the same name. Please rename the cells or - when using multiple datasets with same cell names - use a prefix for each dataset.")
  
  # drop fully zero genes
  if(drop_zero_genes){
    keep = rowSums(countstable > 0) > 0
    countstable = countstable[keep, ]
  }
  
  # drop fully zero genes
  if(drop_zero_genes){
    keep = rowSumscounts(sce) > 0
    sce = sce[keep, ]
  }
  
  
  return(sce)
}

#' Creates a matrix object from kallisto.
#' 
#' @param kallisto_dir A directory containing the kallisto output (one subdirectory per sample). A vector or named vector can be given in order to load several data directories. If a named vector is given, the cell names will be prefixed with the names.
#' @param drop_zero_genes drop gene rows with zero counts
#' @param use_sparse_matrix use a sparse matrix for counts
#' @return a SingleCellExperiment
#' #' @seealso create_sce_from_featurecounts,create_sce_from_10xData
#' @examples
#' sce = create_sce_from_featurecounts(featurecounts_file)
create_sce_from_kallisto = function(kallisto_dir,drop_zero_genes=T,use_sparse_matrix=F){
  require(scater)
  require(Matrix)
  
  # read featureCounts output
  nsets = length(kallisto_dir)
  full_data = vector("list", nsets)
  for(i in seq_len(nsets)){
    if(!file.exists(kallisto_dir[[i]])) stop(paste("Specified kallisto directory ",kallisto_dir[[i]],"does not exist!"))
    
    subdirs = list.dirs(path=kallisto_dir[[i]],full.names=T)
    kallisto_samples = c()
    kallisto_subdirs = c()
    for(l in subdirs){
      if(l==kallisto_dir[[i]]){next}
      
      kallisto_h5 = file.path(l,'abundance.h5')
      kallisto_tsv = file.path(l,'abundance.tsv')
      kallisto_run_info = file.path(l,'run_info.json')
      if(!file.exists(kallisto_h5) | !file.exists(kallisto_tsv) | !file.exists(kallisto_run_info)){next}
      
      kallisto_samples = c(kallisto_samples,basename(l))
      kallisto_subdirs = c(kallisto_subdirs,l)
    }
    
    if(length(kallisto_samples)==0 | length(kallisto_subdirs)==0)stop(paste("No kallisto data found at directory",kallisto_dir[[i]],"!"))
    
    kallisto_sce = readKallistoResults(samples = kallisto_samples,directories = kallisto_subdirs,logExprsOffset = 1)
  }

  
 
  # clear countstable and call gc to free memory
  gc(verbose=F)
}


#' Parses plate information from cell names. Names must correspond to Sample_PlateNumber_RowCol where 'Sample' is the actual name of the Sample, 'PlateNumber' is a digit, 'Row' is a letter between A-H and 'Col' is a digit between 1-12 (96 well plate). For 384 well plates, 'Row' is between A-P and 'Col' between is between 1-24. This information will be added to the colData of the SingleCellExperiment.
#' 
#' @param sce A SingleCellExperiment.
#' @return A SingleCellExperiment with additional information about Sample, PlateCoords, PlateNumber, Row and Col.
#' @examples
#' sce = parse_plate_information(sce)
parse_plate_information = function(sce){
  sample_info = as.data.frame(colData(sce))
  samples_with_plate_info = grepl("_\\d\\d\\d*_[A-Z]\\d\\d$",rownames(sample_info))

  # when we have plate_info, add PlateNumber, Row and Col information, otherwise just add the column Sample
  if(any(samples_with_plate_info)){
    library(tidyr)
    if(!all(samples_with_plate_info)) stop("Some samples do not have correct plate information (Name_PlateNumber_RowCol)!")
  
    sample_info$Row.names = rownames(sample_info)
    sample_info = separate(sample_info,Row.names,c("Sample","PlateNumberTxt","PlateCoords"),sep=c(-7,-4),remove=F)
    sample_info$Sample = gsub("(^_)|(_$)","",sample_info$Sample)
    sample_info$PlateCoords = gsub("(^_)|(_$)","",sample_info$PlateCoords)
    sample_info$PlateNumberTxt = as.integer(gsub("(^_)|(_$)","",sample_info$PlateNumberTxt))
    sample_info$PlateNumber = factor(as.integer(sample_info$PlateNumberTxt),ordered=T,levels=sort(unique(as.integer(sample_info$PlateNumberTxt))))
    sample_info$PlateNumberTxt = NULL
    sample_info$Row = gsub("[0-9]+$","",sample_info$PlateCoords)
    
    if(any(!sample_info$Row %in% c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"))){
      stop("Some samples have invalid row information!")
    }
    
    if(any(!sample_info$Col>24)){
      stop("Some samples have invalid column information!")
    }
    
    # 384 vs 96er plate format
    if("I" %in% sample_info$Row){
      sample_info$Row = factor(sample_info$Row,levels=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"),ordered=T)
    }else{
      sample_info$Row = factor(sample_info$Row,levels=c("A","B","C","D","E","F","G","H"),ordered=T)
    }
    sample_info$Col = as.integer(gsub("^[A-Z]","",sample_info$PlateCoords))
    if(13 %in%  sample_info$Col){
      sample_info$Col = factor(sample_info$Col,levels=1:24,ordered=T)
    }else{
      sample_info$Col = factor(sample_info$Col,levels=1:12,ordered=T)
    }
    
    sample_info$Row.names = NULL
    sample_info$PlateCoords = NULL
  }
  
  colData(sce) = DataFrame(sample_info,row.names = rownames(sample_info))
  
  return(sce)
}

#' Assigns group and control information. Cells can belong to different groups (e.g. tissues) or to positive, negative and bulk controls. This information can be parsed from the name or provided by the user.
#' 
#' @param sce A SingleCellExperiment.
#' @param group_information A named vector assigning the cells to a specific group or to 'PosCtrl', 'NegCtrl' and 'BulkCtrl'.
#' @return A SingleCellExperiment with additional information about Group, ControlType, IsControl and ControlLabel.
#' @examples
#' sce = assign_group_information(sce)
assign_group_information = function(sce,group_information=NULL){
  sample_info = as.data.frame(colData(sce))
  
  # column sample gives the cell name without plate information (if available); otherwise row names are sufficient 
  if("Sample" %in% colnames(sample_info)){ sample_info$Cell__names = sample_info$Sample
  }else{sample_info$Cell__names = rownames(sample_info)}
  
  # provided by user
  if(!is.null(group_information)){
    if(any(!sample_info$Cell__names %in% names(group_information))){
      stop("Some cell names are not in the user-provided group information!")
    }
    sample_info$Group = group_information[sample_info$Cell__names]
    sample_info$ControlType = factor(ifelse(!sample_info$Group %in% c("NegCtrl","PosCtrl","BulkCtrl"),"Sample",sample_info$Group),
                                     levels=c("Sample","BulkCtrl","PosCtrl","NegCtrl"))
    sample_info$IsControl = !sample_info$ControlType=="Sample"
  }else{
    # Is a sample a control? Find out based on sample name
    sample_info$ControlType = factor(ifelse(grepl("^NC",sample_info$Cell__names),"NegCtrl",
                                          ifelse(grepl("^(RNA|UHR)",sample_info$Cell__names),"PosCtrl",
                                                 ifelse(grepl("(^B-|\\d\\d$)",sample_info$Cell__names),"BulkCtrl","Sample"))),
                                          levels=c("Sample","BulkCtrl","PosCtrl","NegCtrl"))
    sample_info$IsControl = !sample_info$ControlType=="Sample"
    sample_info$Group = ifelse(sample_info$IsControl,as.character(sample_info$ControlType),as.character(sample_info$Cell__names))
    sample_info$Group = gsub("\\.\\d$","",sample_info$Group)
  }
  
  # add ControlLabel
  sample_info$ControlLabel = factor(ifelse(sample_info$ControlType=="Sample","S",
                                      ifelse(sample_info$ControlType=="BulkCtrl","B",
                                             ifelse(sample_info$ControlType=="PosCtrl","P","N"))),
                               levels=c("S","B","P","N"))
  
  # finalise group
  group_levels = unique(sample_info$Group)
  group_levels = c(sort(group_levels[!group_levels %in% c("BulkCtrl","PosCtrl","NegCtrl")]),c("BulkCtrl","PosCtrl","NegCtrl"))
  sample_info$Group = factor(sample_info$Group,levels=group_levels)

  sample_info$Cell__names= NULL
  colData(sce) = DataFrame(sample_info,row.names = rownames(sample_info))
  
  return(sce)
}

#' Fetches sample and library information from the LabWeb. Uses the full cell names as well as a bfxid for queries.
#' 
#' Requires RMySQL as well, an installed mysql client and an access to the database. Otherwise, the respective columns will contain NA.
#' 
#' @param sce A SingleCellExperiment.
#' @param bfxid The bfx id to which all cells must be assigned.
#' @param db_conf_file A MySQL configuration file with 'user', 'pass', 'db' and 'host'
#' @return A SingleCellExperiment with additional information about samples and libraries
#' @examples
#' sce = fetch_sample_library_info_from_labweb(sce,'bfx123')
fetch_sample_library_info_from_labweb = function(sce,bfxid=NULL,db_conf_file="/group/sequencing/Bfx/scripts/andreas/perl/db_sampledb1_config.txt"){
  sample_info = as.data.frame(colData(sce))
  sample_info$Cell__names = rownames(sample_info)

  if(is.null(bfxid)){
    stop("The bfxid is missing for LabWeb query!")
  }
  bfxid = gsub("^bfx","",bfxid)
  
  # query LabWeb database
  if("RMySQL" %in% rownames(installed.packages())){
    library(RMySQL)
    
    # parse and check conf file
    if(!file.exists(db_conf_file)){
      stop(paste("MySQL configuration file",db_conf_file,"does not exist!"))
    }
    db_conf = readLines(db_conf_file) %>% grep ("=",.,value=T) %>% strsplit(.,split="=")
    db_conf = setNames(trimws(sapply(db_conf, `[[`, 2)),trimws(sapply(db_conf, `[[`, 1)))
    if(any(!c('dbi','host','user','pass','db') %in% names(db_conf))){
      stop(paste("MySQL configuration file",db_conf_file," is missing 'dbi', 'host', 'user','pass' or 'db' entry!"))
    }
    
    # connect to db and query
    if('port' %in% db_conf){
      mydb = dbConnect(MySQL(), user=db_conf[['user']],password=db_conf[['pass']], dbname=db_conf[['db']], host=db_conf[['host']],port=db_conf[['port']])
    }else{
      mydb = dbConnect(MySQL(), user=db_conf[['user']],password=db_conf[['pass']], dbname=db_conf[['db']], host=db_conf[['host']])
    }
    
    sqlQuery = paste0("SELECT Libraries.ID AS 'LibraryID',Libraries.NAME AS 'LibraryName',Libraries.CONCENTRATION AS 'LibraryConcentration',",
              "Samples.ID AS 'SampleID',Samples.NAME AS 'SampleName',Samples.CONC AS 'SampleConcentration' ",
                    "FROM BfxProjects_Libraries ",
                    "JOIN Libraries ON BfxProjects_Libraries.LIBRARY_ID = Libraries.ID ",
                    "JOIN Samples ON Libraries.SAMPLE_ID = Samples.ID ",
                    "LEFT JOIN TissueSamples ON Samples.TISSUE_SAMPLE_ID = TissueSamples.ID ",
                    "LEFT JOIN Tissues ON TissueSamples.TISSUE_ID = Tissues.ID ",
                    "WHERE BfxProjects_Libraries.BFXPROJECT_ID = ",bfxid)
    sample_library_info = dbGetQuery(mydb, sqlQuery)
    dbDisconnect(mydb)
    
    if(nrow(sample_library_info)==0){
      stop("Fetching sample/library information from the LabWeb was not successfull!")
    }
    
  }else{
    warning("RMySQL not installed. All columns will be set to NA.")
  }
  
  # match db sample names to cell names
  # a) chromium: many cells belong to one sample, in this case
  if(any(grepl("[ACGT]{10,}",sample_info$Cell__names))){
    
  }
  
  # b) Smartseq2: one cell belongs to one sample
  
  
  # check if we have information for all samples
  if(all(sample_info$Cell__names %in% rownames(sample_library_concs))){
    sample_library_concs$LibraryConcentration = as.numeric(sample_library_concs$LibraryConcentration)
    sample_info$LibraryConcentration = sample_library_concs[sample_info$Cell__names,"LibraryConcentration"]
  
    sample_library_concs$SampleConcentration = as.numeric(sample_library_concs$SampleConcentration)
    sample_info$SampleConcentration = sample_library_concs[sample_info$Cell__names,"SampleConcentration"]
  }else{
    warning("Cannot find all samples in LabWeb database. Sample/library concentrations will be set NA.")
    sample_info$SampleConcentration = NA
    sample_library_concs$LibraryConcentration = NA
  }
  
  sample_info$Cell__names = NULL
  
  colData(sce) = DataFrame(sample_info,row.names=rownames(sample_info))
  
  return(sce)
}

#' Adds additional gene annotation fetched from Ensembl or from a user-provided file. Uses the rownames as query.
#' 
#' 
#' @param sce A SingleCellExperiment.
#' @param annotation Either a ensembl version number (when fetching from Ensembl) or a user-provided tab-separated file which must contain at least the columns 'feature_id' (the gene id) and 'feature_symbol' (the gene symbol).
#' @param species A species name supported by Ensembl (when fetching from Ensembl).
#' @return A SingleCellExperiment with additional information about genes
#' @examples
#' sce = add_gene_annotation(sce,species='mus_musculus',annotation='87)
#' sce = add_gene_annotation(sce,annotation='myannots.csv')
add_gene_annotation = function(sce,species=NULL,annotation=NULL){
  row_info = as.data.frame(rowData(sce))
  row_info$orig__row__num = 1:nrow(row_info)
  rownames(row_info) = rownames(sce)
  
  # user-provided annotation
  if(!is.null(annotation) && !grepl("^\\d+$",annotation)){
    if(!file.exists(annotation)){
      stop(paste("The user-provided annotation file",annotation,"does not exists!"))
    }
    
    # annotation file: tab-separated with row names (gene ids) and column names (gene attributes)
    # at least two defined columns: feature_id and feature_symbol
    annotation_table = read.table(annotation,header=T,comment.char = "#",check.names = F,sep="\t")
    
    # need the columns feature_id and feature_symbol
    if(sum(c("feature_id","feature_symbol") %in% colnames(annotation_table))!=2)stop("The annotation table must contain the columns 'feature_id' and 'feature_symbol'!")
    rownames(annotation_table) = annotation_table$feature_id
    
    # add to rowData
    row_info = merge(row_info,annotation_table,by="row.names",all.x=T)
    rownames(row_info) = row_info$Row.names
    row_info = row_info[order(row_info$orig__row__num),]
    row_info$orig__row__num = NULL
    row_info$Row.names = NULL
    if(any(rownames(row_info)!=rownames(sce))) stop("Annotation table could not be merged to sce object!")
    rowData(sce) = DataFrame(row_info,row.names=rownames(row_info))
  }else if(!is.null(annotation) && grepl("^\\d+$",annotation)){
  # Ensembl annotation
    require(biomaRt)
    
    # check version 
    most_recent_annotation = as.integer(gsub(".*\\s(\\d+)$","\\1",listEnsembl()[1,2]))
    ensembl_archives = as.data.frame(listEnsemblArchives())
    ensembl_archives$version = gsub("Ensembl ","",ensembl_archives$version)
    if(!annotation %in% ensembl_archives$version)stop(paste("Ensembl version '",annotation,"' not available on server!"))
    
    # check host, if not most recent then use url otherwise NA
    ensembl_host = ifelse(annotation==most_recent_annotation,NA,tolower(ensembl_archives[ensembl_archives$version==annotation,"url"]))
    ensembl_mart = useEnsembl(biomart="ensembl",version=annotation)
    ensembl_datasets = listDatasets(ensembl_mart)$dataset
    ensembl_dataset = paste(gsub("([a-z])[a-z]+_([a-z]+)","\\1\\2",tolower(species)),"gene","ensembl",sep="_")
    if(!ensembl_dataset %in% ensembl_datasets)stop(paste("Ensembl dataset '",dataset,"' not available on server!"))
    
    # get BM
    if(!is.na(ensembl_host)){
      sce = getBMFeatureAnnos(sce,    
                              filters = "ensembl_gene_id", 
                              attributes = c("ensembl_gene_id","external_gene_name","gene_biotype","chromosome_name","start_position","end_position"), 
                              feature_symbol = "external_gene_name",feature_id = "ensembl_gene_id",biomart = "ENSEMBL_MART_ENSEMBL",dataset = ensembl_dataset,
                              host = ensembl_host)
    }else{
      sce = getBMFeatureAnnos(sce,    
                              filters = "ensembl_gene_id", 
                              attributes = c("ensembl_gene_id","external_gene_name","gene_biotype","chromosome_name","start_position","end_position"), 
                              feature_symbol = "external_gene_name",feature_id = "ensembl_gene_id",biomart = "ENSEMBL_MART_ENSEMBL",dataset = ensembl_dataset)
    }
  } else{
  # no annotation
    row_info$feature_symbol = rownames(row_info)
    row_info$feature_id = rownames(row_info)
    rowData(sce) = DataFrame(row_info,row.names=rownames(row_info))
  }

  # make gene names unique and rename rows using column feature_symbol
  row_info = as.data.frame(rowData(sce))
  row_info$feature_symbol = ifelse(is.na(row_info$feature_symbol) | nchar(row_info$feature_symbol)==0,row_info$feature_id,row_info$feature_symbol)
  row_info$unique_gene_name = make.unique(row_info$feature_symbol)
  rownames(row_info) = row_info$unique_gene_name
  rowData(sce) = DataFrame(row_info,row.names = rownames(row_info))
  rownames(sce) = rownames(row_info)

  return(sce)
}
