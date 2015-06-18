library(rhdf5)
library(biom)

# This generates the matrix columns-wise
generate_matrix <- function(x){
  indptr  = x$sample$matrix$indptr+1
  indices = x$sample$matrix$indices+1
  data    = x$sample$matrix$data
  nr = length(x$observation$ids)
  
  counts = sapply(2:length(indptr),function(i){
    x = rep(0,nr)
    seq = indptr[i-1]:(indptr[i]-1)
    x[indices[seq]] = data[seq]
    x
  })
  rownames(counts) = x$observation$ids
  colnames(counts) = x$sample$ids
  # I wish this next line wasn't necessary
  lapply(1:nrow(counts),function(i){
    counts[i,]
  })
}
generate_metadata <- function(x){
  metadata = x$metadata
  metadata = lapply(1:length(x$ids),function(i){
    id_metadata = lapply(metadata,function(j){
      if(length(dim(j))>1){ as.vector(j[,i,drop=FALSE]) }
      else{ j[i] }
    })
    list(id = x$ids[i],metadata=id_metadata)
  })
  return(metadata)
}
namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}

read_hdf5_biom<-function(file_input){
  x = h5read(file_input,"/",read.attributes = TRUE)
  data = generate_matrix(x)
  rows = generate_metadata(x$observation)
  columns = generate_metadata(x$sample)
  shape = c(length(data),length(data[[1]])) # dim(data)
  
  # Experimental -- need to actually load these from file
  id = attr(x,"id")
  vs = attr(x,"format-version")
  format = sprintf("Biological Observation Matrix %s.%s",vs[1],vs[2])
  format_url = attr(x,"format-url")
  type = "OTU table"
  #type=attr(x,"type")
  generated_by = attr(x,"generated-by")
  date = attr(x,"creation-date")
  matrix_type = "dense"
  matrix_element_type = "int"
  
  namedList(id,format,format_url,type,generated_by,date,matrix_type,matrix_element_type,
            rows,columns,shape,data)
}