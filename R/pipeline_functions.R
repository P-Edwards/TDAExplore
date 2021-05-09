# R Package TDAExplore

# Copyright (C) 2021 Parker Edwards

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 3 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
# USA


#' Sample patches from an image and compute landscapes from them
#'
#' The main pipeline function of TDAExplore for programmatic access. Most users will only be interested in the first 4 parameters.
#' A few parameters are undocumented and it is inadvisable to change them.
#' @param image_name File path (absolute or relative to current directory) to the image. Supported file formats are same as OpenImageR: .png, .jpeg
#' .jpg, and .tiff
#' @param number_of_patches Number of patches to sample from image
#' @param pixel_radius_for_patches All sampled patches will be of this radius, default is 50 pixels.
#' @param number_of_cores The number of CPU cores to parallelize over when computing landscapes for this image. Default is 1 (no parallelization).
#' @param patch_center_image File path (absolute or relative to current directory) to an image to use for selecting the centers of patches. Useful 
#' if you have a masked version of the image where you want patches to be centered. Defaults to FALSE, in which case the image at image_name is
#' used for patch centers.
#' @param return_patches If set to TRUE, output is a list image_landscape_results where image_landscape_results$data is the same data matrix as the default output,
#' image_landscape_results$patches is a list of patch objects (in the same order as the rows of the data matrix) that describe the center and position in the
#' original image for each patch. Default is FALSE, in which case output is described in Value below.
#' @param remove_patch_data If set to FALSE and return_patches is TRUE, the patch details also include the raw subset of the image comprising each patch. 
#' Useful for some imaging tasks. WARNING: Memory costs can quickly get out of hand when setting this to TRUE. Default is FALSE, in which case such data
#' is not retained.
#' @param return_samples If set to TRUE and return_patches is TRUE, will also append the samples to each patch entry. Default is FALSE.
#' Useful for some imaging tasks. WARNING: Memory costs can quickly get out of hand when setting this to TRUE. Default is FALSE, in which case such data
#' is not retained.
#' @param noise_threshold A number between 0 or 1 representing a percentage. The bottom noise_threshold pixels in the image by intensity will be set to 0. 
#' Used to eliminate small intensity background noise. Default is .05, i.e. 5%. 
#' @param number_of_layers For TDA experts: the maximum landscape Lambda_{number_of_layers} to compute. Default is 50.
#' @return Returns a sparse matrix image_landscape_results with number_of_patches rows. image_landscape_results\\[i,\\] is the landscape vector computed for the
#' i'th sampled patch. 
#'
#' Several base R matrix operations are not implemented for SparseM matrices. In this case you may need to use as.matrix(image_landscape_results) instead.
#' Be warned, however, that the memory consumption will generally be much larger.
#'
#' If the matrix has too many entries for SparseM, then instead a regular matrix is output.
#' @export
#' @examples
#' computed_landscapes <- patch_landscapes_from_image("test/triangles.jpg",10)
#' # Plot the average landscape vector
#' plot(colMeans(as.matrix(computed_landscapes)))
#' 
#' computed_landscapes <- patch_landscapes_from_image("test/triangles.jpg",10,return_patches=TRUE)
#' # Where was the center of the 5th patch?
#' computed_landscapes$patches[[5]]$center
patch_landscapes_from_image <- function(image_name,
                                        number_of_patches,
                                        pixel_radius_for_patches=50,
                                        number_of_cores=1,
                                        patch_center_image=FALSE,
                                        return_patches=FALSE,
                                        remove_patch_data=TRUE,
                                        return_samples=FALSE,
                                        lower_threshold=0,
                                        upper_threshold=1,
                                        proportion_of_patch_sparse=.025,
                                        noise_threshold=.05,
                                        number_of_layers=50) {
  # Derived point samples sizes from parameters
  delta <- 2.5
  spacing_for_disc_points <- delta/2
  
  remember <- getOption("warn")
  image_data <- OpenImageR::readImage(image_name) 
  if(length(dim(image_data)) > 2) { 
    image_data <- image_data[,,1]
  }
  image_data <- low_pixel_threshold(image_data,noise_threshold)
  image_length <- nrow(image_data)
  image_width <- ncol(image_data)
 
  
# generate centers
  if(!(patch_center_image==FALSE)) { 
    center_image <- OpenImageR::readImage(patch_center_image)
    if(length(dim(center_image)) > 2) { 
      center_image <- center_image[,,1]
    }
    center_image <- low_pixel_threshold(center_image,noise_threshold)
    image_cdf <- create_cdf_vector(center_image)
    sampling_dims <- dim(center_image)    
  } else {     
    image_cdf <- create_cdf_vector(image_data)
    sampling_dims <- dim(image_data)
  }
  
  center_points <- sample_centers_using_cdf_sparse(image_cdf,number_of_patches,sampling_dims,pixel_radius_for_patches/4)
  on_pixels <- which(image_data!=0,arr.ind=TRUE)
  centroid <- colSums(on_pixels)/nrow(on_pixels)
  
  if(upper_threshold - lower_threshold < 1) { 
    distances_from_centroid <- rdist::cdist(on_pixels,t(centroid))
    estimated_radius <- max(distances_from_centroid)
    center_distances <- rdist::cdist(t(center_points),t(centroid))/estimated_radius
    valid_center_indices <- which((center_distances>= lower_threshold)&(center_distances<=upper_threshold),arr.ind=TRUE)
    center_points <- center_points[,valid_center_indices]
  }

  # Forces sampling up to required number of patches, gradually 
  # relaxing sampling restriction
  k <- 3
  while(ncol(center_points)<number_of_patches) { 
    new_points <- sample_centers_using_cdf_sparse(image_cdf,number_of_patches,sampling_dims,pixel_radius_for_patches/(2**k))
    if(upper_threshold - lower_threshold < 1) { 
      center_distances <- rdist::cdist(t(new_points),t(centroid))/estimated_radius
      valid_center_indices <- which((center_distances>= lower_threshold)&(center_distances<=upper_threshold),arr.ind=TRUE)
      new_points <- new_points[,valid_center_indices]
    }
    center_points <- cbind(center_points,new_points)  
    k <- k+1
  }
  
  center_points <- center_points[,1:number_of_patches,drop=FALSE]
  # generate patches
  mask <- create_disc_mask(pixel_radius_for_patches)
  patches <- list()
  for(i in 1:ncol(center_points)) { 
    patches[[i]] <- get_patch_from_point_and_mask(image_data,center_points[,i],mask)
  }
  

  # generate sparse subsamples of each patch
  samples<-list()
  exclusion_disc <- create_disc_mask(delta)
  for (i in 1:number_of_patches){
    samples[[i]] <- sample_patch_using_top_intensities_sparse(patches[[i]],proportion_of_patch_sparse,exclusion_disc)
    if(remove_patch_data == TRUE) {
      patches[[i]]$data <- NULL
    } 
  }

    # add circular boundaries
  for(i in 1:number_of_patches) {
    samples[[i]] <- patch_sample_with_circ_boundary(samples[[i]],pixel_radius_for_patches,spacing_for_disc_points)
  }
  
  rval <- landscapes_from_samples(samples,pixel_radius_for_patches,number_of_cores=number_of_cores,number_of_layers=number_of_layers)  
  
  if(return_patches==TRUE) { 
    if(return_samples==TRUE) { 
      for(i in 1:number_of_patches) { 
        patches[[i]]$sample <- samples[[i]]
      }
    }
    return(list("data"=rval,"patches"=patches,"centroid"=centroid,"centers"=center_points))
  } else { 
    return(rval)
  }
}

#' Train and cross-validate a linear SVM classifier from image landscapes
#'
#' @param parameters File path to csv file designating which folders to analyze. A template parameters.csv is included at https://github.com/P-Edwards/TDAExplore-ML. 
#' You can alternatively provide parameters via this function's arguments directly. If you provide both a parameters file and function arguments, the function 
#' arguments will overwrite the arguments from the csv file.
#' @param number_of_cores The number of CPU cores for parallelization. Parallelization is over images, so allocating more cores than images provides no benefit.
#' @param experiment_name Prefix that will have additional information appended to form an informative file name for saving. Default is nothing.
#' @param image_directories Character vector of paths to directories holding image files to process. 
#' @param patch_center_image_directories Character vector of the same length as image_directories with paths to directories holding image files. The files in 
#' patch_center_image_directories[i] will be matched, in alphabetical order, to the files in image_directories[i], and used for selecting the centers of patches.
#' This is useful when e.g. patch_center_image_directories[i] contains masks for parts of the images in image_directories[i] to focus on. patch_center_image_directories[i]
#' may also be set to FALSE instead of a path, in which case no masked images will be used for center selection with image_directories[i].
#' @param directory_classes Character vector of same length as image directories. directory_classes[i] specifies the class of the images in image_directories[i].
#' @param data_results_directory File path (relative to current working directory) where data results can be saved. This function will only produce a file path
#' for saving, the user must save the results themselves if desired. ./tda_explores_results is the default.
#' @param radius_of_patches Pixel radius of patches. Default is 50 pixels.
#' @param patch_ratio The number of patches sampled per image will be patch_ratio*(PIXEL AREA OF IMAGE)/(PIXEL AREA OF SINGLE PATCH). Default is 2.
#' @param svm If set to TRUE, trains and tests SVM. If set to FALSE, only computes landscapes for each image. Default is FALSE.
#' @param multisvm If set to TRUE, trains and tests multi-class SVM. If set to FALSE, only computes landscapes for each image. 
#' @param randforest If set to TRUE, trains and tests random forest model. Supports multiclass. Could take prohitively long if --pca is not TRUE. Default is FALSE.
#' @param number_of_folds The number of folds to use in SVM cross-validation. The default is 5.
#' @param pca If set to TRUE, transforms landscapes after computation by projecting onto first 50 PC's, then scaling. Default is FALSE.
#' @param verbose If set to TRUE, outputs some progress information using print. Default is FALSE.
#' @param lower Experimental, leave default. 
#' @param benchmark If set to TRUE, uses a fork cluster type which is more easily trackable by benchmarking software. Does nothing on Windows. Default is FALSE.
#' @param lower Experimental, leave default. 
#' @param upper Experimental, leave default.
#' 
#' @return Returns a list, ml_results, whose entries contain the computation results. 
#' \itemize{
#'  \item{**ml_results$data_results_directory**}{ See arguments, stored for convenience.}
#'  \item{**ml_results$data_name_stem**}{ See arguments, stored for convenience}
#'  \item{**ml_results$summaries**}{ Matrix of computed landscapes for patches.}
#'  \item{**ml_results$patch_types**}{ Vector designating which class of image each landscape (row in the summaries matrix) comes from}
#'  \item{**ml_results$svm**}{ List containing results from each of the 5 folds of SVM training-testing computation}
#' }
#' @export
#' @examples
#' computation_results <- TDAExplore("parameters.csv",number_of_cores=5,verbose=TRUE)
TDAExplore <- function(parameters=FALSE,
                       number_of_cores=2,
                       experiment_name=FALSE,
                       image_directories=FALSE,
                       patch_center_image_directories=FALSE,
                       directory_classes=FALSE,
                       data_results_directory=FALSE,
                       radius_of_patches=FALSE,
                       patch_ratio=2,
                       svm=FALSE,
                       multisvm=FALSE,
                       randforest=FALSE,                       
                       number_of_folds=5,
                       pca=FALSE,                  
                       verbose=FALSE,
                       proportion=.025,
                       benchmark=FALSE,
                       lower=0,
                       upper=1) { 
  data_parameters <- utils::read.csv(parameters,stringsAsFactors=FALSE)
  provided_parameters <- colnames(data_parameters)
  ml_results <- list()

  if(image_directories!=FALSE) { 
    image_type_names <- image_directories
  } else if("image_directories" %in% provided_parameters) { 
    image_type_names <- data_parameters[,"image_directories"]
  } else { 
    stop("Image directories not specified, and there is no default.")
  }
  # Default behavior: Every directory is in different class
  if(directory_classes!=FALSE) { 
    class_names <- directory_classes
  } else if("directory_classes" %in% provided_parameters) { 
    class_names <- data_parameters[,"directory_classes"]
    # Have to be very careful with the factor call here: 
    # Default is levels=sort(unique(class_names)), which 
    # screws up ordering assumptions level if user
    # inputs classes in non-alphabetical order
    class_names <- factor(class_names,levels=unique(class_names),labels=unique(class_names))
  } else { 
    warning("No classes for directory specified, using directory names as a default.")
    class_names <- image_type_names
  }
  if(patch_center_image_directories!=FALSE) { 
    patch_center_names <- patch_center_image_directories
    patch_center_flag <- TRUE
  } else if("patch_center_images_directories" %in% provided_parameters) { 
    patch_center_names <- data_parameters[,"patch_center_images_directories"]
    patch_center_flag <- TRUE
  } else { 
    patch_center_flag <- FALSE 
  }


  image_file_names_by_directory <- list()
  patch_center_file_names_by_directory <- list()

  total_number_of_files <- 0
  for(i in 1:length(image_type_names)) { 
    image_file_names_by_directory[[i]] <- list.files(image_type_names[i],full.names=TRUE)
    if(patch_center_flag) { 
      patch_center_file_names_by_directory[[i]] <- list.files(levels(patch_center_names)[patch_center_names[i]],full.names=TRUE)
    } else { 
      # instead sets up a vector of appropriate length that is all FALSE
      patch_center_file_names_by_directory[[i]] <- as.vector(image_file_names_by_directory[[i]],mode="logical")
      patch_center_file_names_by_directory[[i]][1:length(patch_center_file_names_by_directory[[i]])] <- FALSE
    }
    total_number_of_files <- total_number_of_files + length(image_file_names_by_directory[[i]])
  }

  # Default behavior: Pixel radius is 2.5% of image length.
  if(radius_of_patches==FALSE) { 
    if("radius_of_patches" %in% provided_parameters) { 
        radius_of_patches <- data_parameters[1,"radius_of_patches"]
    } else { 
        warning("No pixel radius for patches specified. Using default, which is 2.5% of image length.")      
        first_image <- OpenImageR::readImage(unlist(image_file_names_by_directory)[2])
        radius_of_patches <- floor(0.025*dim(first_image)[1])
    }
  } 

  remember <- getOption("warn")
  options(warn=-1)
  image_for_dimensions <- OpenImageR::readImage(image_file_names_by_directory[[1]][1])
  options(warn=remember)

  pixel_area_of_images <- nrow(image_for_dimensions)*ncol(image_for_dimensions)
  patches_per_image <- max(floor(patch_ratio*pixel_area_of_images/(pi*(radius_of_patches)**2)),2)

  if(experiment_name==FALSE) {
    if("experiment_name" %in% provided_parameters) {
      experiment_name <- data_parameters[1,"experiment_name"]
    } else { 
      experiment_name <- floor(ruinf(1,min=0,max=1e6))
    }
  }

  if(data_results_directory==FALSE) { 
    if("data_results_directory" %in% provided_parameters) { 
      data_results_directory <- file.path(data_parameters[1,"data_results_directory"])
    } else { 
      data_results_directory <- "tda_explore_results"  
    }
  }
  if(!dir.exists(data_results_directory)) { 
    warning(paste("Results directory didn't exist, creating ",data_results_directory))
    dir.create(data_results_directory)
  }
  ml_results$data_results_directory <- data_results_directory
  current_time <- format(Sys.time(), "%b-%d-%M",digits=3)
  data_name_stem <- paste(experiment_name,"patches_",patch_ratio,"_radius_",radius_of_patches,"_",current_time,sep = "")
  if(upper - lower < 1) { 
    data_name_stem <- paste(data_name_stem,"_radialthresh")
  }
  ml_results$data_name_stem <- data_name_stem


  i <- 1
  type_vector <- vector("character",length=total_number_of_files*patches_per_image)
  offset <- 0
  for(image_list in image_file_names_by_directory) { 
    type_vector[(1+offset):(offset+length(image_list)*patches_per_image)] <- class_names[[i]]
    offset <- offset + length(image_list)*patches_per_image
    i <- i+1 
  }


  image_name_iterator <- iterators::iter(cbind(unlist(image_file_names_by_directory),unlist(patch_center_file_names_by_directory)),by="row")

  if(benchmark==FALSE) {
    cl <- parallel::makeCluster(number_of_cores,outfile="")
    doParallel::registerDoParallel(cl)
    } else { 
      # Forking using DoMC is picked up by SLURM
      # Only available on Linux
      doMC::registerDoMC(cores=number_of_cores)
  }
  `%dopar%` <- foreach::`%dopar%`
  unscrambled_data <- foreach::foreach(image_name=image_name_iterator,.packages=c("TDAExplore")) %dopar% { 
    options(warn=-1)
    if(verbose) {
      print(paste("Started image ", image_name[1]))    
      if(image_name[2]!=FALSE) {
        print(paste("With image for patch centers: ",patch_center_image))
      }
    }
    patch_summaries <- TDAExplore::patch_landscapes_from_image(image_name[1],patches_per_image,patch_center_image=image_name[2],pixel_radius_for_patches = radius_of_patches,proportion_of_patch_sparse=proportion,lower_threshold=lower,upper_threshold=upper) 
    if(verbose) {
      print(paste("Finished image ", image_name[1]))    
    }
    patch_summaries
  }
  if(benchmark==FALSE) {
    parallel::stopCluster(cl)
  }

  if(length(type_vector)>=2147483647) { 
    warning("There are more patches in your data set than the maximum supported by sparse matrices (2147483647 patches is the max). Continuing, but memory costs will be significantly higher.")
    unscrambled_data <- SparseM::as.matrix(do.call(rbind,unscrambled_data))
  } else { 
    unscrambled_data <- as.matrix.csr(do.call(rbind,unscrambled_data))
  }

  if(pca!=FALSE) { 
    ml_results$landscapes_svd <- RSpectra::svds(SparseM::as.matrix(unscrambled_data),k=1000,nu=0,nv=1000,center=FALSE,scale=FALSE,tol=NULL)
    unscrambled_data <- SparseM::as.matrix(unscrambled_data)%*%(ml_results$landscapes_svd$v)
    unscrambled_data <- scale(unscrambled_data)
  }

  image_file_names <- unlist(image_file_names_by_directory)

  
  ml_results$summaries <- unscrambled_data
  ml_results$patch_types <- type_vector
  ml_results$class_names <- class_names
  ml_results$patches_per_image <- patches_per_image
  ml_results$patch_ratio <- patch_ratio
  ml_results$image_file_names <- image_file_names
  ml_results$radius <- radius_of_patches

  if(svm!=FALSE) {
    if(verbose) {
      print("Starting per-landscape SVM")
    }
    if(pca!=FALSE) { 
      # Sets run type to L1 regularized
      svm_run_type <- 5
    } else { 
      # If not projected, use primal L2 regularized + loss
      svm_run_type <- 2
    }
    ml_results$svm <- list()

    # Cross validation on PCA-rotated patch landscapes
    number_of_validation_steps <- number_of_folds
    
    #randomly shuffle the data
    shuffled_order <- sample(length(image_file_names))
    shuffled_image_names <- image_file_names[shuffled_order]

    shuffled_patch_indices <- vector(mode="double",length=length(type_vector))
    image_folds <- cut(seq(1,length(shuffled_image_names)),breaks=number_of_validation_steps,labels=FALSE)
    folds <- vector(mode="double",length=length(type_vector))
    for(i in 1:length(shuffled_order)) { 
      image_index <- shuffled_order[i]
      shuffled_patch_indices[((i-1)*patches_per_image+1):(i*patches_per_image)] <- ((image_index-1)*patches_per_image+1):(image_index*patches_per_image)
      folds[((i-1)*patches_per_image+1):(i*patches_per_image)] <- image_folds[i]
    }

    # First: Use 20% of the data to tune SVM parameters.  

    shuffled_pca <- unscrambled_data[shuffled_patch_indices,]
    shuffled_types <- type_vector[shuffled_patch_indices]
    trainIndexes <- which(folds==1,arr.ind=TRUE)

    # Quick heuristic parameter tuning
    cost <- LiblineaR::heuristicC(SparseM::as.matrix(shuffled_pca[trainIndexes,]))



    # Now do cross validation with the remaining data and selected parameters
    reduced_data <- shuffled_pca
    reduced_types <- shuffled_types

    shuffled_pca <- NULL
    shuffled_types <- NULL
    gc()

    patch_accuracies <- vector("double",number_of_validation_steps)
    image_accuracies <- vector("double",number_of_validation_steps)

    SVM_file_path <- file.path(data_results_directory,paste(data_name_stem,"_SVM_cross_validation.csv",sep=""))
    SVM_avg_per_patch_file_path <- file.path(data_results_directory,paste(data_name_stem,"_patch_feature_average_image_SVM_cross_validation.csv",sep=""))
    proportion_vector <- vector("double",length=max(reduced_types))
    for(i in 1:max(reduced_types)) { 
      proportion_vector[i] <- sum(reduced_types==i)/length(reduced_types)
    }


    reduced_types_unique <- unique(type_vector)
    for(i in 1:number_of_validation_steps) { 
      if(verbose) {
        print(paste("Starting cross validation fold number ",i))
      }
      testIndexes <- which(folds==i,arr.ind=TRUE)
      trainIndexes <- -testIndexes
      weights <- vector(mode="double",length=length(unique(reduced_types[trainIndexes])))
      names(weights) <- unique(reduced_types[trainIndexes]) 
      for(j in 1:length(weights)) { 
        weights[j] <- sum(reduced_types[trainIndexes]==names(weights)[j])
      }
      max_number <- max(weights)
      for(j in 1:length(weights)) { 
        weights[j] <- max_number/weights[j]
      }
      reduced_types_training <- reduced_types[trainIndexes]
      train_data <- reduced_data[trainIndexes,]
      if(reduced_types_training[1]!=reduced_types_unique[1]) { 
        switch_index <- min(which(reduced_types_training==reduced_types_unique[1]))
        reduced_types_training[1] <- reduced_types_unique[1]
        reduced_types_training[switch_index] <- reduced_types_unique[2]
        temporary_row <- train_data[switch_index,]
        train_data[switch_index,] <- train_data[1,]
        train_data[1,] <- temporary_row
      }
      svm_model <- LiblineaR::LiblineaR(data=train_data,target=factor(reduced_types_training),wi=weights,cost=cost,type=svm_run_type)

      ml_results$svm[[i]] <- list()
      ml_results$svm[[i]]$testing_data <- reduced_data[testIndexes,]
      ml_results$svm[[i]]$testing_labels <- reduced_types[testIndexes]
      
      prediction_values <- predict(svm_model,reduced_data[testIndexes,])
      actual_names <- levels(class_names)
      predicted_names <- levels(class_names)
      for(j in 1:length(actual_names)) { 
        actual_names[j] <- paste(actual_names[j],"_actual")
        predicted_names[j] <- paste(predicted_names[j],"_predicted")
      }
      pred_table <- table(factor(prediction_values$predictions,labels=predicted_names,levels=unique(reduced_types)),factor(reduced_types[testIndexes],labels=actual_names,levels=unique(reduced_types)))
      # write.table(pred_table,file=SVM_file_path,sep=",",dec=".",append=TRUE)
      patch_accuracies[i] <- sum(diag(pred_table))/sum(pred_table)
      
      # Classification with average transformed image data
      # i.e. 
      # (1) Transform all patch landscapes into single numbers using model
      # (2) Average all numbers across each image
      # (3) Train and test an SVM model for classifying images
      model_matrix <- svm_model$W
      weight_vector <- model_matrix[1,(1:(dim(model_matrix)[2]-1))]
      bias_term <- model_matrix[1,dim(model_matrix)[2]]
      transformed_data <- as.matrix.csr(reduced_data%*%(weight_vector) + bias_term,nrow=nrow(reduced_data),ncol=1)
      averaged_data <- average_vectors_for_images(transformed_data,patches_per_image,shuffled_image_names,reduced_types,1,1)
      transformed_data <- averaged_data$image_weights
      transformed_types <- averaged_data$image_types
      
      costIndexes <- which(image_folds==1,arr.ind=TRUE)
      testIndexes <- which(image_folds==i,arr.ind=TRUE)
      trainIndexes <- -testIndexes

      # Quick heuristic parameter tuning
      image_cost <- LiblineaR::heuristicC(SparseM::as.matrix(transformed_data[costIndexes,]))

      train_data <- as.matrix.csr(transformed_data[trainIndexes,],ncol=1)
      test_data <- as.matrix.csr(transformed_data[testIndexes,],ncol=1)

      weights <- vector(mode="double",length=length(unique(transformed_types[trainIndexes])))
      names(weights) <- unique(transformed_types[trainIndexes]) 
      for(j in 1:length(weights)) { 
        weights[j] <- sum(transformed_types[trainIndexes]==names(weights)[j])
      }
      max_number <- max(weights)
      for(j in 1:length(weights)) { 
        weights[j] <- max_number/weights[j]
      }

      transformed_types_training <- transformed_types[trainIndexes]
      if(transformed_types_training[1]!=reduced_types_unique[1]) { 
        switch_index <- min(which(transformed_types_training==reduced_types_unique[1]))
        transformed_types_training[1] <- reduced_types_unique[1]
        transformed_types_training[switch_index] <- reduced_types_unique[2]
        temporary_row <- train_data[switch_index,]
        train_data[switch_index,] <- train_data[1,]
        train_data[1,] <- temporary_row
      }
      image_svm_model <- LiblineaR::LiblineaR(data=train_data,target=factor(transformed_types_training),cost=image_cost,wi=weights,type=svm_run_type)
      prediction_values <- predict(image_svm_model,test_data)
      actual_names <- levels(class_names)
      predicted_names <- levels(class_names)
      for(j in 1:length(actual_names)) { 
        actual_names[j] <- paste(actual_names[j],"_actual")
        predicted_names[j] <- paste(predicted_names[j],"_predicted")
      }
      pred_table <- table(factor(prediction_values$predictions,labels=predicted_names,levels=unique(transformed_types)),factor(transformed_types[testIndexes],labels=actual_names,levels=unique(transformed_types)))
      image_accuracies[i] <- sum(diag(pred_table))/sum(pred_table)
      
      ml_results$svm[[i]]$svm_model <- svm_model
      ml_results$svm[[i]]$patch_svm_model <- svm_model
      ml_results$svm[[i]]$image_svm_model <- image_svm_model
      ml_results$svm[[i]]$svm_image_training_indices <- shuffled_order[trainIndexes]
      ml_results$svm[[i]]$svm_patch_accuracies <- patch_accuracies
      ml_results$svm[[i]]$svm_image_accuracies <- image_accuracies
      ml_results$svm[[i]]$svm_cost <- cost 
    }

    if(verbose) {
      print("Radius")
      print(radius_of_patches)
      print("Average patch accuracies")
      print(mean(patch_accuracies))
      print("Average image accuracies")
      print(mean(image_accuracies))
    }
    
  }
  if(multisvm!=FALSE) {
    if(verbose) {
      print("Starting per-landscape SVM")
    }

    svm_run_type <- 4
    ml_results$multisvm <- list()

    # Cross validation on PCA-rotated patch landscapes
    number_of_validation_steps <- number_of_folds
    
    #randomly shuffle the data
    shuffled_order <- sample(length(image_file_names))
    shuffled_image_names <- image_file_names[shuffled_order]

    shuffled_patch_indices <- vector(mode="double",length=length(type_vector))
    image_folds <- cut(seq(1,length(shuffled_image_names)),breaks=number_of_validation_steps,labels=FALSE)
    folds <- vector(mode="double",length=length(type_vector))
    for(i in 1:length(shuffled_order)) { 
      image_index <- shuffled_order[i]
      shuffled_patch_indices[((i-1)*patches_per_image+1):(i*patches_per_image)] <- ((image_index-1)*patches_per_image+1):(image_index*patches_per_image)
      folds[((i-1)*patches_per_image+1):(i*patches_per_image)] <- image_folds[i]
    }

    # First: Use 20% of the data to tune SVM parameters.  

    shuffled_pca <- unscrambled_data[shuffled_patch_indices,]
    shuffled_types <- type_vector[shuffled_patch_indices]
    trainIndexes <- which(folds==1,arr.ind=TRUE)

    # Quick heuristic parameter tuning
    cost <- LiblineaR::heuristicC(SparseM::as.matrix(shuffled_pca[trainIndexes,]))



    # Now do cross validation with the remaining data and selected parameters
    reduced_data <- shuffled_pca
    reduced_types <- shuffled_types
    reduced_types_unique <- unique(type_vector)

    shuffled_pca <- NULL
    shuffled_types <- NULL
    gc()

    patch_accuracies <- vector("double",number_of_validation_steps)
    image_accuracies <- vector("double",number_of_validation_steps)

    
    for(i in 1:number_of_validation_steps) { 
      if(verbose) {
        print(paste("Starting cross validation fold number ",i))
      }
      testIndexes <- which(folds==i,arr.ind=TRUE)
      trainIndexes <- -testIndexes
      weights <- vector(mode="double",length=length(unique(reduced_types[trainIndexes])))
      names(weights) <- levels(class_names)          
      for(j in 1:length(weights)) { 
        weights[j] <- sum(reduced_types[trainIndexes]==unique(type_vector)[j])
      }
      max_number <- max(weights)
      for(j in 1:length(weights)) { 
        weights[j] <- max_number/weights[j]
      }
      reduced_types_training <- reduced_types[trainIndexes]
      train_data <- reduced_data[trainIndexes,]
      svm_model <- LiblineaR::LiblineaR(data=train_data,target=factor(reduced_types_training,levels=unique(type_vector),labels=levels(class_names)),wi=weights,cost=cost,type=svm_run_type)

      
      ml_results$multisvm[[i]] <- list()
      ml_results$multisvm[[i]]$testing_data <- reduced_data[testIndexes,]
      ml_results$multisvm[[i]]$testing_labels <- reduced_types[testIndexes]
      

      prediction_values <- predict(svm_model,SparseM::as.matrix(reduced_data[testIndexes,]))$predictions
      
      actual_names <- levels(class_names)
      predicted_names <- levels(class_names)
      for(j in 1:length(actual_names)) { 
        actual_names[j] <- paste(actual_names[j],"_actual")
        predicted_names[j] <- paste(predicted_names[j],"_predicted")
      }
      pred_table <- table(prediction_values,factor(reduced_types[testIndexes],labels=actual_names,levels=unique(type_vector)))
      patch_accuracies[i] <- sum(diag(pred_table))/sum(pred_table)
      
      # Classifiction with average probability vectors: 
      # Predict class probabilities for each testing patch, then 
      # take an average over all patches in each testing image. 
      # Highest average probability wins the image.
      
      patch_probability_predictions <- predict(svm_model,SparseM::as.matrix(reduced_data),decisionValues=TRUE)$decisionValues

      
      number_of_classes <- ncol(patch_probability_predictions)   

      named_types <- factor(reduced_types[testIndexes],levels = unique(type_vector),labels=levels(class_names))

      image_probability_predictions <- average_vectors_for_images(patch_probability_predictions,patches_per_image,levels(class_names),reduced_types,1,number_of_classes)      
      transformed_data <- image_probability_predictions$image_weights
      transformed_types <- image_probability_predictions$image_types
      
      costIndexes <- which(image_folds==1,arr.ind=TRUE)
      testIndexes <- which(image_folds==i,arr.ind=TRUE)
      trainIndexes <- -testIndexes

      # Quick heuristic parameter tuning
      image_cost <- LiblineaR::heuristicC(SparseM::as.matrix(transformed_data[costIndexes,]))

      train_data <- as.matrix.csr(transformed_data[trainIndexes,])
      test_data <- as.matrix.csr(transformed_data[testIndexes,])

      weights <- vector(mode="double",length=length(unique(transformed_types[trainIndexes])))
      names(weights) <- levels(class_names)
      for(j in 1:length(weights)) { 
        weights[j] <- sum(transformed_types[trainIndexes]==unique(type_vector)[j])
      }
      max_number <- max(weights)
      for(j in 1:length(weights)) { 
        weights[j] <- max_number/weights[j]
      }

      image_svm_model <- LiblineaR::LiblineaR(data=train_data,target=factor(transformed_types[trainIndexes],levels=unique(type_vector),labels=levels(class_names)),cost=image_cost,wi=weights,type=svm_run_type)
      prediction_values <- predict(image_svm_model,test_data)$predictions
                      
      

      image_pred_table <- table(prediction_values,factor(transformed_types[testIndexes],levels=unique(type_vector),labels=actual_names))
      image_accuracies[i] <- sum(diag(image_pred_table))/sum(image_pred_table)

      
      ml_results$multisvm[[i]]$svm_model <- svm_model
      ml_results$multisvm[[i]]$patch_svm_model <- svm_model      
      ml_results$multisvm[[i]]$svm_image_training_indices <- shuffled_order[trainIndexes]
      ml_results$multisvm[[i]]$svm_patch_accuracies <- patch_accuracies
      ml_results$multisvm[[i]]$svm_image_accuracies <- image_accuracies
    }

    if(verbose) {
      print("Radius")
      print(radius_of_patches)
      print("Average patch accuracies")
      print(mean(patch_accuracies))
      print("Average image accuracies")
      print(mean(image_accuracies))
    }
    
  }
  if(randforest!=FALSE) {
    if(verbose) {
      print("Starting per-landscape random forest")
    }
    ml_results$rand_forest <- list()

    # Cross validation on PCA-rotated patch landscapes
    number_of_validation_steps <- number_of_folds
    
    #randomly shuffle the data
    shuffled_order <- sample(length(image_file_names))
    shuffled_image_names <- image_file_names[shuffled_order]

    shuffled_patch_indices <- vector(mode="double",length=length(type_vector))
    image_folds <- cut(seq(1,length(shuffled_image_names)),breaks=number_of_validation_steps,labels=FALSE)
    folds <- vector(mode="double",length=length(type_vector))
    for(i in 1:length(shuffled_order)) { 
      image_index <- shuffled_order[i]
      shuffled_patch_indices[((i-1)*patches_per_image+1):(i*patches_per_image)] <- ((image_index-1)*patches_per_image+1):(image_index*patches_per_image)
      folds[((i-1)*patches_per_image+1):(i*patches_per_image)] <- image_folds[i]
    }

    # Now do cross validation with the remaining data and selected parameters
    reduced_data <- unscrambled_data[shuffled_patch_indices,]
    reduced_types <- type_vector[shuffled_patch_indices]


    patch_accuracies <- vector("double",number_of_validation_steps)
    image_accuracies <- vector("double",number_of_validation_steps)

    reduced_types_unique <- unique(type_vector)
    for(i in 1:number_of_validation_steps) { 
      if(verbose) {
        print(paste("Starting random forest cross validation fold number ",i))
      }
      testIndexes <- which(folds==i,arr.ind=TRUE)
      trainIndexes <- -testIndexes
      reduced_types_training <- reduced_types[trainIndexes]
      train_data <- SparseM::as.matrix(reduced_data[trainIndexes,])


      rf_model <- randomForest::randomForest(train_data,y=factor(reduced_types_training,levels=unique(type_vector),labels=levels(class_names)))

      ml_results$rand_forest[[i]] <- list()
      ml_results$rand_forest[[i]]$testing_data <- SparseM::as.matrix(reduced_data[testIndexes,])
      ml_results$rand_forest[[i]]$testing_labels <- reduced_types[testIndexes]
      
      prediction_values <- predict(rf_model,SparseM::as.matrix(reduced_data[testIndexes,]))
      actual_names <- levels(class_names)
      predicted_names <- levels(class_names)
      for(j in 1:length(actual_names)) { 
        actual_names[j] <- paste(actual_names[j],"_actual")
        predicted_names[j] <- paste(predicted_names[j],"_predicted")
      }
      pred_table <- table(prediction_values,factor(reduced_types[testIndexes],labels=actual_names,levels=unique(type_vector)))
      patch_accuracies[i] <- sum(diag(pred_table))/sum(pred_table)
      
      # Classifiction with average probability vectors: 
      # Predict class probabilities for each testing patch, then 
      # take an average over all patches in each testing image. 
      # Highest average probability wins the image.
      
      patch_probability_predictions <- predict(rf_model,SparseM::as.matrix(reduced_data[testIndexes,]),type="prob")
      number_of_classes <- ncol(patch_probability_predictions)   

      named_types <- factor(reduced_types[testIndexes],levels = unique(type_vector),labels=levels(class_names))

      image_probability_predictions <- average_vectors_for_images(patch_probability_predictions,patches_per_image,levels(class_names),as.character(named_types),1,number_of_classes)
      prediction_values <- vector("character",length=nrow(image_probability_predictions$image_weights))
      for(j in 1:length(prediction_values)) { 
        prediction_values[j] <- class_names[which.max(image_probability_predictions$image_weights[j,])]
      }
      prediction_values <- factor(prediction_values,levels=unique(type_vector),labels=predicted_names)

      image_pred_table <- table(prediction_values,factor(image_probability_predictions$image_types,levels=levels(class_names),labels=actual_names))
      image_accuracies[i] <- sum(diag(image_pred_table))/sum(image_pred_table)
      
      ml_results$rand_forest[[i]]$image_predictions <- prediction_values
      ml_results$rand_forest[[i]]$ground_values <- image_probability_predictions$image_types
      ml_results$rand_forest[[i]]$rf_model <- rf_model      
      ml_results$rand_forest[[i]]$rf_image_training_indices <- shuffled_order[trainIndexes]
      ml_results$rand_forest[[i]]$rf_patch_accuracies <- patch_accuracies
      ml_results$rand_forest[[i]]$rf_image_accuracies <- image_accuracies      
    }
    
    if(verbose) {
      print("Radius")
      print(radius_of_patches)
      print("Average RF patch accuracies")
      print(mean(patch_accuracies))
      print("Average RF image accuracies")
      print(mean(image_accuracies))
    }    

  }

  return(ml_results)
}


#' Convolve an image using patch scores
#'
#' Combines landscape data, data for positioning patches on an image, and a scoring function for landscapes
#' to produce a mask for the image.
#' A few parameters are undocumented and it is inadvisable to change them.
#' @param image_name File path (absolute or relative to current directory) to the image. Supported file formats are same as OpenImageR: .png, .jpeg
#' .jpg, and .tiff
#' @param landscape_data Landscapes for image patches e.g. like those returned from patch_landscapes_for_image(return_patches=TRUE)$data
#' @param patches Patches e.g. like those return from patch_landscapes_for_image(return_patches=TRUE)$patches
#' @param radius Radius of patches in pixels
#' @param transformation_function A function which takes as input a matrix and outputs a vector with scores, one for each row in the matrix
#' @param interval_rep If TRUE, computes the interval representation described below. Default is FALSE.
#' @return Returns a list return_list with sufficient information to plot the image and its mask.
#' return_list$image is a matrix with pixel intensities as values
#' return_list$mask is a matrix of the same dimensions with score intensities as values
#' return_list$interval_rep is a vector of length 50 which compresses the mask into a radial representation
#' return_list$weights is a vector of the patch scores
#' @export
weight_image_using_landscapes_and_transform <- function(image_name,landscape_data,patches,radius,transformation_function=identity,min_weight=0,max_weight=1,interval_rep=FALSE) { 
  remember <- getOption("warn")
  options(warn=-1)
  image_data <- OpenImageR::readImage(image_name)
  options(warn=remember)
  if(length(dim(image_data)) > 2) { 
    image_data <- image_data[,,1]
  }
  image_data <- low_pixel_threshold(image_data,.05)
  patch_weights <- transformation_function(landscape_data)
  unmod_min_weight <- min(patch_weights)
  unmod_max_weight <- max(patch_weights)


  # Scale weights to positive
  unmod_weights <- patch_weights
  patch_weights <- patch_weights - min_weight

  mask <- create_disc_mask(radius)
  image_weight_mask <- as.matrix(image_data)
  image_weight_mask[,] <- 0
  image_average_mask <- as.matrix(image_weight_mask)
  image_average_mask[,] <- 1
  all_ones_submatrix <- as.matrix(mask)
  all_ones_submatrix[,] <- 1
  for(i in 1:length(patches)) { 
    current_patch <- patches[[i]]
    image_average_mask[current_patch$xmin:current_patch$xmax,current_patch$ymin:current_patch$ymax] <- 0
  }  
  for(i in 1:length(patches)) { 
    current_patch <- patches[[i]]
    weight_patch <- all_ones_submatrix*patch_weights[i]
    patch_dims <- dim(image_weight_mask[current_patch$xmin:current_patch$xmax,current_patch$ymin:current_patch$ymax])
    image_weight_mask[current_patch$xmin:current_patch$xmax,current_patch$ymin:current_patch$ymax] <- image_weight_mask[current_patch$xmin:current_patch$xmax,current_patch$ymin:current_patch$ymax] + weight_patch[1:patch_dims[1],1:patch_dims[2]]
    image_average_mask[current_patch$xmin:current_patch$xmax,current_patch$ymin:current_patch$ymax] <- image_average_mask[current_patch$xmin:current_patch$xmax,current_patch$ymin:current_patch$ymax] + (max_weight - min_weight)*all_ones_submatrix[1:patch_dims[1],1:patch_dims[2]]
  }
  total_weight_mask <- image_weight_mask/image_average_mask
  # Set "off cell" points to 0
  total_weight_mask[which(image_data==0,arr.ind=TRUE)] <- 0

  

  if(interval_rep) { 
    # Compute centroid from base image

    # Each 2-entry row in this matrix collects the 
    # (first_index,second_index) of active pixels
    on_pixels <- which(image_data!=0,arr.ind=TRUE)
    centroid <- colSums(on_pixels)/nrow(on_pixels)


    # Compute maximum radius from base image + centroid
    
    distances_from_centroid <- rdist::cdist(on_pixels,t(centroid))
    estimated_radius <- max(distances_from_centroid)

    # Compute "percentage mask": "Radius mask"/(max radius)
    single_dist_from_centroid <- function(i,j) { 
      diff <- centroid - c(i,j)
      return(sqrt(diff%*%diff))
    }
    percentage_mask <- matrix(mapply(single_dist_from_centroid,row(image_data),col(image_data)),nrow=nrow(image_data),ncol=ncol(image_data))/estimated_radius


    # Compute "aggregated interval"
    number_of_intervals <- 50
    aggregated_output <- vector("double",number_of_intervals)
    for(i in 1:number_of_intervals) { 
      lower_threshold <- (i-1)/number_of_intervals
      upper_threshold <- i/number_of_intervals
      active_pixels_within_threshold <- which((image_data>0)&(percentage_mask>=lower_threshold)&(percentage_mask<upper_threshold),arr.ind=TRUE)
      aggregated_output[i] <- sum(total_weight_mask[active_pixels_within_threshold])/length(active_pixels_within_threshold)
    }
  } else { 
    aggregated_output <- 0
  }

  return(list("image"=rotate(image_data),"mask"=rotate(total_weight_mask),"min"=unmod_min_weight,"max"=unmod_max_weight,"weights"=unmod_weights,interval_rep=aggregated_output))
}


image_demonstrator <- function(image_name,
                               number_of_patches,
                               pixel_radius_for_patches=50,
                               number_of_cores=1,
                               proportion_of_patch_sparse=.025,
                               noise_threshold=.05,
                               number_of_layers=50) { 

  # Generate patches and landscapes with all info attached
  base_information <- patch_landscapes_from_image(image_name,number_of_patches,pixel_radius_for_patches,number_of_cores,FALSE,TRUE,FALSE,TRUE,0,1,proportion_of_patch_sparse,noise_threshold,number_of_layers)

  # Plot original image with transparent discs over patch locations
  original_image <- OpenImageR::readImage(image_name)
  if(length(dim(original_image))>2) {
    original_image <- original_image[,,1]
  }
  mask <- matrix(FALSE,nrow=nrow(original_image),ncol=ncol(original_image))
  disc <- create_disc_mask(pixel_radius_for_patches)
  for(i in 1:length(base_information$patches)) { 
    this_patch <- base_information$patches[[i]]
    mask[(this_patch$xmin:this_patch$xmax),(this_patch$ymin:this_patch$ymax)] <- disc[1:nrow(this_patch$data),1:ncol(this_patch$data)]
  }
  par(mar=c(0,0,0,0)); 
  image(rotate(original_image),col=gray.colors(70000,start = 0.0),zlim=c(0,max(original_image)),axes=FALSE) 
  image(rotate(mask),zlim=c(1e-3,1),col=colorspace::divergingx_hcl(2,palette = "RdBu",alpha=.60,rev=TRUE),add=TRUE)

  # Plot patch samples and PD's
  for(i in 1:length(base_information$patches)) { 
    this_patch_sample <- base_information$patches[[i]]$sample
    plot(this_patch_sample[1,],this_patch_sample[2,])
  }
}