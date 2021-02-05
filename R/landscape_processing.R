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

landscape_to_vector <- function(input_landscape,number_of_layers=50) { 
  landscape_0 <- vector(mode="double",length=(dim(input_landscape$degree_0)[2]*number_of_layers))
  landscape_0_layers <- min(number_of_layers,dim(input_landscape$degree_0)[1])
  landscape_0_flattened_length <- landscape_0_layers*dim(input_landscape$degree_0)[2]
  
  landscape_1 <- vector(mode="double",length=dim(input_landscape$degree_1)[2]*number_of_layers)
  landscape_1_layers <- min(number_of_layers,dim(input_landscape$degree_1)[1])
  landscape_1_flattened_length <- landscape_1_layers*dim(input_landscape$degree_1)[2]
  
  landscape_0[1:landscape_0_flattened_length] <- as.vector(t(input_landscape$degree_0[1:landscape_0_layers,,2]))
  landscape_1[1:landscape_1_flattened_length] <- as.vector(t(input_landscape$degree_1[1:landscape_1_layers,,2]))
  return(c(landscape_0,landscape_1))
}

landscapes_from_samples <- function(samples,pixel_radius_for_patches,number_of_cores=1,diag_function=TDA::alphaComplexDiag,number_of_layers=50) { 

  if(number_of_cores > 1) {
    
    cl <- parallel::makeCluster(number_of_cores) 
    doParallel::registerDoParallel(cl)
    sample_iterator <- iterators::iter(samples)
    `%dopar%` <- foreach::`%dopar%`
    persistence_diagrams <-  foreach::foreach(sample=sample_iterator,.packages=c("TDAExplore","TDA")) %dopar% {
      output_diagram <- diag_function(t(sample))
      output_diagram
    }
    parallel::stopCluster(cl)
  } else { 
      persistence_diagrams <-  list() 
      for(i in 1:length(samples)) {
        persistence_diagrams[[i]] <- diag_function(t(samples[[i]]))
      }
    }
  
  #separate between 0 and 1 degree homology
  persistence_diagrams_degree_0<-list()
  persistence_diagrams_degree_1<-list()
  for(i in 1:length(samples)) {
    current_diagram <- persistence_diagrams[[i]]$diagram
    current_diagram <- current_diagram[which(current_diagram[,1]==0),]
    # No longer need 3rd column which contains degree info
    # Take square roots to get in expected range. Note: Works because alpha complex 
    # doesn't take negative filtration values.
    current_diagram <- sqrt(current_diagram[,c(2,3)])
    persistence_diagrams_degree_0[[i]] <- current_diagram
    
  }
  for(i in 1:length(samples)){
    current_diagram <- persistence_diagrams[[i]]$diagram
    current_diagram <- current_diagram[which(current_diagram[,1]==1),,drop=FALSE]
    current_diagram <- sqrt(current_diagram[,c(2,3),drop=FALSE])
    persistence_diagrams_degree_1[[i]] <- current_diagram
  }
  persistence_diagrams <- NULL
  gc()

  #get landscapes
  first_landscape_degree_0 <- landscape(persistence_diagrams_degree_0[[1]],dx=pixel_radius_for_patches/200,max_x=pixel_radius_for_patches)$getInternal()
  first_landscape_degree_1 <- landscape(persistence_diagrams_degree_1[[1]],dx=pixel_radius_for_patches/200,max_x=pixel_radius_for_patches)$getInternal()
  first_landscape <- list("degree_0"=first_landscape_degree_0,"degree_1"=first_landscape_degree_1)
  first_vector <- landscape_to_vector(first_landscape,number_of_layers)
  
  rval <- matrix(0L,nrow=length(samples),ncol=length(first_vector))
  rval[1,] <- first_vector
  if(length(samples) > 1) { 
    for(i in 2:length(samples)) {
      landscape_degree_0 <- landscape(persistence_diagrams_degree_0[[i]],dx=pixel_radius_for_patches/200,max_x=pixel_radius_for_patches)$getInternal()
      landscape_degree_1 <- landscape(persistence_diagrams_degree_1[[i]],dx=pixel_radius_for_patches/200,max_x=pixel_radius_for_patches)$getInternal()
      current_landscape <- list("degree_0"=landscape_degree_0,"degree_1"=landscape_degree_1)
      landscape_vector <- landscape_to_vector(current_landscape,number_of_layers)
      rval[i,] <- landscape_vector
    }
  }
  if(length(rval)>=2147483647) { 
    return(rval)
  }
  return(SparseM::as.matrix.csr(rval))
}

average_vectors_for_images <- function(data_matrix,number_of_patches,names_list,types_list,number_along_each_dimension=3,number_of_dimensions = 10) { 
  
  if(!is.null(dim(data_matrix))) {
    weights <- SparseM::as.matrix.csr(data_matrix[,1:number_of_dimensions])
    number_of_images <- nrow(data_matrix)/number_of_patches
  } else {  
    weights <- matrix(data_matrix,nrow=length(data_matrix),ncol=1)
    number_of_images <- length(weights)/number_of_patches
  }
  if(number_of_dimensions == 1) { 
    averaging_function <- mean
  } else { 
    averaging_function <- colMeans
  }

  image_weights <- matrix(0L,nrow=number_of_images,ncol=number_of_dimensions)
  image_types <- vector("double",length=number_of_images)
  for(i in 1:number_of_images) { 
    if(number_of_patches > 1) {
      image_weights[i,] <- averaging_function(as.matrix(weights[(1+number_of_patches*(i-1)):(number_of_patches*i),]))[1:number_of_dimensions]    
    } else { 
      image_weights[i,] <- weights[i,]
    }
    image_types[i] <- types_list[1+(i-1)*number_of_patches]
  }

  output_names <- list()
  for(i in 1:number_of_dimensions) { 
    dimension_order <- order(image_weights[,i])
    dimension_list <- list() 
    step_size <- number_of_images/number_along_each_dimension
    for(j in 0:(number_along_each_dimension-1)) { 
      dimension_list[[j+1]] <- names_list[dimension_order[j*step_size + 1]]
    }
    dimension_list[[number_along_each_dimension]] <- names_list[dimension_order[number_of_images]]
    output_names[[i]] <- dimension_list
  }
  return(list("names"=output_names,"image_weights"=image_weights,"image_types"=image_types))
}