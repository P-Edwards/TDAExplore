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
patch_landscapes_from_image <- function(image_name,number_of_patches,pixel_radius_for_patches=50,number_of_cores=1,patch_center_image=FALSE,return_patches=FALSE,remove_patch_data=TRUE,lower_threshold=0,upper_threshold=1,proportion_of_patch_sparse=.025,noise_threshold=.05,number_of_layers=50) {
  # Derived point samples sizes from parameters
  delta <- 2.5
  spacing_for_disc_points <- delta/2
  
  image_data <- OpenImageR::readImage(image_name) 
  if(length(dim(image_data)) > 2) { 
    image_data <- image_data[,,1]
  }
  image_data <- low_pixel_threshold(image_data,noise_threshold)
  image_length <- nrow(image_data)
  image_width <- ncol(image_data)
 
  
# generate centers
  if(!(patch_center_image==FALSE)) { 
    print(paste("Started image ", image_name))    
    print(paste("With image for patch centers: ",patch_center_image))
    center_image <- OpenImageR::readImage(patch_center_image)
    center_image <- low_pixel_threshold(center_image,noise_threshold)
    image_cdf <- create_cdf_vector(center_image)
    sampling_dims <- dim(center_image)    
  } else { 
    print(paste("Started image ", image_name))    
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
  average_intensities <- vector("double",ncol(center_points))
  for(i in 1:ncol(center_points)) { 
    patches[[i]] <- get_patch_from_point_and_mask(image_data,center_points[,i],mask)
    average_intensities[i] <- sum(patches[[i]]$data)/sum(which(patches[[i]]$data>0))
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
  print(paste("Finished image",image_name))
  
  if(return_patches==TRUE) { 
    return(list("data"=rval,"patches"=patches,"centroid"=centroid,"centers"=center_points))
  } else { 
    return(rval)
  }
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
#' @return Returns a list return_list with sufficient information to plot the image and its mask.
#' return_list$image is a matrix with pixel intensities as values
#' return_list$mask is a matrix of the same dimensions with score intensities as values
#' return_list$interval_rep is a vector of length 50 which compresses the mask into a radial representation
#' return_list$weights is a vector of the patch scores
#' @export
weight_image_using_landscapes_and_transform <- function(image_name,landscape_data,patches,radius,transformation_function=identity,min_weight=0,max_weight=1) { 
  image_data <- OpenImageR::readImage(image_name)
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

  return(list("image"=rotate(image_data),"mask"=rotate(total_weight_mask),"min"=unmod_min_weight,"max"=unmod_max_weight,"weights"=unmod_weights,interval_rep=aggregated_output))
}




