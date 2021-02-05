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

# Turn the lowest threshold_percentage pixels (by percentage of intensity range) to 0
low_pixel_threshold <- function(image,threshold_percentage) { 
  max <- range(image)[2]; 
  intensity_threshold <- max*threshold_percentage; 
  image[ image < intensity_threshold] <- 0.0; 
  return(image); 
}

create_cdf_vector <- function(image){ 
  flattened_image <- as.vector(image); 
  cdf_vector <- 1:length(flattened_image); 
  sum <- 0; 
  for(i in 1:length(flattened_image)) { 
    sum <- sum + flattened_image[i];
    cdf_vector[i] <- sum; 
  }
  return(cdf_vector);
}

create_disc_mask <- function(radius) { 
  output <- matrix(nrow=(2*radius+1),ncol=(2*radius+1))
  center_point <- c(ceiling((2*radius+1)/2),ceiling((2*radius+1)/2))
  mask_function <- function(i,j) { 
    if(stats::dist(rbind(c(i,j),center_point))[1] > radius) { 
      return(0.0); 
    } else {
      return(1.0);
    }
  }
  rval <- matrix(mapply(mask_function,row(output),col(output)),nrow=nrow(output),ncol=ncol(output))
  rval <- as.vector(rval,mode="logical")
  dim(rval) <- c((2*radius+1),(2*radius+1)) 
  return(rval)
}


get_patch_from_point_and_mask <- function(image,point,mask) { 
  radius <- floor(dim(mask)[1]/2)
  
  patch_lower_x <- max(point[1]-radius+1,1)
  patch_upper_x <- min(point[1]+radius,dim(image)[1])
  patch_lower_y <- max(point[2]-radius+1,1)
  patch_upper_y <- min(point[2]+radius,dim(image)[2])

  mask_stop_x <- patch_upper_x - patch_lower_x + 1
  mask_stop_y <- patch_upper_y - patch_lower_y + 1

  rectangular_patch <- image[patch_lower_x:patch_upper_x,patch_lower_y:patch_upper_y]
  width <- dim(rectangular_patch)[1]
  len <- dim(rectangular_patch)[2]
  mask_start_x <- mask_stop_x- width + 1
  mask_start_y <- mask_stop_y- len + 1
  prod_mask <- as.double(mask[mask_start_x:mask_stop_x,mask_start_y:mask_stop_y]) 
  dim(prod_mask) <- dim(mask[mask_start_x:mask_stop_x,mask_start_y:mask_stop_y])
  
  patch <- list("data" = rectangular_patch*prod_mask, "center" = point,"xmin"=patch_lower_x,"xmax"=patch_upper_x,"ymin"=patch_lower_y,"ymax"=patch_upper_y);
  return(patch); 
}


# ceil(2pi/arccos(1 - eps/2)) for number of roots of unity to get an epsilon-spacing cover of circle
# adding a circular boundary to a point cloud sample of a patch
patch_sample_with_circ_boundary <- function(sample, radius, epsilon){
  n <- ceiling(2*pi/acos(1 - (epsilon^2/(2*radius^2))))
  length <- dim(sample)[1]
  width <- dim(sample)[2]
  center <- c(radius,radius)
  x_coordinates<- vector(mode='numeric', length=n)
  for (i in 1:n){
    x_coordinates[i]<-radius*cos(2*pi*i/n)
  }
  x_coordinates <- x_coordinates+center[1]
  
  y_coordinates<- vector(mode='numeric', length=n)
  for (i in 1:n){
    y_coordinates[i]<-radius*sin(2*pi*i/n)
  }
  y_coordinates <- y_coordinates+center[2]

  circle_boundary <- rbind(x_coordinates, y_coordinates)
  sample_with_circular_boundary<-cbind(sample,circle_boundary)
  return(sample_with_circular_boundary);
}

rotate <- function(x) t(apply(x,2,rev))