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


# Added by Peter Bubenik May 3, 2020; updated with C++ versions Parker Edwards May 9
sample_patch_using_top_intensities_sparse <- function(patch,proportion_of_patch_sparse=.025,exclusion_disc_data=create_disc_mask(2.5)) { 
  nonzero_pixels <- sum(patch$data != 0)
  number_of_pixels_to_sample <- ceiling( nonzero_pixels * proportion_of_patch_sparse )
  pixels_in_order <- order(patch$data,decreasing=TRUE)[1:nonzero_pixels]
  return(sparse_patch_top_intensities(pixels_in_order,number_of_pixels_to_sample,nrow(patch$data),ncol(patch$data),exclusion_disc_data))
}


sample_centers_using_cdf_sparse <- function(image_cdf,number_of_pixels_to_sample,dimensions,delta) { 
  exclusion_disc_data <- create_disc_mask(delta)
  number_of_pixels_currently_sampled <- 0
  sampled_pixels_xy <- list()
  potential_samples <- vector(mode="integer",length=3*number_of_pixels_to_sample)
  random_samples <- stats::runif(3*number_of_pixels_to_sample,0.0,image_cdf[length(image_cdf-1)])
  potential_samples <- findInterval(random_samples,image_cdf,left.open = TRUE);
  return(sparse_patch_top_intensities(potential_samples,number_of_pixels_to_sample,dimensions[1],dimensions[2],exclusion_disc_data))
}