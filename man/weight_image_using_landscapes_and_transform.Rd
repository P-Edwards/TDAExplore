% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline_functions.R
\name{weight_image_using_landscapes_and_transform}
\alias{weight_image_using_landscapes_and_transform}
\title{Convolve an image using patch scores}
\usage{
weight_image_using_landscapes_and_transform(
  image_name,
  landscape_data,
  patches,
  radius,
  transformation_function = identity,
  min_weight = 0,
  max_weight = 1,
  interval_rep = FALSE
)
}
\arguments{
\item{image_name}{File path (absolute or relative to current directory) to the image. Supported file formats are same as OpenImageR: .png, .jpeg
.jpg, and .tiff}

\item{landscape_data}{Landscapes for image patches e.g. like those returned from patch_landscapes_for_image(return_patches=TRUE)$data}

\item{patches}{Patches e.g. like those return from patch_landscapes_for_image(return_patches=TRUE)$patches}

\item{radius}{Radius of patches in pixels}

\item{transformation_function}{A function which takes as input a matrix and outputs a vector with scores, one for each row in the matrix}

\item{interval_rep}{If TRUE, computes the interval representation described below. Default is FALSE.}
}
\value{
Returns a list return_list with sufficient information to plot the image and its mask.
return_list$image is a matrix with pixel intensities as values
return_list$mask is a matrix of the same dimensions with score intensities as values
return_list$interval_rep is a vector of length 50 which compresses the mask into a radial representation
return_list$weights is a vector of the patch scores
}
\description{
Combines landscape data, data for positioning patches on an image, and a scoring function for landscapes
to produce a mask for the image.
A few parameters are undocumented and it is inadvisable to change them.
}
