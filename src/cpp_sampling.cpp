// R Package exPLore

// Copyright (C) 2021 Parker Edwards

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 3 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
// USA

#include <Rcpp.h>
#include <math.h>
#include <algorithm>
#include <cstdio>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
IntegerMatrix sparse_patch_top_intensities(IntegerVector pixels_in_order,int number_of_pixels_to_sample,int data_rows,int data_cols, LogicalMatrix exclusion_disc_data) { 

  LogicalMatrix proximity_mask(data_rows,data_cols);

  int number_of_pixels_currently_sampled = 0;
  IntegerMatrix sampled_pixels_xy(2,number_of_pixels_to_sample);
  
  for(int index=0;index<pixels_in_order.size();index++) {
    if(number_of_pixels_currently_sampled >= number_of_pixels_to_sample) { 
      break;
    }
    // Note: Not shifting +1 at end because C++ indices start from 0,
    // but pixels_in_order (having been returned from R and so indexed from 1) 
    // still need to be shifted down 1 before computation
    int index_1 = ((pixels_in_order[index] - 1) % data_rows);
    int index_2 = ((pixels_in_order[index] - 1)/data_rows);
    // If proximity mask is 1, then skip the point. Otherwise...
    if(proximity_mask(index_1,index_2) == FALSE) { 
      // Collect sample point
      sampled_pixels_xy(0,number_of_pixels_currently_sampled) = index_1;
      sampled_pixels_xy(1,number_of_pixels_currently_sampled) = index_2;
      number_of_pixels_currently_sampled += 1;
      // Add in disc mask
      int length = data_rows;
      int width  = data_cols;      
      double disc_length = (double) exclusion_disc_data.nrow();
      double disc_width = (double) exclusion_disc_data.ncol();
      int middle_point_1 = (int) ceil(disc_length/2);
      int middle_point_2 = (int) ceil(disc_width/2);
      int image_row_minimum = max(index_1 - middle_point_1,0);
      int image_row_maximum = min(index_1+(((int) disc_length)-middle_point_1),length);
      int image_col_minimum = max(index_2 - middle_point_2,0);
      int image_col_maximum = min(index_2+(((int) disc_width)-middle_point_2),width);
      int disc_row_minimum = middle_point_1-(index_1-image_row_minimum);
      int disc_row_maximum = middle_point_1+(image_row_maximum-index_1);
      int disc_col_minimum = middle_point_2-(index_2-image_col_minimum);
      int disc_col_maximum = middle_point_2+(image_col_maximum-index_2);
      for(int i = disc_row_minimum;i<disc_row_maximum;i++) { 
      	for(int j = disc_col_minimum;j<disc_col_maximum;j++) { 
      		proximity_mask(image_row_minimum+i,image_col_minimum+j) = max(exclusion_disc_data(i,j),proximity_mask(image_row_minimum+i,image_col_minimum+j));
      	}
      }
    }
  }
  return sampled_pixels_xy(_,Range(0,max(number_of_pixels_currently_sampled-1,0)));
}
