/*
   MLEM List Mode Reconstruction algorithm with resolution kernel/Time Of Flight Kernel
   13 June  2016 V1.0 P. Solevi and J. Gillam
   Modified by J. Renner in March 2020
   --------------------------------------------------------------------
   This code was originally developed within the framework of the project TOF-Brain-PET

   For any question about the code and its use please contact:
   -- paola.solevi@ovgu.de
   -- paola.solevi@gmail.com
   -- josh.e.renner@gmail.com
*/

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>

using namespace std;

// Constants
const float s2ps = 1000; // our times are in ns
const float pi = 3.141592653589793;
const float c_here = 0.3;

//--------- FUNCTIONS ------------------
float siddon(bool back_project, float * image, float projection_sum, float det1_X, float det1_Y, float det1_Z, float det2_X, float det2_Y, float det2_Z, float time_diff, float ox, float oy, float oz, float vx, float vy, float vz, int nx, int ny, int nz, bool TOF, float TOF_resolution);
float ToFFunction(float dist, float deltaT, float TOF_resolution);

//--------------------------------------
struct weightedVoxel {
     int index; //index of voxel
     float w;   //weight of voxel (length of intersection)
     float dt;  //time difference for voxel
};

extern "C"
float * MLEM_TOF_Reco(int niterations, bool TOF, float TOF_resolution,
	float FOV_XY, float FOV_Z, int NXY, int NZ,
	int ncoinc, float * LOR_X1, float * LOR_Y1, float * LOR_Z1, float * LOR_T1,
	float * LOR_X2, float * LOR_Y2, float * LOR_Z2, float * LOR_T2,
  const char * outfile_prefix, int out_niter) {

	FILE* outfile;
	char outfile_name[100];

	int nvoxels = NXY * NXY * NZ;

	float  ox = -FOV_XY / 2.;  		// location of the first voxel plane in units of measurement - for the Siddon kernel
	float  oy = -FOV_XY / 2.;
	float  oz = -FOV_Z  / 2.;

	int nx = NXY;			// FOV in voxels
	int ny = NXY;
	int nz = NZ;

	float dxy = float(NXY) / FOV_XY;	// voxel/mm
	float dz  = float(NZ) / FOV_Z;

	float vx = 1. / dxy;			// mm/voxel
	float vy = 1. / dxy;
	float vz = 1. / dz;

	float  d1X, d1Y, d1Z, time1;
	float  d2X, d2Y, d2Z, time2;

	// ====== IMAGES ===================================
	float* NEW_IMAGE    = (float*)malloc(nvoxels*sizeof(float));	//IMAGE (k)
	float* IMAGE	      = (float*)malloc(nvoxels*sizeof(float));	//IMAGE (k-1)
	for(int iV = 0; iV < nvoxels; iV++) { IMAGE[iV] = 1.;       }      // Unitary image (k)
	for(int iV = 0; iV < nvoxels; iV++) { NEW_IMAGE[iV] = 0.;   }  // Zero image (k+1)

	// ====== ITERATIONS ==============================
	float time_diff;

	// START ITERATING
	for(int iter = 0; iter < niterations; iter++) {

		float P = 0;

    // Read the coincidences LOR by LOR.
    for(int k_lor = 0; k_lor < ncoinc; k_lor++) {

      if ( k_lor % 10000000 == 0 ) cout << "Iter= " << iter << " Processing " << k_lor << " coincidences " << endl;

			float projection = 0;

			// SAMPLING THE LOR WITHING THE DETECTION VOLUME NLINE TIMES OR USE THE GLOBAL POS PROVIDED BY GATE
			d1X = LOR_X1[k_lor]; d1Y = LOR_Y1[k_lor]; d1Z = LOR_Z1[k_lor]; time1 = LOR_T1[k_lor];
			d2X = LOR_X2[k_lor]; d2Y = LOR_Y2[k_lor]; d2Z = LOR_Z2[k_lor]; time2 = LOR_T2[k_lor];

			// all time variables are in ps, to allow floating precision along with the code.
			time_diff = (time1*s2ps - time2*s2ps); // always negative

    	projection = siddon(false, IMAGE, 0., d1X, d1Y, d1Z, d2X, d2Y, d2Z, time_diff,
																ox, oy, oz, vx, vy, vz, nx, ny, nz, TOF, TOF_resolution);

	    //if (isnan(projection))
	    //    {projection = 0;}
	    //if (isinf(projection))
	    //    {projection = 0;}

	    if (projection > 0) {
        siddon(true, NEW_IMAGE, projection, d1X, d1Y, d1Z, d2X, d2Y, d2Z, time_diff,
									ox, oy, oz, vx, vy, vz, nx, ny, nz, TOF, TOF_resolution);
	    }

			P += projection;
		}

		for(int iV = 0; iV < nvoxels; iV++) {

			// Apply the sensitivity matrix.
			//if(SENS[iV] > 0) NEW_IMAGE[iV] = (NEW_IMAGE[iV] * IMAGE[iV]) / SENS[iV];
			//else NEW_IMAGE[iV] = 0.;

			NEW_IMAGE[iV] = (NEW_IMAGE[iV] * IMAGE[iV]);
			IMAGE[iV] = NEW_IMAGE[iV];
		}

		// Save the first, last, and every out_niter to file.
		if (out_niter > 0 && ((iter % out_niter) == 0 || iter == niterations)) {
			sprintf(outfile_name,"%s%d.raw",outfile_prefix,iter);
			outfile = fopen(outfile_name,"wb");
	 		fwrite(NEW_IMAGE,sizeof(float),nvoxels,outfile);
	  }

		for(int iV = 0; iV < nvoxels; iV++) {
			NEW_IMAGE[iV] = 0.;
		}

		cout << "Completed iteration : " << iter << " reading total coincidences: " << ncoinc << endl;
	}

	delete NEW_IMAGE;
	return IMAGE;
}


//-------------------------------------- Siddon Kernel ------------------------------------------------
float ToFFunction(float dist, float deltaT, float TOF_resolution) {

  float ToFweight = 0;
  float deltaT_line = ( dist + (deltaT*c_here/2.0) );
  float line_resolution_SIGM = TOF_resolution*c_here/2.;

  const float std_trunk = 3.;     //very important to cut the gaussian and give the elements some compact support

  if (abs( deltaT_line ) < (std_trunk*line_resolution_SIGM) ) {
    ToFweight = (1.0/sqrt(2.0*pi*line_resolution_SIGM*line_resolution_SIGM)) * exp( - ((deltaT_line * deltaT_line)/(2*line_resolution_SIGM*line_resolution_SIGM)));\
  }
  else {
    ToFweight = 0;
  }

   //cout << dist << " " << deltaT << " " << deltaT_line << " < " << (time_resolution_SIGM*c_here/2.) << " " << ToFweight  << endl;

   return ToFweight;
}

/**
 Forward or back project using the Siddon algorithm.
 Note that parameter image is treated differently depending on whether
 back_project is true or false.

 If back_project is true, compute the back projection of this line of
 response and store in the array pointed to by input parameter image. In this
 case bp_projection_sum is the corresponding forward projection projection
 sum computed in advance, and the return value will always be 0.

 If back_project is false, compute and return the total forward projection integral of
 the image pointed to by input parameter image for this line of response. In
 this case, parameter bp_projection_sum is ignored.

*/
float siddon(bool back_project, float * image, float bp_projection_sum,
		float det1_X, float det1_Y, float det1_Z, float det2_X, float det2_Y,
		float det2_Z, float time_diff, float ox, float oy, float oz, float vx,
		float vy, float vz, int nx, int ny, int nz, bool TOF, float TOF_resolution){

  float accumulator = 0.0;

  float raylength;
  float SA = 1;
  float X1[3];
  float X2[3];
  X1[0] = det1_X;
  X1[1] = det1_Y;
  X1[2] = det1_Z;
  X2[0] = det2_X;
  X2[1] = det2_Y;
  X2[2] = det2_Z;

  //calculation of distance between the 2 detectors
  raylength = sqrt((X2[0]-X1[0])*(X2[0]-X1[0])+(X2[1]-X1[1])*(X2[1]-X1[1])+(X2[2]-X1[2])*(X2[2]-X1[2]));

  float alpha_x_a = (ox - X1[0]) / (X2[0] - X1[0]);
  float alpha_y_a = (oy - X1[1]) / (X2[1] - X1[1]);
  float alpha_z_a = (oz - X1[2]) / (X2[2] - X1[2]);

  float alpha_x_b = (ox + vx * nx - X1[0]) / (X2[0] - X1[0]);
  float alpha_y_b = (oy + vy * ny - X1[1]) / (X2[1] - X1[1]);
  float alpha_z_b = (oz + vz * nz - X1[2]) / (X2[2] - X1[2]);


  // CALCULATE ALPHA_MIN AND ALPHA_MAX
     //alpha_min is the first touch of the ray with the volume if the ray
     //comes from outside, otherwise it is zero; alpha_max is the last touch
     //of the ray with the volume.
     //in case the ray came from outside, we will also find out which side
     //was it first (this is of essential importance for the continuation
     //of the algorithm)
  float temp = 0;
  float alpha_min = fmax(fmax(temp, fminf(alpha_x_a, alpha_x_b)), fmax(fminf(alpha_y_a, alpha_y_b), fminf(alpha_z_a, alpha_z_b)));

  /* RAY_FIRST_HIT
    saves which plane the ray hit first. if it
    came from inside then it will be -1*/
  int ray_first_hit=-1;
  if (alpha_min > 0)
  {
    // Ray came from outside
    if (alpha_min == fminf(alpha_x_a, alpha_x_b)) {
      ray_first_hit = 1;
    }
    else if (alpha_min == fminf(alpha_y_a, alpha_y_b)) {
      ray_first_hit = 2;
    }
    else if (alpha_min == fminf(alpha_z_a, alpha_z_b)) {
      ray_first_hit = 3;
    }
  }
  else {
    //Ray came from inside
    ray_first_hit = -1;
  }

  //Updates for alpha_x and alpha_y
  float delta_alpha_x = vx / fabs(X2[0]-X1[0]);
  float delta_alpha_y = vy / fabs(X2[1]-X1[1]);
  float delta_alpha_z = vz / fabs(X2[2]-X1[2]);

  int di, dj, dk;
  //signum function for (X2[0]-X1[0])
  if ( X2[0] > X1[0] ) {
    di = 1;
  }
  else if ( X1[0] == X2[0] ) {
    di = 0;
  }
  else {
    di = -1;
  }

  //signum function for (X2[1]-X1[1])
  if ( X2[1] > X1[1] ) {
    dj = 1;
  }
  else if ( X1[1] == X2[1] ) {
    dj = 0;
  }
  else {
    dj = -1;
  }

  //signum function for (X2[2]-X1[2])
  if ( X2[2] > X1[2] ) {
    dk = 1;
  }
  else if ( X1[2] == X2[2] ) {
    dk = 0;
  }
  else {
    dk = -1;
  }

  /*  CALCULATE I_MIN, J_MIN, K_MIN AND ALPHA_X, ALPHA_Y, ALPHA_Z
     I_MIN, J_MIN, K_MIN
         the min values define the plane where the ray first hits the volume
         this means that the one of the min values corresponds to the plane
         where the ray ENTERS the volume, and the other two min values are
         hits INSIDE the volume.

         we find the index i in the x-direction by setting
         x(alpha_min) = ox + (i_min - 1)*dx
         which is equivalent to
         i_min = (x(alpha_min) - ox)/dx + 1

         unfortunately, there is a numerical problem when we calculate the
         intersection with the first plane... one of the coordinates should
         be exactly a grid point - but, it turns out that it can be slightly
         off grid. like 1.0000000000001 or something like that. then, our
         methods for ceil and floor fail. thats why we need to know which of
         the coordinates is the first one and we will not ceil nor floor it
         but just round it (since it is soooo close to its expected value
         anyway)

     ALPHA_X, ALPHA_Y, ALPHA_Z
         these are the alpha values that correspond to the first hit of
         the ray _inside_ the volume. for example, if the ray entered
         from the x-direction and i_min=1, then the alpha_x value would
         correspond to the x-plane with index 2. the values alpha_y and
         alpha_z would correspond to the respective plane with index
         j_min and k_min.*/

  float i_min = (X1[0] + alpha_min * (X2[0]-X1[0]) - ox) / vx;
  float j_min = (X1[1] + alpha_min * (X2[1]-X1[1]) - oy) / vy;
  float k_min = (X1[2] + alpha_min * (X2[2]-X1[2]) - oz) / vz;

  float alpha_x, alpha_y, alpha_z;

  if (alpha_min > 0) {

    // Ray came from outside
    if (ray_first_hit == 1) {

      //Now we know that ray entered from the X-direction
      int i_temp=(int)floorf(i_min + 0.5);//Rounding of i_min
      i_min=i_temp;

      if (X1[1] < X2[1]) {
        j_min = ceilf(j_min);
      }
      else if (X1[1] > X2[1]) {
        j_min = floorf(j_min);
      }
      else {
        j_min = numeric_limits<float>::max();
      }

      if (X1[2] < X2[2]) {
        k_min = ceilf(k_min);
      }
      else if (X1[2] > X2[2]) {
        k_min = floorf(k_min);
      }
      else {
        k_min = numeric_limits<float>::max(); // numeric_limits<float>::max();
      }

      alpha_x = (ox + (i_min + di) * vx - X1[0]) / (X2[0]-X1[0]);
      alpha_y = (oy + (j_min     ) * vy - X1[1]) / (X2[1]-X1[1]);
      alpha_z = (oz + (k_min     ) * vz - X1[2]) / (X2[2]-X1[2]);
    }
    else if (ray_first_hit == 2) {

      // Now we know that ray entered from y-direction
      int j_temp=(int)floorf(j_min + 0.5);    //Rounding of j_min
      j_min=j_temp;

      if (X1[0] < X2[0]) {
        i_min = ceilf(i_min);
      }
      else if (X1[0] > X2[0]) {
        i_min = floorf(i_min);
      }
      else {
				i_min = numeric_limits<float>::max(); //numeric_limits<float>::max();
      }

      if (X1[2] < X2[2]) {
        k_min = ceilf(k_min);
      }
      else if (X1[2] > X2[2]) {
        k_min = floorf(k_min);
      }
      else {
        k_min = numeric_limits<float>::max(); //numeric_limits<float>::max();
      }

      alpha_x = (ox +  i_min       * vx - X1[0]) / (X2[0]-X1[0]);
      alpha_y = (oy + (j_min + dj) * vy - X1[1]) / (X2[1]-X1[1]);
      alpha_z = (oz +  k_min       * vz - X1[2]) / (X2[2]-X1[2]);
    }
    else if (ray_first_hit == 3) {

      // Now we know that ray entered from z direction
      int k_temp=(int)floorf(k_min + 0.5);    //Rounding of z_min
      k_min = k_temp;

      if (X1[0] < X2[0]) {
        i_min = ceilf(i_min);
      }
      else if (X1[0] > X2[0]) {
        i_min = floorf(i_min);
      }
      else {
      	i_min = numeric_limits<float>::max(); //numeric_limits<float>::max();
      }

      if (X1[1] < X2[1]) {
        j_min = ceilf(j_min);
      }
      else if (X1[1] > X2[1]) {
        j_min = floorf(j_min);
      }
      else {
        j_min = numeric_limits<float>::max(); //numeric_limits<float>::max();
      }

      alpha_x = (ox +  i_min        * vx - X1[0]) / (X2[0]-X1[0]);
      alpha_y = (oy +  j_min        * vy - X1[1]) / (X2[1]-X1[1]);
      alpha_z = (oz + (k_min + dk ) * vz - X1[2]) / (X2[2]-X1[2]);
    }
  }
  else {

    //Ray came from inside
    if (X1[0]<X2[0]) {
      i_min = ceilf(i_min);
    }
    else if (X1[0]>X2[0]) {
      i_min = floorf(i_min);
    }
    else {
      i_min = numeric_limits<float>::max(); //numeric_limits<float>::max();
    }

    if (X1[1]<X2[1]) {
      j_min = ceilf(j_min);
    }
    else if (X1[1]>X2[1]) {
      j_min = floorf(j_min);
    }
    else {
      j_min = numeric_limits<float>::max();//numeric_limits<float>::max();
    }

    if (X1[2] < X2[2]) {
      k_min = ceilf(k_min);
    }
    else if (X1[2] > X2[2]) {
      k_min = floorf(k_min);
    }
    else {
      k_min = numeric_limits<float>::max(); //numeric_limits<float>::max();
    }

    alpha_x = (ox + i_min * vx - X1[0]) / (X2[0]-X1[0]);
    alpha_y = (oy + j_min * vy - X1[1]) / (X2[1]-X1[1]);
    alpha_z = (oz + k_min * vz - X1[2]) / (X2[2]-X1[2]);
  }

  float alpha_cur = alpha_min;

  //Calculate indices of starting voxel with a centered alpha
  float alpha_c = (alpha_min + fminf(fminf(alpha_x, alpha_y), alpha_z))/2;
  int i = (int)(floorf((X1[0] + alpha_c * (X2[0] - X1[0]) - ox)/vx));
  int j = (int)(floorf((X1[1] + alpha_c * (X2[1] - X1[1]) - oy)/vy));
  int k = (int)(floorf((X1[2] + alpha_c * (X2[2] - X1[2]) - oz)/vz));

  int finished = 0; //used as boolean as it is not clear if boolean are suported
  float dist ;
  int index ;
  float alpha_dt0 = 0.5;
	weightedVoxel w_voxel;
  float dist_from_center;

  //NOTE TO SELF: an alpha of 0.5 indicates that the point is directly between the two interaction points -- deltaT = 0: this should incorportate ToF nicely!!
	while (finished==0 && i >= 0 && i < nx && j >= 0 && j < ny && k >=0 && k < nz) {

    if (alpha_x <= min(alpha_y, alpha_z)) {

      //We have an x-crossing
      if (alpha_x > 1) {
        alpha_x = 1;
        finished = 1;
      }

      dist = (alpha_x - alpha_cur)*raylength;
      index = (int) (i + j*nx + k*nx*ny);
      dist_from_center = (dist/2.0)  + ((alpha_dt0 - alpha_x)*raylength); //

      alpha_cur = alpha_x;
      alpha_x = alpha_x + delta_alpha_x;
      i += di;

    }
    else if (alpha_y <= min(alpha_x, alpha_z)) {

      //We have a y crossing
      if (alpha_y > 1) {
          alpha_y = 1;
          finished = 1;
      }
      dist = (alpha_y - alpha_cur)*raylength;
      index = (int) (i + j*nx + k*nx*ny);
      dist_from_center = (dist/2.0)  + ((alpha_dt0 - alpha_y)*raylength); //

      alpha_cur = alpha_y;
      alpha_y = alpha_y + delta_alpha_y;
      j += dj;

    }
    else if (alpha_z <= min(alpha_x, alpha_y)) {

      //We have a z crossing
      if (alpha_z > 1) {
          alpha_z = 1;
          finished = 1;
      }
      dist = (alpha_z - alpha_cur)*raylength;
      index = (int) (i + j*nx + k*nx*ny);
      dist_from_center = (dist/2.0)  + ((alpha_dt0 - alpha_z)*raylength); //

      alpha_cur = alpha_z;
      alpha_z = alpha_z + delta_alpha_z;
      k += dk;
    }

    w_voxel.index=index;
    w_voxel.w=dist;
    if(TOF) {
      w_voxel.dt = ToFFunction(dist_from_center, time_diff, TOF_resolution);
    }
    else {
      w_voxel.dt = 1.0;
    }

	  // each element of weighted voxel, so each voxel with its intersection
	  // length will directly be used for the calculation of different types of
	  // the algorithm. the algorithm that should be applied is defined by
	  // the variable calc_type. following values are possible:
	  // SA = 2*atan2(100, fmax((raylength/2.0)+dist_from_center, (raylength/2.0)-dist_from_center));

		// Add to the back projection if selected.
		if(back_project) {
			image[w_voxel.index]+= w_voxel.dt*w_voxel.w / bp_projection_sum;
		}
		// Otherwise, add to the projection sum.
		else {
			accumulator += SA*w_voxel.dt*w_voxel.w*image[w_voxel.index];
		}

	}

	// Return the projection sum (0 if back_project is true).
	return accumulator;
}
