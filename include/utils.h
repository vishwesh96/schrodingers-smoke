#ifndef UTILS_H
#define UTILS_H
#include <complex>
#include <Comp2.h>
#include <GridPoint.h>
#include <fftw3.h>
#include <fstream>
#include <sstream>
#include <string>
#include <limits.h>
#include <algorithm>

#define NUM_VELOCITY_ITER 10

extern int NX;
extern int NY;
extern int NZ;
extern double LX;
extern double LY;
extern double LZ;
extern double H;
extern double DT;

typedef std::vector< std::vector < std::vector< GridPoint> > > grid;
double max_density = INT_MIN;
unsigned int NUM_TIME_STEPS;

struct Volume {
	 std::vector< std::vector < std::vector< bool> > > volume;
	 Eigen::Vector3d velocity;
	 std::vector<double> color;
	 double density;
};

bool is_zero(double p) {
	return(p<EPSILON && p>-EPSILON);
}

void scale_ifft(fftw_complex * out) {
	double scale = NX*NY*NZ;
	for (int i = 0;i< NX*NY*NZ; i++){
		out[i][0] /= scale;
		out[i][1] /= scale;
	}
}


void fftshift(grid & points) {
	
	for (unsigned int i = 0 ; i < points.size()/2 ; i++ ) {
		for(unsigned int j = 0; j< points[i].size()/2; j++) {
			for(unsigned int k = 0; k < points[i][j].size()/2; k++) {
				std::swap(points[i][j][k], points[i+NX/2][j+NY/2][k+NZ/2]);
				std::swap(points[i+NX/2][j][k], points[i][j+NY/2][k+NZ/2]);
				std::swap(points[i+NX/2][j+NY/2][k], points[i][j][k+NZ/2]);
				std::swap(points[i][j+NY/2][k], points[i+NX/2][j][k+NZ/2]);
			}
		}
	}

}


void print(grid &points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j< points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				// Comp2 psi = points[i][j][k].get_psi();
				printf("%d %d %d : %f\n",i,j,k,points[i][j][k].get_density());
			}
		}
	}
}

//Assuming only 1 length
void maintain_bound(Eigen::Vector3d &y) {
	double X = LX * NX, Y = LY * NY, Z = LZ * NZ;
	double epsilon = 0.000001;
	if(y(0)>=X)
		y(0) = y(0) - X + epsilon;
	if(y(1)>=Y) 
		y(1) = y(1) - Y + epsilon;
	if(y(2)>=Z) 
		y(2) = y(2) - Z + epsilon;
	if(y(0)<0) 
		y(0) = y(0) + X - epsilon;
	if(y(1)<0) 
		y(1) = y(1) + Y - epsilon;
	if(y(2)<0) 
		y(2) = y(2) + Z - epsilon;
}

std::vector<int> find_grid_cell(Eigen::Vector3d &y) {
	std::vector<int> n(3);
	n[0] = y(0)/LX;
	n[1] = y(1)/LY;
	n[2] = y(2)/LZ;
	return n;
}

bool is_zero_color(std::vector<double> color) {
	return (color[0]*color[0] + color[1]*color[1] + color[2]*color[2]) < EPSILON;
}

void resize_grid(grid & points) {
	points.resize(NX);
	for(int i=0;i<NX;i++) {
		points[i].resize(NY);
		for(int j=0;j<NY;j++) {
			points[i][j].resize(NZ);
		}
	}
}

void resize_volume(Volume & vol) {
	vol.volume.resize(NX);
	for(int i=0;i<NX;i++) {
		vol.volume[i].resize(NY);
		for(int j=0;j<NY;j++) {
			vol.volume[i][j].resize(NZ);
		}
	}
}

void density_to_vol(grid & points, int time_step, std::string & filename) {
	int num_channels_density = 1;
	int num_channels_color = 3;
	int nx = NX;
	int ny = NY;
	int nz = NZ;
	float xmin = 0;
	float ymin = 0;
	float zmin = 0;
	float xmax = LX * NX;
	float ymax = LY * NY;
	float zmax = LZ * NZ;
	int file_format = 3;
	int data_type = 1;
	std::stringstream ss_density,ss_color;
	ss_density<<"./vol/"<<filename<<"/density/"<<time_step<<".vol";
	ss_color<<"./vol/"<<filename<<"/color/"<<time_step<<".vol";

	std::ofstream file_density,file_color;
	file_density.open(ss_density.str(), std::ios::out | std::ios::binary);
	file_color.open(ss_color.str(), std::ios::out | std::ios::binary);

	file_density.write("VOL",3);
	file_density.write((char*)&file_format, 1);
	file_density.write((char*)&data_type, 4);
	file_density.write((char*)&nx, 4);
	file_density.write((char*)&ny, 4);
	file_density.write((char*)&nz, 4);
	file_density.write((char*)&num_channels_density, 4);
	file_density.write((char*)&xmin, 4);
	file_density.write((char*)&ymin, 4);
	file_density.write((char*)&zmin, 4);
	file_density.write((char*)&xmax, 4);
	file_density.write((char*)&ymax, 4);
	file_density.write((char*)&zmax, 4);


	file_color.write("VOL",3);
	file_color.write((char*)&file_format, 1);
	file_color.write((char*)&data_type, 4);
	file_color.write((char*)&nx, 4);
	file_color.write((char*)&ny, 4);
	file_color.write((char*)&nz, 4);
	file_color.write((char*)&num_channels_color, 4);
	file_color.write((char*)&xmin, 4);
	file_color.write((char*)&ymin, 4);
	file_color.write((char*)&zmin, 4);
	file_color.write((char*)&xmax, 4);
	file_color.write((char*)&ymax, 4);
	file_color.write((char*)&zmax, 4);

	for(int i =0 ;i<NZ;i++) {
		for(int j = 0;j < NY; j++) {
			for(int k = 0; k < NX; k++) {
				float tmp = points[k][j][i].get_density();
				tmp = tmp/(max_density);
				file_density.write((char*)&tmp,4);
			}
		}
	}


	for(int i =0 ;i<NZ;i++) {
		for(int j = 0;j < NY; j++) {
			for(int k = 0; k < NX; k++) {
				std::vector<double> tmp = points[k][j][i].get_color();
				for(int l=0;l<num_channels_color;l++) {
					float tmp_channel = tmp[l];
					file_color.write((char*)&tmp_channel,4);
				}
			}
		}
	}

	file_density.close();
	file_color.close();
}

#endif