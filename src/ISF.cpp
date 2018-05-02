#include <fftw3.h>
#include <GridPoint.h>
#include <fstream>
#include <sstream>
#include <string>
#include <limits.h>
#include <algorithm>

#define NX 128
#define NY 64
#define NZ 64
#define LX 0.078125
#define LY 0.078125
#define LZ 0.078125
#define NUM_TIME_STEPS 10
#define DT 1.0/24
#define H 0.1
#define EPSILON 0.01

typedef std::vector< std::vector < std::vector< GridPoint> > > grid;

double min_density = INT_MAX;
double max_density = INT_MIN;

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

void print(grid &points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				// Comp2 psi = points[i][j][k].get_psi();
				printf("%d %d %d : %f\n",i,j,k,points[i][j][k].get_density());
			}
		}
	}
}

void convert_div(fftw_complex *in, grid  & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				in[k + NZ *(j + NY * i)][0] = points[i][j][k].get_div();
				in[k + NZ *(j + NY * i)][1] = 0;
			}
		}
	}
}

void convert_psi1(fftw_complex * in, grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				in[k + NZ *(j + NY * i)][0] = psi.get_z1().real();
				in[k + NZ *(j + NY * i)][1] = psi.get_z1().imag();
			}
		}
	}
}

void convert_psi2(fftw_complex * in, grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				in[k + NZ *(j + NY * i)][0] = psi.get_z2().real();
				in[k + NZ *(j + NY * i)][1] = psi.get_z2().imag();
			}
		}
	}
}

void reconvert_psi1(fftw_complex * out, grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				psi.set_z1(Complex(out[k + NZ *(j + NY * i)][0], out[k + NZ *(j + NY * i)][1]));
				points[i][j][k].set_psi(psi);
			}
		}
	}
}

void reconvert_psi2(fftw_complex * out, grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				psi.set_z2(Complex(out[k + NZ *(j + NY * i)][0], out[k + NZ *(j + NY * i)][1]));
				points[i][j][k].set_psi(psi);
			}
		}
	}
}

void modify_div(fftw_complex * out, fftw_complex * in, grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				double lambda = (-4.0)*( pow(sin(M_PI*i/NX),2)/(LX*LX) + pow(sin(M_PI*j/NY),2)/(LY*LY) + pow(sin(M_PI*k/NZ),2)/(LZ*LZ) );
			
				if(!is_zero(lambda)){
					in[k + NZ *(j + NY * i)][0] = out[k + NZ *(j + NY * i)][0]/lambda;
					in[k + NZ *(j + NY * i)][1] = out[k + NZ *(j + NY * i)][1]/lambda;
				} else {
					in[k + NZ *(j + NY * i)][0] = 0.0;
					in[k + NZ *(j + NY * i)][1] = 0.0;
				}
			}
		}
	}
}

void reconvert_pressure(fftw_complex * out, grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				double q = out[k + NZ *(j + NY * i)][0];
				if(!is_zero(out[k + NZ *(j + NY * i)][1])) {
					printf("Pressure is not real\n");
					exit(0);
				}
				Complex factor(cos(q),-sin(q));
				psi.set_z1(psi.get_z1() * factor);
				psi.set_z2(psi.get_z2() * factor);
				points[i][j][k].set_pressure(q);
				points[i][j][k].set_psi(psi);
			}
		}
	}
}

void update_velocity(grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				Comp2 psi_x = points[(i+1)%NX][j][k].get_psi();
				Comp2 psi_y = points[i][(j+1)%NY][k].get_psi();
				Comp2 psi_z = points[i][j][(k+1)%NZ].get_psi();
				Eigen::Vector3d v = H*Eigen::Vector3d(std::arg(psi.inner_product(psi_x))/LX,std::arg(psi.inner_product(psi_y))/LY,std::arg(psi.inner_product(psi_z))/LZ);
				points[i][j][k].set_v(v);
			}
		}
	}
}

//Assuming only 1 length
void maintain_bound(Eigen::Vector3d & y) {
	double X = LX * NX, Y = LY * NY, Z = LZ * NZ;
	if(y(0)>X) 
		y(0) -= X;
	if(y(1)>Y) 
		y(1) -= Y;
	if(y(2)>Z) 
		y(2) -= Z;
	if(y(0)<0) 
		y(0) += X;
	if(y(1)<0) 
		y(1) += Y;
	if(y(2)<0) 
		y(2) += Z;
}

std::vector<int> find_grid_cell(Eigen::Vector3d y) {
	std::vector<int> n(3);
	n[0] = y(0)/LX;
	n[1] = y(1)/LY;
	n[2] = y(2)/LZ;
	return n;
}

double interpolate_density(grid & points,Eigen::Vector3d y){
	maintain_bound(y);
	std::vector<int> n = find_grid_cell(y);

	double d1 = points[n[0]][n[1]][n[2]].get_density();
	double d2 = points[(n[0]+1)%NX][n[1]][n[2]].get_density();
	double d3 = points[(n[0]+1)%NX][(n[1]+1)%NY][n[2]].get_density();
	double d4 = points[n[0]][(n[1]+1)%NY][n[2]].get_density();
	double d5 = points[n[0]][n[1]][(n[2]+1)%NZ].get_density();
	double d6 = points[(n[0]+1)%NX][n[1]][(n[2]+1)%NZ].get_density();
	double d7 = points[(n[0]+1)%NX][(n[1]+1)%NY][(n[2]+1)%NZ].get_density();
	double d8 = points[n[0]][(n[1]+1)%NY][(n[2]+1)%NZ].get_density();

	double z1 = n[0]*LX;
	double z2 = n[1]*LY;
	double z3 = n[2]*LZ;
	double z4 = (n[0]+1)*LX;
	double z5 = (n[1]+1)*LY;
	double z6 = (n[2]+1)*LZ;

	double y0 = y(0); 
	double y1 = y(1); 
	double y2 = y(2);

	double density = ( (z6-y2) * (d1*(z4-y0)*(z5-y1) + d2*(y0-z1)*(z5-y1) + d3*(y0-z1)*(y1-z2) + d4*(z4-y0)*(y1-z2)) + (y2-z3)*(d5*(z4-y0)*(z5-y1) + d6*(y0-z1)*(z5-y1) + d7*(y0-z1)*(y1-z2) + d8*(z4-y0)*(y1-z2)) )/(LX*LY*LZ);
	return density;
}

void advect_density(grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Eigen::Vector3d x(i*LX,j*LY,k*LZ);
				Eigen::Vector3d v = points[i][j][k].get_v();
				Eigen::Vector3d y = x - v*DT;
				double density = interpolate_density(points,y);
				min_density = std::min(density, min_density);
				max_density = std::max(density, max_density);
				points[i][j][k].set_density(density);
			}
		}
	}
}

void mult_exp(grid & points) {
for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				// double lambda = -4.0*M_PI*M_PI*( pow( ((i-64.0)/(10.0)),2) + pow( ((j-32.0)/(5.0)),2) + pow(((k-32.0)/(5.0)),2) );
				double lambda = -4.0*M_PI*M_PI*( pow( (i/(LX*NX)),2) + pow( (j/(LY*NY)),2) + pow((k/(LZ*NZ)),2) );
				double theta = 0.5 * lambda * DT * H;
				Complex factor(cos(theta),sin(theta));
				Comp2 psi = points[i][j][k].get_psi();
				psi.set_z1(psi.get_z1() * factor);
				psi.set_z2(psi.get_z2() * factor);
				points[i][j][k].set_psi(psi);
			}
		}
	}
}

void normalize(grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				points[i][j][k].normalize();
			}
		}
	}
}

void schrodinger(fftw_complex * in, fftw_complex * out, fftw_plan fp, fftw_plan bp, grid & points) {
	convert_psi1(in,points);
	fftw_execute(fp);
	reconvert_psi1(out,points);
	convert_psi2(in,points);
	fftw_execute(fp);
	reconvert_psi2(out,points);

	mult_exp(points);

	convert_psi1(in,points);
	fftw_execute(bp);
	scale_ifft(out);
	reconvert_psi1(out,points);
	convert_psi2(in,points);
	fftw_execute(bp);
	scale_ifft(out);
	reconvert_psi2(out,points);
}

void pressure_project(fftw_complex * in, fftw_complex * out, fftw_plan fp, fftw_plan bp, grid & points) {
	
	// Set up divergence for each grid point
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j < points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				double nxp = std::arg(points[i][j][k].get_psi().inner_product(points[(i+1)%NX][j][k].get_psi())); 
				double nxn = std::arg(points[i][j][k].get_psi().inner_product(points[(NX + i-1)%NX][j][k].get_psi()));
				double nyp = std::arg(points[i][j][k].get_psi().inner_product(points[i][(j+1)%NY][k].get_psi()));
				double nyn = std::arg(points[i][j][k].get_psi().inner_product(points[i][(NY+j-1)%NY][k].get_psi()));
				double nzp = std::arg(points[i][j][k].get_psi().inner_product(points[i][j][(k+1)%NZ].get_psi()));
				double nzn = std::arg(points[i][j][k].get_psi().inner_product(points[i][j][(NZ + k-1)%NZ].get_psi()));

				double sum = (nxp+nxn)/(LX*LX) + (nyp+nyn)/(LY*LY) + (nzp+nzn)/(LZ*LZ);
				points[i][j][k].set_div(sum);
			} 
		}
	}

	
	// FFT of div
	convert_div(in, points); //can use real dft for speedup
	fftw_execute(fp);
	modify_div(out,in,points);
	fftw_execute(bp);
	scale_ifft(out);
	reconvert_pressure(out,points);
}

void isf(grid & points) {

	fftw_complex * in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NX * NY * NZ);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NX * NY * NZ);
	fftw_plan fp = fftw_plan_dft_3d(NX,NY,NZ,in,out,FFTW_FORWARD,FFTW_MEASURE);
	fftw_plan bp = fftw_plan_dft_3d(NX,NY,NZ,in,out,FFTW_BACKWARD,FFTW_MEASURE);

	for ( unsigned int t = 0; t < NUM_TIME_STEPS ; t++) {
		min_density = INT_MAX;
		max_density = INT_MIN;
		schrodinger(in,out,fp,bp,points);
		normalize(points);
		pressure_project(in,out,fp,bp,points);
		update_velocity(points);
		advect_density(points);
	}
	fftw_destroy_plan(fp);
	fftw_destroy_plan(bp);
	fftw_free(in);
	fftw_free(out);

}

void initialize_filament(grid & points, Eigen::Vector3d center, Eigen::Vector3d normal, double r, double t) {
for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Eigen::Vector3d p(i*LX,j*LY,k*LZ);
				double d = (p-center).dot(normal.normalized());
				double theta = 0.0;
				/*if(d>-t/2.0 && d<t/2.0 && (p-center).squaredNorm() - pow(d,2) < pow(r,2)) {
					theta = M_PI * (1.0 + 2.0*d/t);
				}*/
				if ( (p-center).squaredNorm() - d*d < r*r ) {
					if ( d > 0 && d <= t/2.0) {
						theta = -M_PI * (2.0*d/t - 1.0);	
					} else if (d <= 0 && d>= -t/2.0){
						theta = -M_PI * (2.0*d/t + 1.0);
					}
				}
				Comp2 psi = points[i][j][k].get_psi();
				psi.set_z1(Complex(cos(theta),sin(theta)));
				psi.set_z2(Complex(EPSILON,0.0));
				points[i][j][k].set_psi(psi);
			}
		}
	}

	fftw_complex * in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NX * NY * NZ);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NX * NY * NZ);
	fftw_plan fp = fftw_plan_dft_3d(NX,NY,NZ,in,out,FFTW_FORWARD,FFTW_MEASURE);
	fftw_plan bp = fftw_plan_dft_3d(NX,NY,NZ,in,out,FFTW_BACKWARD,FFTW_MEASURE);
	
	
	normalize(points);
	pressure_project(in,out,fp,bp,points);
	fftw_destroy_plan(fp);
	fftw_destroy_plan(bp);
	fftw_free(in);
	fftw_free(out);
}

void initialize_density(grid & points) {
for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Eigen::Vector3d vorticity;
				vorticity(0) = (points[i][j][k].get_v()(2) - points[i][(j-1+NY)%NY][k].get_v()(2))/LY - (points[i][j][k].get_v()(1) - points[i][j][(k-1+NZ)%NZ].get_v()(1))/LZ;
				vorticity(1) = (points[i][j][k].get_v()(0) - points[i][j][(k-1+NZ)%NZ].get_v()(0))/LZ - (points[i][j][k].get_v()(2) - points[(i-1+NX)%NX][j][k].get_v()(2))/LX;
				vorticity(2) = (points[i][j][k].get_v()(1) - points[(i-1+NX)%NX][j][k].get_v()(1))/LX - (points[i][j][k].get_v()(0) - points[i][(j-1+NY)%NY][k].get_v()(0))/LY;
				double density = vorticity.norm();
				min_density = std::min(density, min_density);
				max_density = std::max(density, max_density);
				points[i][j][k].set_density(density);
			}
		}
	}
}

//add points1 filament to points filament
void add_filament(grid & points, grid &points1) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				Comp2 psi1 = points1[i][j][k].get_psi();
				psi.set_z1(psi.get_z1() * psi1.get_z1());
				psi.set_z2(psi.get_z2() * psi1.get_z2());
				points[i][j][k].set_psi(psi);
			}
		}
	}
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



void density_to_vol(grid & points) {
	int num_channels = 1;
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

	std::ofstream file;
	file.open("/render/density.vol", std::ios::out | std::ios::binary);
	file.write("VOL",3);
	file.write((char*)&file_format, 1);
	file.write((char*)&data_type, 4);
	file.write((char*)&nx, 4);
	file.write((char*)&ny, 4);
	file.write((char*)&nz, 4);
	file.write((char*)&num_channels, 4);
	file.write((char*)&xmin, 4);
	file.write((char*)&ymin, 4);
	file.write((char*)&zmin, 4);
	file.write((char*)&xmax, 4);
	file.write((char*)&ymax, 4);
	file.write((char*)&zmax, 4);

	for(int i =0 ;i<NX;i++)
	{
		for(int j = 0;j < NY; j++)
		{
			for(int k = 0; k < NZ; k++)
			{
				float tmp = points[i][j][k].get_density();
				tmp = tmp/(max_density-min_density);
				file.write((char*)&tmp,4);
			}
		}
	}

	file.close();
}


int main() {
	grid points;
	Eigen::Vector3d center(5.0,2.5,2.5);
	Eigen::Vector3d normal(-1.0,0.0,0.0);
	double r = 1.5;
	double t = 5*LX;
	resize_grid(points);
	initialize_filament(points,center,normal,r,t);
	update_velocity(points);
	initialize_density(points);
	isf(points);
	density_to_vol(points);
}
