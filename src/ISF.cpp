#include <fftw3.h>
#include <GridPoint.h>

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
				Comp2 psi = points[i][j][k].get_psi();
				printf("%d %d %d : %f\n",i,j,k,points[i][j][k].get_div());
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
	double z4 = ((n[0]+1)%NX)*LX;
	double z5 = ((n[1]+1)%NY)*LY;
	double z6 = ((n[2]+1)%NZ)*LZ;

	double y0 = y(0); 
	double y1 = y(1); 
	double y2 = y(2);

	double density = ( (z6-y2) * (d1*(z4-y0)*(z5-y1) + d2*(z1-y0)*(z5-y1) + d3*(z1-y0)*(z2-y1) + d4*(z4-y0)*(z2-y1)) + (z3-y2)*(d5*(z4-y0)*(z5-y1) + d6*(z1-y0)*(z5-y1) + d7*(z1-y0)*(z2-y1) + d8*(z4-y0)*(z2-y1)) )/(LX*LY*LZ);
	return density;
}

void advect_density(grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Eigen::Vector3d x(i*LX,j*LY,k*LZ);
				Eigen::Vector3d v = points[i][j][k].get_v();
				Eigen::Vector3d y = x - v*DT;
				points[i][j][k].set_density(interpolate_density(points,y));
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

int main() {
	grid points;
	Eigen::Vector3d center(5.0,2.5,2.5);
	Eigen::Vector3d normal(-1.0,0.0,0.0);
	double r = 1.5;
	double t = 5*LX;
	resize_grid(points);
	initialize_filament(points,center,normal,r,t);
	initialize_density(points);
	isf(points);
}
