#include <fftw3.h>
#include <vector>
#include <GridPoint.h>
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <math.h>

#define NX 128
#define NY 64
#define NZ 64
#define LX 0.078125
#define LY 0.078125
#define LZ 0.078125
#define NUM_TIME_STEPS 20
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

void print_psi(grid &points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				printf("%f\n",psi.get_z1().real());
			}
		}
	}
}

void print_rho(grid &points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				// printf("%d %d %d : %f +i%f\n",i,j,k,psi.get_z1().real(),psi.get_z1().imag());
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
				double theta = out[k + NZ *(j + NY * i)][0];
				Complex real_factor(exp(out[k + NZ *(j + NY * i)][1]),0.0);
				Complex factor(cos(theta),-sin(theta));
				psi.set_z1(psi.get_z1() * factor * real_factor);
				psi.set_z2(psi.get_z2() * factor * real_factor);
				points[i][j][k].set_psi(psi);
			}
		}
	}
}

void mult_exp(grid & points) {
for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				double lambda = -4.0*M_PI*M_PI*( pow( ((i-64.0)/(10.0)),2) + pow( ((j-32.0)/(5.0)),2) + pow(((k-32.0)/(5.0)),2) );
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
				points[i][j][k].update_div(sum);
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
	isf(points);
	print_psi(points);
}
