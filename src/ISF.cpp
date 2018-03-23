#include <fftw3.h>
#include <Comp2.h>

#define NX 128
#define NY 128
#define NZ 128
#define LX 1.0
#define LY 1.0
#define LZ 1.0
#define NUM_TIME_STEPS 10
#define DT 0.1
#define H 0.1
#define PI 3.14

typedef vector< vector < vector< GridPoint> > > grid;

void isf(grid & points) {

	fftw_complex * in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NX * NY * NZ);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NX * NY * NZ);
	fftw_plan fp = fftw_plan_dft_3d(NX,NY,NZ,in,out,FFTW_FORWARD,FFTW_MEASURE);
	fftw_plan bp = fftw_plan_dft_3d(NX,NY,NZ,in,out,FFTW_BACKWARD,FFTW_MEASURE);

	for ( int t = 0; t < NUM_TIME_STEPS ; t++) {
		schrodinger(in,out,fp,bp,points);
		normalize(points);
		pressure_project(points);
	}
	fftw_destroy_plan(fp);
	fftw_destroy_plan(bp);
	fftw_free(in);
	fftw_free(out);
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
	reconvert_psi1(out,points);
	convert_psi2(in,points);
	fftw_execute(bp);
	reconvert_psi2(out,points);
}

void pressure_project(fftw_complex * in, fftw_complex * out, fftw_plan fp, fftw_plan bp, grid & points) {
	
	// Set up divergence for each grid point
	for (int i = 0 ; i < points.size() ; i++ ) {
		for(int j = 0; j<points[i].size(); j++) {
			for(int k = 0; k < points[i][j].size(); k++) {
				double V =  8 * LX * LY * LZ;
				double nxp = std::arg(points[i][j][k].get_psi().inner_product(points[(i+1)%NX][j][k])); 
				double nxn = std::arg(points[i][j][k].get_psi().inner_product(points[(NX + i-1)%NX][j][k]));
				double nyp = std::arg(points[i][j][k].get_psi().inner_product(points[i][(j+1)%NY][k]));
				double nyn = std::arg(points[i][j][k].get_psi().inner_product(points[i][(NY+j-1)%NY][k]));
				double nzp = std::arg(points[i][j][k].get_psi().inner_product(points[i][j][(k+1)%NK]));
				double nzn = std::arg(points[i][j][k].get_psi().inner_product(points[i][j][(NK + k-1)%NK]));

				double sum = (nxp*4.0*LY*LZ)/LX + (nxn*4.0*LY*LZ)/LX + (nyp*4.0*LZ*LX)/LY + (nyn*4.0*LZ*LX)/LY + (nzp*4.0*LX*LY)/LZ + (nzn*4.0*LX*LY)/LZ;

				points[i][j][k].update_div(sum/V);
			} 
		}
	}

	// FFT of div
	convert_div(in, points);



	
}

void convert_div(fftw_complex *in, grid & points) {
	for (int i = 0 ; i < points.size() ; i++ ) {
		for(int j = 0; j<points[i].size(); j++) {
			for(int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				in[k + NZ *(j + NY * i)][0] = psi.get_z1().real();
			}
		}
	}
}

void convert_psi1(fftw_complex * in, grid points) {
	for (int i = 0 ; i < points.size() ; i++ ) {
		for(int j = 0; j<points[i].size(); j++) {
			for(int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				in[k + NZ *(j + NY * i)][0] = psi.get_z1().real();
				in[k + NZ *(j + NY * i)][1] = psi.get_z1().img();
			}
		}
	}
}

void convert_psi2(fftw_complex * in, grid points) {
	for (int i = 0 ; i < points.size() ; i++ ) {
		for(int j = 0; j<points[i].size(); j++) {
			for(int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				in[k + NZ *(j + NY * i)][0] = psi.get_z2().real();
				in[k + NZ *(j + NY * i)][1] = psi.get_z2().img();
			}
		}
	}
}

void reconvert_psi1(fftw_complex * out, grid & points) {
	for (int i = 0 ; i < points.size() ; i++ ) {
		for(int j = 0; j<points[i].size(); j++) {
			for(int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				psi.set_z1(Complex(out[k + NZ *(j + NY * i)][0], out[k + NZ *(j + NY * i)][1]));
			}
		}
	}
}

void reconvert_psi2(fftw_complex * out, grid & points) {
	for (int i = 0 ; i < points.size() ; i++ ) {
		for(int j = 0; j<points[i].size(); j++) {
			for(int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				psi.set_z2(Complex(out[k + NZ *(j + NY * i)][0], out[k + NZ *(j + NY * i)][1]));
			}
		}
	}
}

void mult_exp(grid & points) {
for (int i = 0 ; i < points.size() ; i++ ) {
		for(int j = 0; j<points[i].size(); j++) {
			for(int k = 0; k < points[i][j].size(); k++) {
				double lambda = -pow((2*PI),2)*( pow((i/LX*NX),2) + pow((j/LY*NY),2) + pow((k/LZ*NZ),2));
				double theta = 0.5 * lambda * DT * H;
				Complex factor(cos(theta),sin(theta));
				points[i][j][k].set_psi(factor * points[i][j][k].get_psi());
			}
		}
	}
}

void normatlize(grid & points) {
	for (int i = 0 ; i < points.size() ; i++ ) {
		for(int j = 0; j<points[i].size(); j++) {
			for(int k = 0; k < points[i][j].size(); k++) {
				points[i][j][k].normatlize();
			}
		}
	}
}



int main() {
	grid points;



}