#include <utils.h>

int NX = 128;
int NY = 64;
int NZ = 64;
double LX = 0.078125;
double LY = 0.078125;
double LZ = 0.078125;
double H = 0.1;
double DT = 1.0/9.0;

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

void update_velocity(grid & points, const std::vector<Volume> & vols = std::vector<Volume>()) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Comp2 psi = points[i][j][k].get_psi();
				Comp2 psi_x = points[(i+1)%NX][j][k].get_psi();
				Comp2 psi_y = points[i][(j+1)%NY][k].get_psi();
				Comp2 psi_z = points[i][j][(k+1)%NZ].get_psi();
				Eigen::Vector3d v = H*Eigen::Vector3d(std::arg(psi.inner_product(psi_x))/LX,std::arg(psi.inner_product(psi_y))/LY,std::arg(psi.inner_product(psi_z))/LZ);
				points[i][j][k].set_v(v);

				for(unsigned int l=0;l<vols.size();l++) {
					if(vols[l].volume[i][j][k]) {
						points[i][j][k].set_v(vols[l].velocity);
					}
				}
			}
		}
	}
}

double interpolate_density(grid & points,Eigen::Vector3d &y){
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

//Trilinear interpolation
std::vector<double> interpolate_color1(grid & points,Eigen::Vector3d &y){
	maintain_bound(y);
	std::vector<int> n = find_grid_cell(y);

	std::vector<double> c1 = points[n[0]][n[1]][n[2]].get_color();
	std::vector<double> c2 = points[(n[0]+1)%NX][n[1]][n[2]].get_color();
	std::vector<double> c3 = points[(n[0]+1)%NX][(n[1]+1)%NY][n[2]].get_color();
	std::vector<double> c4 = points[n[0]][(n[1]+1)%NY][n[2]].get_color();
	std::vector<double> c5 = points[n[0]][n[1]][(n[2]+1)%NZ].get_color();
	std::vector<double> c6 = points[(n[0]+1)%NX][n[1]][(n[2]+1)%NZ].get_color();
	std::vector<double> c7 = points[(n[0]+1)%NX][(n[1]+1)%NY][(n[2]+1)%NZ].get_color();
	std::vector<double> c8 = points[n[0]][(n[1]+1)%NY][(n[2]+1)%NZ].get_color();

	double z1 = n[0]*LX;
	double z2 = n[1]*LY;
	double z3 = n[2]*LZ;
	double z4 = (n[0]+1)*LX;
	double z5 = (n[1]+1)*LY;
	double z6 = (n[2]+1)*LZ;

	double y0 = y(0); 
	double y1 = y(1); 
	double y2 = y(2);

	std::vector<double> color(3);
	for(int i=0;i<3;i++){
		color[i] = ( (z6-y2) * (c1[i]*(z4-y0)*(z5-y1) + c2[i]*(y0-z1)*(z5-y1) + c3[i]*(y0-z1)*(y1-z2) + c4[i]*(z4-y0)*(y1-z2)) + (y2-z3)*(c5[i]*(z4-y0)*(z5-y1) + c6[i]*(y0-z1)*(z5-y1) + c7[i]*(y0-z1)*(y1-z2) + c8[i]*(z4-y0)*(y1-z2)) )/(LX*LY*LZ);
	}
	return color;
}

//Nearest neighbour interpolation
std::vector<double> interpolate_color2(grid & points,Eigen::Vector3d &y){
	maintain_bound(y);
	std::vector<int> n = find_grid_cell(y);

	double min_dis = INT_MAX;
	std::vector<double> c(3);
	c[0] = c[1] = c[2] = 0.0;
	for (int i = 0;i <2; i++){
		for(int j = 0;j < 2;j++) {
			for(int k = 0;k <2;k++) {
				Eigen::Vector3d corner_loc((n[0]+i)*LX, (n[1]+j)*LY, LZ*(n[2]+k));
				if ((corner_loc-y).norm()<min_dis && !is_zero_color(points[(n[0]+i)%NX][(n[1]+j)%NY][(n[2]+k)%NZ].get_color())) {
					min_dis = (corner_loc - y).norm();
					c = points[(n[0]+i)%NX][(n[1]+j)%NY][(n[2]+k)%NZ].get_color();
				}
			}
		}
	}
	return c;

}

void advect_density_and_color(grid & points) {
	grid old_points = points;
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Eigen::Vector3d x(i*LX,j*LY,k*LZ);
				Eigen::Vector3d v = points[i][j][k].get_v();
				Eigen::Vector3d y = x - v*DT;
				double density = interpolate_density(old_points,y);
				std::vector<double> color = interpolate_color2(old_points,y);
				max_density = std::max(density, max_density);
				points[i][j][k].set_density(density);
				points[i][j][k].set_color(color);
			}
		}
	}
}

void mult_exp(grid & points) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
			for(unsigned int j = 0; j<points[i].size(); j++) {
				for(unsigned int k = 0; k < points[i][j].size(); k++) {
					double lambda = -4.0*M_PI*M_PI*( pow( ((i-64.0)/(10.0)),2) + pow( ((j-32.0)/(5.0)),2) + pow(((k-32.0)/(5.0)),2) );
					//double lambda = -4.0*M_PI*M_PI*( pow( (i/(LX*NX)),2) + pow( (j/(LY*NY)),2) + pow((k/(LZ*NZ)),2) );
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

	fftshift(points);
	mult_exp(points);
	fftshift(points);

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

void pressure_project(grid& points) {
	fftw_complex * in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NX * NY * NZ);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NX * NY * NZ);
	fftw_plan fp = fftw_plan_dft_3d(NX,NY,NZ,in,out,FFTW_FORWARD,FFTW_MEASURE);
	fftw_plan bp = fftw_plan_dft_3d(NX,NY,NZ,in,out,FFTW_BACKWARD,FFTW_MEASURE);
	pressure_project(in,out,fp,bp,points);
	fftw_destroy_plan(fp);
	fftw_destroy_plan(bp);
	fftw_free(in);
	fftw_free(out);
}



void initialize_filament(grid & points, Eigen::Vector3d center, Eigen::Vector3d normal, double r, double t, std::vector<double> color) {
for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Eigen::Vector3d p(i*LX,j*LY,k*LZ);
				double d = (p-center).dot(normal.normalized());
				double theta = 0.0;
				if(d>-t/2.0 && d<t/2.0 && (p-center).squaredNorm() - pow(d,2) < pow(r,2)) {
					theta = M_PI * (1.0 + 2.0*d/t);
					points[i][j][k].set_color(color);
				}
				// if ( (p-center).squaredNorm() - d*d < r*r ) {
				// 	if ( d > 0 && d <= t/2.0) {
				// 		theta = -M_PI * (2.0*d/t - 1.0);	
				// 	} else if (d <= 0 && d>= -t/2.0){
				// 		theta = -M_PI * (2.0*d/t + 1.0);
				// 	}
				// 	points[i][j][k].set_color(color);
				// }
				points[i][j][k].set_psi(Comp2(Complex(cos(theta),sin(theta)),Complex(EPSILON,0.0)));
			}
		}
	}
	normalize(points);
	pressure_project(points);
}

void add_filament(grid& points, Eigen::Vector3d center, Eigen::Vector3d normal, double r, double t, std::vector<double> color) {
for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				Eigen::Vector3d p(i*LX,j*LY,k*LZ);
				double d = (p-center).dot(normal.normalized());
				double theta = 0.0;
				/*if(d>-t/2.0 && d<t/2.0 && (p-center).squaredNorm() - pow(d,2) < pow(r,2)) {
					theta = M_PI * (1.0 + 2.0*d/t);
					points[i][j][k].set_color(color);
				}*/
				if ( (p-center).squaredNorm() - d*d < r*r ) {
					if ( d > 0 && d <= t/2.0) {
						theta = -M_PI * (2.0*d/t - 1.0);	
					} else if (d <= 0 && d>= -t/2.0){
						theta = -M_PI * (2.0*d/t + 1.0);
					}
					points[i][j][k].set_color(color);
				}
				Comp2 psi = points[i][j][k].get_psi();
				points[i][j][k].set_psi(Comp2(psi.get_z1()*Complex(cos(theta),sin(theta)),psi.get_z2()*Complex(EPSILON,0.0)));
			}
		}
	}
	normalize(points);
	pressure_project(points);
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
					max_density = std::max(density, max_density);
					points[i][j][k].set_density(density);
				}
			}
		}
}

Volume cylinder_volume(Eigen::Vector3d center, double radius, double thickness, Eigen::Vector3d velocity, std::vector<double> color, double density) {
	Volume vol;
	resize_volume(vol);
	for (unsigned int i = 0 ; i < vol.volume.size() ; i++ ) {
		for(unsigned int j = 0; j<vol.volume[i].size(); j++) {
			for(unsigned int k = 0; k < vol.volume[i][j].size(); k++) {
				Eigen::Vector3d pos(i*LX,j*LY,k*LZ);
				Eigen::Vector2d pos1(pos(1),pos(2));
				Eigen::Vector2d center1(center(1),center(2));
				vol.volume[i][j][k] = (((pos1-center1).norm()<radius) && (pos(0) - center(0))*(pos(0) - center(0)) < thickness*thickness/4.0);
			}
		}
	}
	vol.velocity = velocity;
	vol.color = color;
	vol.density = density;
	return vol;
}

Volume cylinder_shell_volume(Eigen::Vector3d center, double outer_radius, double inner_radius, double thickness, Eigen::Vector3d velocity, std::vector<double> color, double density) {
	Volume vol;
	resize_volume(vol);
	for (unsigned int i = 0 ; i < vol.volume.size() ; i++ ) {
		for(unsigned int j = 0; j<vol.volume[i].size(); j++) {
			for(unsigned int k = 0; k < vol.volume[i][j].size(); k++) {
				Eigen::Vector3d pos(i*LX,j*LY,k*LZ);
				Eigen::Vector2d pos1(pos(1),pos(2));
				Eigen::Vector2d center1(center(1),center(2));
				vol.volume[i][j][k] = ((pos1-center1).norm()<outer_radius) && ((pos1-center1).norm()>inner_radius) && (pos(0) - center(0))*(pos(0) - center(0)) < thickness*thickness/4.0;
			}
		}
	}
	vol.velocity = velocity;
	vol.color = color;
	vol.density = density;
	return vol;
}

Volume sphere_volume(Eigen::Vector3d center, double radius, Eigen::Vector3d velocity, std::vector<double> color, double density) {
	Volume vol;
	resize_volume(vol);
	for (unsigned int i = 0 ; i < vol.volume.size() ; i++ ) {
		for(unsigned int j = 0; j<vol.volume[i].size(); j++) {
			for(unsigned int k = 0; k < vol.volume[i][j].size(); k++) {
				Eigen::Vector3d pos(i*LX,j*LY,k*LZ);
				vol.volume[i][j][k] = (pos-center).norm()<radius;
			}
		}
	}
	vol.velocity = velocity;
	vol.color = color;
	vol.density = density;
	return vol;
}

void initialize_density_volume(grid &points,const std::vector<Volume> & vols) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
		for(unsigned int j = 0; j<points[i].size(); j++) {
			for(unsigned int k = 0; k < points[i][j].size(); k++) {
				for(unsigned int l=0; l <vols.size(); l++) {
					max_density = std::max(max_density,vols[l].density);
					if(vols[l].volume[i][j][k]){
						points[i][j][k].set_density(vols[l].density);
						points[i][j][k].set_color(vols[l].color);
					}
				}
			}
		}
	}
}	

void constraint_projection(grid &points, const std::vector<Volume> &vols, double t) {
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
			for(unsigned int j = 0; j<points[i].size(); j++) {
				for(unsigned int k = 0; k < points[i][j].size(); k++) {
					for(unsigned int l = 0; l<vols.size(); l++) {
						if(vols[l].volume[i][j][k]) {
							Eigen::Vector3d kw = vols[l].velocity/H;
							Eigen::Vector3d pos(i*LX,j*LY,k*LZ);
							double theta = kw.dot(pos) - vols[l].velocity.squaredNorm()*t/(2.0*H);
							Complex phi(cos(theta),sin(theta));
							Comp2 psi = points[i][j][k].get_psi();
							points[i][j][k].set_psi(Comp2(std::abs(psi.get_z1())*phi,std::abs(psi.get_z2())*phi));
						}
					}
				}
			}
		}
	pressure_project(points);
}

void initialize_psi_velocity(grid & points, std::vector<Volume> & vols) {
	normalize(points);
	for(int i=0;i<NUM_VELOCITY_ITER;i++){
			constraint_projection(points,vols,0.0);	
	}
}

void constraint_projection(grid &points, Volume &vol, double t) {
	Eigen::Vector3d kw = vol.velocity/H;
	for (unsigned int i = 0 ; i < points.size() ; i++ ) {
			for(unsigned int j = 0; j<points[i].size(); j++) {
				for(unsigned int k = 0; k < points[i][j].size(); k++) {
					if(vol.volume[i][j][k]) {
						Eigen::Vector3d pos(i*LX,j*LY,k*LZ);
						double theta = kw.dot(pos) - vol.velocity.squaredNorm()*t/(2.0*H);
						Complex phi(cos(theta),sin(theta));
						Comp2 psi = points[i][j][k].get_psi();
						points[i][j][k].set_psi(Comp2(std::abs(psi.get_z1())*phi,std::abs(psi.get_z2())*phi));
					}
				}
			}
		}
	pressure_project(points);
}

void isf(grid & points, std::string & filename, bool constraint = false, const std::vector<Volume> &vols = std::vector<Volume>()) {

	fftw_complex * in = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NX * NY * NZ);
	fftw_complex * out = (fftw_complex *) fftw_malloc(sizeof(fftw_complex) * NX * NY * NZ);
	fftw_plan fp = fftw_plan_dft_3d(NX,NY,NZ,in,out,FFTW_FORWARD,FFTW_MEASURE);
	fftw_plan bp = fftw_plan_dft_3d(NX,NY,NZ,in,out,FFTW_BACKWARD,FFTW_MEASURE);

	for ( unsigned int t = 0; t < NUM_TIME_STEPS ; t++) {
		fprintf(stderr,"Iteration Number : %d\n",t);
		density_to_vol(points, t, filename);
		max_density = INT_MIN;
		schrodinger(in,out,fp,bp,points);
		normalize(points);
		pressure_project(in,out,fp,bp,points);
		if(constraint)
			constraint_projection(points,vols,(t+1)*DT);
		update_velocity(points,vols);
		advect_density_and_color(points);
		if(constraint)
			initialize_density_volume(points,vols);
	}
	fftw_destroy_plan(fp);
	fftw_destroy_plan(bp);
	fftw_free(in);
	fftw_free(out);

}

void render_jet_obstacle(std::string & filename, bool is_obstacle = true) {
	
	H = 0.1;
	DT = 1.0/7.0;

	grid points;

	resize_grid(points);
	std::vector<Volume> vols;
	std::vector<double> color1(3);
	std::vector<double> color2(3);
	double density1,density2;
	double thickness1;
	double radius1,radius2,radius3;
	
	Eigen::Vector3d velocity1(-50.0,0.0,0.0);
	Eigen::Vector3d velocity2(0.0,0.0,0.0);

	Eigen::Vector3d center1(7.1,2.5,2.5);
	Eigen::Vector3d center2(6.0,2.5,2.5);

	radius1 = 0.60;
	radius2 = 0.3;
	radius3 = 0.5;	
	thickness1 = 1.0;	

	color1[0] = 0.3;
	color1[1] = 0.5;
	color1[2] = 0.7;

	color2[0] = 0.0;
	color2[1] = 0.0;
	color2[2] = 0.0;

	density1 = 0.1;
	density2 = 0.0;

	vols.push_back(cylinder_shell_volume(center1,radius1,radius2,thickness1,velocity1,color1,density1));
	if (is_obstacle)
		vols.push_back(sphere_volume(center2,radius3,velocity2,color2,density2));
	initialize_density_volume(points,vols);
	initialize_psi_velocity(points,vols);
	isf(points,filename,true,vols);
}

void render_sphere_collision(std::string & filename) {

	// H = 0.1 Dt=1/6

	H = 0.1;
	DT = 1.0/15.0;

	grid points;
	Eigen::Vector3d center1(5.4,2.5,2.5);
	Eigen::Vector3d center2(2.8,2.5,2.5);
	double radius1 = 1.0, radius2 = 1.0;	
	Eigen::Vector3d velocity1(-20.0,0.0,0.0);
	Eigen::Vector3d velocity2(20.0,0.0,0.0);

	resize_grid(points);
	std::vector<Volume> vols;
	std::vector<double> color1(3);
	std::vector<double> color2(3);
	double density1,density2;

	color1[0] = 0.99;
	color1[1] = 0.3;
	color1[2] = 0.3;

	color2[0] = 0.3;
	color2[1] = 0.3;
	color2[2] = 0.99;

	density1 = 0.99;
	density2 = 0.99;

	vols.push_back(sphere_volume(center1,radius1,velocity1,color1,density1));
	vols.push_back(sphere_volume(center2,radius2,velocity2,color2,density2));
	
	initialize_density_volume(points,vols);
	initialize_psi_velocity(points,vols);
	
	isf(points,filename);
}

void render_filament(std::string & filename) {
	H = 0.02;
	DT = 1.0/12.0;
	grid points;
	Eigen::Vector3d center1(5.0,2.5,2.5);
	Eigen::Vector3d normal1(-1.0,0.0,0.0);
	double r1 = 1.2;
	double t1 = 10*LX;

	std::vector<double> color1(3);
	color1[0] = 0.7;
	color1[1] = 0.3;
	color1[2] = 0.3;

	resize_grid(points);
	initialize_filament(points,center1,normal1,r1,t1,color1);
	update_velocity(points);
	initialize_density(points);
	isf(points, filename);
}

int main(int argc, char* argv[]) {
	std::string filename = argv[1];
	NUM_TIME_STEPS = std::atoi(argv[2]);

	if (filename == "sphere_collision") {
		render_sphere_collision(filename);
	} else if(filename == "jet_obstacle"){
		render_jet_obstacle(filename);
	} else if(filename == "jet") {
		render_jet_obstacle(filename, false);
	} else if(filename == "filament") {
		render_filament(filename);
	} else {
		printf("Wrong file name");
	}
}
