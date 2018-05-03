#include <complex>
#include <GridPoint.h>

Comp2 GridPoint::get_psi() const {
	return this->psi;
}

double GridPoint::get_density() const {
	return this->rho;
}

double GridPoint::get_pressure() const {
	return this->q;
}

double GridPoint::get_div() const {
	return this->div;
}

Eigen::Vector3d GridPoint::get_v() const {
	return this->v;
}

std::vector<double> GridPoint::get_color() const {
	return this->color;
}

void GridPoint::set_psi(Comp2 psi)  {
	this->psi = psi;
}

void GridPoint::set_density(double rho) {
	this->rho = rho;
}

void GridPoint::set_pressure(double q) {
	this->q = q;
}

void GridPoint::set_div(double div) {
	this->div = div;
}

void GridPoint::set_v(Eigen::Vector3d v) {
	this->v = v;
}

void GridPoint::set_color(std::vector<double> color) {
	this->color = color;
}

void GridPoint::normalize() {
	double norm = psi.norm();
	psi.set_z1(psi.get_z1()/norm);
	psi.set_z2(psi.get_z2()/norm);
}