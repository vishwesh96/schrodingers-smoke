#include <complex>
#include <GridPoint.h>

double GridPoint::get_density() const {
	return psi.norm2();
}

double GridPoint::get_pressure() const {
	return q;
}

double GridPoint::get_div() const {
	return div;
}

Comp2 GridPoint::get_psi() const {
	return psi;
}

void GridPoint::set_psi(Comp2 psi)  {
	this->psi = psi;
}

void GridPoint::update_psi() {

}

void GridPoint::update_pressure() {

}

void GridPoint::update_div(double div) {
	this->div = div;
}

void GridPoint::normalize() {
	double norm = psi.norm();
	psi.set_z1(psi.get_z1()/norm);
	psi.set_z2(psi.get_z2()/norm);
}