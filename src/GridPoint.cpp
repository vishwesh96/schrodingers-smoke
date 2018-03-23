#include <complex>
#include <GridPoint.h>

double GridPoint::get_density() const {
	return psi.norm2();
}

double GridPoint::get_pressure() const {
	return q;
}

void GridPoint::update_psi() {

}

void GridPoint::update_pressure() {

}

void GridPoint::update_div() {

}