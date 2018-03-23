#include <complex>
#include <math.h>

typedef std::complex<double> Complex;

class Comp2 {
	private:
		 Complex z1,z2;

	public:
		
		Comp2(Complex z1, Complex z2): z1(z1), z2(z2) {};
		Comp2() {};

		Complex get_z1() const {
			return z1;
		}

		Complex get_z2() const {
			return z2;
		}

		void set_z1(Complex z) {
			z1 = z;
		}

		void set_z2(Complex z) {
			z2 = z;
		}

		Complex inner_product(Comp2 z) const {
			return z1*std::conj((z.get_z1())) + z2*std::conj(z.get_z2()) ;
		}

		double norm() const {
			return sqrt( (z1 * std::conj(z1) + z2 * std::conj(z2)).real() );
		}

		double norm2() const {
			return ( (z1 * std::conj(z1) + z2 * std::conj(z2)).real() );	
		}


};