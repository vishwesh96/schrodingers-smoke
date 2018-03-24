#include <complex>
#include <Comp2.h>

class GridPoint {

	private:
		Comp2 psi; 
		double q;	// pressure
		double div;

	public:
		GridPoint() {};
		double get_density() const;

		double get_pressure() const;

		double get_div() const;

		void normalize();	

		void update_psi();

		void update_pressure();

		void update_div(double div);

		Comp2 get_psi() const;
		void set_psi(Comp2) ;

};