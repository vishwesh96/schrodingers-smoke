#include <complex>
#include <Comp2.h>

#define EPSILON 0.01

class GridPoint {

	private:
		Comp2 psi; 
		double q;	// pressure
		double div;
		double rho;
		Eigen::Vector3d v;
		std::vector<double> color;
	public:
		GridPoint() {
			this->rho = 0.0;
			this->psi = Comp2(Complex(1.0,0.0),Complex(EPSILON,0.0));
			std::vector<double> c;
			c.push_back(0.0);
			c.push_back(0.0);
			c.push_back(0.0);
			this->color = c;
		};

		Comp2 get_psi() const;

		double get_density() const;

		double get_pressure() const;

		double get_div() const;

		Eigen::Vector3d get_v() const;

		std::vector<double> get_color() const;
		
		void set_psi(Comp2) ;

		void set_density(double rho);

		void set_pressure(double q);

		void set_div(double div);

		void set_v(Eigen::Vector3d v);

		void set_color(std::vector<double> color);

		void normalize();	

};