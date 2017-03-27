#ifndef _Element1D_hpp
#define _Element1D_hpp

class Element1D
{
	public:
		Element1D(double a, double b);
		~Element1D() = default;
		double jacobian();
		double f(double x_hat);
	private:
		double a,b;
};

Element1D::Element1D(double a = 0.0, double b = 1.0) : a(a), b(b) {}

inline double Element1D::f(double x_hat)
{
	return (b - a)*x_hat + a;
}

inline double Element1D::jacobian()
{
	return b - a;
}

#endif