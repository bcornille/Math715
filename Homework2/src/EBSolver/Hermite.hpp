#include <cmath>
#include <array>

#ifndef _Hermite_hpp
#define _Hermite_hpp

class Hermite
{
	public:
		Hermite();
		~Hermite() = default;
		std::array<double, 4> eval(double x);
		std::array<double, 4> eval_d(double x);
		std::array<double, 4> eval_dd(double x);
	private:
		double h00(double x);
		double h10(double x);
		double h01(double x);
		double h11(double x);
		double h00_d(double x);
		double h10_d(double x);
		double h01_d(double x);
		double h11_d(double x);
		double h00_dd(double x);
		double h10_dd(double x);
		double h01_dd(double x);
		double h11_dd(double x);
};

inline double Hermite::h00(double x)
{
	return 2.0*pow(x, 3) - 3.0*pow(x, 2) + 1;
}

inline double Hermite::h10(double x)
{
	return pow(x, 3) - 2.0*pow(x, 2) + x;
}

inline double Hermite::h01(double x)
{
	return 3.0*pow(x, 2) - 2.0*pow(x, 3);
}

inline double Hermite::h11(double x)
{
	return pow(x, 3) - pow(x, 2);
}

inline double Hermite::h00_d(double x)
{
	return 6.0*pow(x, 2) - 6.0*x;
}

inline double Hermite::h10_d(double x)
{
	return 3.0*pow(x, 2) - 4.0*x + 1;
}

inline double Hermite::h01_d(double x)
{
	return 6.0*x - 6.0*pow(x, 2);
}

inline double Hermite::h11_d(double x)
{
	return 3.0*pow(x, 2) - 2.0*x;
}

inline double Hermite::h00_dd(double x)
{
	return 12.0*x - 6.0;
}

inline double Hermite::h10_dd(double x)
{
	return 6.0*x - 4.0;
}

inline double Hermite::h01_dd(double x)
{
	return 6.0 - 12.0*x;
}

inline double Hermite::h11_dd(double x)
{
	return 6.0*x - 2.0;
}

inline std::array<double, 4> Hermite::eval(double x)
{
	return {h00(x), h10(x), h01(x), h11(x)};
}

inline std::array<double, 4> Hermite::eval_d(double x)
{
	return {h00_d(x), h10_d(x), h01_d(x), h11_d(x)};
}

inline std::array<double, 4> Hermite::eval_dd(double x)
{
	return {h00_dd(x), h10_dd(x), h01_dd(x), h11_dd(x)};
}

#endif // _Hermite_hpp
