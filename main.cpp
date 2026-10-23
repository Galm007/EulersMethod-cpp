#include <cmath>
#include <fstream>
#include <functional>
#include <matplot/matplot.h>

using namespace matplot;

float func1(float t, float x) {
	return t + x;
}

float func2(float t, float x) {
	return (1.0 / 3.0) * x * (8.0 - x);
}

float ode_ivp_1_exact_soln(float t) {
	return -t - 1.0 + 2.0 * std::exp(t);
}

float ode_ivp_2_exact_soln(float t) {
	return 8.0 / (1.0 + 7.0 * std::exp(t * -8.0 / 3.0));
}

std::vector<double> EulersMethod(
	double h,
	double from_t,
	double to_t,
	double x0,
	std::function<double(double, double)> dxdt
) {
	std::vector<double> x_values = {};
	x_values.reserve(1 + (to_t - from_t) / h);
	x_values.push_back(x0);

	int n = std::lround((to_t - from_t) / h);
	double x = x0;
	for (int i = 0; i < n; ++i) {
		double t = from_t + i * h;
		x += h * dxdt(t, x);
		x_values.push_back(x);
	}

	return x_values;
}

int main() {
	// Note: if the number of subintervals is n, then the number of
	// elements in the list must be n + 1.
	// That is why I add +1 to all the subinterval amounts.

	// ODE-IVP 1
	float from_t = 0.0, to_t = 1.0, x0 = 1.0;

	std::vector<double> t_space_1_h01      = linspace(from_t, to_t, 1 + (to_t - from_t) / 0.1);
	std::vector<double> t_space_1_h001     = linspace(from_t, to_t, 1 + (to_t - from_t) / 0.01);
	std::vector<double> t_space_1_h0001    = linspace(from_t, to_t, 1 + (to_t - from_t) / 0.001);
	std::vector<double> func1_values_h01   = EulersMethod(0.1,   from_t, to_t, x0, func1);
	std::vector<double> func1_values_h001  = EulersMethod(0.01,  from_t, to_t, x0, func1);
	std::vector<double> func1_values_h0001 = EulersMethod(0.001, from_t, to_t, x0, func1);

	// Solution
	std::vector<double> t_space_1_soln = linspace(from_t, to_t, 100);
	std::vector<double> func1_values_soln = transform(t_space_1_soln, ode_ivp_1_exact_soln);

	// ODE-IVP 2
	from_t = 0.0, to_t = 5.0, x0 = 1.0;

	std::vector<double> t_space_2_n5     = linspace(from_t, to_t, 1 + 5);
	std::vector<double> t_space_2_n10    = linspace(from_t, to_t, 1 + 10);
	std::vector<double> t_space_2_n20    = linspace(from_t, to_t, 1 + 20);
	std::vector<double> func2_values_n5  = EulersMethod((to_t - from_t) / 5.0,  from_t, to_t, x0, func2);
	std::vector<double> func2_values_n10 = EulersMethod((to_t - from_t) / 10.0, from_t, to_t, x0, func2);
	std::vector<double> func2_values_n20 = EulersMethod((to_t - from_t) / 20.0, from_t, to_t, x0, func2);

	// Solution
	std::vector<double> t_space_2_soln = linspace(from_t, to_t, 100);
	std::vector<double> func2_values_soln = transform(t_space_2_soln, ode_ivp_2_exact_soln);

	// ----- Plot graphs -----

	// ODE-IVP 1
	axes_handle window1 = axes({0.15f, 0.60f, 0.8f, 0.30f});
	title(window1, "dy/dx = t + x");
	xlabel(window1, "t");
	ylabel(window1, "x");
	legend(window1, {"h = 0.1", "h = 0.01", "h = 0.001", "exact"})
		->location(legend::general_alignment::topleft);
	plot(
		window1,
		t_space_1_h01,
		func1_values_h01,
		t_space_1_h001,
		func1_values_h001,
		t_space_1_h0001,
		func1_values_h0001,
		t_space_1_soln,
		func1_values_soln
	);

	// ODE-IVP 2
	axes_handle window2 = axes({0.15f, 0.10f, 0.8f, 0.30f});
	title(window2, "dy/dx = (1/3) * x * (8 - x)");
	xlabel(window2, "t");
	ylabel(window2, "x");
	legend(window2, {"n = 5", "n = 10", "n = 20", "exact"})
		->location(legend::general_alignment::topleft);
	plot(
		window2,
      		t_space_2_n5,
      		func2_values_n5,
      		t_space_2_n10,
      		func2_values_n10,
      		t_space_2_n20,
      		func2_values_n20,
		t_space_2_soln,
		func2_values_soln
      	);

	save("plots.svg");

	// Create table of errors for ODE-IVP 1 

	std::ofstream table_file;
	table_file.open("ode-ivp-1-table.txt");
	if (table_file.is_open()) {
		double h01 = func1_values_h01.back();
		double h001 = func1_values_h001.back();
		double h0001 = func1_values_h0001.back();
		double exact = func1_values_soln.back();

		table_file << "--- Approximations of ODE-IVP 1 solution at time t=1 ---\n";
		table_file << "h = 0.1   | value: " << h01   << " | " << "error: " << exact - h01   << "\n";
		table_file << "h = 0.01  | value: " << h001  << " | " << "error: " << exact - h001  << "\n";
		table_file << "h = 0.001 | value: " << h0001 << " | " << "error: " << exact - h0001 << "\n";
		table_file << "Exact     | value: " << exact << "\n";
		table_file.close();
	} else {
		std::cout << "Unable to open file ode-ivp-1-table.txt!" << std::endl;
	}
	return 0;
}
