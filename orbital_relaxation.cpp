

#include <iostream>
#include <cmath>
#include <string>
#include <tgmath.h>
#include <fstream>
#include <vector>
#include <complex>
#include <stdio.h>
#include <gsl/gsl_sf_dawson.h>

using namespace std;


double pi = 3.14159;
double hbar = 6.626 * pow(10, -34) / (2 * pi);
double e = 1.602 * pow(10.0, -19);
double m = 9.11 * pow(10, -31);
double m_perp = .198 * m;
double rho = 2330; 
double v_t = 5420; 
double v_l = 9330; 
double T = 1.0;

double Rashba = 45;
double Dressel = 0;

double xi_d = 5.0 * e; //dilation deformation 5 eV
double xi_u = 8.77 * e; //uniaxial sheer deformation 8.77 eV



double size(double E) {
     return sqrt(2) * hbar / sqrt(E * m_perp);
}

double erfi(double x) {
    return (2.0 / sqrt(pi)) * exp(x * x) * gsl_sf_dawson(x);
}

double A_int(int n, double delta, double v) {
    double q = delta / (hbar * v);
    double L = size(delta);
    double z0 = .25 * L;
    double long a = .25 * (L * L - z0*z0) * q * q;
    if (n == 0) {
        return 0.5 * sqrt(pi) * erfi(sqrt(a)) / sqrt(a);
    }
    if (n == 2) {
        return 0.5 * (-1.0 * sqrt(pi) * erfi(sqrt(a)) / (2 * pow(a, 3.0 / 2.0)) + exp(a) / a);
    }
    if (n == 4) {
        return 0.5 * (3.0 * sqrt(pi) * erfi(sqrt(a)) / (4 * pow(a, 5.0 / 2.0)) + exp(a) * (2 * a - 3.0) / (2 * a * a));
    }
    if (n == 6) {
        return 0.5 * ((4.0 * a * a - 10.0 * a + 15.0) * exp(a) / (4 * pow(a, 3.0)) - (15.0 * sqrt(pi) * erfi(sqrt(a))) / (8 * pow(a, 7.0 / 2.0)));
    }
    else return 0;
}

double long_comp(double E) {
    double delta = E;
    double L = size(E);
    double q = delta / (hbar * v_l);
    double prefactor = exp(- L * L * q *q / 4.0)* pow(delta, 4.0) / (pow(hbar, 4.0) * pow(v_l, 7.0));
    double A_term = (xi_d * xi_d * (A_int(0, delta, v_l) - A_int(2, delta, v_l))
        + 2.0 * xi_u * xi_d * (A_int(2, delta, v_l) - A_int(4, delta, v_l))
        + xi_u * xi_u * (A_int(4, delta, v_l) - A_int(6, delta, v_l)));
    return prefactor * A_term / tanh(delta / (2 * 1.38 * pow(10, -23) * T));

}

double trans_comp(double E) {
    double delta = E;
    double L = size(E);
    double q = delta / (hbar * v_t);
    double prefactor = exp(-L * L * q * q / 4.0) * pow(delta, 4.0) / (pow(hbar, 4.0) * pow(v_t, 7.0));
    double A_term = xi_u * xi_u * (A_int(2, delta, v_t) - 2.0 * A_int(4, delta, v_t) + A_int(6, delta, v_t));
    return prefactor * A_term / tanh( delta / (2 * 1.38 * pow(10, -23) * T));

}

double approx(double E) {
    return 2.0 * xi_u * xi_u * pow(E, 4.0) / (105 * pow(v_t, 7.0) * pow(hbar, 4.0) * pi * rho * m_perp
        * tanh(E / (2 * 1.38 * pow(10, -23) * T)));
}

int main() {
    double data_points = 200;
    double beginning = 0.05;
    double ending = 10.0;
    
    double prefactor = 1.0 / (16.0 * pi * rho * m_perp);
    ofstream myfile;
    myfile.open("relaxation.txt");
    for (double i = beginning * data_points; i <= ending * data_points; i++) {
        double E_mev = i / data_points;
        double E = E_mev * e / 1000.0;
        myfile << E_mev << " ";
        myfile << 1.0 / (prefactor*(long_comp(E) + trans_comp(E))) << " " << 1.0 / approx(E) << endl;
    }
    myfile.close();
    return 0;
}

