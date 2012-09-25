#include "MetricMatrixCube2.h"
#include <dolfin.h>
#include <iostream>

namespace magnumfe {

  // returns the biggest of three double values
  double MetricMatrixCube2::max3(double a, double b, double c)
  {
    double result = a;
    if (result < b) result = b;
    if (result < c) result = c;
    return result;
  }

  // returns the biggest of three double values
  double MetricMatrixCube2::min3(double a, double b, double c)
  {
    double result = a;
    if (result > b) result = b;
    if (result > c) result = c;
    return result;
  }

  // alias Power for mathematica C export
  double MetricMatrixCube2::Power(double x, int y)
  {
    if (y < 0) return 1.0 / Power(x, -y);

    double result = 1.0;
    for (int i=0; i<y; ++i) {
      result *= x;
    }
    return result;
    //return std::pow(x, y);
  }

  void MetricMatrixCube2::eval(dolfin::Array<double>& values, const dolfin::Array<double>& pos) const
  {
    // constats
    //const double D = min3(size[0], size[1], size[2]) * 0.5;
    const double D = min3(size[0], size[1], size[2]);

    // calculate indices
    const int ix = (1+coord)%3;
    const int iy = (2+coord)%3;
    const int iz = (0+coord)%3;

    // point mirror position if z is negative
    const double x = (pos[iz] > 0.0) ? pos[ix] : -pos[ix];
    const double y = (pos[iz] > 0.0) ? pos[iy] : -pos[iy];
    const double z = (pos[iz] > 0.0) ? pos[iz] : -pos[iz];

    const double Lx = size[ix];
    const double Ly = size[iy];
    const double Lz = size[iz];

    double &xx = values[ix + 3*ix];
    double &xy = values[ix + 3*iy];
    double &xz = values[ix + 3*iz];
    double &yx = values[iy + 3*ix];
    double &yy = values[iy + 3*iy];
    double &yz = values[iy + 3*iz];
    double &zx = values[iz + 3*ix];
    double &zy = values[iz + 3*iy];
    double &zz = values[iz + 3*iz];

    // custom 
    const double X = Lz - Lx;
    const double Y = Lz - Ly;
    const double Z = max3(X, Y, 0.0);

    xx      = (D*(Lz - Z)*(D*Ly + (Lz - z)*(Ly - Lz + Z))*(Power(Lx - Lz + z,4) + (Power(x,2)*Power(Power(D,2)*Lx + Power(Lz - z,2)*(Lx - Lz + Z) + D*Lx*(Lz - 2*z + Z),2))/(Power(D,2)*Power(Lz - Z,2))))/ (Power(D + Lz - z,2)*Power(Lx - Lz + z,3)*(Ly - Lz + z)*(D*Lx + (Lz - z)*(Lx - Lz + Z)));
    xy = yx = (x*y*(Power(D,2)*Lx + Power(Lz - z,2)*(Lx - Lz + Z) + D*Lx*(Lz - 2*z + Z))*(Power(D,2)*Ly + Power(Lz - z,2)*(Ly - Lz + Z) + D*Ly*(Lz - 2*z + Z)))/(D*Power(D + Lz - z,2)*Power(Lx - Lz + z,2)*Power(Ly - Lz + z,2)*(Lz - Z));
    xz = zx = (x*(D*Ly + (Lz - z)*(Ly - Lz + Z))*(Power(D,2)*Lx + Power(Lz - z,2)*(Lx - Lz + Z) + D*Lx*(Lz - 2*z + Z)))/(D*(D + Lz - z)*Power(Lx - Lz + z,2)*(Ly - Lz + z)*(Lz - Z));
    yy      = (D*(Lz - Z)*(D*Lx + (Lz - z)*(Lx - Lz + Z))*(Power(Ly - Lz + z,4) + (Power(y,2)*Power(Power(D,2)*Ly + Power(Lz - z,2)*(Ly - Lz + Z) + D*Ly*(Lz - 2*z + Z),2))/(Power(D,2)*Power(Lz - Z,2))))/ (Power(D + Lz - z,2)*(Lx - Lz + z)*Power(Ly - Lz + z,3)*(D*Ly + (Lz - z)*(Ly - Lz + Z)));
    yz = zy = (y*(D*Lx + (Lz - z)*(Lx - Lz + Z))*(Power(D,2)*Ly + Power(Lz - z,2)*(Ly - Lz + Z) + D*Ly*(Lz - 2*z + Z)))/(D*(D + Lz - z)*(Lx - Lz + z)*Power(Ly - Lz + z,2)*(Lz - Z));
    zz      = ((D*Lx + (Lz - z)*(Lx - Lz + Z))*(D*Ly + (Lz - z)*(Ly - Lz + Z)))/(D*(Lx - Lz + z)*(Ly - Lz + z)*(Lz - Z));
  }

};
