#ifndef _METRIC_MATRIX_CUBE_
#define _METRIC_MATRIX_CUBE_

#include <dolfin.h>

namespace magnumfe {

  class MetricMatrixCube: public dolfin::Expression
  {
    protected:
      const dolfin::Array<double> size;
      const int  coord;

    public:
      MetricMatrixCube(const dolfin::Array<double>& size, int coord): dolfin::Expression(3, 3), size(size), coord(coord) {}
      ~MetricMatrixCube(){}

      // returns the biggest of three double values
      static double max3(double a, double b, double c);

      // returns the biggest of three double values
      static double min3(double a, double b, double c);

      // alias Power for mathematica C export
      static double Power(double x, int y);

      // implement eval
      void eval(dolfin::Array<double>& values, const dolfin::Array<double>& pos) const;
  };

};

#endif
