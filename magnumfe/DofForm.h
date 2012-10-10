namespace magnumfe {
  class DofForm {
    public:
      virtual unsigned int rank();

      virtual unsigned int num_coefficients();

      // A: output 
      virtual eval(double* A, const double * const * w);
  }
}
