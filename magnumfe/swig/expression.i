namespace dolfin
{
  class Expression// : public GenericFunction
  {
  public:

    Expression();
    Expression(uint dim);
    Expression(uint dim0, uint dim1);
    Expression(std::vector<uint> value_shape);
    Expression(const Expression& expression);

    virtual ~Expression();

    virtual void eval(Array<double>& values,
                      const Array<double>& x,
                      const ufc::cell& cell) const;

    virtual void eval(Array<double>& values, const Array<double>& x) const;

    virtual uint value_rank() const;

    virtual uint value_dimension(uint i) const;

    virtual void restrict(double* w,
                          const FiniteElement& element,
                          const Cell& dolfin_cell,
                          const ufc::cell& ufc_cell) const;

    virtual void compute_vertex_values(Array<double>& vertex_values,
                                       const Mesh& mesh) const;

  };

  class GenericFunction //: public ufc::function, public Variable
  {
  public:

    GenericFunction();
    virtual ~GenericFunction();
    virtual uint value_rank() const = 0;
    virtual uint value_dimension(uint i) const = 0;
    virtual void eval(Array<double>& values, const Array<double>& x,
                      const ufc::cell& cell) const;
    virtual void eval(Array<double>& values, const Array<double>& x) const;
    virtual void restrict(double* w,
                          const FiniteElement& element,
                          const Cell& dolfin_cell,
                          const ufc::cell& ufc_cell) const = 0;
    virtual void compute_vertex_values(Array<double>& vertex_values,
                                       const Mesh& mesh) const = 0;
    virtual void gather() const {}
    double operator() (double x);
    double operator() (double x, double y);
    double operator() (double x, double y, double z);
    double operator() (const Point& p);
    void operator() (Array<double>& values, double x);
    void operator() (Array<double>& values, double x, double y);
    void operator() (Array<double>& values, double x, double y, double z);
    void operator() (Array<double>& values, const Point& p);
    uint value_size() const;
    virtual void evaluate(double* values,
                          const double* coordinates,
                          const ufc::cell& cell) const;
  };
}
