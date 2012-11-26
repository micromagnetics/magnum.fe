import unittest
from dolfin import *
from magnumfe import *
import numpy

mesh = UnitCubeMesh(1, 1, 1)
VS = FunctionSpace(mesh, "CG", 1)
VV = VectorFunctionSpace(mesh, "CG", 1)
fa = interpolate(Expression(('x[0]*x[1]+1.0', 'x[1]*x[2]+1.0', '1.0')), VV)
fb = interpolate(Expression(('x[1]*x[2]+1.0', '1.0', 'x[0]*x[1]+1.0')), VV)

class DofMapMeshTest(unittest.TestCase):

  def test_trans_scalar_product_matrix(self):
    mesh = UnitCubeMesh(1, 1, 1)

    scalar_product = TransScalarProductMatrix(VS, VV, fa)
    A = DofAssembler.assemble(scalar_product)

    c = Function(VS).vector()
    A.transpmult(fb.vector(), c)
    fc = Function(VS, c)

    # check scalar product on nodes of the unit cube
    self.assertSPEqualAtPoint(fa, fb, fc, (0.0, 0.0, 0.0))
    self.assertSPEqualAtPoint(fa, fb, fc, (0.0, 0.0, 1.0))
    self.assertSPEqualAtPoint(fa, fb, fc, (0.0, 1.0, 0.0))
    self.assertSPEqualAtPoint(fa, fb, fc, (1.0, 1.0, 1.0))

  def test_scalar_product_matrix(self):
    mesh = UnitCubeMesh(1, 1, 1)

    scalar_product = ScalarProductMatrix(VS, VV, fa)
    A = DofAssembler.assemble(scalar_product)

    c = A*fb.vector()
    fc = Function(VS, c)

    # check scalar product on nodes of the unit cube
    self.assertSPEqualAtPoint(fa, fb, fc, (0.0, 0.0, 0.0))
    self.assertSPEqualAtPoint(fa, fb, fc, (0.0, 0.0, 1.0))
    self.assertSPEqualAtPoint(fa, fb, fc, (0.0, 1.0, 0.0))
    self.assertSPEqualAtPoint(fa, fb, fc, (1.0, 1.0, 1.0))

  def test_normalized_vector(self):
    normalized_vector = NormalizedVector(VV, fa)
    c = DofAssembler.assemble(normalized_vector)
    fc = Function(VV, c)

    # check norm on nodes of the unit cube
    self.assertNormalizedAtPoint(fa, fc, (0.0, 0.0, 0.0))
    self.assertNormalizedAtPoint(fa, fc, (0.0, 0.0, 1.0))
    self.assertNormalizedAtPoint(fa, fc, (0.0, 1.0, 0.0))
    self.assertNormalizedAtPoint(fa, fc, (1.0, 1.0, 1.0))

  def assertSPEqualAtPoint(self, f1, f2, f12, point):
    v1  = numpy.zeros((3,), dtype="d")
    v2  = numpy.zeros((3,), dtype="d")
    v12 = numpy.zeros((1,), dtype="d")
    f1.eval(v1, numpy.array(point))
    f2.eval(v2, numpy.array(point))
    f12.eval(v12, numpy.array(point))
    sp = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
    self.assertAlmostEqual(v12[0], sp)

  def assertNormalizedAtPoint(self, f1, f2, point):
    v1  = numpy.zeros((3,), dtype="d")
    v2  = numpy.zeros((3,), dtype="d")
    f1.eval(v1, numpy.array(point))
    f2.eval(v2, numpy.array(point))
    
    #norm1 = v1[0]**2 + v1[1]**2 + v1[2]**2
    norm2 = v2[0]**2 + v2[1]**2 + v2[2]**2

    self.assertAlmostEqual(norm2, 1.0)
    self.assertAlmostEqual(v1[0] / v1[1], v2[0] / v2[1])
    self.assertAlmostEqual(v1[1] / v1[2], v2[1] / v2[2])

if __name__ == '__main__':
    unittest.main()
