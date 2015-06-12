import unittest
from magnumfe import *

set_log_active(False)

class CacheTest(unittest.TestCase):

  def test_initial_update(self):
    mesh = UnitCubeMesh(1,1,1)
    state = State(mesh)

    cache = Cache()
    self.assertTrue(cache.requires_update(state))

  def test_change_state(self):
    mesh = UnitCubeMesh(1,1,1)
    state1 = State(mesh)
    state2 = State(mesh)

    cache = Cache()
    count = 0

    if cache.requires_update(state1): count += 1
    if cache.requires_update(state1): count += 1
    self.assertEqual(1, count)

    if cache.requires_update(state2): count += 1
    self.assertEqual(2, count)

  def test_update_required(self):
    mesh = UnitCubeMesh(2, 2, 2)
    state = State(mesh, m = Constant((1.0, 0.0, 0.0)), j = Constant((0.0, 0.0, 0.0)))

    cache = Cache("m", "t")
    count = 0

    if cache.requires_update(state): count += 1
    self.assertEqual(1, count)
    if cache.requires_update(state): count += 1
    self.assertEqual(1, count)

    state.t = 1.0

    if cache.requires_update(state): count += 1
    self.assertEqual(2, count)
    if cache.requires_update(state): count += 1
    self.assertEqual(2, count)

    state.m = Constant((0.0, 1.0, 0.0))

    if cache.requires_update(state): count += 1
    self.assertEqual(3, count)
    if cache.requires_update(state): count += 1
    self.assertEqual(3, count)

    state.j = Constant((1.0, 0.0, 0.0))

    if cache.requires_update(state): count += 1
    self.assertEqual(3, count)

if __name__ == '__main__':
    unittest.main()
