import unittest

testmodules = [
    'demag_field_fk_test',
    'demag_field_st_test',
    'dof_map_test',
    'llg_test',
    'material_test',
    'mesher_test',
    'state_test',
    'wrapped_mesh_test'
    ]

suite = unittest.TestSuite()

for t in testmodules:
  try:
    # If the module defines a suite() function, call it to get the suite.
    mod = __import__(t, globals(), locals(), ['suite'])
    suitefn = getattr(mod, 'suite')
    suite.addTest(suitefn())
  except (ImportError, AttributeError):
    subsuite = unittest.TestSuite()
    # else, just load all the test cases from the module.
    subsuite.addTest(unittest.defaultTestLoader.loadTestsFromName(t))
    suite.addTest(subsuite)

unittest.TextTestRunner().run(suite)
