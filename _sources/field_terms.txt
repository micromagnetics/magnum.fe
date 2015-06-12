.. module:: magnumfe

Field Terms
===========

Each contribution to the effective field is represented by a subclass of :class:`LLGTerm`.
A subclass has at least to implement the :code:`form_term_rhs` method.

Available field classes are :class:`DemagField`, :class:`ExchangeField`, :class:`ExternalField`, :class:`OerstedField`, :class:`SpinCurrent`, :class:`UniaxialAnisotropyField`.

:class:`LLGTerm`
++++++++++++++++

.. autoclass:: LLGTerm
   :members:

:class:`DemagField`
+++++++++++++++++++

.. autoclass:: DemagField
   :members:

:class:`ExchangeField`
++++++++++++++++++++++

.. autoclass:: ExchangeField
   :members:

:class:`ExternalField`
++++++++++++++++++++++

.. autoclass:: ExternalField
   :members:

:class:`OerstedField`
+++++++++++++++++++++

.. autoclass:: OerstedField
   :members:

:class:`SpinTorque`
+++++++++++++++++++++++++++++++++++

.. autoclass:: SpinTorque
   :members:

:class:`UniaxialAnisotropyField`
++++++++++++++++++++++++++++++++

.. autoclass:: UniaxialAnisotropyField
   :members:
