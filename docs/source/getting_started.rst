Getting Started
===============


.. code-block:: python

  >>> from pyma import pyhma
  >>> from pyhma import ReadVASP, Compute
  >>> rv = ReadVASP(*filenames)
  >>> rv.read()
  >>> c = Compute()
  >>> c.compute()
  >>> data = c.get_statistics()
  >>> print(data_['uahc'] , data_['uahm'] , data_['pahc'] , data_['pahm'])
  >>> print('Uah_conv (eV/atom):', data['uahc'])
  >>> print('Uah_hma  (eV/atom):', data['uahm'])
  >>> print('Pah_conv (GPa):', data['pahc'])
  >>> print('Pah_hma  (GPa):', data['pahm'])




.. math::
   
