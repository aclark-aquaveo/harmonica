Using the Python API
=====================================

This tutorial demonstrates the use of the hamonica Python API. Amplitude, phase,
and speed tidal harmonics are extracted for a small set of points from the TPXO8, ADCIRC,
and LeProvost tidal databases. The tutorial also provides an example of obtaining the 
amplitude, frequency, speed, earth tidal reduction factor, equilibrium argument, and
nodal factor of a constituent at a specified time.

To complete this tutorial you must have the harmonica package installed. See :ref:`Installation` for
instructions on installing the harmonica package to a Python environment.

We will be using the script located at :file:`{harmonica root}/tutorials/python_api/test.py`.
The code is also provided below.

.. literalinclude:: ../../tutorials/python_api/test.py
   :linenos:

Executing :file:`test.py` will fetch the required tidal resources from the Internet if
they do not already exist in the harmonica data directory. Edit the config variables in
:file:`__init__.py` to change the data directory or specify a preexisting data directory.

After executing the script, output from the tidal harmonic extraction can be viewed in
:file:`tidal_test.out` located in the Python API tutorial directory.

The following lines set up the point locations where tidal harmonic components will be
extracted. Locations should be specified as tuples of latitude and longitude degrees.

.. literalinclude:: ../../tutorials/python_api/test.py
   :lines: 7-12
   :lineno-match:

In addition to the points of interest, create a list of constituents to query for.

.. literalinclude:: ../../tutorials/python_api/test.py
   :lines: 14-14
   :lineno-match:

Next, construct the tidal constituent extractor interface. This example starts out by
using the 'leprovost' model. 'tpxo8' is the default model used when none is specified.

.. literalinclude:: ../../tutorials/python_api/test.py
   :lines: 15-15
   :lineno-match:

The next section of code demonstrates using the interface to get amplitudes, frequencies,
speeds, earth tidal reduction factors, equilibrium arguments, and nodal factors for
specified constituents at a specified time. The tidal model specified at construction has
no effect on this functionality.

.. literalinclude:: ../../tutorials/python_api/test.py
   :lines: 17-18
   :lineno-match:

The last block uses the ADCIRC, LeProvost, and TPXO interfaces to extract tidal harmonic
constituents for a list of locations and constituents. The optional 'model' argument can
be specified to switch between tidal models. Note that when getting the TPXO components,
we pass in a third positional argument (kwarg is 'positive_ph'). If True, all output
phases will be positive.

.. literalinclude:: ../../tutorials/python_api/test.py
   :lines: 25-42
   :lineno-match:

The LeProvost model is freely distributed. FES2014 is an updated version of the model that
significantly increases the grid resolution and number of supported constituents. The
data files for FES2014 cannot be openly distributed, but local copies of the files
can be used if they exist. See :file:`resource.py` for the expected filenames.

.. literalinclude:: ../../tutorials/python_api/test.py
   :lines: 44-48
   :lineno-match: