Using the Python API
=====================================

This tutorial demonstrates the use of the hamonica Python API. Amplitude, phase,
and speed tidal harmonics are extracted for a small set of points from the TPXO8, ADCIRC
(Northwest Atlantic and Northeast Pacific), and LeProvost tidal databases. The tutorial
also provides an example of obtaining the frequency, earth tidal reduction factor, amplitude,
nodal factor, and equilibrium argument of a constituent at a specified time using the ADCIRC
and LeProvost tidal databases.

To complete this tutorial you must have the harmonica package installed. See :ref:`Installation` for
instructions on installing the harmonica package to a Python environment.

We will be using the script located at :file:`{harmonica root}/tutorials/python_api/test.py`.
The code is also provided below.

.. literalinclude:: ../../tutorials/python_api/test.py
   :linenos:

Executing :file:`test.py` will fetch the required tidal resources from the Internet if
they do not already exist in the working directory. To specify existing resources use the
following command line arguments:

-l path     Path to the LeProvost \*.legi files
-a file     Path and filename of the ADCIRC Northwest Atlantic executable
-p file     Path and filename of the ADCIRC Northeast Pacific executable

You may also specify a temporary working directory that will be used by the ADCIRC
database extractors (default is current working directory).

-w path     Path ADCIRC extractors will use as a temporary working directory

After executing the script, output from the tidal harmonic extraction can be viewed in
:file:`tidal_test.out` located in the Python API tutorial directory.

The following lines set up the point locations we will be extracting tidal harmonic
components for. Locations should be specified as tuples of latitude and longitude
degrees. Note that we have split locations in the Atlantic and Pacific. LeProvost and
TPXO support all ocean locations, but the ADCIRC database is restricted to one or the other.

.. literalinclude:: ../../tutorials/python_api/test.py
   :lines: 19-33
   :lineno-match:

The next section of code sets up the tidal harmonic extraction interfaces. For the ADCIRC
extractors, the tidal region must be specified at construction.

.. literalinclude:: ../../tutorials/python_api/test.py
   :lines: 36-43
   :lineno-match:

The ADCIRC and LeProvost extractors can also provide frequencies, earth tidal reduction factors, amplitudes, nodal factors,
and equilibrium arguments for specified constituents at a specified time. The next section of code demonstrates this.

.. literalinclude:: ../../tutorials/python_api/test.py
   :lines: 48-51
   :lineno-match:

The last block uses the ADCIRC, LeProvost, and TPXO interfaces to extract tidal harmonic constituents for a list of
locations and constituents.

.. literalinclude:: ../../tutorials/python_api/test.py
   :lines: 62-85
   :lineno-match: