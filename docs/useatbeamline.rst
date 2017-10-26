==================================
Use at beamline (windows or leroy)
==================================

The main use for this code is to enable the conversion from the tif files that
the SITK writes and various datafiles written by LabView into fits files.
At the same time, it can provide a quicklook of the data to diagnose
problems before the next frame is taken.

The software is installed on both the windows computer and on leroy.

Data directories
================
When starting to take data, select (or make) a directory on the windows copmuter
*inside* the Dropbox folder. This way, data gets synced to leroy automatically
and you cna work with Python on either the windows computer or leroy as you like.
Additionally, we can use Dropbox to immidiately transfer the data to other machines.
Hopefully, leroy will be nfs mounted at some point, too...

Conversion and Analysis
=======================

Start the IPython command prompt on windows or on leroy by typing
``ipython --matplotlib`` in the terminal. All the following code then can be executed in
IPython.
First, import the module with the code for the X-ray polarimetry beamline::

  >>> import xpolbeamline

Then, create an object that remembers the path to the raw data and the statsfile
so you don't have to type it every time. File path can be absolute or relative to
the current working directory.

.. warning::

   In Python, ``\`` is an escape character for strings, for example ``\n`` mean "newline"
   or ``\t`` means TAB. On Windows ``\`` is also used to separate directories.
   Thus, when you type the pathnames into the Python command line, you need to
   escape the ``\`` with another ``\`` like this: ``C:\\Users\\name``.

So, let's make an object that remembers our directories::

  >>> myrun = xpolbeamline.UI(inpath='C:\\Users\\name\\Dropbox\\xpoldata\\3Oct17',
  ...             statspath='C:\\Users\\name\\Dropbox\\xpolstats\\stats_10_03_17.txt')

In general, white spaces and indentation at the beginning of code lines are
important in Python, but here we just break the line inside of the ``( ... )``
because the pathnames just happen to be too long to fit, so that's OK.

After LabView has written a the tif and txt file to that directory, we can just do:

.. doctest-skip::

  >>> myrun.convert_display_newest()

to convert the most recently modified tif file in that directory to a fits image and
event list and display diagnostic quick view plots.

The plots are interactive, there are little buttons on the top that allow you to
zoom in and out, pan, or save the image to a file.

Using IPython
=============

Both `Python <http://www.python.org>`_ and `IPython <https://ipython.org/>`_ have
extensive documentation. Essentially, `IPython <https://ipython.org/>`_ is just
an enhanced shell for Python. Anything that can be done in plain Python, can be done
in `IPython <https://ipython.org/>`_, too, but the shell is a little nicer and supports
additional features, e.g. in plain Python plotting windows as blocking, while
`IPython <https://ipython.org/>`_ has the ``--matplotlib`` mode where plotting
windows do not block the console.

Here are a few quick tricks:

  - IPython autocompletes with the TAB key: ``myrun.<TAB>``
    will list all members of the ``run`` object, in case you do not remember the exact name
    of a function.

  - In the shell, use the "up arrow" and "down arrow" keys to find the most recently
    used commands. Edit and execute again as needed.
  
  - Quick access to documentation: ``myrun.convert_display?``
    will print the documentation for the method in the terminal (typically
    exit by pressing "q" like in the "less" program on the unix terminal).


Reference / API
===============
.. automodapi:: xpolbeamline.ui
