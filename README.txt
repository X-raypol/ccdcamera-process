README

The SITKConverter class is made to be run in the polarimetry lab at MIT.  For the converter aspect to run, specific Tiff and Txt files from a labview script are necessary. However, the analysis can be run with fits image files from almost anywhere.

The three included scripts cover most of the class can do, but feel free to make custom scripts as most methods are self contained.

Converter_Analysis.py is the main script.  It must be run from the command line with ipython initialized with the matplotlib extension in order for the graphs to show up. Once run it will prompt the user for a fits or tiff file. Tiff files will be analyzed, displayed and converted into both fits image and fits table files, fits image files will be analyzed and then displayed, and fits tables will simply be displayed.

There are two lines commented out in the script representing alternate hot pixel removal options, one makes a hot pixel list based on the given file, and the other removes hot pixels based on a given hot pixel text file as opposed to the simple removal by occurrence currently being implemented.

Drag_drop_converter.py is a secondary converter if the quick view graphs are not necessary. In order to be run there needs to be a batch file made in the same directory with the command: “<path to python executable> <path to this file> %*” Once that is made simply drag the tiff file onto the batch file, and it will make the fits img and table files.

Finally add_info_to_header.py is a fallback for if something breaks and the header is not properly formatted. This can be run from any python executable, whether it is IDLE, ipython, or another, and the prompts, once run, are pretty self explanatory.