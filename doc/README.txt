To build the documentation for harmonica execute the following command in a console:

make <format>

Format is typically "html", but to view all available formats execute the make command with no arguments.

After executing the make command, the output files will be in the build directory. "Index.html"
(or "Index" with other format extension) is the root of the documentation pages.

Edit the *.rst reStructuredText files in the source directory to change the layout/content/format
of the pages. Edit the conf.py file in the source directory to change the settings used by Sphinx
to generate the documentation.