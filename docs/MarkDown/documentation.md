# Doxygen documentation
The code as been documented using javadoc style comments.
Doxygen is used to produce a html interactive documentation based on those.
It also crates a $\LaTeX$ based pdf.

## Requirements
* `doxygen` to compile the documentation. if not available in the system the option to compile the documentation will be turned off
* `graphviz` to create dependency/hierarchy graph for each classes.
* some kind of latex compiler is also necessary, usually provided with doxygen under linux.

## Configuration
All doxygen configurations can be found at [Doxyfile](../Doxygen/Doxyfile).

## Building the documentation
Follow the instruction in [Compilation and execution](building.md) and then run: 
```bash
make docs
```

## Reading the documentation
The results can be found in `docs/Doxygen/html/index.html` for the html version
and in `docs/Doxygen/latex/refman.pdf` for the pdf version.
Simply open them with your favorite file.