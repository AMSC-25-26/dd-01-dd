# Visualization
This project saves the result of the solvers as `.vtk` files.
The program [ParaView](https://www.paraview.org/) is suggested to open and visualize these files.
Additional python script for automatic visualization is provided

## How to Visualize the solution
Open [ParaView](https://www.paraview.org/) and import the solution by clicking on: 
`File > Open` and then select the `.vtk` file to visualize, normally it will be in the `output` directory.

To correctly visualize the VTK files either create your own filter
or use the provided `src/visualization_pipeline.py` script.

## Using the custom script
Check that the [ParaView](https://www.paraview.org/) Python shell is enabled by clicking on: `View`
and ensure that `Python Shell` is ticked.

If the shell is enabled you should see a `Run Script` button on the bottom right of the program.
To apply the script click it and chose the file `src/visualization_pipeline.py`.