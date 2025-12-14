from paraview.simple import *


data_source = GetActiveSource()

if data_source is None:
    raise RuntimeError("No data source selected. Please load your VTK file and select it in the Pipeline Browser before running this script.")


global_data = MergeBlocks(Input=data_source)


layout = GetLayout()

chart_view = CreateView('XYChartView')
AssignViewToLayout(view=chart_view, layout=layout, hint=0)

display = Show(global_data, chart_view)

display.UseIndexForXAxis = 0 

display.XArrayName = 'Points_X'

display.SeriesVisibility = ['u'] 
display.SeriesColor = ['u', '0', '0', '0'] 

chart_view.ChartTitle = "Global Solution"
chart_view.LeftAxisTitle = "Solution"
chart_view.BottomAxisTitle = "X Coordinate"

Render()