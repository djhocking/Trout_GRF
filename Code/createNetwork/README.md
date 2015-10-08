# Establish Network Structure

This section incorporates the observed sampling locations into the stream network structure. Scripts are set up to work with the spatial layers of the custom National Hydrography Dataset High Resolution Delineation - Version 1 (NHDHRDV1), and should be transferable to Version 2 when it is released. Sample points are identified by a unique ID specified by the user. Confluence nodes are defined by the FEATUREID of the downstream segment. The next section describes the data layers required to establish the stream network structure for the GMRF model.



## Data Requirements

Establishing the network structure for the GMRF model requires the use of existing spatial layers. Watershed hydrography and observed sample locations are needed as input to the process. The current versions rely on the the NHDHRD dataset, though other products may be substituted with some adjustment to the scripts.

#### 1. Sample Locations Shapefile
A point layer developed entirely by the user represents the sample locations. This layer is identified in scripts as `pointsFile`. Each site should have a unique ID. The points may or may not be snapped to the watershed flowlines. The most accurate snapping method is described in the optional "Optional Preprocessing" section below.

#### 2. Watershed Boundary Shapefile
A polygon layer represents the outline of the watershed being evaluated. This layer is identified in scripts as `processingBoundary`. This is typically developed from the NHD High Resolution Watershed Boundary Dataset (WBD). A single Hydrologic Unit Code (HUC) outline or a custom combination of multiple HUCs may be used for this layer. Whatever the area, it should have one and only one outflow location.

#### 3. Watershed Flowlines Shapefile
A vector layer represents the flowlines of the watershed being evaluated. This layer is identified in scripts as `flowlinesFile`. These flowlines must be selected from the NHDHRD dataset for the script to function properly. If a different hydrography layer is desired, the scripts will need to be adapated to the network structure specifications. These flowlines are typically selected from the full-range flowline layer by using the Watershed Boundary layer. A caution: the two datasets differ slightly so the flowlines should be visually inspected along the boundary to ensure that the network within the desired watershed is properly selected.

#### 4. Delineated Catchments RData File
The delineated catchments for the NHDHRD are defined in an RData file, specifically as a list object containing a sub-list for each catchment. Each sub-list identifies all of the FEATUREIDs of the catchments in the upstream network from the parent catchment.



## Processing Steps

The steps in this section create the input family table for the GMRF model.

### Network Spatial Processing

The `observedSitesNetworkProcessing.py` script is executed in ArcPython. All input layers must be defined in the same projected coordinate system. The sites that fall inside of the watershed are selected from the input dataset. These points are snapped to the stream network (within the specified distance). For each sample location, the position along the stream segment ("FMEAS") and the latitude and longitude coordinates are calulated. Spatial coordinates are also calculated for the upstream point of each stream segment, which represent confluence nodes. These values are output in two separate tables: `occupancySites.dbf` for the sample locations and `flowlinesTable.dbf` for the confluence nodes.

#### User Inputs Description

|      Variable           |    Layer Type     |                                    Description                                                             |
|      :-----:            |   :----------:    |                                   :-----------                                                             |
| `baseDirectory`         | Character String  | The path to the directory where all of the output files will be saved                                      |
| `pointsFile`            | Point Shapefile   | The path to the point shapefile of observed sample locations                                               |
| `uniqueIDField`         | Character String  | The field representing the unique identifiers of the observed sample locations                             |
| `flowlinesFile`         | Line Shapefile    | The path to the flowlines of the subject watershed                                                         |       
| `processingBoundary`    | Polygon Shapefile | A polygon representing the outline of the subject watershed                                                |
| `snapDistanceM`         | Character String  | The maximum distance, in meters, from the flowlines where the points will be snapped to the stream network |



### Create the Family Table
The "createNetworkFamily.R" script uses the output from the spatial processing to create the family table for input to the GMRF model. The RData file containing the dealineated catchments lists is also used in this process. The output of this script is the `family.csv` file.

#### User Inputs Description
|      Variable           |        Type       |                                    Description                                    |
|      :-----:            |       :----:      |                                   :-----------                                    |
| `baseDirectory`         | Character String  | The path to the directory where all of the output files will be saved             |
| `mouthID`               | Numeric           | The FEATUREID of the farthest downstream segment                                  |
| `delineationFilePath`   | Character String  | The filepath to the RData file with the list of delineated catchments from NHDHRD |



## Optional Preprocessing

### Map Sample Locations to Stream Network 

The `mapSampleSitesToNetwork.py` script snaps the observed sample points to the appropriate stream network. The various flowline versions from the NHDHRD are utilized to ensure the most accurate snapping. All high resolution streams as well as potentially unidentified streams (DEM-based channels) are accounted for. Points are identified based on which network classification they are snapped to. The `mapSampleSitesToNetwork - Description` PDF file provides an in depth description of the processing completed by this script. 

The script is executed in Arc Python and assumes all input spatial layers are in the same projected coordinate system. Two point layers are output by the script, one where all of the points assigned to the high resolution flowlines (`snappedPointsDetailed_[datasetID]`) and one with all of the flowlines assigned to the truncated flowlines (`snappedPointsTruncated_[datasetID]`). The latter of these two layers is intended to be used as the `pointsFile` variable in the "Network Spatial Processing". For more on these two flowline versions consult the `mapSampleSitesToNetwork - Description` PDF file or the NHDHRD documentation. A table describing the breakdown of sample point assignment to flowlines is also output.


### Specify inputs
|      Variable           |    Layer Type    |                                    Description                                                                                       |
|      :-----:            |   :----------:   |                                   :-----------                                                                                       |
| `networkGrid_DEM`       | Raster           | An NHDHRD product - the stream raster based on the DEM layer and a minimum flow accumulation threshold area                          |
| `networkGrid_Detailed`  | Raster           | An NHDHRD product - the high resolution stream raster                                                                                |
| `networkGrid_Truncated` | Raster           | An NHDHRD product - the truncated streams raster, based on masking the high resolution layer to the threshold                        | 
| `points`                | Point Shapefile  | The locations of the sample sites containing a unique ID for each site                                                               |
| `datasetID`             | Character String | A unique name used to identify the specific dataset being processed                                                                  |
| `bufferInMeters`        | Character String | The maximum distance over which points will be snapped to the stream network                                                         |
| `gisDirectory`          | Character String | The path to the directory where the spatial files will be saved. A geodatabase ('hydrography.gdb') will be created in this directory if it does not exist |