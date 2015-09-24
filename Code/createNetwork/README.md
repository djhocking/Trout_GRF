# Establish Network Structure

This section incorporates the observed sampling locations into the stream network structure. Scripts are set up to work with the spatial layers of the custom-made delineation of high resolution catchments. 



## mapFishSitesToNetwork.py

- This script is used to snap points to to the appropriate stream network.
- This script snaps the sampling locations to the nearest stream. Points are identified based on which network classification they are snapped to. All points are included in the output, but only the points falling on the "Truncated Flowlines" are used in the network analysis. The script is executed in Arc Python.
- Only points within the defined distance (`bufferInMeters`) will be snapped and used. 
- The script assumes all input spatial layers are in the same projected coordinate system (Albers).
- The mapFishSitesToNetwork.pdf file provides an in depth description of the processing in this script.

### Specify inputs

|      Variable           |    Layer Type    |                                    Description                                                                                       |
|      :-----:            |   :----------:   |                                   :-----------                                                                                       |
| `networkGrid_DEM`       | Raster           | An NHDHRD product - the stream raster based on the DEM layer and a minimum flow accumulation threshold.                              |
| `networkGrid_Detailed`  | Raster           | An NHDHRD product - the high resolution stream raster.                                                                               |
| `networkGrid_Truncated` | Raster           | An NHDHRD product - the truncated streams raster, based on masking the high resolution layer to the threshold.                       | 
| `points`                | Point Shapefile  | The locations of the sample sites containing a unique ID for each site.                                                              |
| `datasetID`             | Character String | A unique name used to identify the specific dataset being processed.                                                                 |
| `bufferInMeters`        | Character String | The maximum distance over which points will be snapped to the stream network.                                                        |
| `workingDirectory`      | Character String | The path to the directory where the spatial files will be saved. A geodatabase ('hydrography.gdb') will be created in this directory |



## observedSitesNetworkProcessing.py
