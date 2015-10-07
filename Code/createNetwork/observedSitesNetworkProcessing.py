# ==============
# Specify inputs
# ==============

# Base Directory
baseDirectory = "C:/KPONEIL/GitHub/projects/Trout_GRF/Data/gisFiles/modelVersions/SusquehannaWest_Threepass"

# Sample Locations Shapefile
pointsFile = "C:/KPONEIL/GitHub/projects/Trout_GRF/Data/gisFiles/hydrography.gdb/snappedPointsTruncated_threepass"

# Unique ID
uniqueIDField = "GIS_Key"

# Flowlines Shapefile
flowlinesFile = baseDirectory + "/spatialModel_Inputs.gdb/flowlines"

# Processing Boundary Polygon
processingBoundary = baseDirectory + "/spatialModel_Inputs.gdb/processingBoundary"

# Snapping Distance (in meters)
snapDistanceM = "100"


# ===================
# File pre-processing
# ===================

# Create version geodatabase
workingDirectory = baseDirectory + "/spatialModel_Processing.gdb"
if not arcpy.Exists(workingDirectory): arcpy.CreateFileGDB_management (baseDirectory, "spatialModel_Processing", "CURRENT")

# Map Definitions
mxd = arcpy.mapping.MapDocument("CURRENT")
df = arcpy.mapping.ListDataFrames(mxd)[0]

# Copy flowlines file for editing
flowlines = arcpy.FeatureClassToFeatureClass_conversion(flowlinesFile, 
														workingDirectory, 
														"flowlinesEditing")


# Select points in the watershed
# ------------------------------
arcpy.MakeFeatureLayer_management(pointsFile,    "pointsFileLyr")

arcpy.SelectLayerByLocation_management ("pointsFileLyr", 
											"INTERSECT", 
											processingBoundary, 
											"", 
											"NEW_SELECTION")
											
points = arcpy.FeatureClassToFeatureClass_conversion("pointsFileLyr", 
														workingDirectory, 
														"networkPoints")

arcpy.Delete_management("pointsFileLyr")


# Calculate fields for processing
# -------------------------------
# Rename unique IDs
arcpy.AddField_management(points, "nodeID", "TEXT")
arcpy.CalculateField_management (points, "nodeID", "!" + uniqueIDField + "!", "PYTHON_9.3")

arcpy.AddField_management(flowlines, "nodeID", "TEXT")
arcpy.CalculateField_management (flowlines, "nodeID", """"N_" + str(!FEATUREID!)""", "PYTHON_9.3")


# Create field for creating routes (ensures correct direction)
arcpy.AddField_management(flowlines, "fromMeas", "SHORT")
arcpy.CalculateField_management (flowlines, "fromMeas", 0, "PYTHON_9.3")


# ==================
# Network Processing
# ==================

# Snap the points to the flowlines
arcpy.Snap_edit(points,
					[[flowlines, "EDGE", snapDistanceM + " Meters"]])


# Determine positions
# -------------------
# Create routes
routes = arcpy.CreateRoutes_lr(flowlines,
									"FEATUREID",
									workingDirectory + "/flowRoutes",
									"TWO_FIELDS",
									"fromMeas", 
									"LengthKM", 
									"", "", "",
									"IGNORE",
									"INDEX")


# Find the distance of the points along the routes
occupancySites = arcpy.LocateFeaturesAlongRoutes_lr(points,
														routes,
														"FEATUREID",
														snapDistanceM + " METERS", 
														baseDirectory + "/pointLocations.dbf",
														"FEATUREID POINT FMEAS TMEAS",
														"FIRST",
														"DISTANCE",
														"ZERO",
														"FIELDS",
														"M_DIRECTON")


# Calculate flowline node locations
# ---------------------------------

# Project to a geographic coordinate system to calculate lat/lon
spatialRef_flow = arcpy.Describe(flowlines).spatialReference

flowlinesGCS = arcpy.Project_management(flowlines,
											workingDirectory + "/flowlinesGCS",
											spatialRef_flow.GCS)

# Calculate lat/lon fields (at upstream point of line)											
arcpy.AddField_management(flowlinesGCS, "NodeLon", "DOUBLE")
arcpy.AddField_management(flowlinesGCS, "NodeLat",  "DOUBLE")			
arcpy.CalculateField_management (flowlinesGCS, "NodeLon", "!Shape!.positionAlongLine(1.0,True).firstPoint.X", "PYTHON_9.3")
arcpy.CalculateField_management (flowlinesGCS, "NodeLat", "!Shape!.positionAlongLine(1.0,True).firstPoint.Y", "PYTHON_9.3")

# Add lat/lon fields to the main flowlines table
arcpy.MakeTableView_management(flowlinesGCS, "flowlinesGCSTable")
arcpy.MakeTableView_management(flowlines, "flowlinesTable")	

arcpy.JoinField_management("flowlinesTable", "FEATUREID", "flowlinesGCSTable", "FEATUREID", ["NodeLat","NodeLon"])

arcpy.TableToTable_conversion("flowlinesTable", 
								baseDirectory,
								"flowlinesTable.dbf")


								
# Point Node locations
# --------------------

# Project to a geographic coordinate system to calculate lat/lon
spatialRef_points = arcpy.Describe(points).spatialReference

pointsGCS = arcpy.Project_management(points,
										workingDirectory + "/pointsGCS",
										spatialRef_points.GCS)

# Calculate lat/lon fields										
arcpy.AddField_management(pointsGCS, "NodeLon", "DOUBLE")
arcpy.AddField_management(pointsGCS, "NodeLat",  "DOUBLE")			
arcpy.CalculateField_management (pointsGCS, "NodeLon", "!shape.centroid.X!", "PYTHON_9.3")
arcpy.CalculateField_management (pointsGCS, "NodeLat", "!shape.centroid.Y!", "PYTHON_9.3")

# Add lat/lon fields to the main points table
arcpy.MakeTableView_management(pointsGCS, "pointsGCSTable")	
arcpy.MakeTableView_management(occupancySites, "occupancyTable")	

arcpy.JoinField_management("occupancyTable", "nodeID", "pointsGCSTable", "nodeID", ["NodeLat","NodeLon"])

arcpy.TableToTable_conversion("occupancyTable", 
								baseDirectory,
								"occupancySites.dbf")

							
								

# Delete Interim Files
# --------------------
arcpy.Delete_management("occupancyTable")
arcpy.Delete_management("pointsGCSTable")
arcpy.Delete_management("flowlinesGCSTable")
#arcpy.Delete_management(flowlines)								
#arcpy.Delete_management(flowlinesGCS)							
#arcpy.Delete_management(routes)
#arcpy.Delete_management(points)
#arcpy.Delete_management(pointsGCS)
#arcpy.Delete_management(occupancySites)

