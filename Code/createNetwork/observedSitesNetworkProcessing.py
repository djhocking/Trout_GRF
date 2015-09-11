# Specify inputs
# --------------

versionID = 

#baseDirectory = "C:/KPONEIL/streamNetwork/UpperWhiteRiver"
baseDirectory = "C:/KPONEIL/streamNetwork/PA/networkAnalysis/SusquehannaWest_Threepass"

#pointsFile = "C:/KPONEIL/streamNetwork/occupancy/networkAnalysis.gdb/snappedPointsDetailed_VTFS"

pointsFile = "C:/KPONEIL/streamNetwork/PA/occupancy/spatial/hydrography.gdb/snappedPointsDetailed_threepass"

# File pre-processing
# -------------------

# Create version geodatabase
workingDirectory = baseDirectory + "/spatialModel.gdb"
if not arcpy.Exists(workingDirectory): arcpy.CreateFileGDB_management (baseDirectory, "spatialModel", "CURRENT")

flowlinesSource = workingDirectory + "/flowlines"
processingBoundary = workingDirectory + "/processingBoundary"


# Map Definitions
mxd = arcpy.mapping.MapDocument("CURRENT")
df = arcpy.mapping.ListDataFrames(mxd)[0]


flowlines = arcpy.FeatureClassToFeatureClass_conversion(flowlinesSource, 
														workingDirectory, 
														"flowlinesEditing")


# Process points
# --------------
arcpy.MakeFeatureLayer_management(pointsFile,    "pointsFileLyr")



arcpy.SelectLayerByLocation_management ("pointsFileLyr", 
											"INTERSECT", 
											processingBoundary, 
											"", 
											"NEW_SELECTION")

arcpy.SelectLayerByAttribute_management ("pointsFileLyr", 
											"SUBSET_SELECTION", 
											""" LocationClass =  'Truncated Network'""")	
											
points = arcpy.FeatureClassToFeatureClass_conversion("pointsFileLyr", 
														workingDirectory, 
														"networkPoints")

arcpy.Delete_management("pointsFileLyr")


# Calculate unique IDs
arcpy.AddField_management(points, "nodeID", "TEXT")
arcpy.CalculateField_management (points, "nodeID", "!location_id!", "PYTHON_9.3")

arcpy.AddField_management(flowlines, "nodeID", "TEXT")
arcpy.CalculateField_management (flowlines, "nodeID", """"N_" + str(!FEATUREID!)""", "PYTHON_9.3")

# Create field for creating routes
arcpy.AddField_management(flowlines, "fromMeas", "SHORT")
arcpy.CalculateField_management (flowlines, "fromMeas", 0, "PYTHON_9.3")


# Snap the points to the flowlines
arcpy.Snap_edit(points,
					[[flowlines, "EDGE", " 100 Meters"]])


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
														"100 METERS", 
														baseDirectory + "/pointLocations.dbf",
														"FEATUREID POINT FMEAS TMEAS",
														"FIRST",
														"DISTANCE",
														"ZERO",
														"FIELDS",
														"M_DIRECTON")


# Flowline Node locations
# -----------------------
spatialRef_flow = arcpy.Describe(flowlines).spatialReference

flowlinesGCS = arcpy.Project_management(flowlines,
											workingDirectory + "/flowlinesGCS",
											spatialRef_flow.GCS)
			
arcpy.AddField_management(flowlinesGCS, "NodeLon", "DOUBLE")
arcpy.AddField_management(flowlinesGCS, "NodeLat",  "DOUBLE")			
arcpy.CalculateField_management (flowlinesGCS, "NodeLon", "!Shape!.positionAlongLine(1.0,True).firstPoint.X", "PYTHON_9.3")
arcpy.CalculateField_management (flowlinesGCS, "NodeLat", "!Shape!.positionAlongLine(1.0,True).firstPoint.Y", "PYTHON_9.3")


arcpy.MakeTableView_management(flowlinesGCS, "flowlinesGCSTable")
arcpy.MakeTableView_management(flowlines, "flowlinesTable")	


arcpy.JoinField_management("flowlinesTable", "FEATUREID", "flowlinesGCSTable", "FEATUREID", ["NodeLat","NodeLon"])

arcpy.TableToTable_conversion("flowlinesTable", 
								baseDirectory,
								"flowlinesTable.dbf")


								
# Point Node locations
# --------------------
spatialRef_points = arcpy.Describe(points).spatialReference

pointsGCS = arcpy.Project_management(points,
										workingDirectory + "/pointsGCS",
										spatialRef_points.GCS)
			
arcpy.AddField_management(pointsGCS, "NodeLon", "DOUBLE")
arcpy.AddField_management(pointsGCS, "NodeLat",  "DOUBLE")			
arcpy.CalculateField_management (pointsGCS, "NodeLon", "!shape.centroid.X!", "PYTHON_9.3")
arcpy.CalculateField_management (pointsGCS, "NodeLat", "!shape.centroid.Y!", "PYTHON_9.3")

arcpy.MakeTableView_management(pointsGCS, "pointsGCSTable")	
arcpy.MakeTableView_management(occupancySites, "occupancyTable")	

arcpy.JoinField_management("occupancyTable", "location_i", "pointsGCSTable", "location_id", ["NodeLat","NodeLon"])

arcpy.TableToTable_conversion("occupancyTable", 
								baseDirectory,
								"occupancySites.dbf")

							
								

# Delete Interim Files
# --------------------
arcpy.Delete_management(flowlines)								
arcpy.Delete_management(flowlinesGCS)							
arcpy.Delete_management(routes)
#arcpy.Delete_management(points)
arcpy.Delete_management(pointsGCS)
arcpy.Delete_management("pointsGCSTable")
arcpy.Delete_management("occupancyTable")
arcpy.Delete_management("flowlinesGCSTable")
arcpy.Delete_management(occupancySites)

