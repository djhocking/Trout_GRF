rm(list = ls())

# ==============
# Load libraries
# ==============
library(foreign)
library(dplyr)
library(reshape2)


# ==============
# Specify Inputs
# ==============
# Base Directory
baseDirectory <- "C:/KPONEIL/GitHub/projects/Trout_GRF/Data/gisFiles/modelVersions/SusquehannaWest_Threepass"

# The FEATUREID of the most downstream segment
mouthID <- 922301

# Catchment Delineation File
delineationFilePath <- "F:/KPONEIL/SourceData/streamStructure/northeastHRD/NortheastHRD_delineatedCatchments.RData"


# =========
# Read data
# =========
# Read in the flowlines table
flowlines <- read.dbf (file.path(baseDirectory, "flowlinesTable.dbf"), as.is = T)

# Read in the sites table with network locations
sites <- read.dbf (file.path(baseDirectory, "occupancySites.dbf"), as.is = T)

# Read in the delineated catchments file
load(delineationFilePath)


# ============================
# Eliminate non-network values
# ============================
flow1 <- nrow(flowlines)
site1 <- nrow(sites)

network <- delineatedCatchments[ names(delineatedCatchments) == mouthID ][[1]]

flowlines <- filter(flowlines, FEATUREID %in% network)
sites <- filter(sites, FEATUREID %in% network)%>%
            mutate(fmeasKM = FMEAS)

flow2 <- nrow(flowlines)
site2 <- nrow(sites)

print(paste0(flow2 - flow1, " segments and ", site2 - site1, " sites removed."))


# ================
# Generate vectors
# ================

# Confluence nodes
confluences <- flowlines[ ,c('nodeID', 'NodeLon', 'NodeLat')]

confluences <- cbind(confluences, colsplit(confluences$nodeID, "_", c("type", "ID")))

# Site nodes
points <- sites[ ,c('nodeID', 'NodeLon', 'NodeLat')]

points$type <- "P"
points$ID <- points$nodeID

# Join
family <- rbind(confluences, points)
colnames(family)[1] <- 'child_b'

family[,c('parent_raw', 'parent_b', 'dist_b')] <- NA

# Check for unique child_b values
if(length(unique(family$child_b)) != nrow(family)){stop("Duplicate values for child_b exist.")}


# =====================================
# Evaluate network and populate vectors
# =====================================
for( i in 1:nrow(family) ){
  
  print(i)
  
  # ID of current node
  childID <- family$child_b[i]
  
  
  # First node in network
  # ---------------------
  if(childID ==  paste0("N_", mouthID)){
    
    family$parent_raw[i] <- NA
    family$parent_b[i]   <- NA
    family$dist_b[i]     <- NA
  }

  # Node is a confluence
  # --------------------
  if( family$type[i] == "N" & childID != paste0("N_", mouthID) ){
    
    # FEATUREID
    fid <- family$ID[i]
    
    # Segment
    segments <- filter(flowlines, FEATUREID == fid) 
    
    # Sample points on segment
    segPoints <- filter(sites, FEATUREID == fid)

    if( nrow(segPoints) == 0){
      
      # Case 1: child_b is a confluence and there are no sample points on the downstream segment
      #   The parent node is the nodeID of the downstream confluence, which is named by the next downstream segment (NextDownID)
      #   The distance is the length of the segment
      
      parentNode <- paste0("N_", segments$NextDownID)
      
      family$parent_raw[i] <- parentNode
      family$parent_b[i]   <- which(family$child_b == parentNode)
      family$dist_b[i]     <- segments$LengthKM
      
      rm(parentNode)
      
    }else{
      
      # Case 2: child_b is a confluence and there are 1 or more sample points on the downstream segment
      #   The parent node is the nodeID of the closest sample point, which is the sample point on the segment with the minimum FMEAS value.
      #   The distance is the length between the point and the upstream confluence(FMEAS)
      
      parentNode <- filter(segPoints, fmeasKM == min(fmeasKM))%>%
                          select(nodeID, fmeasKM)

      family$parent_raw[i] <- parentNode$nodeID
      family$parent_b[i]   <- which(family$child_b == parentNode$nodeID)
      family$dist_b[i]     <- parentNode$fmeasKM
      
      rm(parentNode)
    }
    
    rm(list = c("fid", "segPoints", "segments"))
  }
  
  # Node is a sample point
  # ----------------------
  if(family$type[i] == "P"){  
  
    #FEATUREID
    fid <- sites$FEATUREID[sites$nodeID == childID]
  
    # Distance to child node from upstream confluence
    childFmeasKM <- sites$fmeasKM[sites$nodeID == childID]

    #Segment
    segments <- filter(flowlines, FEATUREID == fid)   
    
    # Next downstream point on segment
    parentPoint <- filter(sites, FEATUREID == fid & fmeasKM > childFmeasKM )%>%
                          filter(fmeasKM == min(fmeasKM))
    
    if( nrow(parentPoint) == 0){
      
      #Case 3: child_b is a sample point with no downstream sample points on the same segment
      #   The parent node is the nodeID of the downstream confluence, which is named by the next downstream segment (NextDownID)
      #   The distance is the difference between the full length of the current segment and the distance to the point from upstream confluence(FMEAS)
      
      parentNode <- paste0("N_", segments$NextDownID)
      
      family$parent_raw[i] <- parentNode
      family$parent_b[i]   <- which(family$child_b == parentNode)
      family$dist_b[i]     <- segments$LengthKM - childFmeasKM
      
      rm(parentNode)
      
    }else{
      
      # Case 4: child_b is a sample point with at least 1 downstream sample point on the same segment
      #   The parent node is the nodeID of the closest sample point, which is the next downstream sample point on the segment (determined by the minimum FMEAS value of the downstream points)
      #   The distance is the difference between the FMEAS values of the 2 points
      
      family$parent_raw[i] <- parentPoint$nodeID
      family$parent_b[i]   <- which(family$child_b == parentPoint$nodeID)
      family$dist_b[i]     <- parentPoint$fmeasKM - childFmeasKM
    }

    rm(list = c("fid", "childFmeasKM", "segments", "parentPoint"))
    
  }
}


# ===================
# Prep output objects
# ===================

# Sample point vector
b_p <- which(family$type == "P")

# child_b, parent_b, and dist_b vectors
family <- select(family, c(child_b, parent_b, dist_b, NodeLat, NodeLon))

#save(family, b_p, file = file.path(baseDirectory, "sample.RData"))


# Save family
write.csv(family, 
          file = file.path(baseDirectory, "family.csv"), 
          row.names = F)



# ===============
# Manual checking
# ===============
# Check the parent node of a specific node:

#nom <- "N_889464"
#nom <- family$child_b[family$parent_b[which(family$child_b == nom)]]; 
#print(nom)

