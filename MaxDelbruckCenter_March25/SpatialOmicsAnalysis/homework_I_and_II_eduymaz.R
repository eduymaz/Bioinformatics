##############################################################
BiocManager::install("rhdf5")
BiocManager::install("RBioFormats") 
install.packages("Seurat")
devtools::install_github("dmcable/spacexr")
BiocManager::install("ComplexHeatmap")
BiocManager::install("HDF5Array")
devtools::install_github("BIMSBbioinfo/HDF5DataFrame")
devtools::install_github("BIMSBbioinfo/ImageArray")

if (!require("devtools", quietly = TRUE))
  install.packages("devtools")
devtools::install_github("BIMSBbioinfo/VoltRon")

#devtools::install_github("bnprks/BPCells/r@v0.3.0")
#remotes::install_github("bnprks/BPCells/r")
install.packages('BPCells', repos = c('https://bnprks.r-universe.dev', 'https://cloud.r-project.org'))

devtools::install_github("vitessce/vitessceR")
BiocManager::install("basilisk")
install.packages("ggnewscale")
devtools::install_github('immunogenomics/presto')

##############################################################
library(rhdf5)
library(RBioFormats)
library(spacexr)
library(Seurat)
library(ComplexHeatmap)
library(HDF5DataFrame)
library(ImageArray)
library(VoltRon)
library(BPCells)
library(vitessceR)
library(basilisk)
library(ggnewscale)
library(presto)


##############################################################

setwd("/Users/eduymaz/VScode/COMPGEN25/module2")

## module2_hw2_ElifDuymaz.R
options(java.parameters = "-Xmx8g")
library(VoltRon)
library(dplyr)

# Step 1: Upload images
Xen_R1_image <- importImageData(
  "workshop/data/BreastCancer/Xenium_R1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_highres.tif",
  sample_name = "XeniumHE",
  channel_names = "H&E", 
  tile.size = 100
)
Visium_image <- importImageData(
  "workshop/data/BreastCancer/Visium/spatial/tissue_hires_image.png",
  sample_name = "VisiumHE",
  channel_names = "H&E", 
  tile.size = 100
)

# Step 2: Visualize the images
vrSpatialPlot(Xen_R1_image, channel = "H&E")
vrSpatialPlot(Visium_image, channel = "H&E")

# Step 3: process the image using magick (negate and rotate)
img_path <- "workshop/data/BreastCancer/Visium/spatial/tissue_hires_image.png"
img <- image_read(img_path)

#with negative
img_neg <- image_negate(img)

# Rotate 90 degrees
img_rot <- image_rotate(img_neg, 90)

# Save processed image
image_write(img_rot, "processed_visium_image.png")

# Transfer processed image to VoltRon
Visium_image <- importImageData(
  "processed_visium_image.png",
  sample_name = "VisiumHE_processed",
  channel_names = "H&E", 
  tile.size = 100
)

# Check processed image
vrSpatialPlot(Visium_image, channel = "H&E")


img_processed <- image_read("processed_visium_image.png")

# Scaling
img_scaled <- image_scale(img_processed, "521.28x576!")

# Save and transfer to VoltRon
image_write(img_scaled, "scaled_visium_image.png")
Visium_image <- importImageData(
  "scaled_visium_image.png",
  sample_name = "VisiumHE_scaled",
  channel_names = "H&E", 
  tile.size = 100
)

# Scale Xenium image
img_xen <- image_read("workshop/data/BreastCancer/Xenium_R1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_highres.tif")
img_xen_scaled <- image_scale(img_xen, "796.86x579.96!")
image_write(img_xen_scaled, "scaled_xenium_image.tif")
Xen_R1_image <- importImageData(
  "scaled_xenium_image.tif",
  sample_name = "XeniumHE_scaled",
  channel_names = "H&E", 
  tile.size = 100
)

# Step 4: Use magick for scaling

# Scale Visium image
img_processed <- image_read("processed_visium_image.png")
img_scaled <- image_scale(img_processed, "521.28x576!")
image_write(img_scaled, "scaled_visium_image.png")

Visium_image <- importImageData(
  "scaled_visium_image.png",
  sample_name = "VisiumHE_scaled",
  channel_names = "H&E", 
  tile.size = 100
)

# Scale Xenium image
img_xen <- image_read("workshop/data/BreastCancer/Xenium_R1/Xenium_FFPE_Human_Breast_Cancer_Rep1_he_image_highres.tif")
img_xen_scaled <- image_scale(img_xen, "796.86x579.96!")
image_write(img_xen_scaled, "scaled_xenium_image.tif")

Xen_R1_image <- importImageData(
  "scaled_xenium_image.tif",
  sample_name = "XeniumHE_scaled",
  channel_names = "H&E", 
  tile.size = 100
)


# Step 5: Shiny
registration_result <- registerSpatialData(
  object_list = list(Visium_image, Xen_R1_image)
)



#

# module2_hw1_ElifDuymaz.R

Visium_data <- importVisium("workshop/data/BreastCancer/Visium/",
                            sample_name = "Visium_BC")


#Xenium_data <- loadVoltRon("workshop/data/ondisk/Xen_R1/")

#install.packages("imager")
#library(EBImage)

head(vrFeatures(Visium_data))
length(vrFeatures(Visium_data))

# normalize and select features
Visium_data <- normalizeData(Visium_data)
Visium_data <- getFeatures(Visium_data, n = 3000)

# selected features
head(vrFeatureData(Visium_data))
selected_features <- getVariableFeatures(Visium_data)
head(selected_features, 20)

### Dimensional Reduction ####
####

# embedding
Visium_data <- getPCA(Visium_data, features = selected_features, dims = 30)
Visium_data <- getUMAP(Visium_data, dims = 1:30)
vrEmbeddingNames(Visium_data)

# embedding visualization
vrEmbeddingPlot(Visium_data, embedding = "umap")


####
### Clustering ####
####

# graph for neighbors
Visium_data <- getProfileNeighbors(Visium_data, dims = 1:30, k = 10, method = "SNN")
vrGraphNames(Visium_data)

# clustering
Visium_data <- getClusters(Visium_data, resolution = 0.5, label = "Clusters", graph = "SNN")

####
### Visualization ####
####

# embedding
vrEmbeddingPlot(Visium_data, embedding = "umap", group.by = "Clusters")


# new version hw 1 : 

Visium_data <- importVisium("workshop/data/BreastCancer/Visium/",
                            sample_name = "Visium_BC")


# Checking the content of the data
head(vrFeatures(Visium_data))
length(vrFeatures(Visium_data))


### **Feature Selection & Normalization** ###

# Normalize data
Visium_data <- normalizeData(Visium_data)

# Selecting the 3000 features with the most variability
Visium_data <- getFeatures(Visium_data, n = 3000)

# Checking selected features
head(vrFeatureData(Visium_data))
selected_features <- getVariableFeatures(Visium_data)
head(selected_features, 20)


### **Dimensional Reduction ** ###

# PCA ve UMAP ile gömme işlemi (embedding)
Visium_data <- getPCA(Visium_data, features = selected_features, dims = 30)
Visium_data <- getUMAP(Visium_data, dims = 1:30)
vrEmbeddingNames(Visium_data)

# PCA ve UMAP visualization
vrEmbeddingPlot(Visium_data, embedding = "umap")

### **Clustering** ###

# Neighborhood graphing
Visium_data <- getProfileNeighbors(Visium_data, dims = 1:30, k = 10, method = "SNN")
vrGraphNames(Visium_data)

# Clustering process
Visium_data <- getClusters(Visium_data, resolution = 0.5, label = "Clusters", graph = "SNN")

# Visualizing the clustering result on UMAP
vrEmbeddingPlot(Visium_data, embedding = "umap", group.by = "Clusters")

### **Spatial Plot** ###

# Visium spatial plot
vrSpatialPlot(Visium_data, group.by = "Clusters")
vrSpatialPlot(Visium_data, group.by = "Clusters", pt.size = 0.18)

####

vrSpatialPlot(Visium_data, group.by = "Clusters", crop = TRUE, alpha = 1)

vrHeatmapPlot(Visium_data, features = vrFeatures(Visium_data), group.by = "Clusters",
              show_row_names = TRUE, show_heatmap_legend = TRUE)

# ComplexHeatmap
ht_opt$message = FALSE

# Heatmap
vrHeatmapPlot(Visium_data, features = vrFeatures(Visium_data), group.by = "Clusters",
              show_row_names = TRUE, show_heatmap_legend = TRUE)




