# Description

**_NetworkAnalystR_** is the underlying R package synchronized with NetworkAnalyst web server. It is designed for systems interpretation of gene list data through comprehensive collections of knowledge-based networks - protein-protein interactions, gene regulatory networks, co-expression networks, as well as chemical-, disease- or drug- associated networks. The R package is composed of R functions necessary for the web-server to perform network creation, trimming and analysis.

Following installation and loading of NetworkAnalystR, users will be able to reproduce web server results from their local computers using the R command history downloaded from NetworkAnalystR. Running the R functions will allow more flexibility and reproducibility.

Note - NetworkAnalystR is still under development - we cannot guarantee full functionality

# Installation

## 1. Install package dependencies

To use NetworkAnalystR, make sure your R version is >4.0.3 and install all package dependencies. Ensure that you are able to download packages from Bioconductor. To install package dependencies, use the pacman R package. Note that some of these packages may require additional library dependencies that need to be installed prior to their own successful installation.

```{r eval=FALSE}
install.packages("pacman")

library(pacman)

pacman::p_load(igraph, RColorBrewer, qs, rjson, RSQLite)
```

## 2. Install the package

NetworkAnalystR is freely available from GitHub. The package documentation, including the vignettes for each module and user manual is available within the downloaded R package file. If all package dependencies were installed, you will be able to install the NetworkAnalystR. 

Install the package directly from github using the _devtools_ package. Open R and enter:

```{r eval=FALSE}
# Step 1: Install devtools
install.packages('devtools')
library(devtools)

# Step 2: Install NetworkAnalystR WITHOUT documentation
devtools::install_github("xia-lab/NetworkAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual", "--no-build-vignettes"))

# Step 2: Install NetworkAnalystR WITH documentation
devtools::install_github("xia-lab/NetworkAnalystR", build = TRUE, build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
```

## Tips for using the NetworkAnalystR package

1. The first function that you will use in every module is the `Init.Data` function, which initiates R objects that stores user's data, parameters for further processing and analysis.
2. The NetworkAnalystR package will output data files/tables/analysis/networks outputs in your current working directory.
3. Every function must be executed in sequence as it is shown on the R Command history, please do not skip any commands as this can result in errors downstream.
4. Main functions in NetworkAnalystR are documented. Use the ?Function format to open its documentation. For instance, use ?NetworkAnalystR::ReadTabExpression to find out more about this function.
5. It is recommended to set the working folder to an empty folder because numerous files will be generated through the process.
6. R package can not be used for visual analytics as they are hosted on the website.
7. R package is derived from R scripts used for powering web server. The values returned may not be useful in the context of local usage. For gene list named “datalist1”. use this function to access: dataSet <- readDataset(fileName)

# Examples

## 1. Create interaction network from a list of genes

### 1.1 Load NetworkAnalystR library and initialize R objects
```{r eval=FALSE}
library(NetworkAnalystR)

#boolean FALSE indicates local mode (vs web mode);
dataSets <- Init.Data(FALSE);
```
#### 1.2 Read example gene list
```{r eval=FALSE}
genelist <- "#Entrez  logFC
4495  61.12
4496  51.06
4499  23.79
6354  21.04
6369  19.76
4494  16.24
4501  14.76
11026  14.04
199675  12.65
4316  12.04
771  8.19
6346  7.07
6367  6.97
5473  6.76
2357  5.71
5265  5.65
1462  5.27
2358  4.92
22918  4.58
56729  4.39
8710  4.25
3563  4.21
57214  4.13
3290  4.08
6363  4.08
27063  4.05"
dataSets <- MapListIds("datalist1", genelist, "hsa", "entrez");
```
Take a look at the mapped genelist by using readDataset("datalist1") function.
```{r eval=FALSE}
dataSet <- readDataset("datalist1") 
print(head(dataSet$prot.mat))
#      [,1]
#4495 61.12
#4496 51.06
#4499 23.79
#6354 21.04
#6369 19.76
#4494 16.24
```

### 1.3 Search Protein-protein interaction database
This step uses the proteins from uploaded the gene list as seeds to search for direct interacting partners in selected database. The same seeds can be used to query multiple databases, resulting in integrated networks.
```{r eval=FALSE}
analSet <- SearchNetDB("ppi", "hsa_innate", TRUE, 900, 1);
##Check the result 
##The resulting network is saved in as edge list analSet$ppi.net
View(analSet$ppi.net)
```
### 1.4 Convert edge list to igraph object.
The edge list can contain multiple subnetworks that are not interconnected. The individual subnetworks are stored in ```analSet$ppi.comps```. By default, the largest subnetwork “subnetwork1” is used for visualization in the web. You can access through ```analSet$ppi.comps[["subnetwork1"]]```.

```{r eval=FALSE}
analSet <- CreateGraph();
##Check the size of resulting subnetworks
print(analSet$net.stats);
#            Node Edge Query
#subnetwork1   72   76     9
#subnetwork2   14   13     2
#subnetwork3   11   10     1
#subnetwork4    5    4     1
```
### 1.5 Check network topology
This is to assess the degree and betweenness distribution of nodes within a specific subnetwork.
```{r eval=FALSE}
PlotDegreeHistogram("degree_distribution_0_" , "subnetwork1");
PlotBetweennessHistogram("betweenness_distribution_0_" , "subnetwork1");
```
Check your working directory for png images named ``degree_distribution_0_dpi72.png`` and ``betweenness_distribution_0_dpi72``, open them. <br />
![Degree distribution](https://dev.expressanalyst.ca/ExpressAnalyst/resources/images/RTutorial/degree_distribution_0_dpi72.png) <br />
![Betweenness distribution](https://dev.expressanalyst.ca/ExpressAnalyst/resources/images/RTutorial/betweenness_distribution_0_dpi72.png) <br />
You can see that there is a large number of nodes with one connection from the degree distribution and only few nodes with large betweenness. Overall, it looks like a small-world network.

### 1.6 Convert igraph object to json file for interactive network visualization
```{r eval=FALSE}
analSet <- PrepareNetwork("subnetwork1", "networkanalyst_0.json");
#The output is saved in "networkanalyst_0.json"
#it can be accessed as R object
View(analSet[[filenm]]);
#You can plot the overall graph (containing all subnetworks)
plot(analSet$overall.graph, layout=layout_with_fr, vertex.size=5, vertex.label=NA)
```

![Default network](https://dev.expressanalyst.ca/ExpressAnalyst/resources/images/RTutorial/defaultnet.png)

### 1.7 Perform PCSF on the network (trimming)
If the network is too large, you can also perform PCSF or minimum connected network. You need to convert igraph object to json file again after.
```{r eval=FALSE}
analSet <- ComputePCSFNet(); 
analSet <- PrepareNetwork("subnetwork1", "networkanalyst_0.json");
#You can plot the overall graph (containing all subnetworks);
plot(analSet$overall.graph, layout=layout_with_fr, vertex.size=5, vertex.label=NA)
```
![Default network](https://dev.expressanalyst.ca/ExpressAnalyst/resources/images/RTutorial/pcsfnet.png)
