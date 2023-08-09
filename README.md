# R-Tree Implementation Documentation

## Introduction

The code implements an R-Tree data structure, a balanced tree structure designed for indexing and querying multidimensional data. It optimizes search and nearest neighbor queries in multi-dimensional spaces by partitioning the space into smaller rectangles, grouping nearby points together.

<img src = "https://github.com/Shubham-me/R_TREES/assets/87982899/8571dad6-15cd-461c-8a52-7b9ca443a2cd" width = "700" height = "250">


## Node Class

The `node` class represents internal and leaf nodes of the R-Tree. It has attributes:

- `type`: Indicates internal (0) or leaf (1) node.
- `numberChildren`: Number of children in the node.
- `Rectangle`: Bounding rectangle of the node in D-dimensional space.
- `children`: Array of pointers to child nodes.

At the leaf node level its children are just points in D-dimension, while for the internal nodes the children are again Minimum Bounding Rectangles.
### Overlapping Check

The `overlap` function checks if two MBRs (Minimum Bounding Rectangles) overlap in D-dimensional space.
If the two MBRs are isolated in any dimension, i.e if $[l_{1i},u_{1i}] \space [l_{2i},u_{2i}]$ are exclusive then MBRs don't overlap. 

### Search Operation

The `search` function traverses the R-Tree to find nodes intersecting a query MBR. Found nodes are stored in `nodeList`.
The overlapping function is used to find the query children nodes that overlap with query MBR and are recursively called for search.

### Area and Area Increase Calculations

The `area` function calculates the MBR's area, while `areaIncrease` calculates area increase when adding a point or MBR.

### Split Operation

The `split` function splits a node into sub-MBRs to maintain balance. It redistributes children between nodes when the MBR reaches its maximum limit. It is called during insert operation.

### Insert Operation

The `insert` function inserts a point. For leaves, it inserts directly. For internals, it chooses child with minimum area increase. In case the MBR is filled to max limit the split function is called and node is splitted into two nodes and recursively the parent MBR max limit condition is handled. In case the root splits, a new root is made.

### Nearest Neighbor Search

The `Nearest_Neighbour` function finds nearest neighbors for a query point, pruning the search space efficiently. 
The metrics used for pruning are `MinDistance` and `MinmaxDistance` which are defined as below:

<img src = "https://github.com/Shubham-me/R_TREES/assets/87982899/bcf3cf71-9c86-4f0c-b3cf-55e967c0a376" width = "400" height = "200"><br>

`MinDistance` between MBR Q and any point P is the distance of closest point in MBR to P.  
<img src = "https://github.com/Shubham-me/R_TREES/assets/87982899/05a59846-c269-4859-8e25-7aade838ff99" width = "400" height = "200"><br>

The `MinMaxDistance` guarantees there is an object within the MBR at a distance less than or equal to this distance, it is defined as below.  
<img src = "https://github.com/Shubham-me/R_TREES/assets/87982899/1ef8a2e3-a26c-49de-bb3d-ccc381a3ff42" width = "400" height = "300"><br>

While Querying nerest neighbours of any point P, For MBR Q and R if MinDistance of Q is greater than MinMaxDistance of R then we can prune the search in MBR Q and we have a conform nearer point in R.

Algorithm : For the children MBRs we compute the MinMaxDistance for each and take the minimum of them.
The MBRs with MinDistance > min(MinMaxDistance) can be excluded from search.

## Main Functions

### Insertion of Points

`insertInTree` inserts points using `insertMain` while considering new root nodes.

### Generating Points

`generatePoint` generates a random D-dimensional point, while `generatePoints` creates multiple points.

### Nearest Neighbor Search and Comparison

The main function sets up an R-Tree, inserts and queries points using R-Tree and brute-force approaches, then compares results.

## Usage and Example

The code demonstrates R-Tree usage for multidimensional nearest neighbor search. It inserts points, generates query points, and finds neighbors using R-Tree and brute-force.

## Conclusion

The R-Tree optimizes spatial queries through its hierarchical structure. This implementation showcases insertion, search, and nearest neighbor queries, making it useful for spatial databases, GIS, and computer graphics.

## References
1) R-TREES. A DYNAMIC INDEX STRUCTURE FOR SPATIAL SEARCHING	- Antonin Guttman (University of California Berkeley)
2) Nearest Neighbour Queries - Nick Roussopoulos, Stephen Kelley, Frederic Vincent (University of Maryland)
