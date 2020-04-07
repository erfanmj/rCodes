library(igraph)
library(dplyr)
library(eulerian)
library(BiocGenerics)
library(parallel)
library(base)
library(stats)
library(readxl)

# 1- Importing Original data ;
Combined <- read_excel("D:/1- IUST/3- Article/3- Data/sampleData.xlsx", 
                       col_types = c("numeric", "numeric", "numeric"))


#2- Making the original data as an directed graph
Graph <- graph.data.frame(Combined, directed = TRUE)

#3- Set vertex namse equal to labels for better understanding visually
V(Graph)$name <- as.numeric(V(Graph)$name)
V(Graph)$name <- sort(V(Graph)$name)
V(Graph)$label <- V(Graph)$name


#4- Drawing the graph for better working
tkplot(Graph, canvas.width = 700,
       canvas.height = 700,
       vertex.color = "yellow", edge.label = E(Graph)$weight)


#5- Calculating Edge Betweenness then sorting by edgeBetweenness
GraphWithEB <- transform(
  edgeID = 1:gsize(Graph),
  Combined,
  edgeBetweenness = edge_betweenness(Graph, directed = TRUE))

GraphWithEB <- GraphWithEB[order(GraphWithEB$edgeBetweenness,decreasing = TRUE),]

#6- LOCAL COMMUNITY DETECTION:
#FUNCTIONS: 
# This function finds the common neighbors of some selected vertices
# Arguments: v: the set of vertices g: is your graph
# first: all neighbors are selected by their names 
# if a neighbor is occured more than once and wasn't in the original vertice 
# set is considered as the common neighbor of the set of vertices 
# function returns vertice IDs if it finds any, otherwise null 
commonNeighbors <- function(v,g)
{
  allNeigh <- list()
  for (i in v)
  {
    allNeigh <- append(allNeigh,c(neighbors(as.undirected(g),paste(i,sep=""))$name))
  }
  allNeigh <- cbind(allNeigh)
  allNeigh <- table(as.numeric(allNeigh))
  allNeigh <- as.data.frame(allNeigh)
  colnames(allNeigh) <- c('vertexID','freq')
  allNeigh <- allNeigh %>% dplyr::filter (freq > 1 & !vertexID %in% v )
  allNeigh <- as.matrix(allNeigh)
  if(length(allNeigh > 0))
  {
    allNeigh <- as.data.frame(allNeigh)
    return(as.matrix(allNeigh$vertexID))
  }
  else
  {
    return(NULL)
  }
}

insideOfCommEdgeIDs <- function(graph,vertices)
{
  out <- matrix()
  if(length(vertices) > 1)
  {
    for(i in vertices)
    {
      for (j in vertices)
      {
        if (are_adjacent(graph,i,j))
        {
          out <- rbind(out,get.edge.ids(graph,c(i,j),error = FALSE,directed = TRUE))
        }
      }
    }
    out <-out[!is.na(out)]
    return(out)
  }
  return(NULL)
}

## This Function Returns the edgeIDs of a vertex which are out of the community
## adjMat should be GraphWithEB
## v is the vertice you want to get its edges
## e is the set of edges which are in the community (we expel them)
outOfCommEdgeIDs <- function(adjMat,v,e)
{
  output <- matrix(ncol = 5)
  colnames(output) <- c("from","to","weight","edgeID","edgeBetweenness")
  output2 <- matrix(ncol = 5)
  colnames(output2) <- c("from","to","weight","edgeID","edgeBetweenness")
  for (i in v)
  {
    output <- rbind(output,adjMat %>% dplyr::filter(from == i | to == i))
  }
  output2 <- output %>% dplyr::filter(!edgeID %in% e)
  return(output2$edgeID)
}

## This Function Returns summation of edge weights of a vertex which are out of the community
## adjMat should be GraphWithEB
## v is the vertex you want to get its edge weights
## e is the set of edges which are in the community (we expel them)
mExternal <- function(adjMat,v,e)
{
  output <- matrix(ncol = 5)
  colnames(output) <- c("from","to","weight","edgeID","edgeBetweenness")
  output2 <- matrix(ncol = 5)
  colnames(output2) <- c("from","to","weight","edgeID","edgeBetweenness")
  for(i in v)
  {
    output <- rbind(output,adjMat %>% dplyr::filter(from == i | to == i))
  }
  output2 <- output %>% dplyr::filter(!edgeID %in% e)
  return(sum(output2$weight,na.rm = TRUE))
}

mInternal <- function(graph,adjMat,vertexList)
{
  output <- matrix(ncol = 5)
  colnames(output) <- c("from","to","weight","edgeID","edgeBetweenness")
  edgeIDs <- matrix()
  edgeIDs <- insideOfCommEdgeIDs(graph,vertexList)
  for(i in edgeIDs)
  {
    output <- rbind(output,adjMat %>% dplyr::filter(edgeID == i))
  }
  return(sum(output$weight,na.rm = TRUE))
}

secondPartWFM <- function(subgraphnel,graph,vertices,adjMat) {
  out <- tryCatch(
    {
      eulerian::eulerian(subgraphnel)
      edges <- insideOfCommEdgeIDs(graph,vertices)
      filtered <- adjMat %>% dplyr::filter(edgeID %in% edges)
      return(sum(filtered$weight)) 
    },
    error=function(cond) {
      return(0)
    }
  )    
  return(out)
}

WFM <- function(graph,adjMat,vertices)
{
  inerEdges <- insideOfCommEdgeIDs(graph,vertices)
  outerEdges <- matrix()
  outerEdges <- outOfCommEdgeIDs(adjMat,vertices,inerEdges)
  mint <- mInternal(graph,adjMat,vertices)
  mext <- mExternal(adjMat,vertices,outerEdges)
  firstPartWFM <- mint/mext
  induced <- as_graphnel(induced_subgraph(graph,vertices))
  secondPartWFM <- secondPartWFM(induced,graph,vertices,adjMat)
  return(firstPartWFM + secondPartWFM)
}

localCommunities <- function(graph,adjMat)
{
  LC <- matrix()
  for (i in adjMat$edgeID)
  {
    LCV <- matrix()
    LCV <- ends(graph,i,names = TRUE)
    LCV <- t(LCV)
    if(length(LCV) == 0)
    {
      rm(LCV)
      next
    }
    CN <- matrix()
    insideCommEdges <- matrix()
    CN <- commonNeighbors(LCV,graph)
    if(is.null(CN))
    {
      rm(CN)
      rm(insideCommEdges)
      rm(LCV)
      next
    }
    insideCommEdges <- insideOfCommEdgeIDs(graph,LCV)
    mInternalLCV <- mInternal(graph,adjMat,LCV)
    for(i in CN)
    {
      if(WFM(graph,adjMat,rbind(LCV,i)) >= WFM(graph,adjMat,LCV))
      {
        LCV <- rbind(LCV,i)
      }
    }
    LCV <- cbind(LCV)
    LCV <- table(LCV)
    LCV <- as.data.frame(LCV)
    LC <- cbind(LC,LCV)
    rm(LCV)
  }
  return(LC)
}














