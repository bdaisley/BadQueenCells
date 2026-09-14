library(dplyr)
library(readr)
library(stringr)
library(igraph)
library(ggraph)
library(tidyr)
library(ggraph)
library(scatterpie)
library(scales)
library(V8)
library(jsonlite)

goeBURST_layout <- function(
    graph,
    viva_file = "vivagraph_julbra.js",
    distance_attribute = "distance",
    max_steps = 100000L,
    final_drag = 0.015
) {
  if (!igraph::is_igraph(graph)) {
    stop("graph must be an igraph object.")
  }
  
  if (!file.exists(viva_file)) {
    stop("VivaGraphJS file not found: ", viva_file)
  }
  
  if (is.null(igraph::V(graph)$name)) {
    igraph::V(graph)$name <- as.character(
      seq_len(igraph::vcount(graph))
    )
  }
  
  node.names <- as.character(
    igraph::V(graph)$name
  )
  
  edges <- igraph::as_data_frame(
    graph,
    what = "edges"
  )
  
  if (!distance_attribute %in% names(edges)) {
    stop(
      "The edge attribute '",
      distance_attribute,
      "' was not found."
    )
  }
  
  edges$distance <- as.numeric(
    edges[[distance_attribute]]
  )
  
  if (
    anyNA(edges$distance) ||
    any(!is.finite(edges$distance)) ||
    any(edges$distance <= 0)
  ) {
    stop("All edge distances must be positive finite numbers.")
  }
  
  graph.input <- list(
    nodes = data.frame(
      id = node.names,
      stringsAsFactors = FALSE
    ),
    links = data.frame(
      source = as.character(edges$from),
      target = as.character(edges$to),
      distance = edges$distance,
      stringsAsFactors = FALSE
    )
  )
  
  graph.json <- jsonlite::toJSON(
    graph.input,
    dataframe = "rows",
    auto_unbox = TRUE,
    digits = NA
  )
  
  context <- V8::v8()
  
  # Make the browser-oriented library available in V8
  context$eval(
    "var window = this; var global = this; var self = this;"
  )
  
  context$source(viva_file)
  
  context$eval(
    '
function runPhyloVizLayout(graphJSON, maximumSteps, finalDrag) {
    var input = JSON.parse(graphJSON);

    var graph = Viva.Graph.graph();

    /*
     * Preserve input order exactly, since PHYLOViZ inserts nodes
     * and links in the order supplied.
     */
    for (var i = 0; i < input.nodes.length; i++) {
        graph.addNode(
            input.nodes[i].id,
            input.nodes[i]
        );
    }

    for (var j = 0; j < input.links.length; j++) {
        var edge = input.links[j];

        graph.addLink(
            edge.source,
            edge.target,
            {
                connectionStrength: Number(edge.distance),
                value: Number(edge.distance)
            }
        );
    }

    var idealSpringLength = 7;

    var layout = Viva.Graph.Layout.forceDirected(
        graph,
        {
            springLength: idealSpringLength,
            springCoeff: 0.0001,
            dragCoeff: 0.001,
            gravity: -10,
            theta: 0.8,

            /*
             * Exact PHYLOViZ spring transformation:
             * visual spring length is proportional to allelic distance.
             */
            springTransform: function(link, spring) {
                spring.length =
                    idealSpringLength *
                    link.data.connectionStrength;
            }
        }
    );

    /*
     * PHYLOViZ selects the first node encountered with the
     * highest degree and pins it at the origin.
     */
    var maximumDegree = 0;
    var topNode = null;

    graph.forEachNode(function(node) {
        var degree = node.links ? node.links.length : 0;

        if (maximumDegree < degree) {
            maximumDegree = degree;
            topNode = node;
        }
    });

    /*
     * Exact initial placements used by PHYLOViZ.
     */
    if (input.nodes.length > 0) {
        layout.setNodePosition(
            input.nodes[0].id,
            20,
            -20
        );
    }

    if (topNode !== null) {
        layout.setNodePosition(
            topNode.id,
            0,
            0
        );

        layout.pinNode(
            topNode,
            true
        );
    }

    /*
     * PHYLOViZ performs one initial layout step per node.
     */
    for (var step = 0; step < input.nodes.length; step++) {
        layout.step();
    }

    /*
     * The default PHYLOViZ interface then sets its drag slider
     * to 15, corresponding to a final coefficient of 0.015.
     */
    layout.simulator.dragCoeff(
        Number(finalDrag)
    );

    /*
     * The PHYLOViZ renderer continues until layout.step()
     * reports that the simulation is stable.
     */
    var stepsTaken = 0;
    var stable = false;

    while (
        !stable &&
        stepsTaken < Number(maximumSteps)
    ) {
        stable = layout.step();
        stepsTaken++;
    }

    var coordinates = [];

    graph.forEachNode(function(node) {
        var position =
            layout.getNodePosition(node.id);

        coordinates.push({
            name: String(node.id),
            x: position.x,
            y: position.y
        });
    });

    return JSON.stringify({
        nodes: coordinates,
        topNode:
            topNode === null
                ? null
                : String(topNode.id),
        maximumDegree: maximumDegree,
        stepsTaken: stepsTaken,
        stable: stable
    });
}
'
  )
  
  output.json <- context$call(
    "runPhyloVizLayout",
    graph.json,
    as.integer(max_steps),
    as.numeric(final_drag)
  )
  
  output <- jsonlite::fromJSON(
    output.json,
    simplifyDataFrame = TRUE
  )
  
  coordinates <- output$nodes[
    match(node.names, output$nodes$name),
    ,
    drop = FALSE
  ]
  
  if (anyNA(coordinates$name)) {
    stop("Not all graph nodes received coordinates.")
  }
  
  rownames(coordinates) <- coordinates$name
  
  list(
    coordinates = coordinates,
    top_node = output$topNode,
    maximum_degree = output$maximumDegree,
    steps_taken = output$stepsTaken,
    stable = output$stable
  )
}
