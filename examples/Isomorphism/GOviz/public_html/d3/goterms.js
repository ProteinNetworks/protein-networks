var width = document.getElementById('networks').clientWidth,
    height = document.getElementById('networks').clientHeight;
var width2 = document.getElementById('networkfocus').clientWidth,
    height2 = document.getElementById('networkfocus').clientHeight;



var radius = 5;
var color = d3.scale.category20();

var force = d3.layout.force()
              .charge(-100)
              .linkDistance(20)
              .gravity(1)
              .size([width,height]);

var svg = d3.select("#networks").append("svg")
    .attr("width", width)
    .attr("height", height)
    .style("margin", "auto")
    .style("display", "block");

var svg2 = d3.select("#networkfocus").append("svg")
    .attr("width", width2)
    .attr("height", height2)
    .style("margin", "auto")
    .style("display", "block");

var globalgraph; // for use in generating the popups

d3.json("d3/GOtermsD3.json", function(error, graph) {
  if (error) throw error;

  globalgraph = graph;

  force
      .nodes(graph.nodes)
      .links(graph.links)
      .start();

  var link = svg.selectAll(".link")
      .data(graph.links)
    .enter().append("line")
      .attr("class", "link")
      .style("stroke-width", 1);

  var node = svg.selectAll(".node")
      .data(graph.nodes)
    .enter().append("circle")
    .attr("class", "node")
    .on("dblclick", click);

  function click(d) {
    force.stop();
    var isoGroup = d.group;
    console.log("Isomorphism Class:", isoGroup);
    svg2.selectAll("*").remove();
    renderGraph(isoGroup);
//    node.style("fill", function(e) {
  //     if (e.group === isoGroup) {
    //    return color(d.group);
    //   }
   //   else {
   //     return color(d3.rgb("white"));
  //   }
//  });
  };
  


  node.attr("r", radius)
      .style("fill", function(d) { return color(d.group); });
      //.call(force.drag);
  resize();
  d3.select(window).on("resize", resize);

  function resize() {
    width = document.getElementById('networks').clientWidth,
    height = document.getElementById('networks').clientHeight;
    svg.attr("width", width).attr("height", height);
    force.size([width, height]).resume();
  }

  force.on("tick", function() {
    link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node.attr("cx", function(d) { return d.x; }) //return Math.max(radius, Math.min(width - radius, d.x));})
        .attr("cy", function(d) { return d.y; }); //return Math.max(radius, Math.min(height - radius, d.y));});
  });
});


function renderGraph(d) {
// Given an isomorphism group label, identify the nodes and links corresponding to that group.
// Then render.
  var subnodes = globalgraph.nodes.filter(function(entry) { return entry.group === d;});

  var sublinks = globalgraph.links.filter(function(entry) {
    return entry.source.group === d;      
  });        



  var force2 = d3.layout.force()
              .charge(-100)
              .linkDistance(20)
              .size([width2,height2]);

  force2
      .nodes(subnodes)
      .links(sublinks)
      .start();

  var link2 = svg2.selectAll(".link")
      .data(sublinks)
    .enter().append("line")
      .attr("class", "link")
      .style("stroke-width", 1);

  var node2 = svg2.selectAll(".node")
      .data(subnodes)
    .enter().append("circle")
    .attr("class", "node");

  node2.attr("r", radius)
      .style("fill", function(d) { return color(d.group); })
      .call(force2.drag);
  resize();
  d3.select(window).on("resize", resize);

  function resize() {
    width2 = document.getElementById('networkfocus').clientWidth,
    height2 = document.getElementById('networkfocus').clientHeight;
    svg2.attr("width", width2).attr("height", height2);
    force2.size([width2, height2]).resume();
  }

  force2.on("tick", function() {
    link2.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });

    node2.attr("cx", function(d) { return d.x; })
        .attr("cy", function(d) { return d.y; });
  });
};
