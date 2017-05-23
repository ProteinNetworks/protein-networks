var width = document.getElementById('networks').clientWidth,
    height = document.getElementById('networks').clientHeight;
var width2 = document.getElementById('networkfocus').clientWidth,
    height2 = document.getElementById('networkfocus').clientHeight;
var width3 = document.getElementById('networktable').clientWidth,
    height3 = document.getElementById('networktable').clientHeight;
var width4 = document.getElementById('networktable2').clientWidth,
    height4 = document.getElementById('networktable2').clientHeight;



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


var PDBtable = d3.select("#networktable").append("table")
    .attr("width", width3)
    .attr("height", height3)
    .attr("class", "table table-striped")
    .style("margin", "auto");
    //.style("display", "block");

var GOtable = d3.select("#networktable2").append("table")
    .attr("width", width4)
    .attr("height", "auto")
    .attr("class", "table table-striped")
    .style("margin", "auto");
    //.style("display", "block");

var globalgraph; // for use in generating the popups

d3.json("d3/GOtermsD3WithNames.json", function(error, graph) {
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
    PDBtable.selectAll("*").remove();
    GOtable.selectAll("*").remove();
        
    renderGraph(isoGroup);
    renderTable(isoGroup);
  };
  


  node.attr("r", function(d) {
      var numPDBs = graph.data.filter(function(entry) {return entry.group === d.group;})[0].PDBs.length;
      if (numPDBs < 10) {
        return 3;
      }
      else if (numPDBs < 100) {
        return 5;
      }
      else {
        return 7;
      }
      })
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

    node.attr("cx", function(d) { return d.x; })
        .attr("cy", function(d) { return d.y; });
  });
});


function renderTable(isoGroup) {
// Given an isomorphism group label, pull the corresponding "data" field and plot as a table.

   var tableData = globalgraph.data.filter(function(entry) {return entry.group === isoGroup;})[0];


   var  thead = PDBtable.append("thead"),
        tbody = PDBtable.append("tbody");

    // append the header row
    thead.append("tr")
        .append("th")
            .text("PDBs");

    // create a row for each object in the data
    var rows = tbody.selectAll("tr")
        .data(tableData.PDBs)
        .enter()
        .append("td")
        // This html adds a link with the PDB reference as text, and linking to an NGL Viewer instance.
        // The image is pulled from the RCSB and sized to give 5 pdbs per row.
        .html(function(d) {
                return "<a href=http://www.rcsb.org/pdb/ngl/ngl.do?pdbid=" + d + " target='_blank'\>" + d + "</a><br/>"
                         + "<img src=http://www.rcsb.org/pdb/images/" + d + "_asym_r_500.jpg width=" + width/6 + " />";});
   // Now sort the GO data
   var columnOrder = ["GO Label", "GO Name", "Corrected p value", "Unique"];
   var thead = GOtable.append("thead").selectAll("th")
                      .data(columnOrder) //d3.keys(tableData.GOterms[0]))
                      //.data(function(d) {
                     //   return columnOrder.map(function(m) { return d3.keys(tableData.GOterms[0]); });
                     // })
                      .enter().append("th").text(function(d){return d});
  console.log(d3.keys(tableData.GOterms[0]));
// fill the table

// create rows
var tr = GOtable.append("tbody").selectAll("tr")
               .data(tableData.GOterms).enter().append("tr")
// cells
var td = tr.selectAll("td")
  .data(function(d){
        return [d['GO Label'], d['GO Name'], d['Corrected p value'], d['Unique']] //columnOrder.map(function(m) {return d3.values(d)[m]})
  })
  .enter().append("td")
  .text(function(d) {return d})


};



function renderGraph(isoGroup) {
// Given an isomorphism group label, identify the nodes and links corresponding to that group.
// Then render.
  var subnodes = globalgraph.nodes.filter(function(entry) { return entry.group === isoGroup;});

  var sublinks = globalgraph.links.filter(function(entry) {
    return entry.source.group === isoGroup;      
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