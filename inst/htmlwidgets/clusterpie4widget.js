HTMLWidgets.widget({

  name: 'clusterpie4widget',

  type: 'output',

  initialize: function(el,width, height) {
    var mapdiv = document.createElement("div");
    mapdiv.id = "map";
    el.appendChild(mapdiv);

    return {
      // TODO: add instance fields as required
    };

  },

  renderValue: function(el, x,  instance) {


 "use strict";
//    iconField = '5065', //This is the fieldame for marker icon
  var geojson,
    metadata,
    geojsonPath,
    categoryField = 'allele', //This is the fieldname for marker category (used in the pie and legend)
    popupFields = ['popname', 'sampid','allele'], //Popup will display these fields
    tileServer = 'http://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png',
    tileAttribution =  'Map data: <a href="http://openstreetmap.org">OSM</a>',
    rmax = 32, //Maximum radius for cluster pies
    markerclusters = L.markerClusterGroup({
      maxClusterRadius: 1.75*rmax,
      iconCreateFunction: defineClusterIcon, //this is where the magic happens
    showCoverageOnHover: false,
		zoomToBoundsOnClick: true
    }),
    map = L.map('map').setView([40, 120], 2);

  //Add basemap
 // L.tileLayer(tileServer, {attribution: tileAttribution,  maxZoom: 15, minZoom: 2}).addTo(map);
  L.tileLayer('http://mt{s}.google.com/vt/v=w2.106&x={x}&y={y}&z={z}&s=', {subdomains:'0123',attribution:'&copy;Google2012',maxZoom: 10, minZoom: 2}).addTo(map);
 // L.tileLayer('http://c.tiles.mapbox.com/v3/examples.map-szwdot65/{Z}/{X}/{Y}.png').addTo(map);


  //and the empty markercluster layer
  map.addLayer(markerclusters);

  //Ready to go, load the geojson

  if (x.server != null) {
    geojsonPath = 'http://' + x.server + '/' + x.geojsonfile;

  d3.json(geojsonPath, function(error, data) {
      if (!error){
          geojson = data;
          metadata = data.properties;
          var markers = L.geoJson(geojson, {
  		pointToLayer: defineFeature,
			onEachFeature: defineFeaturePopup
          });
          markerclusters.addLayer(markers);
          map.fitBounds(markers.getBounds());
          map.attributionControl.addAttribution(metadata.attribution);
          renderLegend();
      } else {
	  console.log('Could not load data...');
      }
  });
} else {
var data = JSON.parse(x.data);
geojson = data;

metadata = data.properties;
var markers = L.geoJson(geojson, {
    pointToLayer: defineFeature,
		onEachFeature: defineFeaturePopup
    });
markerclusters.addLayer(markers);
map.fitBounds(markers.getBounds());
map.attributionControl.addAttribution(metadata.attribution);
renderLegend();
}

function defineFeature(feature, latlng) {
  	var categoryVal = feature.properties[categoryField],
     	name = feature.properties['name'];

    //iconVal = feature.properties[iconField];
    var myClass = 'marker icon-'+metadata.fields[categoryField].lookup[categoryVal]+' allele-'+metadata.fields[categoryField].lookup[categoryVal];
    var myIcon = L.divIcon({
        className: myClass,
        iconSize:null
    });
    return L.marker(latlng, {icon: myIcon});
}

function defineFeaturePopup(feature, layer) {
  var props = feature.properties,
    fields = metadata.fields,
    popupContent = '';

  popupFields.map( function(key) {
    if (props[key]) {
      var val = props[key],
        label = fields[key].name;
      if (fields[key].lookup) {
        val = fields[key].lookup[val];
      }
      popupContent += '<span class="attribute"><span class="label">'+label+':</span> '+val+'</span>';
    }
  });
  popupContent = '<div class="map-popup">'+popupContent+'</div>';
  layer.bindPopup(popupContent,{offset: L.point(1,-2)});
}

function defineClusterIcon(cluster) {

    var children = cluster.getAllChildMarkers(),
        n = children.length, //Get number of markers in cluster
        strokeWidth = 1, //Set clusterpie stroke width
        r = rmax-2*strokeWidth-(n<10?12:n<100?8:n<1000?4:0), //Calculate clusterpie radius...
        iconDim = (r+strokeWidth)*2, //...and divIcon dimensions (leaflet really want to know the size)
        data = d3.nest() //Build a dataset for the pie chart
          .key(function(d) { return d.feature.properties[categoryField]; })
          .entries(children, d3.map),
        //bake some svg markup
        html = bakeThePie({data: data,
                            valueFunc: function(d){return d.values.length;},
                            strokeWidth: 1,
                            outerRadius: r,
                            innerRadius: r-10,
                            pieClass: 'cluster-pie',
                            pieLabel: n,
                            pieLabelClass: 'marker-cluster-pie-label',
                            pathClassFunc: function(d){return "allele-"+metadata.fields[categoryField].lookup[d.data.key];},
                            pathTitleFunc: function(d){
                              var poplist = d3.nest()
          .key(function(d) { return d.feature.properties["popname"];})
          .entries(d.data.values);
                              var pops = "";
                              for (var k = 0; k < poplist.length; ++k) {
                                var p = poplist[k];
                                pops = pops + p["key"] + "\n";
                              }
                              return   pops + "\n" + metadata.fields[categoryField].lookup[d.data.key]+' ('+d.data.values.length+' chrom'+(d.data.values.length!=1?'s':'')+')';}
                          }),
        //Create a new divIcon and assign the svg markup to the html property
        myIcon = new L.DivIcon({
            html: html,
            className: 'marker-cluster',
            iconSize: new L.Point(iconDim, iconDim)
        });



        for (var item in data.data) {
        console.log(item+":"+data.data[item]);
        }
    return myIcon;
}

/*function that generates a svg markup for the pie chart*/
function bakeThePie(options) {
    /*data and valueFunc are required*/
    if (!options.data || !options.valueFunc) {
        return '';
    }
        // rInner = options.innerRadius?options.innerRadius:r-10, //Default inner radius = r-10
  var data = options.data,
        valueFunc = options.valueFunc,
        r = options.outerRadius?options.outerRadius:28, //Default outer radius = 28px
		rInner=0,
       strokeWidth = options.strokeWidth?options.strokeWidth:1, //Default stroke is 1
        pathClassFunc = options.pathClassFunc?options.pathClassFunc:function(){return '';}, //Class for each path
        pathTitleFunc = options.pathTitleFunc?options.pathTitleFunc:function(){return '';}, //Title for each path
        pieClass = options.pieClass?options.pieClass:'marker-cluster-pie', //Class for the whole pie
        pieLabel = options.pieLabel?options.pieLabel:d3.sum(data,valueFunc), //Label for the whole pie
        pieLabelClass = options.pieLabelClass?options.pieLabelClass:'marker-cluster-pie-label',//Class for the pie label


        origo = (r+strokeWidth), //Center coordinate
        w = origo*2, //width and height of the svg element
        h = w,
        donut = d3.layout.pie(),
        arc = d3.svg.arc().innerRadius(rInner).outerRadius(r);

    //Create an svg element
    var svg = document.createElementNS(d3.ns.prefix.svg, 'svg');
    //Create the pie chart
    var vis = d3.select(svg)
        .data([data])
        .attr('class', pieClass)
        .attr('width', w)
        .attr('height', h);

    var arcs = vis.selectAll('g.arc')
        .data(donut.value(valueFunc))
        .enter().append('svg:g')
        .attr('class', 'arc')
        .attr('transform', 'translate(' + origo + ',' + origo + ')');

    arcs.append('svg:path')
        .attr('class', pathClassFunc)
        .attr('stroke-width', strokeWidth)
        .attr('d', arc)
        .append('svg:title')
          .text(pathTitleFunc);

    vis.append('text')
        .attr('x',origo)
        .attr('y',origo)
        .attr('class', pieLabelClass)
        .attr('text-anchor', 'middle')
        //.attr('dominant-baseline', 'central')
        //*IE doesnt seem to support dominant-baseline, but setting dy to .3em does the trick*/
        .attr('dy','.3em')
        .text(pieLabel);
    //Return the svg-markup rather than the actual element
    return serializeXmlNode(svg);
}

/*Function for generating a legend with the same categories as in the clusterPie*/
function renderLegend() {
    var data = d3.entries(metadata.fields[categoryField].lookup),
      legenddiv = d3.select('body').append('div')
        .attr('id','legend');

    var heading = legenddiv.append('div')
        .classed('legendheading', true)
        .text(metadata.fields['description']);

    var legenditems = legenddiv.selectAll('.legenditem')
        .data(data);

    legenditems
        .enter()
        .append('div')
        .attr('class',function(d){return 'allele-'+d.value;})
        .classed({'legenditem': true})
        .text(function(d){return (d.value+": "+d.key);});
}

/*Helper function*/
function serializeXmlNode(xmlNode) {
    if (typeof window.XMLSerializer != "undefined") {
        return (new window.XMLSerializer()).serializeToString(xmlNode);
    } else if (typeof xmlNode.xml != "undefined") {
        return xmlNode.xml;
    }
    return "";
}


  },

  resize: function(el, width, height, instance) {

  }

});
