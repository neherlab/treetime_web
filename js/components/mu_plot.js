var EventEmitter = require('events').EventEmitter;
var d3 = require('d3');
var d3tip = require('d3-tip');
var Globals = require('./globals.js')
var ANIMATION_DURATION = Globals.ANIMATION_DURATION;


var nodeTooltip = d3tip()
  .direction('n')
  .attr('class', 'd3-tip')
  .offset([-10, 0])
  .html(function(d){
    var string = "";
    string += "<h4>" + d.name +"</h4>";
    return string;
  });

var getLeafNodes = function(leafNodes, obj){
        if(obj.children){
            obj.children.forEach(function(child){getLeafNodes(leafNodes,child)});
            leafNodes.push(obj);
        } else{
            leafNodes.push(obj);
        }
};

var getInternalNodes = function(internalNodes, obj){

    if(obj.children){
        obj.children.forEach(function(child){getInternalNodes(internalNodes,child)});
        internalNodes.push(obj);
    }

};

var _internalStrokeColor = function(){
  return '#bc5a45';
}

var _terminalStrokeColor = function(){
  return '#034f84';
}

var _internalFillColor = function(){
  return "#f18973"
}

var _terminalFillColor = function(){
  return '#92a8d1';
}


var MuPlot = {};

MuPlot.mu = [];

MuPlot.old_state = {};

MuPlot.padding_bottom = 80;

MuPlot.padding_text = 20;

MuPlot.padding_top = 10;

MuPlot.padding_right = 10;

MuPlot.padding_left = 120;

MuPlot.regression = {};

MuPlot.create = function(el, props, state){
    //console.log("CREATING MU PLOT")
    this.width = el.offsetWidth
    this.height = el.offsetHeight
    console.log("MU height: " + this.height);
    this.svg = d3.select(el).append('svg')
      .attr('class', 'd3_mu')
      .attr('width',  this.width)
      .attr('height', this.height);

    this.svg.append('g')
      .attr('class', 'd3_mu_axis')

    this.svg.append('g')
      .attr('class', 'd3_mu_points')

    this.legend = this.svg.append("g")
      .attr("class","d3_mu_legend")
      .attr("transform","translate(50,30)")
      .style("font-size","12px")
      //.call(d3.legend)

    this.svg.call(nodeTooltip);

    var dispatcher = new EventEmitter();
        //this.update(el, props.root, state, dispatcher);
    return dispatcher;

};

MuPlot._set_points_from_root = function(dispatcher){

  //console.log("ROOT updates LH...");
  if (!this.tree) {return ;}
  var tip_lhs = [];
  getLeafNodes(tip_lhs, this.tree);
  this.points = tip_lhs
      .map(function(d){
        return ({
          name: d.name,
          x: d.numdate,
          y: d.xvalue,
          terminal: !d.children,
          fill_color: d.children ? _internalFillColor() : _terminalFillColor(),
          stroke_color: d.children ? _internalStrokeColor() : _terminalStrokeColor()
        });
      }).filter(function(d){
        return (
          (d.numdate != 0.0) && (d.xvalue != 0)
          );
      });

  this._update_lin_regression(dispatcher);

};

MuPlot._update_lin_regression = function(dispatcher){

    //console.log("Updateting linear regression for MU plot...");
    var terminals = this.points.filter(function(d){return d.terminal;});

    var n = terminals.length;

    if (n == 0) {
      this.regression = {};
      return;
    }

    console.log(this.points);
    var sum_x =  d3.sum(terminals.map(function(d){return d.x}));
    var sum_y =  d3.sum(terminals.map(function(d){return d.y}));
    var sum_xy = d3.sum(terminals.map(function(d){return d.x * d.y}));
    var sum_xx = d3.sum(terminals.map(function(d){return d.x * d.x}));
    var sum_yy = d3.sum(terminals.map(function(d){return d.y * d.y}));

    var slope = (n * sum_xy - sum_x * sum_y) / (n*sum_xx - sum_x * sum_x);
    var intercept =  (sum_y - slope * sum_x)/n;
    var r2 = Math.pow((n*sum_xy - sum_x*sum_y)/Math.sqrt((n*sum_xx-sum_x*sum_x)*(n*sum_yy-sum_y*sum_y)),2);

    this.regression = {
      'slope' :  slope,
      'intercept' :  intercept,
      'r2' :  r2
    };
    dispatcher.emit('mol_clock:regression_changed', this.regression);

};

MuPlot._draw_regression = function(el, scales) {

  if (!this.regression) return;

  var max_x = d3.max(scales.x.domain());
  var min_x = d3.min(scales.x.domain());

  var max_y = this.regression.slope * max_x + this.regression.intercept
  var min_y = this.regression.slope * min_x + this.regression.intercept

  var svg = d3.select(el).select('.d3_mu_points').append('line')
    .attr('class', 'd3_mu_regression')
    .attr('x1', function(d){return scales.x(min_x)})
    .attr('y1', function(d){return scales.y(min_y)})
    .attr('x2', function(d){return scales.x(max_x)})
    .attr('y2', function(d){return scales.y(max_y)})
    .style("stroke", "#4D92BF")
    .style("stroke-width", '2px')

};

MuPlot.update = function(el, root, state, dispatcher){

    //console.log("UPDATING MU");

    if (this.width != el.offsetWidth || this.height != el.offsetHeight){

        var g = d3.select(el).select('.d3_mu_axis').selectAll("*");
        g.remove();
        g = d3.select(el).select('.d3_mu_points').selectAll("*");
        g.remove();
        this.svg
          .attr('width',  el.offsetWidth)
          .attr('height', el.offsetHeight);
        var scales = this._scales(el);
        this._draw_axis(el, scales)
        this._draw_points(el, scales, dispatcher)
        this._draw_regression(el, scales)
        this.width = el.offsetWidth;
        this.height = el.offsetHeight;

    }

    if (this.tree != root){
      // update all points
      //console.log("MuPlot detected Tree changes, recreating the plot...")
      this.tree = root;
      this._set_points_from_root(dispatcher);
      this._update_lin_regression(dispatcher);
      var scales = this._scales(el);
      this._draw_points(el, scales, dispatcher)
      this._draw_regression(el, scales)
      this._draw_axis(el, scales)

    }

    if(this.old_state.selected_tip != state.selected_tip){
        //console.log("MU: tip selection changed.")
        var selected_tip = state.selected_tip;

        this.points.map(function(d){
          if (selected_tip && d.name == selected_tip){
              d.selected = true;
          }else{
              d.selected = false;
          }
        });

        this._refresh_selected_tip(el);
    }

    this.old_state = state;
    // selected node

};

MuPlot._refresh_selected_tip = function(el){

    var g = d3.select(el).selectAll('.d3_mu_points');
    var tip = g.selectAll('.d3_mu_point')
    tip.attr("r", this._tipRadius)

};

MuPlot._scales = function(el){


  var width = el.offsetWidth;
  var height = el.offsetHeight;
  var xs = this.points.map(function(d){return d.x});
  var ys = this.points.map(function(d){return d.y});
  var x = d3.scale.linear()
    .domain([d3.min(xs) , d3.max(xs)])
    .range([this.padding_left, width - this.padding_right]);

  var y = d3.scale.linear()
      .domain([d3.max(ys), d3.min(ys)])
      .range([this.padding_top,height-this.padding_bottom])
  return {x: x, y: y};

};

MuPlot._draw_text = function(el, scales){

  var text_x  = scales.x(0.1 * d3.max(scales.x.domain())  +  0.9 * d3.min(scales.x.domain()));
  var text_y  = scales.y(0.9 * d3.max(scales.y.domain()) - 0.1 * d3.min(scales.y.domain()));

  var g = d3.select(el).select('.d3_mu_axis').append("g")
    .attr("class", "d3_mu_legend")
    .attr("x", text_x)
    .attr("y", text_y)
    .style("text-anchor", "left")

  g.append("rect")
    .attr("x", text_x-10)
    .attr("y", text_y-10)
    .attr("width", 150)
    .attr("height", 85)
    .attr("fill", "white")
    .attr("stroke", "black")

  var c1 = g.append("circle")
    .attr("cx", text_x + 1)
    .attr("cy", text_y + 0)
    .attr("r", 6)
    .attr("width", 10)
    .attr("height", 10)
    .style("fill", _internalFillColor)
    .style("stroke", _internalStrokeColor)

  var t1 = g.append("text")
    .attr("x", text_x + 12)
    .attr("y", text_y + 0)
    .attr("dy", 6)
    .text("Internal nodes")

  g.append("circle")
    .attr("cx", text_x + 1)
    .attr("cy", text_y + 20)
    .attr("r", 6)
    .attr("width", 10)
    .attr("height", 10)
    .style("fill", _terminalFillColor)
    .style("stroke", _terminalStrokeColor)

  g.append("text")
    .attr("x", text_x + 12)
    .attr("y", text_y + 20)
    .attr("dy", 6)
    .text("Terminal nodes")

  var superscript = "⁰¹²³⁴⁵⁶⁷⁸⁹";
  var html = "<span> &mu; = " + this.regression.slope.toExponential(3) + " </br> R" + superscript[2] +" = " + Math.round(this.regression.r2 * 1000) / 1000 + "</span>"
  console.log(html)
  g.append('foreignObject')
    .attr("x", text_x + 12)
    .attr("y", text_y + 30)
    .attr('width', 100)
    .attr('height', 50)
    .append("xhtml:body")
    .html(html)

};

MuPlot._tipRadius = function(d){
    return d.selected ? 10.0 : 6.0;
};


MuPlot._draw_axis = function(el, scales){

    var width = el.offsetWidth;
    var height = el.offsetHeight;

    var xAxis = d3.svg.axis()
         .scale(scales.x)
         .orient("bottom")
         .ticks(5)
         .tickSize(-height + this.padding_top, 0, 0)
         .tickFormat(d3.format("d"))

    var yAxis = d3.svg.axis()
        .scale(scales.y)
        .orient("left")
        .ticks(10)
        .tickSize(-width+this.padding_left, 0, 0)

    var svg = d3.select(el).select('.d3_mu_axis')

    svg.append("g")
         .attr("class", "d3_mu_x_axis")
         .attr("transform", "translate(0," + (height -  this.padding_bottom) + ")")
         .call(xAxis)

    var axis_start = scales.x.range()[0];
    var axis_width = scales.x.range()[1]-scales.x.range()[0];
    svg.append("text")      // text label for the x axis
        .attr("x", axis_start + axis_width / 2 )
        .attr("y", height - this.padding_text )
        .style("text-anchor", "middle")
        .text("Sampling date");

    svg.append("g")
        .attr("class", "d3_mu_y_axis")
        .attr("transform", "translate(" + (this.padding_left) + ",0)")
        .call(yAxis);


    svg.append("text")
        .attr("transform", "rotate(-90)")
        .attr("y", this.padding_text)
        .attr("x",  - (height - this.padding_bottom)/2)
        .attr("dy", "1em")
        .style("text-anchor", "middle")
        .text("Distance to root");

    this._draw_text(el, scales)
};

MuPlot._draw_points = function(el, scales, dispatcher){

    //console.log("MU DRAW POINTS...")
    var g = d3.select(el).selectAll('.d3_mu_points');

    var tip = g.selectAll('.d3_mu_point')
        .data(this.points);

    tip.enter()
      .append("circle")
      .attr("class", "d3_mu_point")
      .attr("id", function(d){
        return "NAME" //(d.name).replace(/\//g, "")
      })

    tip
      .attr("cx", function(d){return scales.x(d.x)})
      .attr("cy", function(d){return scales.y(d.y)})
      .attr("r", this._tipRadius)
      .style("fill", function(d){return d.fill_color;})
      .style('stroke',function(d){return d.stroke_color;})
      .style('stroke-width',"1")

      .on('mouseover', function(d) {
          nodeTooltip.show(d);
          dispatcher.emit('point:point_mouseover', d);
      })
      .on('mouseout', function(d) {
          nodeTooltip.hide();
          dispatcher.emit('point:point_mouseout', d);
      })

        var tip = g.selectAll('.d3_mu_point')
        .data(this.points);

};

export default MuPlot;
