TCAD = {
  TWO : {}
};

TCAD.TWO.Styles = {
  DEFAULT : {
    lineWidth : 2,
    strokeStyle : "#ffffff", 
    fillStyle : "#000000"
  },
  
  SERVICE : {
    lineWidth : 0.3,
    strokeStyle : "#ff0000",
    fillStyle : "#FF0000"
  },

  MARK : {
    lineWidth : 2,
    strokeStyle : "#ff0000",
    fillStyle : "#FF0000"
  },

  SNAP : {
    lineWidth : 2,
    strokeStyle : "#00FF00",
    fillStyle : "#00FF00"
  },

  DIM : {
    lineWidth : 1,
    strokeStyle : "#bcffc1",
    fillStyle : "#00FF00"
  },

  BOUNDS : {
    lineWidth : 2,
    strokeStyle : "#fff5c3",
    fillStyle : "#000000"
  },
  
  CONSTRUCTION : {
    lineWidth : 1,
    strokeStyle : "#aaaaaa",
    fillStyle : "#000000"
  }
};

TCAD.TWO.utils = {};

TCAD.TWO.utils.extend = function(func, parent) {
  for(var prop in parent.prototype) {
    if(parent.prototype.hasOwnProperty(prop))
      func.prototype[prop] = parent.prototype[prop];
  }
};

TCAD.TWO.utils.point = function(x, y){ return {x: x, y: y} };

TCAD.TWO.utils.drawPoint = function (ctx, x, y, rad, scale) {
  ctx.beginPath();
  ctx.arc(x, y, rad / scale, 0, 2 * Math.PI, false);
  ctx.fill();
};

TCAD.TWO.utils.setStyle = function(style, ctx, scale) {
  ctx.lineWidth  = style.lineWidth / scale;
  ctx.strokeStyle  = style.strokeStyle;
  ctx.fillStyle  = style.fillStyle;
};

/** @constructor */
TCAD.TWO.Viewer = function(canvas) {
  
  this.canvas = canvas;
  this.params = new TCAD.Parameters();
  this.io = new TCAD.IO(this);
  var viewer = this;
  this.retinaPxielRatio = window.devicePixelRatio > 1 ? window.devicePixelRatio : 1;
  function updateCanvasSize() {
    var canvasWidth = canvas.parentNode.offsetWidth;;
    var canvasHeight = canvas.parentNode.offsetHeight;;

    canvas.width = canvasWidth * viewer.retinaPxielRatio;
    canvas.height = canvasHeight * viewer.retinaPxielRatio;

    canvas.style.width = canvasWidth + "px";
    canvas.style.height = canvasHeight + "px";
  }

  this.onWindowResize = function() {
    updateCanvasSize();
    viewer.refresh();
  };
  updateCanvasSize();
  window.addEventListener( 'resize', this.onWindowResize, false );

  Object.defineProperty(this, "activeLayer", {
    get: viewer.getActiveLayer ,
    set: viewer.setActiveLayer
  });

  this.bus = new TCAD.Bus();
  this.ctx = this.canvas.getContext("2d");
  this._activeLayer = null;
  this.layers = [];
  this._serviceLayers = [];
  this.dimLayer = new TCAD.TWO.Layer("_dim", TCAD.TWO.Styles.DIM);
  this.dimLayers = [this.dimLayer];
  this.bus.defineObservable(this, 'dimScale', 'dimScale', 1);
  this.bus.subscribe('dimScale', function(){ viewer.refresh(); });
  
  this._workspace = [this.dimLayers, this.layers, this._serviceLayers];
  this.toolManager = new TCAD.TWO.ToolManager(this, new TCAD.TWO.PanTool(this));
  this.parametricManager = new TCAD.TWO.ParametricManager(this);

  this.translate = {x : 0.0, y : 0.0};
  this.scale = 1.0;

  this.selected = [];
  this.snapped = [];
  
  this._setupServiceLayer();

  this.historyManager = new TCAD.HistoryManager(this);
  this.refresh();
};

TCAD.TWO.Viewer.prototype.validateGeom = function() {
  for (var i = 0; i < this.layers.length; i++) {
    var objs = this.layers[i].objects;
    for (var j = 0; j < objs.length; j++) {
      if (!objs[j].validate()) {
        return false;        
      }
    }
  }
  return true;
};

TCAD.TWO.Viewer.prototype.addSegment = function(x1, y1, x2, y2, layer) {
  var a = new TCAD.TWO.EndPoint(x1, y1);
  var b = new TCAD.TWO.EndPoint(x2, y2);
  var line = new TCAD.TWO.Segment(a, b);
  layer.objects.push(line);
  line.layer = layer;
  return line;
};

TCAD.TWO.Viewer.prototype.remove = function(obj) {
  if (obj.layer != null) {
    var idx = obj.layer.objects.indexOf(obj);
    if (idx != -1) {
      this.parametricManager.removeConstraintsByObj(obj);
      obj.layer.objects.splice(idx, 1);
    }
  }
};

TCAD.TWO.Viewer.prototype.add = function(obj, layer) {
  layer.objects.push(obj);
  obj.layer = layer;
};

TCAD.TWO.Viewer.prototype.search = function(x, y, buffer, deep, onlyPoints, filter) {

  buffer *= 0.5;
  
  var pickResult = [];
  var aim = new TCAD.Vector(x, y);

  var heroIdx = 0;
  var unreachable = buffer * 2;
  var heroLength = unreachable; // unreachable

  function isFiltered(o) {
    for (var i = 0; i < filter.length; ++i) {
      if (filter[i] === o) return true;
    }
    return false;
  }

  for (var i = 0; i < this.layers.length; i++) {
    var objs = this.layers[i].objects;
    for (var j = 0; j < objs.length; j++) {
      var l = unreachable + 1;
      var before = pickResult.length;
      objs[j].acceptV(true, function(o) {
        if (onlyPoints && o._class !== 'TCAD.TWO.EndPoint') {
          return false;  
        }
        l = o.normalDistance(aim);
        if (l >= 0 && l <= buffer && !isFiltered(o)) {
          pickResult.push(o);
          return false;
        }
        return true;
      });
      var hit = before - pickResult.length != 0;
      if (hit) {
        if (!deep && pickResult.length != 0) return pickResult;
        if (l >= 0 && l < heroLength) {
          heroLength = l;
          heroIdx = pickResult.length - 1;
        }
      }
    }
  }
  if (pickResult.length > 0) {
    var _f = pickResult[0];
    pickResult[0] = pickResult[heroIdx];
    pickResult[heroIdx] = _f;
  }
  return pickResult;
};

TCAD.TWO.Viewer.prototype._setupServiceLayer = function() {
  var layer = new TCAD.TWO.Layer("_service", TCAD.TWO.Styles.SERVICE);
//  layer.objects.push(new TCAD.TWO.CrossHair(0, 0, 20));
  layer.objects.push(new TCAD.TWO.BasisOrigin(null, this));
  layer.objects.push(new TCAD.TWO.Point(0, 0, 2));
  this._serviceLayers.push(layer);

  layer = new TCAD.TWO.Layer("_selection", TCAD.TWO.Styles.DEFAULT);
  layer.objects = this.selected;
  this._serviceLayers.push(layer);
};

TCAD.TWO.Viewer.prototype.refresh = function() {
  var viewer = this;
  window.requestAnimationFrame( function() {
    viewer.repaint();     
  });  
};

TCAD.TWO.Viewer.prototype.repaint = function() {

  var ctx = this.ctx;
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  
  ctx.fillStyle = "#808080";
  ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
  
  //Order is important!
  ctx.transform(1, 0, 0, -1, 0, this.canvas.height );
  ctx.transform(1, 0, 0, 1, this.translate.x , this.translate.y );
  ctx.transform(this.scale, 0, 0, this.scale, 0, 0);

  var prevStyle = null;
  var style;
  for (var w = 0; w < this._workspace.length; w++) {
    var layers = this._workspace[w];
    for (var l = 0; l < layers.length; l++) {
      var layer = layers[l];
      for (var o = 0; o < layer.objects.length; o++) {
        var obj = layer.objects[o];
        style = obj.style != null ? obj.style : layer.style;
        if (style != prevStyle) TCAD.TWO.utils.setStyle(style, ctx, this.scale / this.retinaPxielRatio);
        obj.draw(ctx, this.scale / this.retinaPxielRatio, this);
      }
    }
  }
};

TCAD.TWO.Viewer.prototype.snap = function(x, y, excl) {
  this.cleanSnap();
  var snapTo = this.search(x, y, 20 / this.scale, true, true, excl);
  if (snapTo.length > 0) {
    snapTo = snapTo[0];
    this.mark(snapTo, TCAD.TWO.Styles.SNAP);
    this.snapped.push(snapTo);
    return snapTo;
  }
  return null;
};

TCAD.TWO.Viewer.prototype.cleanSnap = function() {
  while(this.snapped.length > 0) {
    this.snapped.pop().marked = null;
  }
};

TCAD.TWO.Viewer.prototype.showBounds = function(x1, y1, x2, y2, offset) {
  var dx = x2 - x1;
  var dy = y2 - y1;
  if (this.canvas.width > this.canvas.height) {
    this.scale = this.canvas.height / dy;
  } else {
    this.scale = this.canvas.width / dx;
  }
  this.translate.x = -x1 * this.scale;
  this.translate.y = -y1 * this.scale;
};

TCAD.TWO.Viewer.prototype.screenToModel2 = function(x, y, out) {

  out.x = x * this.retinaPxielRatio;
  out.y = this.canvas.height - y * this.retinaPxielRatio;

  out.x -= this.translate.x;
  out.y -= this.translate.y;

  out.x /= this.scale;
  out.y /= this.scale;
};

TCAD.TWO.Viewer.prototype.screenToModel = function(e) {
  return this._screenToModel(e.offsetX, e.offsetY);
};

TCAD.TWO.Viewer.prototype._screenToModel = function(x, y) {
  var out = {x: 0, y: 0};
  this.screenToModel2(x, y, out);
  return out;
};

TCAD.TWO.Viewer.prototype.accept = function(visitor) {
  for (var i = 0; i < this.layers.length; i++) {
    var objs = this.layers[i].objects;
    var result = null;
    for (var j = 0; j < objs.length; j++) {
      if (!objs[j].accept(visitor)) {
        return false;
      }
    }
  }
};

TCAD.TWO.Viewer.prototype.findLayerByName = function(name) {
  for (var i = 0; i < this.layers.length; i++) {
    if (this.layers[i].name == name) {
      return this.layers[i];
    }
  }
  return null;
};

TCAD.TWO.Viewer.prototype.findById = function(id) {
  var result = null;
  this.accept(function(o) {
    if (o.id === id) {
      result = o;
      return false;
    }
    return true;
  });
  return result;
};

TCAD.TWO.Viewer.prototype.select = function(objs, exclusive) {
  if (exclusive) this.deselectAll();
  for (var i = 0; i < objs.length; i++) {
    this.mark(objs[i]);
  }
};

TCAD.TWO.Viewer.prototype.pick = function(e) {
  var m = this.screenToModel(e);
  return this.search(m.x, m.y, 20 / this.scale, true, false, []);
};

TCAD.TWO.Viewer.prototype.mark = function(obj, style) {
  if (style === undefined) {
    style = TCAD.TWO.Styles.MARK;
  }
  obj.marked = style;
  this.selected.push(obj);
};

TCAD.TWO.Viewer.prototype.getActiveLayer = function() {
  var layer = this._activeLayer;
  if (layer == null || layer.readOnly) {
    layer = null;
    for (var i = 0; i < this.layers.length; i++) {
      var l = this.layers[i];
      if (!l.readOnly) {
        layer = l;
        break;
      }
    }
  }
  if (layer == null) {
    layer = new TCAD.TWO.Layer("JustALayer", TCAD.TWO.Styles.DEFAULT);
    this.layers.push(layer);
  }
  return layer;
};

TCAD.TWO.Viewer.prototype.setActiveLayer = function(layer) {
  if (!layer.readOnly) {
    this._activeLayer = layer;
    this.bus.notify("activeLayer");
  }
};

TCAD.TWO.Viewer.prototype.deselect = function(obj) {
  for (var i = 0; i < this.selected.length; i++) {
    if (obj.id == this.selected[i].id) {
      this.selected.splice(i, 1)[0].marked = null;
      break;
    }
  }
};

TCAD.TWO.Viewer.prototype.deselectAll = function() {
  for (var i = 0; i < this.selected.length; i++) {
    this.selected[i].marked = null;
  }
  while(this.selected.length > 0) this.selected.pop();
};

/** @constructor */
TCAD.TWO.Layer = function(name, style) {
  this.name = name;
  this.style = style;
  this.objects = [];
  this.readOnly = false; // This is actually a mark for boundary layers coming from 3D
};

TCAD.TWO.Viewer.prototype.fullHeavyUIRefresh = function() {
  this.refresh();
  this.parametricManager.notify();
};

/** @constructor */
TCAD.TWO.Polygon = function(points) {
  this.points = points;
  this.style = null;
};

TCAD.TWO.Polygon.prototype.draw = function(ctx) {

  if (this.points.length < 3) {
    return;    
  }
  
  ctx.beginPath();
  var first = this.points[0];
  ctx.moveTo(first.x, first.y);
  for (var i = 1; i < this.points.length; i++) {
    var p = this.points[i];
    ctx.lineTo(p.x, p.y);    
  }
  ctx.closePath();
  ctx.stroke(); 
};

/** @constructor */
TCAD.TWO.Polyline = function(points) {
  this.points = points;
  this.style = null;
};

TCAD.TWO.Polyline.prototype.draw = function(ctx) {

  if (this.points.length < 2) {
    return;
  }

  var first = this.points[0];
  ctx.moveTo(first.x, first.y);
  for (var i = 1; i < this.points.length; i++) {
    var p = this.points[i];
    ctx.lineTo(p.x, p.y);
  }
  ctx.stroke();
};

TCAD.TWO.utils.ID_COUNTER = 0;

TCAD.TWO.utils.genID = function() {
  return TCAD.TWO.utils.ID_COUNTER ++;
};

/** @constructor */
TCAD.TWO.SketchObject = function() {
  this.id = TCAD.TWO.utils.genID();
  this.aux = false;
  this.marked = null;
  this.visible = true;
  this.children = [];
  this.linked = [];
  this.layer = null;
};

TCAD.TWO.SketchObject.prototype.accept = function(visitor) {
  return this.acceptV(false, visitor);
};

TCAD.TWO.SketchObject.prototype.acceptV = function(onlyVisible, visitor) {
  if (onlyVisible && !this.visible) return true;
  for (var i = 0; i < this.children.length; i++) {
    var child = this.children[i];
    if (!child.acceptV(onlyVisible, visitor)) {
      return false;
    }
  }
  return visitor(this);
};

TCAD.TWO.SketchObject.prototype.validate = function() {
  return true;
};


TCAD.TWO.SketchObject.prototype.getDefaultTool = function(viewer) {
  return new TCAD.TWO.DragTool(this, viewer);
};

TCAD.TWO.SketchObject.prototype.isAuxOrLinkedTo = function() {
  if (!!this.aux) {
    return true;
  }
  for (var i = 0; i < this.linked.length; ++i) {
    if (!!this.linked[i].aux) {
      return true;
    }
  }
  return false;
};

TCAD.TWO.SketchObject.prototype._translate = function(dx, dy, translated) {
  translated[this.id] = 'x';
  for (var i = 0; i < this.linked.length; ++i) {
    if (translated[this.linked[i].id] != 'x') {
      this.linked[i]._translate(dx, dy, translated);
    }
  }
  this.translateImpl(dx, dy);
};

TCAD.TWO.SketchObject.prototype.translate = function(dx, dy) {
//  this.translateImpl(dx, dy);
  if (this.isAuxOrLinkedTo()) {
    return;
  }
  this._translate(dx, dy, {});
};

TCAD.TWO.SketchObject.prototype.draw = function(ctx, scale, viewer) {
  if (!this.visible) return;
  if (this.marked != null) {
    ctx.save();
    TCAD.TWO.utils.setStyle(this.marked, ctx, scale);
  }
  this.drawImpl(ctx, scale, viewer);
  if (this.marked != null) ctx.restore();
  for (var i = 0; i < this.children.length; i++) {
    this.children[i].draw(ctx, scale);
  }
};

/** @constructor */
TCAD.TWO.Ref = function(value) {
  this.id = TCAD.TWO.utils.genID();
  this.value = value;
};

TCAD.TWO.Ref.prototype.set = function(value) {
  this.value = value;
};

TCAD.TWO.Ref.prototype.get = function() {
  return this.value;
};

/** @constructor */
TCAD.TWO.Param = function(obj, prop) {
  this.id = TCAD.TWO.utils.genID();
  this.obj = obj;
  this.prop = prop;
};

TCAD.TWO.Param.prototype.set = function(value) {
  this.obj[this.prop] = value;
};

TCAD.TWO.Param.prototype.get = function() {
  return this.obj[this.prop];
};

/** @constructor */
TCAD.TWO.EndPoint = function(x, y) {
  TCAD.TWO.SketchObject.call(this);
  this.x = x;
  this.y = y;
  this.parent = null;
  this._x =  new TCAD.TWO.Param(this, 'x');
  this._y =  new TCAD.TWO.Param(this, 'y');
};

TCAD.TWO.utils.extend(TCAD.TWO.EndPoint, TCAD.TWO.SketchObject);

TCAD.TWO.EndPoint.prototype._class = 'TCAD.TWO.EndPoint';

TCAD.TWO.EndPoint.prototype.collectParams = function(params) {
  params.push(this._x);
  params.push(this._y);
};

TCAD.TWO.EndPoint.prototype.normalDistance = function(aim) {
  return aim.minus(new TCAD.Vector(this.x, this.y)).length();
};

TCAD.TWO.EndPoint.prototype.getReferencePoint = function() {
  return this;
};

TCAD.TWO.EndPoint.prototype.translateImpl = function(dx, dy) {
  this.x += dx;
  this.y += dy;
};

TCAD.TWO.EndPoint.prototype.drawImpl = function(ctx, scale) {
  TCAD.TWO.utils.drawPoint(ctx, this.x, this.y, 3, scale)
};

/** @constructor */
TCAD.TWO.Segment = function(a, b) {
  TCAD.TWO.SketchObject.call(this);
  this.a = a;
  this.b = b;
  a.parent = this;
  b.parent = this;
  this.children.push(a, b);
};

TCAD.TWO.utils.extend(TCAD.TWO.Segment, TCAD.TWO.SketchObject);

TCAD.TWO.Segment.prototype._class = 'TCAD.TWO.Segment';
TCAD.TWO.Segment.MIN_LENGTH = 0.1; // 0.08; // if length < 0.08 canvas doesn't even draw a line

TCAD.TWO.Segment.prototype.validate = function() {
  return TCAD.math.distanceAB(this.a, this.b) > TCAD.TWO.Segment.MIN_LENGTH;
};

TCAD.TWO.Segment.prototype.addFixingGeomConstraints = function(constrs) {
  constrs.push(new TCAD.TWO.Constraints.MinLength(this.a, this.b, TCAD.TWO.Segment.MIN_LENGTH));
};

TCAD.TWO.Segment.prototype.collectParams = function(params) {
  this.a.collectParams(params);
  this.b.collectParams(params);
};

TCAD.TWO.Segment.prototype.normalDistance = function(aim) {
  var x = aim.x;
  var y = aim.y;

  var ab = new TCAD.Vector(this.b.x - this.a.x, this.b.y - this.a.y)
  var e = ab.normalize();
  var a = new TCAD.Vector(aim.x - this.a.x, aim.y - this.a.y);
  var b = e.multiply(a.dot(e));
  var n = a.minus(b);

  //Check if vector b lays on the vector ab
  if (b.length() > ab.length()) {
    return -1;
  }

  if (b.dot(ab) < 0) {
    return -1;
  }

  return n.length();
};

TCAD.TWO.Segment.prototype.getReferencePoint = function() {
  return this.a;
};

TCAD.TWO.Segment.prototype.translateImpl = function(dx, dy) {
  this.a.translate(dx, dy);
  this.b.translate(dx, dy);
};

TCAD.TWO.Segment.prototype.drawImpl = function(ctx, scale) {
  ctx.beginPath();
  ctx.moveTo(this.a.x, this.a.y);
  ctx.lineTo(this.b.x, this.b.y);
//  ctx.save();
//  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.stroke();
//  ctx.restore();
};

/** @constructor */
TCAD.TWO.Point = function(x, y, rad) {
  this.x = x;
  this.y = y;
  this.rad = rad;
  this.style = null;
};

TCAD.TWO.Point.prototype.draw = function(ctx, scale) {
  TCAD.TWO.utils.drawPoint(ctx, this.x, this.y, this.rad, scale);
};

/** @constructor */
TCAD.TWO.CrossHair = function(x, y, rad) {
  this.x = x;
  this.y = y;
  this.rad = rad;
  this.style = null;
};

TCAD.TWO.CrossHair.prototype.draw = function(ctx, scale) {
  ctx.beginPath();
  var rad = this.rad / scale;
  ctx.moveTo(this.x - rad, this.y);
  ctx.lineTo(this.x + rad, this.y);
  ctx.closePath();
  ctx.moveTo(this.x, this.y - rad);
  ctx.lineTo(this.x, this.y + rad);
  ctx.closePath();

  ctx.save();
  ctx.setTransform(1, 0, 0, 1, 0, 0);
  ctx.stroke();
  ctx.restore();
};

/** @constructor */
TCAD.TWO.BasisOrigin = function(basis, viewer) {
  this.viewer = viewer;
  this.inverseX = false;
  this.inverseY = false;
  this.lineWidth = 100;
  this.xColor = '#FF0000';
  this.yColor = '#00FF00';
};

TCAD.TWO.BasisOrigin.prototype.draw = function(ctx, scale) {
  ctx.save();
  if (this.inverseX) {
    this.xScale = -1;
    this.xShift = this.lineWidth + 10;
  } else {
    this.xScale = 1;
    this.xShift = 10;
  }
  if (this.inverseY) {
    this.yScale = -1;
    this.yShift = this.viewer.canvas.height - this.lineWidth - 10;
  } else {
    this.yScale = 1;
    this.yShift = this.viewer.canvas.height - 10;
  }

  ctx.setTransform( this.xScale, 0, 0, this.yScale, this.xShift, this.yShift);
  ctx.beginPath();

  ctx.lineWidth  = 1;
  ctx.strokeStyle  = this.yColor;

  var headA = 5;
  var headB = 10;

  ctx.moveTo(0, 0);
  ctx.lineTo(0, - this.lineWidth);
  
  ctx.moveTo(0, - this.lineWidth);
  ctx.lineTo(headA, 0 - this.lineWidth + headB);

  ctx.moveTo(0, - this.lineWidth);
  ctx.lineTo(- headA, - this.lineWidth + headB);
  ctx.closePath();
  ctx.stroke();

  ctx.beginPath();
  ctx.strokeStyle  = this.xColor;
  ctx.moveTo(0, 0);
  ctx.lineTo(this.lineWidth, 0);

  ctx.moveTo(this.lineWidth, 0);
  ctx.lineTo(this.lineWidth - headB, headA);

  ctx.moveTo(this.lineWidth, 0);
  ctx.lineTo(this.lineWidth - headB, - headA);
  ctx.closePath();
  ctx.stroke();

  ctx.restore();
};

/** @constructor */
TCAD.TWO.ToolManager = function(viewer, defaultTool) {
  this.defaultTool = defaultTool;
  this.tool = defaultTool;
  var canvas = viewer.canvas;
  var tm = this;
  canvas.addEventListener('mousemove', function (e) {
    e.preventDefault();
    //e.stopPropagation(); // allow propagation for move in sake of dynamic layout 
    tm.getTool().mousemove(e);
  }, false);
  canvas.addEventListener('mousedown', function (e) {
    e.preventDefault();
    e.stopPropagation();
    tm.getTool().mousedown(e);
  }, false);
  canvas.addEventListener('mouseup', function (e) {
    e.preventDefault();
    e.stopPropagation();
    tm.getTool().mouseup(e);
  }, false);
  canvas.addEventListener('mousewheel', function (e) {
    e.preventDefault();
    e.stopPropagation();
    tm.getTool().mousewheel(e);
  }, false);
  canvas.addEventListener('dblclick', function (e) {
    e.preventDefault();
    e.stopPropagation();
    if (tm.tool.dblclick !== undefined) {
      tm.tool.dblclick(e);
    }
  }, false);

  window.addEventListener("keydown", function (e) {
    tm.getTool().keydown(e);
    if (e.keyCode == 27) {
      tm.releaseControl();
    } else if (e.keyCode == 46 || e.keyCode == 8) {
      var selection = viewer.selected.slice();
      viewer.deselectAll();
      for (var i = 0; i < selection.length; i++) {
        viewer.remove(selection[i]);
      }
      viewer.refresh();
    }
  }, false);
  window.addEventListener("keypress", function (e) {
    tm.getTool().keydown(e);
  }, false);
  window.addEventListener("keyup", function (e) {
    tm.getTool().keydown(e);
  }, false);
};

TCAD.TWO.ToolManager.prototype.takeControl = function(tool) {
  this.tool = tool;
};

TCAD.TWO.ToolManager.prototype.releaseControl = function() {
  this.tool.cleanup();
  this.tool = this.defaultTool;
};

TCAD.TWO.ToolManager.prototype.getTool = function() {
  return this.tool;
};

/** @constructor */
TCAD.TWO.PanTool = function(viewer) {
  this.viewer = viewer;
  this.dragging = false;
  this.x = 0.0;
  this.y = 0.0;
};

TCAD.TWO.PanTool.prototype.keydown = function(e) {};
TCAD.TWO.PanTool.prototype.keypress = function(e) {};
TCAD.TWO.PanTool.prototype.keyup = function(e) {};
TCAD.TWO.PanTool.prototype.cleanup = function(e) {};

TCAD.TWO.PanTool.prototype.mousemove = function(e) {
  if (!this.dragging) {
    return;    
  }
  var dx = e.pageX - this.x;
  var dy = e.pageY - this.y;
  dy *= -1;
  
  this.viewer.translate.x += dx * this.viewer.retinaPxielRatio;
  this.viewer.translate.y += dy * this.viewer.retinaPxielRatio;
  
  this.x = e.pageX;
  this.y = e.pageY;
  this.deselectOnUp = false;
  this.viewer.refresh();
};

TCAD.TWO.PanTool.prototype.mousedown = function(e) {
  if (e.button == 0) {
    var picked = this.viewer.pick(e);
    if (picked.length > 0) {
      if (e.shiftKey) {
        var toSelect = picked[0];
        var ids = this.viewer.selected.map(function(s){return s.id});        
        for (var i = 0; i < picked.length; i++) {
          if (ids.indexOf(picked[i].id) != -1) {
            this.viewer.deselect(picked[i]);
          } else {
            toSelect = picked[i];                
          }
        }
        this.viewer.select([toSelect], false);
        this.deselectOnUp = false;
      } else {
        var toSelect = picked[0];
        if (this.viewer.selected.length === 1) {
          for (var i = 0; i < picked.length - 1; i++) {
            if (picked[i].id == this.viewer.selected[0].id) {
              toSelect = picked[i + 1];
              break;
            }
          }
        }
        this.viewer.select([toSelect], true);
        if (!toSelect.isAuxOrLinkedTo()) {
          var tool = toSelect.getDefaultTool(this.viewer);
          tool.mousedown(e);
          this.viewer.toolManager.takeControl(tool);
        }
      }
      this.viewer.refresh();
      return;
    }
  }

  this.dragging = true;
  this.deselectOnUp = true;
  this.x = e.pageX;
  this.y = e.pageY;
};

TCAD.TWO.PanTool.prototype.mouseup = function(e) {
  this.dragging = false;
  if (this.deselectOnUp) {
    this.viewer.deselectAll();
    this.viewer.refresh();
  }
  this.deselectOnUp = false;
};

TCAD.TWO.PanTool.prototype.mousewheel = function(e) {

  var delta = 0;

  if ( e.wheelDelta ) { // WebKit / Opera / Explorer 9
    delta = e.wheelDelta;
  } else if ( e.detail ) { // Firefox
    delta = - e.detail;
  }

  var before = this.viewer.screenToModel(e);

  var step = 0.05;
  delta = delta < 0 ? 1 - step : 1 + step;
  this.viewer.scale *= delta;

  var after = this.viewer.screenToModel(e);
  
  var dx = after.x - before.x;
  var dy = after.y - before.y;

  this.viewer.translate.x += dx * this.viewer.scale;
  this.viewer.translate.y += dy * this.viewer.scale;

  this.viewer.refresh();
};

/** @constructor */
TCAD.TWO.DragTool = function(obj, viewer) {
  this.obj = obj;
  this.viewer = viewer;
  this._point = {x: 0, y: 0};
  this.origin = {x: 0, y: 0};
  this.lastSolved = {x: 0, y: 0};
  this.ref = this.obj.getReferencePoint();
  this.solver = null;
};
TCAD.TWO.DragTool.prototype.snapshots = []
TCAD.TWO.DragTool.prototype.keydown = function(e) {};
TCAD.TWO.DragTool.prototype.keypress = function(e) {};
TCAD.TWO.DragTool.prototype.keyup = function(e) {};
TCAD.TWO.DragTool.prototype.cleanup = function(e) {};

TCAD.TWO.DragTool.prototype.mousemove = function(e) {
  var x = this._point.x;
  var y = this._point.y;
  this.viewer.screenToModel2(e.offsetX, e.offsetY, this._point);
  var dx = this._point.x - x;
  var dy = this._point.y - y;
  for (var i = 0; i < this.lockedShifts.length; i += 2) {
    this.lockedValues[i] = this._point.x - this.lockedShifts[i];
    this.lockedValues[i + 1] = this._point.y - this.lockedShifts[i + 1];
  }
  this.solver.updateLock(this.lockedValues);
  if (!e.altKey && !e.ctrlKey) {
    var distanceSinceLastSolved = TCAD.math.distance(this.lastSolved.x, this.lastSolved.y, this._point.x, this._point.y)
    var MIN_STEP = 5;
    if (distanceSinceLastSolved > MIN_STEP) {
      var dir = new TCAD.Vector(this._point.x - this.lastSolved.x, this._point.y - this.lastSolved.y);
      dir._normalize();
      dir._multiply(MIN_STEP);
      var increment = dir.copy();
      for (var keepMoving = true; keepMoving; dir._plus(increment)) {
        var emulatedX, emulatedY;
        if (dir.length() >= distanceSinceLastSolved) {
          keepMoving = false;
          emulatedX = this._point.x;
          emulatedY = this._point.y;
        } else {
          emulatedX = this.lastSolved.x + dir.x;
          emulatedY = this.lastSolved.y + dir.y;
        }
        for (var i = 0; i < this.lockedShifts.length; i += 2) {
          this.lockedValues[i] = emulatedX - this.lockedShifts[i];
          this.lockedValues[i + 1] = emulatedY - this.lockedShifts[i + 1];
        }
        this.solver.updateLock(this.lockedValues);
        this.solveRequest(true);
      }
    } else {
      this.solveRequest(true);
    }
    this.lastSolved.x = this._point.x;
    this.lastSolved.y = this._point.y;
  } else {
    this.obj.translate(dx, dy);
  }

  this.viewer.refresh();
};

TCAD.TWO.DragTool.prototype.mousedown = function(e) {
  this.origin.x = e.offsetX;
  this.origin.y = e.offsetY;
  this.viewer.screenToModel2(e.offsetX, e.offsetY, this._point);
  this.lastSolved.x = this._point.x;
  this.lastSolved.y = this._point.y;
  this.prepareSolver([]);
};

TCAD.TWO.DragTool.prototype.mouseup = function(e) {
  this.solveRequest(false);
  this.viewer.refresh();
  this.viewer.toolManager.releaseControl();
  var traveled = TCAD.math.distance(this.origin.x, this.origin.y, e.offsetX, e.offsetY);
  if (traveled >= 10) {
    this.viewer.historyManager.lightCheckpoint(10);
  }
  //this.animateSolution();
};

TCAD.TWO.DragTool.prototype.mousewheel = function(e) {
};

TCAD.TWO.DragTool.prototype.solveRequest = function(rough) {
  this.solver.solve(rough, 1);
  this.solver.sync();
  if (false) {
    var fixingConstrs = [];
    this.viewer.accept(function(obj) {
      if (obj._class === 'TCAD.TWO.Segment') {
        if (!obj.validate()) {
          obj.addFixingGeomConstraints(fixingConstrs);
        }
      }
      return true;
    });
    if (fixingConstrs.length !== 0) {
      this.prepareSolver(fixingConstrs);
      this.solver.solve(rough, 1);
      this.solver.sync();
    }
  }
};

TCAD.TWO.DragTool.prototype.getParamsToLock = function() {
  var params = [];
  this.obj.accept(function(obj) {
    if (obj._class === 'TCAD.TWO.EndPoint') {
      params.push(obj._x);
      params.push(obj._y);
    }
    return true;
  });
  return params;
};

TCAD.TWO.DragTool.prototype.prepareSolver = function(extraConstraints) {
  var locked = this.getParamsToLock();
  this.lockedShifts = [];
  this.lockedValues = [];
  for (var i = 0; i < locked.length; i += 2) {
    this.lockedShifts[i] = this._point.x - locked[i].get();
    this.lockedShifts[i + 1] = this._point.y - locked[i + 1].get();
  }
  this.solver = this.viewer.parametricManager.prepare(locked, extraConstraints);
  //this.enableRecording();
};

TCAD.TWO.DragTool.prototype.enableRecording = function() {
  var solver = this.solver;
  var snapshots = this.snapshots = [];
  optim.DEBUG_HANDLER = function()  {
    snapshots.push([]);
    for (var i = 0; i < solver.solvers.length; i++) {
      var sys = solver.solvers[i].system;
      snapshots[i].push(sys.params.map(function(p) {return p.get()}))
    }
  };
};

TCAD.TWO.DragTool.prototype.animateSolution = function() {
  if (this.snapshots.length === 0) return; 
  var stepNum = 0;
  var scope = this;
  var then = Date.now();
  var speed = 3000; 
  function step() {
    var now = Date.now();
    var elapsed = now - then;

    if (elapsed > speed) {
      for (var i = 0; i < scope.solver.solvers.length; i++) {
        var sys = scope.solver.solvers[i].system;
        if (stepNum >= scope.snapshots[i].length) continue;
        var values = scope.snapshots[i][stepNum];
        for (var k = 0; k < values.length; k++) {
          sys.params[k]._backingParam.set(values[k]);
        }
      }
      stepNum ++;
      
      then = now;
      scope.viewer.repaint();
    }

    if (stepNum < scope.snapshots[0].length) {
      window.requestAnimationFrame(step);
    } 
  }
  window.requestAnimationFrame(step);
};