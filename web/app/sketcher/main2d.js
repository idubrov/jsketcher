import { Viewer } from './canvas.js';
import { WinManager, Window, List, Dock, Terminal, dockBtn, openWin, closeWin, DIRECTIONS } from '../ui';
import $ from 'jquery';
import { io } from './io';
import { LinearDimension } from './shapes/dim';
import { AddPointTool, AddSegmentTool } from './shapes/segment';
import { FilletTool } from './helpers';
import { AddArcTool } from './shapes/arc';
import { EditCircleTool } from './shapes/circle';
import { AddDimTool, AddCircleDimTool, Dimension, HDimension, VDimension } from './shapes/dim';

export const STORAGE_PREFIX = "TCAD.projects.";

/** @constructor */
function App2D() {
  var app = this;

  this.viewer = new Viewer(document.getElementById('viewer'));
  this.winManager = new WinManager();

  this.initSketchManager();
  this._exportWin = new Window($('#exportManager'), app.winManager);

  $('#exportManager li').click(function() {closeWin(app._exportWin);});

  this.constraintFilter = {};
  this.actions = {};

  //For debug view
  this._actionsOrder = [];

  var dockEl = $('#dock');
  var statusEl = $('#status');
  this.dock = new Dock(dockEl, statusEl, App2D.views);
  this.dock.show('Constraints');

  var consoleBtn = dockBtn('Commands', 'list');
  statusEl.append(consoleBtn);
  var commandsWin = new Window($('#commands'), this.winManager);
  commandsWin.tileUpRelative = $('#viewer');
  consoleBtn.click(function() {
    commandsWin.toggle();
  });
  new Terminal(commandsWin, function(command) {
    return "Command " + command + " executed";
  });

  this.winManager.registerResize(dockEl, DIRECTIONS.EAST, function() {$('body').trigger('layout'); });
  $('body').on('layout', this.viewer.onWindowResize);

  this.registerAction = function(id, desc, action) {
    app.actions[id] = {id: id, desc: desc, action: action};
    app._actionsOrder.push(id);
  };

  this.registerAction('new', "Create New Sketch", function () {
    app.newSketch();
  });

  this.registerAction('open', "Open Sketch", function (e) {
    app._sketchesList.refresh();
    openWin(app._sketchesWin, e);
  });

  this.registerAction('clone', "Clone Sketch", function () {
    app.cloneSketch();
  });

  this.registerAction('export', "Export", function (e) {
    openWin(app._exportWin, e);
  });

  this.registerAction('exportSVG', "Export To SVG", function () {
    io.exportTextData(app.viewer.io.svgExport(), app.getSketchId() + ".svg");
  });

  this.registerAction('exportDXF', "Export To DXF", function () {
    io.exportTextData(app.viewer.io.dxfExport(), app.getSketchId() + ".dxf");
  });

  this.registerAction('undo', "Undo", function () {
    app.viewer.historyManager.undo();
  });

  this.registerAction('redo', "Redo", function () {
    app.viewer.historyManager.redo();
  });

  this.registerAction('checkpoint', "Checkpoint", function () {
    app.viewer.historyManager.checkpoint();
  });

  this.registerAction('addPoint', "Add Point", function () {
    app.viewer.toolManager.takeControl(new AddPointTool(app.viewer));
  });

  this.registerAction('addSegment', "Add Segment", function () {
    app.viewer.toolManager.takeControl(new AddSegmentTool(app.viewer, false));
  });

  this.registerAction('addMultiSegment', "Add Multi Segment", function () {
    app.viewer.toolManager.takeControl(new AddSegmentTool(app.viewer, true));
  });

  this.registerAction('addArc', "Add Arc", function () {
    app.viewer.toolManager.takeControl(new AddArcTool(app.viewer));
  });

  this.registerAction('addCircle', "Add Circle", function () {
    app.viewer.toolManager.takeControl(new EditCircleTool(app.viewer));
  });

  this.registerAction('pan', "Pan", function () {
    app.viewer.toolManager.releaseControl();
  });

  this.registerAction('addFillet', "Add Fillet", function () {
    app.viewer.toolManager.takeControl(new FilletTool(app.viewer));
  });

  this.registerAction('addDim', "Add Dimension", function () {
    app.viewer.toolManager.takeControl(new AddDimTool(app.viewer, app.viewer.dimLayer, function(a,b) {return new Dimension(a,b)} ));
  });

  this.registerAction('addHDim', "Add Horizontal Dimension", function () {
    app.viewer.toolManager.takeControl(new AddDimTool(app.viewer, app.viewer.dimLayer, function(a,b) {return new HDimension(a,b)} ));
  });

  this.registerAction('addVDim', "Add Vertical Dimension", function () {
    app.viewer.toolManager.takeControl(new AddDimTool(app.viewer, app.viewer.dimLayer, function(a,b) {return new VDimension(a,b)} ));
  });

  this.registerAction('addCircleDim', "Add Circle Dimension", function () {
    app.viewer.toolManager.takeControl(new AddCircleDimTool(app.viewer, app.viewer.dimLayer));
  });

  this.registerAction('save', "Save", function () {
      var sketchData = app.viewer.io.serializeSketch();
      var sketchId = app.getSketchId();
      localStorage.setItem(app.getSketchId(), sketchData);
      app.viewer.historyManager.checkpoint();
  });

  this.registerAction('coincident', "Coincident", function () {
    app.viewer.parametricManager.coincident(app.viewer.selected);
  });

  this.registerAction('verticalConstraint', "Vertical Constraint", function () {
    app.viewer.parametricManager.vertical(app.viewer.selected);
  });

  this.registerAction('horizontalConstraint', "Horizontal Constraint", function () {
    app.viewer.parametricManager.horizontal(app.viewer.selected);
  });

  this.registerAction('parallelConstraint', "Parallel Constraint", function () {
    app.viewer.parametricManager.parallel(app.viewer.selected);
  });

  this.registerAction('perpendicularConstraint', "Perpendicular Constraint", function () {
    app.viewer.parametricManager.perpendicular(app.viewer.selected);
  });

  this.registerAction('P2LDistanceConstraint', "Distance Between Point and Line", function () {
    app.viewer.parametricManager.p2lDistance(app.viewer.selected, prompt);
  });

  this.registerAction('P2PDistanceConstraint', "Distance Between two Points", function () {
    app.viewer.parametricManager.p2pDistance(app.viewer.selected, prompt);
  });

  this.registerAction('RadiusConstraint', "Radius Constraint", function () {
    app.viewer.parametricManager.radius(app.viewer.selected, prompt);
  });

  this.registerAction('EntityEqualityConstraint', "Radius Equals Constraint", function () {
    app.viewer.parametricManager.entityEquality(app.viewer.selected);
  });

  this.registerAction('tangentConstraint', "Tangent Constraint", function () {
    app.viewer.parametricManager.tangent(app.viewer.selected);
  });

  this.registerAction('lockConstraint', "Lock Constraint", function () {
    app.viewer.parametricManager.lock(app.viewer.selected);
  });

  this.registerAction('pointOnLine', "Point On Line", function () {
    app.viewer.parametricManager.pointOnLine(app.viewer.selected);
  });

  this.registerAction('pointOnArc', "Point On Arc", function () {
    app.viewer.parametricManager.pointOnArc(app.viewer.selected);
  });

  this.registerAction('pointInMiddle', "Point In the Middle", function () {
    app.viewer.parametricManager.pointInMiddle(app.viewer.selected);
  });

  this.registerAction('llAngle', "Angle Between 2 Lines", function () {
    app.viewer.parametricManager.llAngle(app.viewer.selected, prompt);
  });

  this.registerAction('symmetry', "Symmetry", function () {
    app.viewer.parametricManager.symmetry(app.viewer.selected, prompt);
  });
  this.registerAction('lockConvex', "Lock Convexity", function () {
    app.viewer.parametricManager.lockConvex(app.viewer.selected, alert);
  });

  this.registerAction('analyzeConstraint', "Analyze Constraint", function () {
    app.viewer.parametricManager.analyze(alert);
  });

  this.registerAction('solve', "Solve System", function () {
    app.viewer.parametricManager.solve();
    app.viewer.refresh();
  });

  this.registerAction('CLEAN UP', "Clean All Draw", function () {
    app.cleanUpData();
    app.viewer.refresh();
  });

  this.registerAction('fit', "Fit Sketch On Screen", function () {
    app.fit();
    app.viewer.refresh();
  });
};

App2D.views = [
  {
    name: 'Dimensions',
    icon: 'arrows-v'
  },
  {
    name: 'Properties',
    icon: 'sliders'
  },
  {
    name: 'Constraints',
    icon: 'cogs'
  }
];

App2D.bottomViews = [
  {
    name: 'Commands',
    icon: 'desktop'
  }
];

App2D.faBtn = function(iconName) {
  return $('<i>', {'class' : 'fa fa-'+iconName});
};

App2D.prototype.fit = function() {

  var bbox = new io.BBox();

  for (var l = 0; l < this.viewer.layers.length; ++l) {
    var layer = this.viewer.layers[l];
    for (var i = 0; i < layer.objects.length; ++i) {
      var obj = layer.objects[i];
      bbox.check(obj);
    }
  }
  if (!bbox.isValid()) {
    return;
  }
  var bounds = bbox.bbox;
  this.viewer.showBounds(bounds[0], bounds[1], bounds[2], bounds[3]);
  bbox.inc(20 / this.viewer.scale);
  this.viewer.showBounds(bounds[0], bounds[1], bounds[2], bounds[3]);
};

App2D.prototype.cloneSketch = function() {
  var name = prompt("Name for sketch clone");
  if (name != null) {
    if (this.isSketchExists(name)) {
      alert("Sorry, a sketch with the name '" + name + "' already exists. Won't override it.");
      return;
    }
    localStorage.setItem(STORAGE_PREFIX + name, this.viewer.io.serializeSketch())
    this.openSketch(name);
  }
};

App2D.prototype.isSketchExists = function(name) {
  return localStorage.getItem(STORAGE_PREFIX + name) != null;
};

App2D.prototype.openSketch = function(name) {
  var uri = window.location.href.split("#")[0];
  if (name !== "untitled") {
    uri += "#" + name;
  }
  var win = window.open(uri, '_blank');
  win.focus();
};

App2D.prototype.newSketch = function() {
  var name = prompt("Name for sketch");
  if (name != null) {
    if (this.isSketchExists(name)) {
      alert("Sorry, a sketch with the name '" + name + "' already exists. Won't override it.");
      return;
    }
    this.openSketch(name);
  }
};

App2D.prototype.initSketchManager = function(data, ext) {
  this._sketchesWin = new Window($('#sketchManager'), this.winManager);
  var app = this;
  var sketchesList = new List('sketchList', {
    items : function() {
      var theItems = [];
      for (var name in localStorage) {
        if (!localStorage.hasOwnProperty(name)) {
          continue;
        }
        if (name.indexOf(STORAGE_PREFIX) === 0) {
          name = name.substring(STORAGE_PREFIX.length);
        }
        theItems.push({name : name});
      }
      return theItems;
    },

    remove : function(item) {
      if (confirm("Selected sketch will be REMOVED! Are you sure?")) {
        localStorage.removeItem(STORAGE_PREFIX + item.name);
        sketchesList.refresh();
      }
    },

    mouseleave : function(item) {},
    hover : function(item) {},

    click : function(item) {
      app.openSketch(item.name);
    }
  });
  $('#sketchManager').find('.content').append(sketchesList.ul);
  sketchesList.refresh();
  this._sketchesList = sketchesList;
};

App2D.prototype.loadFromLocalStorage = function() {
  var sketchId = this.getSketchId();
  var sketchData = localStorage.getItem(sketchId);
  if (sketchData != null) {
    this.viewer.historyManager.init(sketchData);
    this.viewer.io.loadSketch(sketchData);
  }
  this.viewer.repaint();
};

App2D.prototype.getSketchId = function() {
  var id = window.location.hash.substring(1);
  if (!id) {
    id = "untitled";
  }
  return STORAGE_PREFIX + id;
};

export default App2D;
