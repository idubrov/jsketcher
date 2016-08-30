import { PNLTRI } from 'pnltri';
import { CSG } from 'openjscad-csg';
import { Vector, Matrix } from './math/vector';
import math from './math/math';
import { view, DPR } from './3d/viewer';
import struct from './3d/hashmap';
import graph from './math/graph';
import { craft } from './workbench';
import { Vector3, Object3D, Geometry, ShaderMaterial, SpriteMaterial, Sprite,
  LineBasicMaterial, Line, SphereGeometry, MeshBasicMaterial, Mesh, MeshPhongMaterial,
  FaceColors, DoubleSide, Face3 } from 'three';

const utils = {};

utils.createSquare = function(width) {

  width /= 2;

  return [
    new Vector(-width, -width, 0),
    new Vector( width, -width, 0),
    new Vector( width,  width, 0),
    new Vector(-width,  width, 0)
  ];
};

utils.csgVec = function(v) {
  return new CSG.Vector3D(v.x, v.y, v.z);
};

utils.vec = function(v) {
  return new Vector(v.x, v.y, v.z);
};

utils.createBox = function(width) {
  var square = utils.createSquare(width);
  var rot = math.rotateMatrix(3/4, math.AXIS.Z, math.ORIGIN);
  square.forEach(function(v) { rot._apply(v) } );
  var normal = geom.normalOfCCWSeq(square);
  return geom.extrude(square, normal, normal.multiply(width), 1);
};

utils.createCSGBox = function(width) {
  var csg = CSG.fromPolygons(utils.createBox(width));
  return utils.createSolid(csg);
};

utils.toCsgGroups = function(polygons) {
  var groups = [];
  for (var i = 0; i < polygons.length; i++) {
    var p = polygons[i];
    if (p.holes.length === 0) {
      groups.push( new CSGGroup([new TCAD.SimplePolygon(p.shell, p.normal)], p.normal) );
    } else {
      // TODO: triangulation needed
      groups.push( new TCAD.CSGGroup([new TCAD.SimplePolygon(p.shell, p.normal)], p.normal) );
    }
  }
  return groups;
};

utils.checkPolygon = function(poly) {
  if (poly.length < 3) {
    throw new Error('Polygon should contain at least 3 point');
  }
};

utils.createPoint = function(x, y, z) {
//  var g = new PlaneGeometry(0.05, 0.05);
//  var m = new MeshBasicMaterial({color: 0x0000ff, side: DoubleSide});
//  return new Mesh(g, m);

  var material = new ShaderMaterial({
//    color: 0xff0000,
//    linewidth: 5
    vertexShader :
      'void main() {\n\t' +
      'gl_Position = projectionMatrix * modelViewMatrix * vec4( position, 1.0 );' +
      'gl_PointSize =10.0;\n\t' +
     '\n}',

    fragmentShader :
        'void main() {\n\t' +
        "vec2 coord = gl_PointCoord - vec2(0.5);  //from [0,1] to [-0.5,0.5]\n" +
        "if(length(coord) > 0.5)                  //outside of circle radius?\n" +
        "    discard;\n"+
        "else\n"+
        "    gl_FragColor = vec4( 1.0, 0.0, 0.0, 1.0 );\n"
    +'\n}'
  });

  var geometry = new Geometry();
  geometry.vertices.push(new Vector3(x, y, z));
//  geometry.vertices.push(new Vector3(x+.001, y+.001, z+.001));

//  var line = new PointCloud(geometry, material);
//  line.position.x = x;
//  line.position.y = y;
//  line.position.z = z;
//  return line;

  material = new SpriteMaterial( { color: 0xffffff, fog: false } );
  var sprite = new Sprite( material );
  sprite.position.set( x, y, z );
  return sprite;
};

utils.createLine = function (a, b, color) {
  var material = new LineBasicMaterial({
    color: color,
    linewidth: 1
  });
  var geometry = new Geometry();
  geometry.vertices.push(new Vector3(a.x, a.y, a.z));
  geometry.vertices.push(new Vector3(b.x, b.y, b.z));
  return new Line(geometry, material);
};

utils.createPoint = function (x, y, z) {
  var geometry = new SphereGeometry( 5, 16, 16 );
  var material = new MeshBasicMaterial( {color: 0xff0000} );
  var sphere = new Mesh(geometry, material);
  sphere.position.x = x;
  sphere.position.y = y;
  sphere.position.z = z;
  return sphere;
};

utils.createSolidMaterial = function() {
  return new MeshPhongMaterial({
    vertexColors: FaceColors,
    color: view.FACE_COLOR,
    shininess: 0,
    polygonOffset : true,
    polygonOffsetFactor : 1,
    polygonOffsetUnits : 2,
    side : DoubleSide
  });
};

utils.createSolid = function(csg) {
  var material = utils.createSolidMaterial();
  return new Solid(csg, material);
};

utils.intercept = function(obj, methodName, aspect) {
  var originFunc = obj[methodName];
  obj[methodName] = function() {
    var $this = this;
    aspect(function() {originFunc.apply($this, arguments)}, arguments);
  }
};

utils.createPlane = function(basis, depth) {
  var tu = utils;

  var initWidth = 1;
  var boundingPolygon = [
      new Vector(0,  0, 0),
      new Vector(initWidth,  0, 0),
      new Vector(initWidth, initWidth, 0),
      new Vector(0, initWidth, 0)
    ];
  var shared = tu.createShared();

  var material = tu.createSolidMaterial();
  material.transparent = true;
  material.opacity = 0.5;
  material.side = DoubleSide;

  var tr = new Matrix().setBasis(basis);
  var currentBounds = new BBox();
  var points = boundingPolygon.map(function(p) { p.z = depth; return tr._apply(p); });
  var polygon = new CSG.Polygon(points.map(function(p){return new CSG.Vertex(utils.csgVec(p))}), shared);
  var plane = new Solid(CSG.fromPolygons([polygon]), material, 'PLANE');
  plane.wireframeGroup.visible = false;
  plane.mergeable = false;
  var _3d = tr.invert();

  function setBounds(bbox) {
    var corner = new Vector(bbox.minX, bbox.minY, 0);
    var size = new Vector(bbox.width(), bbox.height(), 1);
    _3d._apply(size);
    _3d._apply(corner);
    plane.mesh.scale.set(size.x, size.y, size.z);
    plane.mesh.position.set(corner.x, corner.y, corner.z);
    currentBounds = bbox;
    var poly = new CSG.Polygon(bbox.toPolygon().map(function(p){return new CSG.Vertex(utils.csgVec( _3d._apply(p) ))}), shared);
    plane.csg = CSG.fromPolygons([poly]);
  }
  var bb = new BBox();
  bb.checkBounds(-400, -400);
  bb.checkBounds( 400,  400);
  setBounds(bb);

  var sketchFace = plane.polyFaces[0];
  tu.intercept(sketchFace, 'syncSketches', function(invocation, args) {
    var geom = args[0];
    invocation(geom);
    var bbox = new BBox();
    var connections = geom.connections.concat(utils.arrFlatten1L(geom.loops));
    for (var i = 0; i < connections.length; ++i) {
      var l = connections[i];
      bbox.checkBounds(l.a.x, l.a.y);
      bbox.checkBounds(l.b.x, l.b.y);
    }
    if (bbox.maxX > currentBounds.maxX || bbox.maxY > currentBounds.maxY || bbox.minX < currentBounds.minX || bbox.minY < currentBounds.minY) {
      bbox.expand(50);
      setBounds(bbox);
    }
  });

  return plane;
};


utils.fixCCW = function(path, normal) {
  var _2DTransformation = new Matrix().setBasis(geom.someBasis(path, normal)).invert();
  var path2D = [];
  for (var i = 0; i < path.length; ++i) {
    path2D[i] = _2DTransformation.apply(path[i]);
  }

  if (!geom.isCCW(path2D)) {
    path = path.slice(0);
    path.reverse();
  }
  return path;
};

const TOLERANCE = 1E-6;

utils.areEqual = function(v1, v2, tolerance) {
  return Math.abs(v1 - v2) < tolerance;
};

utils.areVectorsEqual = function(v1, v2, tolerance) {
  return utils.areEqual(v1.x, v2.x, tolerance) &&
      utils.areEqual(v1.y, v2.y, tolerance) &&
      utils.areEqual(v1.z, v2.z, tolerance);
};

utils.vectorsEqual = function(v1, v2) {
  return utils.areVectorsEqual(v1, v2, TOLERANCE);
};

utils.equal = function(v1, v2) {
  return utils.areEqual(v1, v2, TOLERANCE);
};

utils.strictEqual = function(a, b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
};

utils.isPointInsidePolygon = function( inPt, inPolygon ) {
  var EPSILON = TOLERANCE;

  var polyLen = inPolygon.length;

  // inPt on polygon contour => immediate success    or
  // toggling of inside/outside at every single! intersection point of an edge
  //  with the horizontal line through inPt, left of inPt
  //  not counting lowerY endpoints of edges and whole edges on that line
  var inside = false;
  for( var p = polyLen - 1, q = 0; q < polyLen; p = q ++ ) {
    var edgeLowPt  = inPolygon[ p ];
    var edgeHighPt = inPolygon[ q ];

    var edgeDx = edgeHighPt.x - edgeLowPt.x;
    var edgeDy = edgeHighPt.y - edgeLowPt.y;

    if ( Math.abs(edgeDy) > EPSILON ) {			// not parallel
      if ( edgeDy < 0 ) {
        edgeLowPt  = inPolygon[ q ]; edgeDx = - edgeDx;
        edgeHighPt = inPolygon[ p ]; edgeDy = - edgeDy;
      }
      if ( ( inPt.y < edgeLowPt.y ) || ( inPt.y > edgeHighPt.y ) ) 		continue;

      if ( inPt.y == edgeLowPt.y ) {
        if ( inPt.x == edgeLowPt.x )		return	true;		// inPt is on contour ?
        // continue;				// no intersection or edgeLowPt => doesn't count !!!
      } else {
        var perpEdge = edgeDy * (inPt.x - edgeLowPt.x) - edgeDx * (inPt.y - edgeLowPt.y);
        if ( perpEdge == 0 )				return	true;		// inPt is on contour ?
        if ( perpEdge < 0 ) 				continue;
        inside = ! inside;		// true intersection left of inPt
      }
    } else {		// parallel or colinear
      if ( inPt.y != edgeLowPt.y ) 		continue;			// parallel
      // egde lies on the same horizontal line as inPt
      if ( ( ( edgeHighPt.x <= inPt.x ) && ( inPt.x <= edgeLowPt.x ) ) ||
         ( ( edgeLowPt.x <= inPt.x ) && ( inPt.x <= edgeHighPt.x ) ) )		return	true;	// inPt: Point on contour !
      // continue;
    }
  }

  return	inside;
};

utils.sketchToPolygons = function(geom2) {

  var dict = struct.hashTable.forVector2d();
  var edges = struct.hashTable.forDoubleArray();

  var lines = geom2.connections;

  function edgeKey(a, b) {
    return [a.x, a.y, b.x, b.y];
  }

  var size = 0;
  var points = [];
  function memDir(a, b) {
    var dirs = dict.get(a);
    if (dirs === null) {
      dirs = [];
      dict.put(a, dirs);
      points.push(a);
    }
    dirs.push(b);
  }

  for (var i = 0; i < lines.length; i++) {
    var a = lines[i].a;
    var b = lines[i].b;

    memDir(a, b);
    memDir(b, a);
    edges.put(edgeKey(a, b), lines[i]);
  }

  var graph2 = {

    connections : function(e) {
      var dirs = dict.get(e);
      return dirs === null ? [] : dirs;
    },

    at : function(index) {
      return points[index];
    },

    size : function() {
      return points.length;
    }
  };

  var loops = graph.findAllLoops(graph2, dict.hashCodeF, dict.equalsF);
  var polygons = [];
  for (var li = 0; li < loops.length; ++li) {
    var loop = loops[li];
    if (!geom.isCCW(loop)) loop.reverse();
    var polyPoints = [];
    for (var pi = 0; pi < loop.length; ++pi) {
      var point = loop[pi];
      var next = loop[(pi + 1) % loop.length];

      var edge = edges.get(edgeKey(point, next));
      if (edge === null) {
        edge = edges.get(edgeKey(next, point));
      }
      polyPoints.push(point);
      point.sketchConnectionObject = edge.sketchObject;
    }
    if (polyPoints.length >= 3) {
      polygons.push(polyPoints);
    } else {
      console.warn("Points count < 3!");
    }
  }
  for (var li = 0; li < geom2.loops.length; ++li) {
    var loop = geom2.loops[li];
    var polyPoints = loop.slice(0);
    for (var si = 0; si < polyPoints.length; si++) {
      var conn = polyPoints[si];
      //reuse a point and ignore b point since it's a guaranteed loop
      conn.a.sketchConnectionObject = conn.sketchObject;
      polyPoints[si] = conn.a;
    }
    // we assume that connection object is the same al other the loop. That's why reverse is safe.
    if (!geom.isCCW(polyPoints)) polyPoints.reverse();
    if (polyPoints.length >= 3) {
      polygons.push(polyPoints);
    }
  }
  return polygons;
};

const geom = {};

geom.someBasis2 = function(normal) {
  var x = normal.cross(normal.randomNonParallelVector());
  var y = normal.cross(x).unit();
  return [x, y, normal];
};

geom.someBasis = function(twoPointsOnPlane, normal) {
  var a = twoPointsOnPlane[0];
  var b = twoPointsOnPlane[1];

  var x = b.minus(a).normalize();
  var y = normal.cross(x).normalize();

  return [x, y, normal];
};

geom.normalOfCCWSeq = function(ccwSequence) {
  var a = ccwSequence[0];
  var b = ccwSequence[1];
  var c = ccwSequence[2];

  return b.minus(a).cross(c.minus(a)).normalize();
};

geom.normalOfCCWSeqTHREE = function(ccwSequence) {
  var a = ccwSequence[0];
  var b = ccwSequence[1].clone();
  var c = ccwSequence[2].clone();

  return b.sub(a).cross(c.sub(a)).normalize();
};


// http://en.wikipedia.org/wiki/Shoelace_formula
geom.area = function (contour) {
  var n = contour.length;
  var a = 0.0;
  for ( var p = n - 1, q = 0; q < n; p = q ++ ) {
    a += contour[ p ].x * contour[ q ].y - contour[ q ].x * contour[ p ].y;
  }
  return a * 0.5;
};

geom.isCCW = function(path2D) {
  return geom.area(path2D) >= 0;
};

function BBox() {
  this.minX = Number.MAX_VALUE;
  this.minY = Number.MAX_VALUE;
  this.maxX = -Number.MAX_VALUE;
  this.maxY = -Number.MAX_VALUE;
  this.checkBounds = function(x, y) {
    this.minX = Math.min(this.minX, x);
    this.minY = Math.min(this.minY, y);
    this.maxX = Math.max(this.maxX, x);
    this.maxY = Math.max(this.maxY, y);
  };

  this.center = function() {
    return new Vector(this.minX + (this.maxX - this.minX) / 2, this.minY + (this.maxY - this.minY) / 2, 0)
  };

  this.width = function() {
    return this.maxX - this.minX;
  };

  this.height = function() {
    return this.maxY - this.minY;
  };

  this.expand = function(delta) {
    this.minX -= delta;
    this.minY -= delta;
    this.maxX += delta;
    this.maxY += delta;
  };

  this.toPolygon = function() {
    return [
      new Vector(this.minX, this.minY, 0),
      new Vector(this.maxX, this.minY, 0),
      new Vector(this.maxX, this.maxY, 0),
      new Vector(this.minX, this.maxY, 0)
    ];
  }
};

geom.calculateExtrudedLid = function(sourcePolygon, normal, direction, expansionFactor) {
  var lid = [];
  var length = sourcePolygon.length;
  var work;
  var si;
  if (!!expansionFactor && expansionFactor != 1) {
    var source2d = [];
    work = [];

    var _3dTr = new Matrix().setBasis(geom.someBasis2(new CSG.Vector3D(normal))); // use passed basis
    var _2dTr = _3dTr.invert();
    var sourceBBox = new BBox();
    var workBBox = new BBox();
    for (si = 0; si < length; ++si) {
      var sourcePoint = _2dTr.apply(sourcePolygon[si]);
      source2d[si] = sourcePoint;
      work[si] = sourcePoint.multiply(expansionFactor);
      work[si].z = source2d[si].z = 0;
      sourceBBox.checkBounds(sourcePoint.x, sourcePoint.y);
      workBBox.checkBounds(work[si].x, work[si].y)
    }
    var alignVector = workBBox.center().minus(sourceBBox.center());
    var depth = normal.dot(sourcePolygon[0]);
    for (si = 0; si < length; ++si) {
      work[si] = work[si].minus(alignVector);
      work[si].z = depth;
      work[si] = _3dTr.apply(work[si]);
    }
  } else {
    work = sourcePolygon;
  }

  for (si = 0; si < length; ++si) {
    lid[si] = work[si].plus(direction);
  }

  return lid;
};

geom.extrude = function(source, sourceNormal, target, expansionFactor) {

  var extrudeDistance = target.normalize().dot(sourceNormal);
  if (extrudeDistance == 0) {
    return [];
  }
  var negate = extrudeDistance < 0;

  var poly = [null, null];
  var lid = geom.calculateExtrudedLid(source, sourceNormal, target, expansionFactor);

  var bottom, top;
  if (negate) {
    bottom = lid;
    top = source;
  } else {
    bottom = source;
    top = lid;
  }

  var n = source.length;
  for ( var p = n - 1, i = 0; i < n; p = i ++ ) {
    var shared = utils.createShared();
    shared.__tcad.csgInfo = {derivedFrom:  source[p].sketchConnectionObject};
    var face = new CSG.Polygon([
      new CSG.Vertex(utils.csgVec(bottom[p])),
      new CSG.Vertex(utils.csgVec(bottom[i])),
      new CSG.Vertex(utils.csgVec(top[i])),
      new CSG.Vertex(utils.csgVec(top[p]))
    ], shared);
    poly.push(face);
  }

  var bottomNormal, topNormal;
  if (negate) {
    lid.reverse();
    bottomNormal = sourceNormal;
    topNormal = sourceNormal.negate();
  } else {
    source = source.slice(0);
    source.reverse();
    bottomNormal = sourceNormal.negate();
    topNormal = sourceNormal;
  }

  function vecToVertex(v) {
    return new CSG.Vertex(utils.csgVec(v));
  }

  var sourcePlane = new CSG.Plane(bottomNormal.csg(), bottomNormal.dot(source[0]));
  var lidPlane = new CSG.Plane(topNormal.csg(), topNormal.dot(lid[0]));

  poly[0] = new CSG.Polygon(source.map(vecToVertex), utils.createShared(), sourcePlane);
  poly[1] = new CSG.Polygon(lid.map(vecToVertex), utils.createShared(), lidPlane);
  return poly;
};

geom.SOLID_COUNTER = 0;

geom.triangulate = function(path, normal) {
  var _3dTransformation = new Matrix().setBasis(geom.someBasis2(normal));
  var _2dTransformation = _3dTransformation.invert();
  var i;
  var shell = [];
  for (i = 0; i < path.length; ++i) {
    shell[i] = _2dTransformation.apply(path[i].pos);
  }
  var myTriangulator = new PNLTRI.Triangulator();
  return  myTriangulator.triangulate_polygon( [ shell ] );
//  return Shape.utils.triangulateShape( f2d.shell, f2d.holes );
};

utils.groupCSG = function(csg) {
  var csgPolygons = csg.toPolygons();
  var groups = {};
  for (var i = 0; i < csgPolygons.length; i++) {
    var p = csgPolygons[i];
    var tag = p.shared.getTag();
    if (groups[tag] === undefined) {
      groups[tag] = {
        tag : tag,
        polygons : [],
        shared : p.shared,
        plane : p.plane
      };
    }
    groups[tag].polygons.push(p);
  }
  return groups;
};

utils.SHARED_COUNTER = 0;
utils.createShared = function() {
  var id = utils.SHARED_COUNTER ++;
  var shared = new CSG.Polygon.Shared([id, id, id, id]);
  shared.__tcad = {};
  return shared;
};

utils.isSmoothPiece = function(shared) {
  return shared.__tcad && !!shared.__tcad.csgInfo && !!shared.__tcad.csgInfo.derivedFrom &&
  (shared.__tcad.csgInfo.derivedFrom._class === 'TCAD.TWO.Arc' || shared.__tcad.csgInfo.derivedFrom._class === 'TCAD.TWO.Circle');
};

utils.sameID = function(id1, id2) {
  if (id1 === null || id2 === null) {
    return false;
  }
  return id1 === id2;
};

utils.getDerivedID = function(shared) {
  return shared.__tcad && !!shared.__tcad.csgInfo && !!shared.__tcad.csgInfo.derivedFrom ? shared.__tcad.csgInfo.derivedFrom.id : null;
};

utils.getDerivedFrom = function(shared) {
  return shared.__tcad && !!shared.__tcad.csgInfo && !!shared.__tcad.csgInfo.derivedFrom ? shared.__tcad.csgInfo.derivedFrom : null;
};

/** @constructor */
function Solid(csg, material, type) {
  csg = csg.reTesselated().canonicalized();
  this.tCadType = type || 'SOLID';
  this.csg = csg;

  this.cadGroup = new Object3D();
  this.cadGroup.__tcad_solid = this;

  var geometry = new Geometry();
  geometry.dynamic = true;
  this.mesh = new Mesh(geometry, material);
  this.cadGroup.add(this.mesh);

  this.tCadId = geom.SOLID_COUNTER ++;
  this.faceCounter = 0;

  this.wireframeGroup = new Object3D();
  this.cadGroup.add(this.wireframeGroup);

  this.polyFaces = [];
  this.wires = struct.hashTable.forEdge();
  this.curvedSurfaces = {};
  this.mergeable = true;

  this.setupGeometry();
};

Solid.prototype.setupGeometry = function() {
  function threeV(v) {return new Vector3( v.x, v.y, v.z )}

  var off = 0;
  var groups = utils.groupCSG(this.csg);
  var geom = this.mesh.geometry;
  for (var gIdx in groups)  {
    var group = groups[gIdx];
    if (group.shared.__tcad === undefined) group.shared.__tcad = {};
    var polyFace = new SketchFace(this, group);
    this.polyFaces.push(polyFace);
    for (var p = 0; p < group.polygons.length; ++p) {
      var poly = group.polygons[p];
      var vLength = poly.vertices.length;
      if (vLength < 3) continue;
      var firstVertex = poly.vertices[0];
      geom.vertices.push(threeV(firstVertex.pos));
      geom.vertices.push(threeV(poly.vertices[1].pos));
      var normal = threeV(poly.plane.normal);
      for (var i = 2; i < vLength; i++) {
        geom.vertices.push(threeV(poly.vertices[i].pos));

        var a = off;
        var b = i - 1 + off;
        var c = i + off;
        var face = new Face3(a, b, c);
        polyFace.faces.push(face);
        face.__TCAD_polyFace = polyFace;
        face.normal = normal;
        face.materialIndex = gIdx;
        geom.faces.push(face);
        //face.color.set(new Color().setRGB( Math.random(), Math.random(), Math.random()));
      }
      //TCAD.view.setFaceColor(polyFace, utils.isSmoothPiece(group.shared) ? 0xFF0000 : null);
      off = geom.vertices.length;
    }
    this.collectCurvedSurface(polyFace);
    this.collectWires(polyFace);
  }

  geom.mergeVertices();

  this.processWires();
};

Solid.prototype.vanish = function() {
  this.cadGroup.parent.remove( this.cadGroup );
  this.mesh.material.dispose();
  this.mesh.geometry.dispose();
};

Solid.prototype.collectCurvedSurface = function(face) {
  var derivedFrom = utils.getDerivedFrom(face.csgGroup.shared);
  if (derivedFrom === null || derivedFrom._class !== "TCAD.TWO.Arc" && derivedFrom._class !== "TCAD.TWO.Circle" ) return;
  var surfaces = this.curvedSurfaces[derivedFrom.id];
  if (surfaces === undefined) {
    surfaces = [];
    this.curvedSurfaces[derivedFrom.id] = surfaces;
  }
  surfaces.push(face);
  face.curvedSurfaces = surfaces;
};

Solid.prototype.collectWires = function(face) {

  function contains(planes, plane) {
    for (var j = 0; j < planes.length; j++) {
      if (planes[j].equals(plane)) {
        return true;
      }
    }
    return false;
  }
  var paths = craft.reconstructSketchBounds(this.csg, face, true);
  for (var i = 0; i < paths.length; i++) {
    var path = paths[i];
    var p, q, n = path.vertices.length;
    for (q = 0, p = n - 1; q < n; p = q++) {
      var edge = [path.vertices[p], path.vertices[q]];
      var data = this.wires.get(edge);

      if (data === null) {
        data = {
          sharedPlanes : [face.csgGroup.plane],
          sharedFaces : [face]
        };
        this.wires.put(edge, data);
      } else {
        if (!contains(data.sharedPlanes, face.csgGroup.plane)) {
          data.sharedPlanes.push(face.csgGroup.plane);
        }
        data.sharedFaces.push(face);
      }
    }
  }
};

Solid.SMOOTH_LIMIT = 10 * Math.PI / 180;

Solid.prototype.processWires = function() {
  var solid = this;
  this.wires.entries(function(edge, data) {
    var u = utils;
    if (data.sharedPlanes.length > 1) {
      var plane0 = data.sharedPlanes[0];
      var plane1 = data.sharedPlanes[1];
      var angle = Math.acos(plane0.normal.dot(plane1.normal));
      if (angle < Solid.SMOOTH_LIMIT) {
        return;
      }
    }
    for (var i = 0; i < data.sharedFaces.length; ++i) {
      for (var j = i + 1; j < data.sharedFaces.length; ++j) {
        var face0 = data.sharedFaces[0];
        var face1 = data.sharedFaces[1];
        if (u.sameID(u.getDerivedID(face0.csgGroup.shared), u.getDerivedID(face1.csgGroup.shared))) {
          return;
        }
      }
    }

    solid.addLineToScene(edge[0], edge[1]);
  });
};

Solid.prototype.addLineToScene = function(a, b) {
  var lg = new Geometry();
  lg.vertices.push(a);
  lg.vertices.push(b);
  var line = new Line(lg, SketchFace.prototype.WIREFRAME_MATERIAL);
  this.wireframeGroup.add(line);
};

/** @constructor */
function SketchFace(solid, csgGroup) {
  csgGroup.__face = this;
  if (csgGroup.shared.__tcad.faceId === undefined) {
    this.id = solid.tCadId + ":" + (solid.faceCounter++);
  } else {
    this.id = csgGroup.shared.__tcad.faceId;
  }
  csgGroup.shared.__tcad.faceId = this.id;

  this.solid = solid;
  this.csgGroup = csgGroup;
  this.faces = [];
  this.sketch3DGroup = null;
  this.curvedSurfaces = null;
};

SketchFace.prototype.SKETCH_MATERIAL = new LineBasicMaterial({
  color: 0xFFFFFF, linewidth: 3/DPR});
SketchFace.prototype.WIREFRAME_MATERIAL = new LineBasicMaterial({
  color: 0x2B3856, linewidth: 3/DPR});

SketchFace.prototype.calcBasis = function() {
  var vec = utils.vec;
  var normal = vec(this.csgGroup.plane.normal);
  var alignPlane, x, y;
  if (Math.abs(normal.dot(math.AXIS.Y)) < 0.5) {
    alignPlane = normal.cross(math.AXIS.Y);
  } else {
    alignPlane = normal.cross(math.AXIS.Z);
  }
  y = alignPlane.cross(normal);
  x = y.cross(normal);
  return [x, y, normal];
};

SketchFace.prototype.basis = function() {
  if (!this._basis) {
    this._basis = this.calcBasis();
  }
  return this._basis;
  //return geom.someBasis(this.csgGroup.polygons[0].vertices.map(function (v) {
  //  return vec(v.pos)
  //}), vec(this.csgGroup.plane.normal));
};

SketchFace.prototype.depth = function() {
  return this.csgGroup.plane.w;
};

SketchFace.prototype.syncSketches = function(geom) {
  var i;
  var normal = this.csgGroup.plane.normal;
  var offVector = normal.scale(0); // disable it. use polygon offset feature of material

  if (this.sketch3DGroup != null) {
    for (i = this.sketch3DGroup.children.length - 1; i >= 0; --i) {
      this.sketch3DGroup.remove(this.sketch3DGroup.children[i]);
    }
  } else {
    this.sketch3DGroup = new Object3D();
    this.solid.cadGroup.add(this.sketch3DGroup);
  }

  var basis = this.basis();
  var _3dTransformation = new Matrix().setBasis(basis);
  //we lost depth or z off in 2d sketch, calculate it again
  var depth = this.csgGroup.plane.w;
  var connections = geom.connections.concat(utils.arrFlatten1L(geom.loops));
  for (i = 0; i < connections.length; ++i) {
    var l = connections[i];
    var lg = new Geometry();
    l.a.z = l.b.z = depth;
    var a = _3dTransformation.apply(l.a);
    var b = _3dTransformation.apply(l.b);

    lg.vertices.push(a.plus(offVector).three());
    lg.vertices.push(b.plus(offVector).three());
    var line = new Line(lg, this.SKETCH_MATERIAL);
    this.sketch3DGroup.add(line);
  }
};

let POLYGON_COUNTER = 0;

/** @constructor */
function Polygon(shell, holes, normal) {
  this.id = POLYGON_COUNTER ++;
  if (!holes) {
    holes = [];
  }
  utils.checkPolygon(shell);
  for (var h = 0; h < holes.length; ++h) {
    utils.checkPolygon(holes[h]);
  }

  if (normal === undefined) {
    normal = geom.normalOfCCWSeq(shell);
  } else {
    shell = utils.fixCCW(shell, normal);
    if (holes.length > 0) {
      var neg = normal.negate();
      for (var h = 0; h < holes.length; ++h) {
        holes[h] = utils.fixCCW(holes[h], neg);
      }
    }

  }

  this.normal = normal;
  this.shell = shell;
  this.holes = holes;
};

Polygon.prototype.reverse = function(triangle) {
  var first = triangle[0];
  triangle[0] = triangle[2];
  triangle[2] = first;
};

Polygon.prototype.flip = function() {
  return new Polygon(this.shell, this.holes, this.normal.negate());
};

Polygon.prototype.shift = function(target) {
  var shell = [];
  var i;
  for (i = 0; i < this.shell.length; ++i) {
    shell[i] = this.shell[i].plus(target);
  }
  var holes = [];
  for (var h = 0; h < this.holes.length; ++h) {
    holes[h] = [];
    for (i = 0; i < this.holes[h].length; ++i) {
      holes[h][i] = this.holes[h][i].plus(target);
    }
  }
  return new Polygon(shell, holes, this.normal);
};

Polygon.prototype.get2DTransformation = function() {
  var _3dTransformation = new Matrix().setBasis(geom.someBasis(this.shell, this.normal));
  var _2dTransformation = _3dTransformation.invert();
  return _2dTransformation;
};

Polygon.prototype.to2D = function() {

  var _2dTransformation = this.get2DTransformation();

  var i, h;
  var shell = [];
  var holes = [];
  for (i = 0; i < this.shell.length; ++i) {
    shell[i] = _2dTransformation.apply(this.shell[i]);
  }
  for (h = 0; h < this.holes.length; ++h) {
    holes[h] = [];
    for (i = 0; i < this.holes[h].length; ++i) {
      holes[h][i] = _2dTransformation.apply(this.holes[h][i]);
    }
  }
  return {shell: shell, holes: holes};
};

Polygon.prototype.collectPaths = function(paths) {
  paths.push(this.shell);
  paths.push.apply(paths, this.holes);
};

Polygon.prototype.triangulate = function() {

  function triangulateShape( contour, holes ) {
    var myTriangulator = new PNLTRI.Triangulator();
    return  myTriangulator.triangulate_polygon( [ contour ].concat(holes) );
  }

  var i, h;
  var f2d = this.to2D();

  for (i = 0; i < f2d.shell.length; ++i) {
    f2d.shell[i] = f2d.shell[i].three();
  }
  for (h = 0; h < f2d.holes.length; ++h) {
    for (i = 0; i < f2d.holes[h].length; ++i) {
      f2d.holes[h][i] = f2d.holes[h][i].three();
    }
  }
  return triangulateShape( f2d.shell, f2d.holes );
//  return Shape.utils.triangulateShape( f2d.shell, f2d.holes );
};

Polygon.prototype.eachVertex = function(handler) {
  var i, h;
  for (i = 0; i < this.shell.length; ++i) {
    if (handler(this.shell, i) === true) return;
  }
  for (h = 0; h < this.holes.length; ++h) {
    for (i = 0; i < this.holes[h].length; ++i) {
      if (handler(this.holes[h], i) === true) return;
    }
  }
};

/** @constructor */
function Sketch() {
  this.group = new Object3D();
};

utils.iteratePath = function(path, shift, callback) {
  var p, q, n = path.length;
  for (p = n - 1,q = 0;q < n; p = q++) {
    var ai = (p + shift) % n;
    var bi = (q + shift) % n;
    if (!callback(path[ai], path[bi], ai, bi, q, path)) {
      break
    }
  }
};

utils.addAll = function(arr, arrToAdd) {
  for (var i = 0; i < arrToAdd.length; i++) {
    arr.push(arrToAdd[i]);
  }
};

utils.arrFlatten1L = function(arr) {
  var result = [];
  for (var i = 0; i < arr.length; i++) {
    utils.addAll(result, arr[i]);
  }
  return result;
};

export { utils, geom, Sketch, Solid, BBox, Polygon, TOLERANCE };
