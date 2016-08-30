import { DPR, view } from '../viewer';
import { LineBasicMaterial, Object3D, Vector3, Geometry, Line, PlaneGeometry, MeshLambertMaterial, Mesh, DoubleSide } from 'three';
import math from '../../math/math';
import { geom } from '../../engine';

const wizards = {};

wizards.IMAGINE_MATERIAL = new LineBasicMaterial({
  color: 0xFA8072,
  linewidth: 1/DPR,
  depthWrite: false,
  depthTest: false
});

wizards.BASE_MATERIAL = new LineBasicMaterial({
  color: 0x8B0000,
  linewidth: 3/DPR,
  depthWrite: false,
  depthTest: false
});

wizards.OpWizard = function(viewer) {
  this.previewGroup = new Object3D();
  this.lines = [];
  this.viewer = viewer;
  viewer.scene.add(this.previewGroup);
};

wizards.OpWizard.prototype.setupLine = function(lineId, a, b, material) {
  var line = this.lines[lineId];
  if (line === undefined) {
    var lg = new Geometry();
    lg.vertices.push(new Vector3().copy(a));
    lg.vertices.push(new Vector3().copy(b));
    line = new Line(lg, material);
    line.renderOrder = 1e10;
    this.previewGroup.add(line);
    this.lines[lineId] = line;
  } else {
    line.geometry.vertices[0] = new Vector3().copy(a);
    line.geometry.vertices[1] = new Vector3().copy(b);
    line.geometry.verticesNeedUpdate = true;
  }
};

wizards.OpWizard.prototype.dispose = function() {
  this.viewer.scene.remove(this.previewGroup);
  this.viewer.render();
};

wizards.ExtrudeWizard = function(viewer, polygons) {
  wizards.OpWizard.call(this, viewer);
  this.polygons = polygons;
};

wizards.ExtrudeWizard.prototype = Object.create( wizards.OpWizard.prototype );

wizards.ExtrudeWizard.prototype.update = function(basis, normal, depth, scale, deflection, angle) {
  var linesCounter = 0;
  var target;
  if (deflection != 0) {
    target = normal.copy();
    if (depth < 0) target._negate();
    target = math.rotateMatrix(deflection * Math.PI / 180, basis[0], math.ORIGIN)._apply(target);
    if (angle != 0) {
      target = math.rotateMatrix(angle * Math.PI / 180, basis[2], math.ORIGIN)._apply(target);
    }
    target._multiply(Math.abs(depth));
  } else {
    target = normal.multiply(depth)
  }
  for (var i = 0; i < this.polygons.length; i++) {
    var poly = this.polygons[i];
    var lid = geom.calculateExtrudedLid(poly, normal, target, scale);
    var p, q, n = poly.length;
    for (p = n - 1, q = 0; q < n; p = q++) {
      this.setupLine(linesCounter ++, poly[p], poly[q], wizards.BASE_MATERIAL);
      this.setupLine(linesCounter ++, lid[p], lid[q], wizards.IMAGINE_MATERIAL);
    }
    for (q = 0; q < n; q++) {
      this.setupLine(linesCounter ++, poly[q], lid[q], wizards.IMAGINE_MATERIAL);
    }
  }
  this.operationParams = {
    target : target,
    expansionFactor : scale
  }
};


wizards.PlaneWizard = function(viewer) {
  this.previewGroup = new Object3D();
  this.viewer = viewer;
  viewer.scene.add(this.previewGroup);
  this.previewGroup.add(this.plane = this.createPlane());
  this.viewer.render();
  this.operationParams = {
    basis : math.IDENTITY_BASIS,
    depth : 0
  };
};

wizards.PlaneWizard.prototype.createPlane = function() {
  var geometry = new PlaneGeometry(750,750,1,1,1);
  var material = new MeshLambertMaterial( { color : view.FACE_COLOR, transparent: true, opacity:0.5, side: DoubleSide });
  var plane = new Mesh(geometry, material);
  return plane;
};

wizards.PlaneWizard.prototype.update = function(orientation, w) {
  if (orientation === 'XY') {
    this.plane.rotation.x = 0;
    this.plane.rotation.y = 0;
    this.plane.rotation.z = 0;
    this.plane.position.x = 0;
    this.plane.position.y = 0;
    this.plane.position.z = w;
    this.operationParams.basis = math.IDENTITY_BASIS;
  } else if (orientation === 'XZ') {
    this.plane.rotation.x = Math.PI / 2;
    this.plane.rotation.y = 0;
    this.plane.rotation.z = 0;
    this.plane.position.x = 0;
    this.plane.position.y = w;
    this.plane.position.z = 0;
    this.operationParams.basis = [math.AXIS.X, math.AXIS.Z, math.AXIS.Y];
  } else if (orientation === 'ZY') {
    this.plane.rotation.x = 0;
    this.plane.rotation.y = Math.PI / 2;
    this.plane.rotation.z = 0;
    this.plane.position.x = w;
    this.plane.position.y = 0;
    this.plane.position.z = 0;
    this.operationParams.basis = [math.AXIS.Z, math.AXIS.Y, math.AXIS.X];
  } else {
    throw orientation + " isn't supported yet";
  }
  this.operationParams.depth = w;
  this.viewer.render();
};

wizards.PlaneWizard.prototype.dispose = function() {
  this.viewer.scene.remove(this.previewGroup);
  this.viewer.render();
};

export default wizards;
