import parametric from './constr/solver';
import constraints from './constr/constraints';
import { Vector } from '../math/math';
import { Ref } from './canvas';
import math from '../math/math';

const EQUALS_ELIMINATION_ENABLED = true;
const utils = {};
const Constraints = {};

/** @constructor */
function SubSystem() {
  this.alg = 1;
  this.error = 0;
  this.reduce = false;
  this.constraints = [];
};

/** @constructor */
function ParametricManager(viewer) {
  this.viewer = viewer;
  this.subSystems = [];
  this.listeners = [];
  this.constantTable = {};

  this.viewer.params.define("constantDefinition", null);
  this.viewer.params.subscribe("constantDefinition", "parametricManager", this.rebuildConstantTable, this)();
  this.constantResolver = this.createConstantResolver();
};

ParametricManager.prototype.createConstantResolver = function() {
  var pm = this;
  return function(value) {
    var _value = pm.constantTable[value];
    if (_value !== undefined) {
      value = _value;
    } else if (typeof(value) != 'number') {
      console.error("unable to resolve constant " + value);
    }
    return value;
  }
};

ParametricManager.prototype.notify = function(event) {
  for (var i = 0; i < this.listeners.length; ++i) {
    var l = this.listeners[i];
    l(event);
  }
};

ParametricManager.prototype.rebuildConstantTable = function(constantDefinition) {
  this.constantTable = {};
  if (constantDefinition == null) return;
  var lines = constantDefinition.split('\n');
  var prefix = "(function() { \n";
  for (var i = 0; i < lines.length; i++) {
    var line = lines[i];
    var m = line.match(/^\s*([^\s]+)\s*=(.+)$/);
    if (m != null && m.length == 3) {
      var constant = m[1];
      try {
        var value = eval(prefix + "return " + m[2] + "; \n})()");
        this.constantTable[constant] = value;
        prefix += constant + " = " + value + ";\n"
      } catch(e) {
        console.log(e);
      }
    }
  }
  this.refresh();
};


ParametricManager.prototype.findComponents = function(constr) {
  if (this.subSystems.length === 0) {
    this.subSystems.push(new SubSystem());
  }
  return [0];
};

ParametricManager.prototype.tune = function(subSystem) {

};

ParametricManager.prototype._add = function(constr) {
  var subSystemIds = this.findComponents(constr);
  var subSystem;
  switch (subSystemIds.length) {
    case 0:
      subSystem = new SubSystem();
      this.subSystems.push(subSystem);
      break;
    case 1:
      subSystem = this.subSystems[subSystemIds[0]];
      break;
    default:
      subSystem = this.subSystems[subSystemIds[0]];
      for (var i = 1; i < subSystemIds.length; i++) {
        var toMerge = subSystemIds[i];
        for (var j = 0; j < toMerge.constraints.length; j++) {
          subSystem.push(toMerge.constraints[j]);
        }
      }
    break;
  }
  subSystem.constraints.push(constr);
  return subSystem;
};

ParametricManager.prototype.checkRedundancy = function (subSystem, constr) {
  var solver = this.prepareForSubSystem([], subSystem.constraints);
  if (solver.diagnose().conflict) {
    alert("Most likely this "+constr.NAME + " constraint is CONFLICTING!")
  }
};

ParametricManager.prototype.refresh = function() {
  this.solve();
  this.notify();
  this.viewer.refresh();
};

ParametricManager.prototype.add = function(constr) {
  this.viewer.historyManager.checkpoint();
  var subSystem = this._add(constr);
  this.checkRedundancy(subSystem, constr);
  this.refresh();
};

ParametricManager.prototype.addAll = function(constrs) {
  for (var i = 0; i < constrs.length; i++) {
    var subSystem = this._add(constrs[i]);
    this.checkRedundancy(subSystem, constrs[i]);
  }
  this.refresh();
};

ParametricManager.prototype.remove = function(constr) {
  this.viewer.historyManager.checkpoint();
  for (var j = 0; j < this.subSystems.length; j++) {
    var sub = this.subSystems[j];
    for (var i = 0; i < sub.constraints.length; ++i) {
      var p = sub.constraints[i];
      if (p === constr) {
        sub.constraints.splice(i, 1);
        if (p.NAME === 'coi') {
          this.unlinkObjects(p.a, p.b);
        }
        break;
      }
    }
  }
  this.refresh();
};

ParametricManager.prototype.removeConstraintsByObj = function(obj) {
  var ownedParams = [];
  obj.collectParams(ownedParams);
  this.removeConstraintsByParams(ownedParams);
};

ParametricManager.prototype.removeConstraintsByParams = function(ownedParams) {
  for (var s = 0; s < this.subSystems.length; s++) {
    var toRemove = [];
    var sub = this.subSystems[s];
    for (var i = 0; i < sub.constraints.length; ++i) {
      var sdataArr = sub.constraints[i].getSolveData(this.constantResolver);
      MAIN:
      for (var j = 0; j < sdataArr.length; j++) {
        var sdata = sdataArr[j];
        var params = sdata[1];
        for (var oi = 0; oi < ownedParams.length; ++oi) {
          for (var k = 0; k < params.length; ++k) {
            if (ownedParams[oi].id === params[k].id) {
              toRemove.push(i);
              break MAIN;
            }
          }
        }
      }
    }
    toRemove.sort();

    for (var i = toRemove.length - 1; i >= 0 ; --i) {
      sub.constraints.splice(  toRemove[i], 1);
    }
  }

  this.notify();
};

ParametricManager.prototype.lock = function(objs) {
  var p = this._fetchPoints(objs);
  for (var i = 0; i < p.length; ++i) {
    this._add(new Constraints.Lock(p[i], { x : p[i].x, y : p[i].y} ));
  }
  this.refresh();
};

ParametricManager.prototype.vertical = function(objs) {
  this.add(new Constraints.Vertical(this._fetchLine(objs)));
};

ParametricManager.prototype.horizontal = function(objs) {
  this.add(new Constraints.Horizontal(this._fetchLine(objs)));
};

ParametricManager.prototype.parallel = function(objs) {
  var lines = this._fetchTwoLines(objs);
  this.add(new Constraints.Parallel(lines[0], lines[1]));
};

ParametricManager.prototype.perpendicular = function(objs) {
  var lines = this._fetchTwoLines(objs);
  this.add(new Constraints.Perpendicular(lines[0], lines[1]));
};

ParametricManager.prototype.lockConvex = function(objs, warnCallback) {
  var lines = this._fetchTwoLines(objs);
  var l1 = lines[0];
  var l2 = lines[1];
  var pts =[l1.a, l1.b, l2.a, l2.b];
  function isLinked(p1, p2) {
    for (var i = 0; i < p1.linked.length; ++i) {
      if (p1.linked[i].id === p2.id) {
        return true;
      }
    }
    return false;
  }

  function swap(arr, i1, i2) {
    var _ = arr[i1];
    arr[i1] = arr[i2];
    arr[i2] = _;
  }

  if (isLinked(pts[0], pts[2])) {
    swap(pts, 0, 1);
  } else if (isLinked(pts[0], pts[3])) {
    swap(pts, 0, 1);
    swap(pts, 2, 3);
  } else if (isLinked(pts[1], pts[3])) {
    swap(pts, 2, 3);
  } else if (isLinked(pts[1], pts[2])) {
    //we are good
  } else {
    warnCallback("Lines must be connected");
    return;
  }

  var c = pts[0];
  var a = pts[1];
  var t = pts[3];

  // ||ac x at|| > 0
  var crossNorma = (c.x - a.x) * (t.y - a.y) - (c.y - a.y) * (t.x - a.x);

  if (crossNorma < 0) {
    var _ =  c;
    c = t;
    t = _;
  }

  this.add(new Constraints.LockConvex(c, a, t));
};

ParametricManager.prototype.tangent = function(objs) {
  var al = this._fetchArcCircAndLine(objs);
  var arc  = al[0];
  var line  = al[1];
  this.add(new Constraints.Tangent( arc, line));
};

ParametricManager.prototype.rr = function(arcs) {
  var prev = arcs[0];
  for (var i = 1; i < arcs.length; ++i) {
    this._add(new Constraints.RR(prev, arcs[i]));
    prev = arcs[i];
  }
  this.refresh();
};

ParametricManager.prototype.ll = function(lines) {
  var prev = lines[0];
  for (var i = 1; i < lines.length; ++i) {
    this._add(new Constraints.LL(prev, lines[i]));
    prev = lines[i];
  }
  this.refresh();

};

ParametricManager.prototype.entityEquality = function(objs) {
  var arcs = this._fetch(objs, ['TCAD.TWO.Arc', 'TCAD.TWO.Circle'], 0);
  var lines = this._fetch(objs, ['TCAD.TWO.Segment'], 0);
  if (arcs.length > 0) this.rr(arcs);
  if (lines.length > 0) this.ll(lines);
};

ParametricManager.prototype.p2lDistance = function(objs, promptCallback) {
  var pl = this._fetchPointAndLine(objs);

  var target = pl[0];
  var segment = pl[1];

  var ex = new Vector(-(segment.b.y - segment.a.y), segment.b.x - segment.a.x).normalize();
  var distance = Math.abs(ex.dot(new Vector(segment.a.x - target.x, segment.a.y - target.y)));

  var promptDistance = utils.askNumber(Constraints.P2LDistance.prototype.SettableFields.d, distance.toFixed(2), promptCallback, this.constantResolver);

  if (promptDistance != null) {
    this.add(new Constraints.P2LDistance(target, segment, promptDistance));
  }
};

ParametricManager.prototype.pointInMiddle = function(objs) {
  var pl = this._fetchPointAndLine(objs);
  this.add(new Constraints.PointInMiddle(pl[0], pl[1]));
};

ParametricManager.prototype.symmetry = function(objs) {
  var pl = this._fetchPointAndLine(objs);
  this.add(new Constraints.Symmetry(pl[0], pl[1]));
};

ParametricManager.prototype.pointOnArc = function(objs) {
  var points = this._fetch(objs, ['TCAD.TWO.EndPoint'], 1);
  var arcs = this._fetch(objs, ['TCAD.TWO.Arc', 'TCAD.TWO.Circle'], 1);
  this.add(new Constraints.PointOnArc(points[0], arcs[0]));
};

ParametricManager.prototype.pointOnLine = function(objs) {
  var pl = this._fetchPointAndLine(objs);
  var target = pl[0];
  var segment = pl[1];
  this.add(new Constraints.PointOnLine(target, segment));
};

ParametricManager.prototype.llAngle = function(objs, promptCallback) {
  var lines = this._fetch(objs, 'TCAD.TWO.Segment', 2);
  var l1 = lines[0];
  var l2 = lines[1];

  var points = [l1.a, l1.b, l2.a, l2.b];

  if (l1.b.x < l1.a.x) {
    points[0] = l1.b;
    points[1] = l1.a;
  }

  if (l2.b.x < l2.a.x) {
    points[2] = l2.b;
    points[3] = l2.a;
  }

  var dx1 = points[1].x - points[0].x;
  var dy1 = points[1].y - points[0].y;
  var dx2 = points[3].x - points[2].x;
  var dy2 = points[3].y - points[2].y;

  var angle = Math.atan2(dy2,dx2) - Math.atan2(dy1,dx1);
  angle *= 1 / Math.PI * 180;
  angle = utils.askNumber(Constraints.Angle.prototype.SettableFields.angle, angle.toFixed(2), promptCallback, this.constantResolver);
  if (angle === null) return;
  this.add(new Constraints.Angle(points[0], points[1], points[2], points[3], angle));
};

utils.constRef = function(value) {
  return function() {
    return value;
  };
};

ParametricManager.prototype.p2pDistance = function(objs, promptCallback) {
  var p = this._fetchTwoPoints(objs);
  var distance = new Vector(p[1].x - p[0].x, p[1].y - p[0].y).length();
  var promptDistance = utils.askNumber(Constraints.P2PDistance.prototype.SettableFields.d, distance.toFixed(2), promptCallback, this.constantResolver);

  if (promptDistance != null) {
    this.add(new Constraints.P2PDistance(p[0], p[1], promptDistance));
  }
};

utils.askNumber = function(promptText, initValue, promptCallback, resolver) {
  var promptValueStr = promptCallback(promptText, initValue);
  if (promptValueStr != null) {
    var promptValue = Number(promptValueStr);
    if (promptValue == promptValue) { // check for NaN
      return promptValue;
    } else {
      if (!!resolver) {
        promptValue = resolver(promptValueStr);
        if (promptValue == promptValue) {
          return promptValueStr;
        }
      }
    }
  }
  return null;
};

ParametricManager.prototype.radius = function(objs, promptCallback) {
  var arcs = this._fetchArkCirc(objs, 1);
  var radius = arcs[0].r.get();
  var promptDistance = utils.askNumber(Constraints.Radius.prototype.SettableFields.d, radius.toFixed(2), promptCallback, this.constantResolver);
  if (promptDistance != null) {
    for (var i = 0; i < arcs.length; ++i) {
      this._add(new Constraints.Radius(arcs[i], promptDistance));
    }
    this.refresh();
  }
};

ParametricManager.prototype.linkObjects = function(objs) {
  var i;
  var masterIdx = -1;
  for (i = 0; i < objs.length; ++i) {
    if (ParametricManager.isAux(objs[i])) {
      if (masterIdx !== -1) {
        throw "not allowed to have a coincident constraint between two or more auxiliary objects";
      }
      masterIdx = i;
    }
  }
  if (masterIdx === -1) masterIdx = objs.length - 1;


  for (i = 0; i < objs.length; ++i) {
    if (i === masterIdx) continue;
    objs[i].x = objs[masterIdx].x;
    objs[i].y = objs[masterIdx].y;
    var c = new Constraints.Coincident(objs[i], objs[masterIdx]);
    this._add(c);
  }
  this.notify();
};

ParametricManager.prototype.unlinkObjects = function(a, b) {

  function _unlink(a, b) {
    for (var i = 0; i < a.linked.length; ++i) {
      var obj = a.linked[i];
      if (obj.id === b.id) {
        a.linked.splice(i, 1);
        break;
      }
    }
  }
  _unlink(a, b);
  _unlink(b, a);
};

ParametricManager.prototype.coincident = function(objs) {
  if (objs.length == 0) return;
  this.linkObjects(objs);
  this.solve();
  this.viewer.refresh();
};

ParametricManager.prototype.getSolveData = function() {
  var sdata = [];
  for (var i = 0; i < this.subSystems.length; i++) {
    this.__getSolveData(this.subSystems[i].constraints, sdata);
  }
  return sdata;
};

ParametricManager.prototype.__getSolveData = function(constraints, out) {
  for (var i = 0; i < constraints.length; ++i) {
    var constraint = constraints[i];
    var data = constraint.getSolveData(this.constantResolver);
    for (var j = 0; j < data.length; ++j) {
      data[j].push(constraint.reducible !== undefined);
      out.push(data[j]);
    }
  }
  return out;
};

ParametricManager.prototype.solve = function() {
  var solver = this.prepare([]);
  solver.solve(false);
  solver.sync();
};

ParametricManager.prototype.solveWithLock = function(lock) {
  var solver = this.prepare(lock);
  solver.solve(false);
  solver.sync();
};

ParametricManager.prototype.prepare = function(locked, extraConstraints) {
  return this._prepare(locked, this.subSystems, extraConstraints);
};

ParametricManager.prototype._prepare = function(locked, subSystems, extraConstraints) {
  var solvers = [];
  for (var i = 0; i < subSystems.length; i++) {
    solvers.push(this.prepareForSubSystem(locked, subSystems[i].constraints, extraConstraints));
  }
  if (subSystems.length == 0 && locked.length != 0) {
    solvers.push(this.prepareForSubSystem(locked, [], extraConstraints));
  }
  return {
    solvers : solvers,

    solve : function(rough) {
      for (var i = 0; i < solvers.length; i++) {
        var alg =  i < subSystems.length ? subSystems[i].alg : 1;
        var res = solvers[i].solve(rough, alg);
        if (res.returnCode !== 1) {
          alg = alg == 1 ? 2 : 1;
          //if (solvers[i].solve(rough, alg).returnCode == 1) {
            //subSystems[i].alg = alg;
          //}
        }
      }
    },

    sync : function() {
      for (var i = 0; i < solvers.length; i++) {
        solvers[i].sync();
      }
    },

    updateLock : function(values) {
      for (var i = 0; i < solvers.length; i++) {
        solvers[i].updateLock(values);
      }
    }
  }
};

ParametricManager.isAux = function(obj) {
  while (!!obj) {
    if (!!obj.aux) {
      return true;
    }
    obj = obj.parent;
  }
  return false;
};

ParametricManager.fetchAuxParams = function(system, auxParams, auxDict) {
  for (var i = 0; i < system.length; ++i) {
    for (var p = 0; p < system[i][1].length; ++p) {
      var parameter = system[i][1][p];
      if (parameter.obj !== undefined) {
        if (ParametricManager.isAux(parameter.obj)) {
          if (auxDict[parameter.id] === undefined) {
            auxDict[parameter.id] = parameter;
            auxParams.push(parameter);
          }
        }
      }
    }
  }
};

ParametricManager.__toId = function(v) {
  return v.id;
};

ParametricManager.prototype.prepareForSubSystem = function(locked, subSystemConstraints, extraConstraints) {

  var pdict = {};
  var params;
  var _constrs = [];

  var equalsDict = {};
  var equalsIndex = [];
  var eqcElimination = {};

  var system = [];
  this.__getSolveData(subSystemConstraints, system);
  if (!!extraConstraints) this.__getSolveData(extraConstraints, system);

  var auxParams = [];
  var auxDict = {};

  ParametricManager.fetchAuxParams(system, auxParams, auxDict);

  var links = [];
  if (EQUALS_ELIMINATION_ENABLED) {
    var c, pi, paramToConstraints = {};

    function intersect(array1, array2) {
      if (!array1 || !array2) return false;
      return array1.filter(function(n) {
        return array2.indexOf(n) != -1
      }).length != 0;
    }

    for (i = 0; i < system.length; ++i) {
      c = system[i];
      if (c[3] !== true) {
        var sameParams = {};
        for (pi = 0; pi < c[1].length; pi++) {
          var param = c[1][pi];
          var paramConstrs = paramToConstraints[param.id];
          if (paramConstrs === undefined) {
            paramConstrs = [];
            paramToConstraints[param.id] = paramConstrs;
          }
          paramConstrs.push(i);
        }
      }
    }


    function Link(a, b, constr) {
      this.a = a;
      this.b = b;
      this.constr = constr;
      this.invalid = false;
      this.processed = false;
    }

    for (i = 0; i < system.length; ++i) {
      c = system[i];
      if (c[3] === true) { //Reduce flag
        var cp1 = c[1][0];
        var cp2 = c[1][1];
        links.push(new Link(cp1, cp2, i));
      }
    }
  }

  function shared(param1, param2) {
    if (param1 == param2) return false;
    var assoc0 = paramToConstraints[param1];
    var assoc1 = paramToConstraints[param2];
    return intersect(assoc0, assoc1);
  }

  var linkTuples = [];

  function mergeLinks(startIndex, into) {
    var linkI = links[startIndex];
    if (linkI.processed) return;
    linkI.processed = true;
    into.push(linkI);
    for (var j = startIndex + 1; j < links.length; j++) {
      var linkJ = links[j];
      if (linkI.a.id == linkJ.a.id || linkI.a.id == linkJ.b.id || linkI.b.id == linkJ.a.id || linkI.b.id == linkJ.b.id) {
        mergeLinks(j, into);
      }
    }
  }
  for (i = 0; i < links.length; i++) {
    if (links[i].processed) continue;
    var linkTuple = [];
    linkTuples.push(linkTuple);
    mergeLinks(i, linkTuple)
  }

  function resolveConflicts() {
    for (var i = 0; i < linkTuples.length; i++) {
      var tuple = linkTuples[i];

      for (var j = 0; j < tuple.length; j++) {
        var linkA = tuple[j];
        if (linkA.invalid) continue;
        if (shared(linkA.a.id, linkA.b.id)) {
          linkA.invalid = true;
          continue;
        }
        for (var k = j + 1; k < tuple.length; k++) {
          var linkB = tuple[k];
          if (shared(linkA.a.id, linkB.a.id) || shared(linkA.a.id, linkB.b.id) || shared(linkA.b.id, linkB.a.id) || shared(linkA.b.id, linkB.b.id)) {
            linkB.invalid = true;
          }
        }
      }
    }
  }
  resolveConflicts();

  function _merge(arr1, arr2) {
    for (var i = 0; i < arr2.length; ++i) {
      if (arr1.indexOf(arr2[i]) < 0) {
        arr1.push(arr2[i]);
      }
    }
  }

  function linksToTuples(linkTuples) {
    var tuples = [];
    for (var i = 0; i < linkTuples.length; i++) {
      var linkTuple = linkTuples[i];
      var tuple = [];
      tuples.push(tuple);
      for (var j = 0; j < linkTuple.length; j++) {
        var link = linkTuple[j];
        if (!link.invalid) {
          _merge(tuple, [link.a.id, link.b.id]);
          eqcElimination[link.constr] = true;
          equalsDict[link.a.id] = link.a;
          equalsDict[link.b.id] = link.b;
        }
      }
    }
    return tuples;
  }
  var tuples = linksToTuples(linkTuples);

  var readOnlyParams = auxParams.concat(locked);
  for (var i = 0; i < tuples.length; ++i) {
    var tuple = tuples[i];
    equalsIndex.push(tuple);
    for (var mi = 0; mi < readOnlyParams.length; ++mi) {
      var master = readOnlyParams[mi];
      var masterIdx = tuple.indexOf(master.id);
      if (masterIdx >= 0) {
        var tmp = tuple[0];
        tuple[0] = tuple[masterIdx];
        tuple[masterIdx] = tmp;
        break;
      }
    }
  }

  var equalsElimination = {};
  for (ei = 0; ei < equalsIndex.length; ++ei) {
    var master = equalsIndex[ei][0];
    for (i = 1; i < equalsIndex[ei].length; ++i) {
      equalsElimination[equalsIndex[ei][i]] = master;
    }
  }

  function getParam(p) {
    var master = equalsElimination[p.id];
    if (master !== undefined) {
      p = equalsDict[master];
    }
    var _p = pdict[p.id];
    if (_p === undefined) {
      if (p.__cachedParam__ === undefined) {
        _p = new parametric.Param(p.id, p.get());
        p.__cachedParam__ = _p;
      } else {
        _p = p.__cachedParam__;
        _p.reset(p.get());
      }

      _p._backingParam = p;
      pdict[p.id] = _p;
    }
    return _p;
  }

  var i;
  var p;
  var _p;
  var ei;

  var aux = [];
  for (i = 0; i < system.length; ++i) {


    var sdata = system[i];
    params = [];

    for (p = 0; p < sdata[1].length; ++p) {
      var param = sdata[1][p];
      _p = getParam(param);
      params.push(_p);

      (function () {
        if (auxDict[param.id] !== undefined) {
          aux.push(_p);
          return;
        }
        for (var i = 0; i < equalsIndex.length; ++i) {
          var eqParms = equalsIndex[i];
          if (eqParms.indexOf(param.id) != -1) {
            for (var j = 0; j < eqParms.length; j++) {
              var eqp = eqParms[j];
              if (auxDict[eqp] !== undefined) {
                aux.push(_p);
                return;
              }
            }
          }
        }
      })();
    }
    if (eqcElimination[i] === true) continue;

    var _constr = constraints.create(sdata[0], params, sdata[2]);
    _constrs.push(_constr);
  }

  var _locked = [];
  if (locked !== undefined) {
    for (p = 0; p < locked.length; ++p) {
      _locked[p] = getParam(locked[p]);
    }
  }

  var solver = parametric.prepare(_constrs, _locked, aux);
  function solve(rough, alg) {
    return solver.solveSystem(rough, alg);
  }
  var viewer = this.viewer;
  function sync() {
    for (p in pdict) {
      _p = pdict[p];
      if (!!_p._backingParam.__aux) continue;
      _p._backingParam.set(_p.get());
    }

    //Make sure all coincident constraints are equal
    for (ei = 0; ei < equalsIndex.length; ++ei) {
      var master = equalsDict[ equalsIndex[ei][0]];
      for (i = 1; i < equalsIndex[ei].length; ++i) {
        var slave = equalsDict[equalsIndex[ei][i]];
        slave.set(master.get());
      }
    }
  }

  solver.solve = solve;
  solver.sync = sync;
  return solver;
};

Constraints.ParentsCollector = function() {
  this.parents = [];
  var parents = this.parents;
  var index = {};
  function add(obj) {
    if (index[obj.id] === undefined) {
      index[obj.id] = obj;
      parents.push(obj);
    }
  }
  this.check = function(obj) {
    if (obj.parent !== null) {
      add(obj.parent);
    } else {
      add(obj);
    }
  };
};

Constraints.Factory = {};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.Coincident = function(a, b) {
  this.a = a;
  this.b = b;
  a.linked.push(b);
  b.linked.push(a);
};

Constraints.Coincident.prototype.NAME = 'coi';
Constraints.Coincident.prototype.UI_NAME = 'Coincident';
Constraints.Coincident.prototype.reducible = true;

Constraints.Coincident.prototype.getSolveData = function() {
  return [
    ['equal', [this.a._x, this.b._x], []],
    ['equal', [this.a._y, this.b._y], []]
  ];
};

Constraints.Coincident.prototype.serialize = function() {
  return [this.NAME, [this.a.id, this.b.id]];
};

Constraints.Factory[Constraints.Coincident.prototype.NAME] = function(refs, data) {
  return new Constraints.Coincident(refs(data[0]), refs(data[1]));
};

Constraints.Coincident.prototype.getObjects = function() {
  return [this.a, this.b];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.Lock = function(p, c) {
  this.p = p;
  this.c = c;
};

Constraints.Lock.prototype.NAME = 'lock';
Constraints.Lock.prototype.UI_NAME = 'Lock';

Constraints.Lock.prototype.getSolveData = function() {
  return [
    ['equalsTo', [this.p._x], [this.c.x]],
    ['equalsTo', [this.p._y], [this.c.y]]
  ];
};

Constraints.Lock.prototype.serialize = function() {
  return [this.NAME, [this.p.id, this.c]];
};

Constraints.Factory[Constraints.Lock.prototype.NAME] = function(refs, data) {
  return new Constraints.Lock(refs(data[0]), data[1]);
};


Constraints.Lock.prototype.getObjects = function() {
  return [this.p];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.Parallel = function(l1, l2) {
  this.l1 = l1;
  this.l2 = l2;
};

Constraints.Parallel.prototype.NAME = 'parallel';
Constraints.Parallel.prototype.UI_NAME = 'Parallel';

Constraints.Parallel.prototype.getSolveData = function() {
  var params = [];
  this.l1.collectParams(params);
  this.l2.collectParams(params);
  return [[this.NAME, params, []]];
};

Constraints.Parallel.prototype.serialize = function() {
  return [this.NAME, [this.l1.id, this.l2.id]];
};

Constraints.Factory[Constraints.Parallel.prototype.NAME] = function(refs, data) {
  return new Constraints.Parallel(refs(data[0]), refs(data[1]));
};

Constraints.Parallel.prototype.getObjects = function() {
  return [this.l1, this.l2];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.Perpendicular = function(l1, l2) {
  this.l1 = l1;
  this.l2 = l2;
};

Constraints.Perpendicular.prototype.NAME = 'perpendicular';
Constraints.Perpendicular.prototype.UI_NAME = 'Perpendicular';

Constraints.Perpendicular.prototype.getSolveData = function() {
  var params = [];
  this.l1.collectParams(params);
  this.l2.collectParams(params);
  return [[this.NAME, params, []]];
};

Constraints.Perpendicular.prototype.serialize = function() {
  return [this.NAME, [this.l1.id, this.l2.id]];
};

Constraints.Factory[Constraints.Perpendicular.prototype.NAME] = function(refs, data) {
  return new Constraints.Perpendicular(refs(data[0]), refs(data[1]));
};

Constraints.Perpendicular.prototype.getObjects = function() {
  return [this.l1, this.l2];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.P2LDistance = function(p, l, d) {
  this.p = p;
  this.l = l;
  this.d = d;
};

Constraints.P2LDistance.prototype.NAME = 'P2LDistance';
Constraints.P2LDistance.prototype.UI_NAME = 'Distance P & L';

Constraints.P2LDistance.prototype.getSolveData = function(resolver) {
  var params = [];
  this.p.collectParams(params);
  this.l.collectParams(params);
  return [[this.NAME, params, [resolver(this.d)]]];
};

Constraints.P2LDistance.prototype.serialize = function() {
  return [this.NAME, [this.p.id, this.l.id, this.d]];
};

Constraints.Factory[Constraints.P2LDistance.prototype.NAME] = function(refs, data) {
  return new Constraints.P2LDistance(refs(data[0]), refs(data[1]), data[2]);
};

Constraints.P2LDistance.prototype.getObjects = function() {
  return [this.p, this.l];
};

Constraints.P2LDistance.prototype.SettableFields = {'d' : "Enter the distance"};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.MinLength = function(a, b, min) {
  this.a = a;
  this.b = b;
  this.min = min;
};

Constraints.MinLength.prototype.aux = true;
Constraints.MinLength.prototype.NAME = 'MinLength';
Constraints.MinLength.prototype.UI_NAME = 'MinLength';

Constraints.MinLength.prototype.getSolveData = function() {
  var params = [];
  this.a.collectParams(params);
  this.b.collectParams(params);
  return [[this.NAME, params, [this.min]]];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.P2LDistanceV = function(p, l, d) {
  this.p = p;
  this.l = l;
  this.d = d;
};

Constraints.P2LDistanceV.prototype.aux = true;
Constraints.P2LDistanceV.prototype.NAME = 'P2LDistanceV';
Constraints.P2LDistanceV.prototype.UI_NAME = 'Distance P & L';

Constraints.P2LDistanceV.prototype.getSolveData = function() {
  var params = [];
  this.p.collectParams(params);
  this.l.collectParams(params);
  params.push(this.d);
  return [[this.NAME, params]];
};

// We don't serialize auxiliary constraints
//
//Constraints.P2LDistanceV.prototype.serialize = function() {
//  return [this.NAME, [this.p.id, this.l.id, this.d.id]];
//};
//
//Constraints.Factory[Constraints.P2LDistanceV.prototype.NAME] = function(refs, data) {
//  return new Constraints.P2LDistanceV(refs(data[0]), refs(data[1]), refs(data[2]));
//};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.P2PDistance = function(p1, p2, d) {
  this.p1 = p1;
  this.p2 = p2;
  this.d = d;
};

Constraints.P2PDistance.prototype.NAME = 'P2PDistance';
Constraints.P2PDistance.prototype.UI_NAME = 'Distance Points';

Constraints.P2PDistance.prototype.getSolveData = function(resolver) {
  var params = [];
  this.p1.collectParams(params);
  this.p2.collectParams(params);
  return [[this.NAME, params, [resolver(this.d)]]];
};

Constraints.P2PDistance.prototype.serialize = function() {
  return [this.NAME, [this.p1.id, this.p2.id, this.d]];
};

Constraints.Factory[Constraints.P2PDistance.prototype.NAME] = function(refs, data) {
  return new Constraints.P2PDistance(refs(data[0]), refs(data[1]), data[2]);
};

Constraints.P2PDistance.prototype.getObjects = function() {
  return [this.p1, this.p2];
};

Constraints.P2PDistance.prototype.SettableFields = {'d' : "Enter the distance"};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.P2PDistanceV = function(p1, p2, d) {
  this.p1 = p1;
  this.p2 = p2;
  this.d = d;
};

Constraints.P2PDistanceV.prototype.aux = true;
Constraints.P2PDistanceV.prototype.NAME = 'P2PDistanceV';
Constraints.P2PDistanceV.prototype.UI_NAME = 'Distance Points';

Constraints.P2PDistanceV.prototype.getSolveData = function() {
  var params = [];
  this.p1.collectParams(params);
  this.p2.collectParams(params);
  params.push(this.d);
  return [[this.NAME, params]];
};

// We don't serialize auxiliary constraints
//
//Constraints.P2PDistanceV.prototype.serialize = function() {
//  return [this.NAME, [this.p1.id, this.p2.id, this.d.id]];
//};
//
//Constraints.Factory[Constraints.P2PDistanceV.prototype.NAME] = function(refs, data) {
//  return new Constraints.P2PDistanceV(refs(data[0]), refs(data[1]), refs(data[2]));
//};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.Radius = function(arc, d) {
  this.arc = arc;
  this.d = d;
};

Constraints.Radius.prototype.NAME = 'Radius';
Constraints.Radius.prototype.UI_NAME = 'Radius Value';


Constraints.Radius.prototype.getSolveData = function(resolver) {
  return [['equalsTo', [this.arc.r], [resolver(this.d)]]];
};

Constraints.Radius.prototype.serialize = function() {
  return [this.NAME, [this.arc.id, this.d]];
};

Constraints.Factory[Constraints.Radius.prototype.NAME] = function(refs, data) {
  return new Constraints.Radius(refs(data[0]), data[1]);
};

Constraints.Radius.prototype.getObjects = function() {
  return [this.arc];
};

Constraints.Radius.prototype.SettableFields = {'d' : "Enter the radius value"};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.RR = function(arc1, arc2) {
  this.arc1 = arc1;
  this.arc2 = arc2;
};

Constraints.RR.prototype.NAME = 'RR';
Constraints.RR.prototype.UI_NAME = 'Radius Equality';
//Constraints.RR.prototype.reducible = true;


Constraints.RR.prototype.getSolveData = function() {
  return [['equal', [this.arc1.r, this.arc2.r], []]];
};

Constraints.RR.prototype.serialize = function() {
  return [this.NAME, [this.arc1.id, this.arc2.id]];
};

Constraints.Factory[Constraints.RR.prototype.NAME] = function(refs, data) {
  return new Constraints.RR(refs(data[0]), refs(data[1]));
};

Constraints.RR.prototype.getObjects = function() {
  return [this.arc1, this.arc2];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.LL = function(line1, line2) {
  this.line1 = line1;
  this.line2 = line2;
  this.length = new Ref(math.distanceAB(line1.a, line1.b));
};

Constraints.LL.prototype.NAME = 'LL';
Constraints.LL.prototype.UI_NAME = 'Lines Equality';

Constraints.LL.prototype.getSolveData = function() {
  var params1 = [];
  var params2 = [];
  this.line1.collectParams(params1);
  this.line2.collectParams(params2);
  params1.push(this.length);
  params2.push(this.length);
  return [
    ['P2PDistanceV', params1, []],
    ['P2PDistanceV', params2, []]
  ];
};

Constraints.LL.prototype.serialize = function() {
  return [this.NAME, [this.line1.id, this.line2.id]];
};

Constraints.Factory[Constraints.LL.prototype.NAME] = function(refs, data) {
  return new Constraints.LL(refs(data[0]), refs(data[1]));
};

Constraints.LL.prototype.getObjects = function() {
  return [this.line1, this.line2];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.Vertical = function(line) {
  this.line = line;
};

Constraints.Vertical.prototype.NAME = 'Vertical';
Constraints.Vertical.prototype.UI_NAME = 'Vertical';
//Constraints.Vertical.prototype.reducible = true;

Constraints.Vertical.prototype.getSolveData = function() {
  return [['equal', [this.line.a._x, this.line.b._x], []]];
};

Constraints.Vertical.prototype.serialize = function() {
  return [this.NAME, [this.line.id]];
};

Constraints.Factory[Constraints.Vertical.prototype.NAME] = function(refs, data) {
  return new Constraints.Vertical(refs(data[0]));
};

Constraints.Vertical.prototype.getObjects = function() {
  return [this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.Horizontal = function(line) {
  this.line = line;
};

Constraints.Horizontal.prototype.NAME = 'Horizontal';
Constraints.Horizontal.prototype.UI_NAME = 'Horizontal';
//Constraints.Horizontal.prototype.reducible = true;

Constraints.Horizontal.prototype.getSolveData = function() {
  return [['equal', [this.line.a._y, this.line.b._y], []]];
};

Constraints.Horizontal.prototype.serialize = function() {
  return [this.NAME, [this.line.id]];
};

Constraints.Factory[Constraints.Horizontal.prototype.NAME] = function(refs, data) {
  return new Constraints.Horizontal(refs(data[0]));
};

Constraints.Horizontal.prototype.getObjects = function() {
  return [this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.Tangent = function(arc, line) {
  this.arc = arc;
  this.line = line;
};

Constraints.Tangent.prototype.NAME = 'Tangent';
Constraints.Tangent.prototype.UI_NAME = 'Tangent';

Constraints.Tangent.prototype.getSolveData = function() {
  var params = [];
  this.arc.c.collectParams(params);
  this.line.collectParams(params);
  params.push(this.arc.r);
  return [['P2LDistanceV', params, []]];
};

Constraints.Tangent.prototype.serialize = function() {
  return [this.NAME, [this.arc.id, this.line.id]];
};

Constraints.Factory[Constraints.Tangent.prototype.NAME] = function(refs, data) {
  return new Constraints.Tangent(refs(data[0]), refs(data[1]));
};

Constraints.Tangent.prototype.getObjects = function() {
  return [this.arc, this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.PointOnLine = function(point, line) {
  this.point = point;
  this.line = line;
};

Constraints.PointOnLine.prototype.NAME = 'PointOnLine';
Constraints.PointOnLine.prototype.UI_NAME = 'Point On Line';

Constraints.PointOnLine.prototype.getSolveData = function() {
  var params = [];
  this.point.collectParams(params);
  this.line.collectParams(params);
  return [['P2LDistance', params, [0]]];
};

Constraints.PointOnLine.prototype.serialize = function() {
  return [this.NAME, [this.point.id, this.line.id]];
};

Constraints.Factory[Constraints.PointOnLine.prototype.NAME] = function(refs, data) {
  return new Constraints.PointOnLine(refs(data[0]), refs(data[1]));
};

Constraints.PointOnLine.prototype.getObjects = function() {
  return [this.point, this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.PointOnArc = function(point, arc) {
  this.point = point;
  this.arc = arc;
};

Constraints.PointOnArc.prototype.NAME = 'PointOnArc';
Constraints.PointOnArc.prototype.UI_NAME = 'Point On Arc';

Constraints.PointOnArc.prototype.getSolveData = function() {
  var params = [];
  this.point.collectParams(params);
  this.arc.c.collectParams(params);
  params.push(this.arc.r);
  return [['P2PDistanceV', params, []]];
};

Constraints.PointOnArc.prototype.serialize = function() {
  return [this.NAME, [this.point.id, this.arc.id]];
};

Constraints.Factory[Constraints.PointOnArc.prototype.NAME] = function(refs, data) {
  return new Constraints.PointOnArc(refs(data[0]), refs(data[1]));
};

Constraints.PointOnArc.prototype.getObjects = function() {
  return [this.point, this.arc];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.PointInMiddle = function(point, line) {
  this.point = point;
  this.line = line;
  this.length = new Ref(math.distanceAB(line.a, line.b) / 2);
};

Constraints.PointInMiddle.prototype.NAME = 'PointInMiddle';
Constraints.PointInMiddle.prototype.UI_NAME = 'Point In the Middle';

Constraints.PointInMiddle.prototype.getSolveData = function() {
  var params1 = [];
  var params2 = [];

  this.line.a.collectParams(params1);
  this.point.collectParams(params1);
  params1.push(this.length);

  this.line.b.collectParams(params2);
  this.point.collectParams(params2);
  params2.push(this.length);

  return [
    ['P2PDistanceV', params1, []],
    ['P2PDistanceV', params2, []]
  ];
};

Constraints.PointInMiddle.prototype.serialize = function() {
  return [this.NAME, [this.point.id, this.line.id]];
};

Constraints.Factory[Constraints.PointInMiddle.prototype.NAME] = function(refs, data) {
  return new Constraints.PointInMiddle(refs(data[0]), refs(data[1]));
};

Constraints.PointInMiddle.prototype.getObjects = function() {
  return [this.point, this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.Symmetry = function(point, line) {
  this.point = point;
  this.line = line;
  this.length = new Ref(math.distanceAB(line.a, line.b) / 2);
};

Constraints.Symmetry.prototype.NAME = 'Symmetry';
Constraints.Symmetry.prototype.UI_NAME = 'Symmetry';

Constraints.Symmetry.prototype.getSolveData = function(resolver) {
  var pointInMiddleData = Constraints.PointInMiddle.prototype.getSolveData.call(this, [resolver]);
  var pointOnLineData = Constraints.PointOnLine.prototype.getSolveData.call(this, [resolver]);
  return pointInMiddleData.concat(pointOnLineData);
};

Constraints.Symmetry.prototype.serialize = function() {
  return [this.NAME, [this.point.id, this.line.id]];
};

Constraints.Factory[Constraints.Symmetry.prototype.NAME] = function(refs, data) {
  return new Constraints.Symmetry(refs(data[0]), refs(data[1]));
};

Constraints.Symmetry.prototype.getObjects = function() {
  return [this.point, this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.Angle = function(p1, p2, p3, p4, angle) {
  this.p1 = p1;
  this.p2 = p2;
  this.p3 = p3;
  this.p4 = p4;
  this._angle = new Ref(0);
  this.angle = angle;
};

Constraints.Angle.prototype.NAME = 'Angle';
Constraints.Angle.prototype.UI_NAME = 'Lines Angle';

Constraints.Angle.prototype.getSolveData = function(resolver) {
  this._angle.set(resolver(this.angle) / 180 * Math.PI);
  var params = [];
  this.p1.collectParams(params);
  this.p2.collectParams(params);
  this.p3.collectParams(params);
  this.p4.collectParams(params);
  params.push(this._angle);
  return [['angleConst', params, []]];
};

Constraints.Angle.prototype.serialize = function() {
  return [this.NAME, [this.p1.id, this.p2.id, this.p3.id, this.p4.id, this.angle]];
};

Constraints.Factory[Constraints.Angle.prototype.NAME] = function(refs, data) {
  return new Constraints.Angle( refs(data[0]), refs(data[1]), refs(data[2]), refs(data[3]), data[4] );
};

Constraints.Angle.prototype.getObjects = function() {
  var collector = new Constraints.ParentsCollector();
  collector.check(this.p1);
  collector.check(this.p2);
  collector.check(this.p3);
  collector.check(this.p4);
  return collector.parents;
};

Constraints.Angle.prototype.SettableFields = {'angle' : "Enter the angle value"};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
Constraints.LockConvex = function(c, a, t) {
  this.c = c;
  this.a = a;
  this.t = t;
};

Constraints.LockConvex.prototype.NAME = 'LockConvex';
Constraints.LockConvex.prototype.UI_NAME = 'Lock Convexity';

Constraints.LockConvex.prototype.getSolveData = function() {
  var params = [];
  this.c.collectParams(params);
  this.a.collectParams(params);
  this.t.collectParams(params);
  return [['LockConvex', params, []]];
};

Constraints.LockConvex.prototype.serialize = function() {
  return [this.NAME, [this.c.id, this.a.id, this.t.id]];
};

Constraints.Factory[Constraints.LockConvex.prototype.NAME] = function(refs, data) {
  return new Constraints.LockConvex(refs(data[0]), refs(data[1]), refs(data[2]));
};

Constraints.LockConvex.prototype.getObjects = function() {
  var collector = new Constraints.ParentsCollector();
  collector.check(this.c);
  collector.check(this.a);
  collector.check(this.t);
  return collector.parents;
};

// FIXME: From fetchers.js



ParametricManager.prototype._fetchTwoPoints = function(objs) {
  var points = [];
  for (var i = 0; i < objs.length; ++i) {
    if (objs[i]._class == 'TCAD.TWO.EndPoint') {
      points.push(objs[i]);
    } else if (objs[i]._class == 'TCAD.TWO.Segment') {
      points.push(objs[i].a);
      points.push(objs[i].b);
    }
  }
  if (points.length < 2) {
    throw "Illegal Argument. Constraint requires 2 points or 1 line."
  }
  return points;
};

ParametricManager.prototype._fetchPoints = function(objs) {
  var points = [];
  for (var i = 0; i < objs.length; ++i) {
    objs[i].accept(function(o) {
      if (o._class === 'TCAD.TWO.EndPoint')  {
        points.push(o);
      }
      return true;
    });
  }
  if (points.length == 0) {
    throw "Illegal Argument. Constraint requires at least 1 point/line/arc/circle."
  }
  return points;
};

ParametricManager.prototype._fetchArkCirc = function(objs, min) {
  var arcs = [];
  for (var i = 0; i < objs.length; ++i) {
    if (objs[i]._class === 'TCAD.TWO.Arc' || objs[i]._class === 'TCAD.TWO.Circle') {
      arcs.push(objs[i]);
    }
  }
  if (arcs.length < min) {
    throw "Illegal Argument. Constraint requires at least " + min + " arcs/circles."
  }
  return arcs;
};

ParametricManager.prototype._fetch = function(objs, types, min) {
  var result = [];
  for (var i = 0; i < objs.length; ++i) {
    if (types.indexOf(objs[i]._class)  > -1 ) {
      result.push(objs[i]);
    }
  }
  if (result.length < min) {
    throw "Illegal Argument. Constraint requires at least " + min + " of " + types;
  }
  return result;
};

ParametricManager.prototype._fetchPointAndLine = function(objs) {

  var point = null;
  var line = null;

  for (var i = 0; i < objs.length; ++i) {
    if (objs[i]._class == 'TCAD.TWO.EndPoint') {
      point = objs[i];
    } else if (objs[i]._class == 'TCAD.TWO.Segment') {
      line = objs[i];
    }
  }
  if (point == null || line == null) {
    throw "Illegal Argument. Constraint requires point and line."
  }

  return [point, line];
};

ParametricManager.prototype._fetchLine = function(objs) {
  for (var i = 0; i < objs.length; ++i) {
    if (objs[i]._class == 'TCAD.TWO.Segment') {
      return objs[i];
    }
  }
  throw "Illegal Argument. Constraint requires a line."
};

ParametricManager.prototype._fetchArcCircAndLine = function(objs) {

  var arc = null;
  var line = null;

  for (var i = 0; i < objs.length; ++i) {
    if (objs[i]._class === 'TCAD.TWO.Arc' || objs[i]._class === 'TCAD.TWO.Circle') {
      arc = objs[i];
    } else if (objs[i]._class == 'TCAD.TWO.Segment') {
      line = objs[i];
    }
  }
  if (arc == null || line == null) {
    throw "Illegal Argument. Constraint requires arc and line."
  }

  return [arc, line];
};

ParametricManager.prototype._fetchTwoLines = function(objs) {
  var lines = [];
  for (var i = 0; i < objs.length; ++i) {
    if (objs[i]._class == 'TCAD.TWO.Segment') {
      lines.push(objs[i]);
    }
  }
  if (lines.length < 2) {
    throw "Illegal Argument. Constraint requires 2 lines."
  }
  return lines;
};


// ------------------------------------------------------------------------------------------------------------------ //

export { Constraints, ParametricManager, utils };
