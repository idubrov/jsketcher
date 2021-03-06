TCAD.TWO.Constraints = {};
TCAD.EQUALS_ELIMINATION_ENABLED = true;

/** @constructor */
TCAD.TWO.SubSystem = function() {
  this.alg = 1;
  this.error = 0;
  this.reduce = false;
  this.constraints = [];
};

/** @constructor */
TCAD.TWO.ParametricManager = function(viewer) {
  this.viewer = viewer;
  this.subSystems = [];
  this.listeners = [];
  this.constantTable = {};
  
  this.viewer.params.define("constantDefinition", null);
  this.viewer.params.subscribe("constantDefinition", "parametricManager", this.rebuildConstantTable, this)();
  this.constantResolver = this.createConstantResolver();
};

TCAD.TWO.ParametricManager.prototype.createConstantResolver = function() {
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

TCAD.TWO.ParametricManager.prototype.notify = function(event) {
  for (var i = 0; i < this.listeners.length; ++i) {
    var l = this.listeners[i];
    l(event);
  }
};

TCAD.TWO.ParametricManager.prototype.rebuildConstantTable = function(constantDefinition) {
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


TCAD.TWO.ParametricManager.prototype.findComponents = function(constr) {
  if (this.subSystems.length === 0) {
    this.subSystems.push(new TCAD.TWO.SubSystem());
  }
  return [0];  
};

TCAD.TWO.ParametricManager.prototype.tune = function(subSystem) {
  
};

TCAD.TWO.ParametricManager.prototype._add = function(constr) {
  var subSystemIds = this.findComponents(constr);
  var subSystem; 
  switch (subSystemIds.length) {
    case 0:
      subSystem = new TCAD.TWO.SubSystem();
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

TCAD.TWO.ParametricManager.prototype.checkRedundancy = function (subSystem, constr) {
  var solver = this.prepareForSubSystem([], subSystem.constraints);
  if (solver.diagnose().conflict) {
    alert("Most likely this "+constr.NAME + " constraint is CONFLICTING!")
  }
};

TCAD.TWO.ParametricManager.prototype.refresh = function() {
  this.solve();
  this.notify();
  this.viewer.refresh();
};

TCAD.TWO.ParametricManager.prototype.add = function(constr) {
  this.viewer.historyManager.checkpoint();
  var subSystem = this._add(constr);
  this.checkRedundancy(subSystem, constr);
  this.refresh();
};

TCAD.TWO.ParametricManager.prototype.addAll = function(constrs) {
  for (var i = 0; i < constrs.length; i++) {
    var subSystem = this._add(constrs[i]);
    this.checkRedundancy(subSystem, constrs[i]);
  }
  this.refresh();
};

TCAD.TWO.ParametricManager.prototype.remove = function(constr) {
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

TCAD.TWO.ParametricManager.prototype.removeConstraintsByObj = function(obj) {
  var ownedParams = [];
  obj.collectParams(ownedParams);
  this.removeConstraintsByParams(ownedParams);
};

TCAD.TWO.ParametricManager.prototype.removeConstraintsByParams = function(ownedParams) {
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

TCAD.TWO.ParametricManager.prototype.lock = function(objs) {
  var p = this._fetchPoints(objs);
  for (var i = 0; i < p.length; ++i) {
    this._add(new TCAD.TWO.Constraints.Lock(p[i], { x : p[i].x, y : p[i].y} ));
  }
  this.refresh();
};

TCAD.TWO.ParametricManager.prototype.vertical = function(objs) {
  this.add(new TCAD.TWO.Constraints.Vertical(this._fetchLine(objs)));
};

TCAD.TWO.ParametricManager.prototype.horizontal = function(objs) {
  this.add(new TCAD.TWO.Constraints.Horizontal(this._fetchLine(objs)));
};

TCAD.TWO.ParametricManager.prototype.parallel = function(objs) {
  var lines = this._fetchTwoLines(objs);
  this.add(new TCAD.TWO.Constraints.Parallel(lines[0], lines[1]));
};

TCAD.TWO.ParametricManager.prototype.perpendicular = function(objs) {
  var lines = this._fetchTwoLines(objs);
  this.add(new TCAD.TWO.Constraints.Perpendicular(lines[0], lines[1]));
};

TCAD.TWO.ParametricManager.prototype.lockConvex = function(objs, warnCallback) {
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
  
  this.add(new TCAD.TWO.Constraints.LockConvex(c, a, t));
};

TCAD.TWO.ParametricManager.prototype.tangent = function(objs) {
  var al = this._fetchArcCircAndLine(objs);
  var arc  = al[0];
  var line  = al[1];
  this.add(new TCAD.TWO.Constraints.Tangent( arc, line));
};

TCAD.TWO.ParametricManager.prototype.rr = function(arcs) {
  var prev = arcs[0];
  for (var i = 1; i < arcs.length; ++i) {
    this._add(new TCAD.TWO.Constraints.RR(prev, arcs[i]));
    prev = arcs[i];
  }
  this.refresh();
};

TCAD.TWO.ParametricManager.prototype.ll = function(lines) {
  var prev = lines[0];
  for (var i = 1; i < lines.length; ++i) {
    this._add(new TCAD.TWO.Constraints.LL(prev, lines[i]));
    prev = lines[i];
  }
  this.refresh();

};

TCAD.TWO.ParametricManager.prototype.entityEquality = function(objs) {
  var arcs = this._fetch(objs, ['TCAD.TWO.Arc', 'TCAD.TWO.Circle'], 0);
  var lines = this._fetch(objs, ['TCAD.TWO.Segment'], 0);
  if (arcs.length > 0) this.rr(arcs);
  if (lines.length > 0) this.ll(lines);
};

TCAD.TWO.ParametricManager.prototype.p2lDistance = function(objs, promptCallback) {
  var pl = this._fetchPointAndLine(objs);

  var target = pl[0];
  var segment = pl[1];
  
  var ex = new TCAD.Vector(-(segment.b.y - segment.a.y), segment.b.x - segment.a.x).normalize();
  var distance = Math.abs(ex.dot(new TCAD.Vector(segment.a.x - target.x, segment.a.y - target.y)));

  var promptDistance = TCAD.TWO.utils.askNumber(TCAD.TWO.Constraints.P2LDistance.prototype.SettableFields.d, distance.toFixed(2), promptCallback, this.constantResolver);

  if (promptDistance != null) {
    this.add(new TCAD.TWO.Constraints.P2LDistance(target, segment, promptDistance));
  }
};

TCAD.TWO.ParametricManager.prototype.pointInMiddle = function(objs) {
  var pl = this._fetchPointAndLine(objs);
  this.add(new TCAD.TWO.Constraints.PointInMiddle(pl[0], pl[1]));
};

TCAD.TWO.ParametricManager.prototype.symmetry = function(objs) {
  var pl = this._fetchPointAndLine(objs);
  this.add(new TCAD.TWO.Constraints.Symmetry(pl[0], pl[1]));
};

TCAD.TWO.ParametricManager.prototype.pointOnArc = function(objs) {
  var points = this._fetch(objs, ['TCAD.TWO.EndPoint'], 1);
  var arcs = this._fetch(objs, ['TCAD.TWO.Arc', 'TCAD.TWO.Circle'], 1);
  this.add(new TCAD.TWO.Constraints.PointOnArc(points[0], arcs[0]));
};

TCAD.TWO.ParametricManager.prototype.pointOnLine = function(objs) {
  var pl = this._fetchPointAndLine(objs);
  var target = pl[0];
  var segment = pl[1];
  this.add(new TCAD.TWO.Constraints.PointOnLine(target, segment));
};

TCAD.TWO.ParametricManager.prototype.llAngle = function(objs, promptCallback) {
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
  angle = TCAD.TWO.utils.askNumber(TCAD.TWO.Constraints.Angle.prototype.SettableFields.angle, angle.toFixed(2), promptCallback, this.constantResolver);
  if (angle === null) return;
  this.add(new TCAD.TWO.Constraints.Angle(points[0], points[1], points[2], points[3], angle));
};

TCAD.TWO.utils.constRef = function(value) {
  return function() {
    return value;    
  };
};

TCAD.TWO.ParametricManager.prototype.p2pDistance = function(objs, promptCallback) {
  var p = this._fetchTwoPoints(objs);
  var distance = new TCAD.Vector(p[1].x - p[0].x, p[1].y - p[0].y).length();
  var promptDistance = TCAD.TWO.utils.askNumber(TCAD.TWO.Constraints.P2PDistance.prototype.SettableFields.d, distance.toFixed(2), promptCallback, this.constantResolver);

  if (promptDistance != null) {
    this.add(new TCAD.TWO.Constraints.P2PDistance(p[0], p[1], promptDistance));
  }
};

TCAD.TWO.utils.askNumber = function(promptText, initValue, promptCallback, resolver) {
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

TCAD.TWO.ParametricManager.prototype.radius = function(objs, promptCallback) {
  var arcs = this._fetchArkCirc(objs, 1);
  var radius = arcs[0].r.get();
  var promptDistance = TCAD.TWO.utils.askNumber(TCAD.TWO.Constraints.Radius.prototype.SettableFields.d, radius.toFixed(2), promptCallback, this.constantResolver);
  if (promptDistance != null) {
    for (var i = 0; i < arcs.length; ++i) {
      this._add(new TCAD.TWO.Constraints.Radius(arcs[i], promptDistance));
    }
    this.refresh();
  }
};

TCAD.TWO.ParametricManager.prototype.linkObjects = function(objs) {
  var i;
  var masterIdx = -1;
  for (i = 0; i < objs.length; ++i) {
    if (TCAD.TWO.ParametricManager.isAux(objs[i])) {
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
    var c = new TCAD.TWO.Constraints.Coincident(objs[i], objs[masterIdx]);
    this._add(c);
  }
  this.notify();
};

TCAD.TWO.ParametricManager.prototype.unlinkObjects = function(a, b) {
  
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

TCAD.TWO.ParametricManager.prototype.coincident = function(objs) {
  if (objs.length == 0) return;
  this.linkObjects(objs);
  this.solve();
  this.viewer.refresh();
};

TCAD.TWO.ParametricManager.prototype.getSolveData = function() {
  var sdata = []; 
  for (var i = 0; i < this.subSystems.length; i++) {
    this.__getSolveData(this.subSystems[i].constraints, sdata);
  }
  return sdata;
};

TCAD.TWO.ParametricManager.prototype.__getSolveData = function(constraints, out) {
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

TCAD.TWO.ParametricManager.prototype.solve = function() {
  var solver = this.prepare([]);
  solver.solve(false);
  solver.sync();
};

TCAD.TWO.ParametricManager.prototype.solveWithLock = function(lock) {
  var solver = this.prepare(lock);
  solver.solve(false);
  solver.sync();
};

TCAD.TWO.ParametricManager.prototype.prepare = function(locked, extraConstraints) {
  return this._prepare(locked, this.subSystems, extraConstraints);
};

TCAD.TWO.ParametricManager.prototype._prepare = function(locked, subSystems, extraConstraints) {
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

TCAD.TWO.ParametricManager.isAux = function(obj) {
  while (!!obj) {
    if (!!obj.aux) {
      return true;
    }
    obj = obj.parent;
  }
  return false;
};

TCAD.TWO.ParametricManager.fetchAuxParams = function(system, auxParams, auxDict) {
  for (var i = 0; i < system.length; ++i) {
    for (var p = 0; p < system[i][1].length; ++p) {
      var parameter = system[i][1][p];
      if (parameter.obj !== undefined) {
        if (TCAD.TWO.ParametricManager.isAux(parameter.obj)) {
          if (auxDict[parameter.id] === undefined) {
            auxDict[parameter.id] = parameter;
            auxParams.push(parameter);
          }
        }
      }
    }
  }
};

TCAD.TWO.ParametricManager.__toId = function(v) {
  return v.id;
};

TCAD.TWO.ParametricManager.prototype.prepareForSubSystem = function(locked, subSystemConstraints, extraConstraints) {

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

  TCAD.TWO.ParametricManager.fetchAuxParams(system, auxParams, auxDict);

  var links = [];
  if (TCAD.EQUALS_ELIMINATION_ENABLED) {
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
        _p = new TCAD.parametric.Param(p.id, p.get());
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

    var _constr = TCAD.constraints.create(sdata[0], params, sdata[2]);
    _constrs.push(_constr);
  }

  var _locked = [];
  if (locked !== undefined) {
    for (p = 0; p < locked.length; ++p) {
      _locked[p] = getParam(locked[p]);
    }
  }
  
  var solver = TCAD.parametric.prepare(_constrs, _locked, aux);
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

TCAD.TWO.Constraints.ParentsCollector = function() {
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

TCAD.TWO.Constraints.Factory = {};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.Coincident = function(a, b) {
  this.a = a;
  this.b = b;
  a.linked.push(b);
  b.linked.push(a);
};

TCAD.TWO.Constraints.Coincident.prototype.NAME = 'coi';
TCAD.TWO.Constraints.Coincident.prototype.UI_NAME = 'Coincident';
TCAD.TWO.Constraints.Coincident.prototype.reducible = true;

TCAD.TWO.Constraints.Coincident.prototype.getSolveData = function() {
  return [
    ['equal', [this.a._x, this.b._x], []],
    ['equal', [this.a._y, this.b._y], []]
  ];
};

TCAD.TWO.Constraints.Coincident.prototype.serialize = function() {
  return [this.NAME, [this.a.id, this.b.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.Coincident.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.Coincident(refs(data[0]), refs(data[1]));  
};

TCAD.TWO.Constraints.Coincident.prototype.getObjects = function() {
  return [this.a, this.b];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.Lock = function(p, c) {
  this.p = p;
  this.c = c;
};

TCAD.TWO.Constraints.Lock.prototype.NAME = 'lock';
TCAD.TWO.Constraints.Lock.prototype.UI_NAME = 'Lock';

TCAD.TWO.Constraints.Lock.prototype.getSolveData = function() {
  return [
    ['equalsTo', [this.p._x], [this.c.x]],
    ['equalsTo', [this.p._y], [this.c.y]]
  ];
};

TCAD.TWO.Constraints.Lock.prototype.serialize = function() {
  return [this.NAME, [this.p.id, this.c]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.Lock.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.Lock(refs(data[0]), data[1]);
};


TCAD.TWO.Constraints.Lock.prototype.getObjects = function() {
  return [this.p];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.Parallel = function(l1, l2) {
  this.l1 = l1;
  this.l2 = l2;
};

TCAD.TWO.Constraints.Parallel.prototype.NAME = 'parallel';
TCAD.TWO.Constraints.Parallel.prototype.UI_NAME = 'Parallel';

TCAD.TWO.Constraints.Parallel.prototype.getSolveData = function() {
  var params = [];
  this.l1.collectParams(params);
  this.l2.collectParams(params);
  return [[this.NAME, params, []]];
};

TCAD.TWO.Constraints.Parallel.prototype.serialize = function() {
  return [this.NAME, [this.l1.id, this.l2.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.Parallel.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.Parallel(refs(data[0]), refs(data[1]));
};

TCAD.TWO.Constraints.Parallel.prototype.getObjects = function() {
  return [this.l1, this.l2];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.Perpendicular = function(l1, l2) {
  this.l1 = l1;
  this.l2 = l2;
};

TCAD.TWO.Constraints.Perpendicular.prototype.NAME = 'perpendicular';
TCAD.TWO.Constraints.Perpendicular.prototype.UI_NAME = 'Perpendicular';

TCAD.TWO.Constraints.Perpendicular.prototype.getSolveData = function() {
  var params = [];
  this.l1.collectParams(params);
  this.l2.collectParams(params);
  return [[this.NAME, params, []]];
};

TCAD.TWO.Constraints.Perpendicular.prototype.serialize = function() {
  return [this.NAME, [this.l1.id, this.l2.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.Perpendicular.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.Perpendicular(refs(data[0]), refs(data[1]));
};

TCAD.TWO.Constraints.Perpendicular.prototype.getObjects = function() {
  return [this.l1, this.l2];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.P2LDistance = function(p, l, d) {
  this.p = p;
  this.l = l;
  this.d = d;
};

TCAD.TWO.Constraints.P2LDistance.prototype.NAME = 'P2LDistance';
TCAD.TWO.Constraints.P2LDistance.prototype.UI_NAME = 'Distance P & L';

TCAD.TWO.Constraints.P2LDistance.prototype.getSolveData = function(resolver) {
  var params = [];
  this.p.collectParams(params);
  this.l.collectParams(params);
  return [[this.NAME, params, [resolver(this.d)]]];
};

TCAD.TWO.Constraints.P2LDistance.prototype.serialize = function() {
  return [this.NAME, [this.p.id, this.l.id, this.d]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.P2LDistance.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.P2LDistance(refs(data[0]), refs(data[1]), data[2]);
};

TCAD.TWO.Constraints.P2LDistance.prototype.getObjects = function() {
  return [this.p, this.l];
};

TCAD.TWO.Constraints.P2LDistance.prototype.SettableFields = {'d' : "Enter the distance"};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.MinLength = function(a, b, min) {
  this.a = a;
  this.b = b;
  this.min = min;
};

TCAD.TWO.Constraints.MinLength.prototype.aux = true;
TCAD.TWO.Constraints.MinLength.prototype.NAME = 'MinLength';
TCAD.TWO.Constraints.MinLength.prototype.UI_NAME = 'MinLength';

TCAD.TWO.Constraints.MinLength.prototype.getSolveData = function() {
  var params = [];
  this.a.collectParams(params);
  this.b.collectParams(params);
  return [[this.NAME, params, [this.min]]];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.P2LDistanceV = function(p, l, d) {
  this.p = p;
  this.l = l;
  this.d = d;
};

TCAD.TWO.Constraints.P2LDistanceV.prototype.aux = true;
TCAD.TWO.Constraints.P2LDistanceV.prototype.NAME = 'P2LDistanceV';
TCAD.TWO.Constraints.P2LDistanceV.prototype.UI_NAME = 'Distance P & L';

TCAD.TWO.Constraints.P2LDistanceV.prototype.getSolveData = function() {
  var params = [];
  this.p.collectParams(params);
  this.l.collectParams(params);
  params.push(this.d);
  return [[this.NAME, params]];
};

// We don't serialize auxiliary constraints
//
//TCAD.TWO.Constraints.P2LDistanceV.prototype.serialize = function() {
//  return [this.NAME, [this.p.id, this.l.id, this.d.id]];
//};
//
//TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.P2LDistanceV.prototype.NAME] = function(refs, data) {
//  return new TCAD.TWO.Constraints.P2LDistanceV(refs(data[0]), refs(data[1]), refs(data[2]));
//};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.P2PDistance = function(p1, p2, d) {
  this.p1 = p1;
  this.p2 = p2;
  this.d = d;
};

TCAD.TWO.Constraints.P2PDistance.prototype.NAME = 'P2PDistance';
TCAD.TWO.Constraints.P2PDistance.prototype.UI_NAME = 'Distance Points';

TCAD.TWO.Constraints.P2PDistance.prototype.getSolveData = function(resolver) {
  var params = [];
  this.p1.collectParams(params);
  this.p2.collectParams(params);
  return [[this.NAME, params, [resolver(this.d)]]];
};

TCAD.TWO.Constraints.P2PDistance.prototype.serialize = function() {
  return [this.NAME, [this.p1.id, this.p2.id, this.d]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.P2PDistance.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.P2PDistance(refs(data[0]), refs(data[1]), data[2]);
};

TCAD.TWO.Constraints.P2PDistance.prototype.getObjects = function() {
  return [this.p1, this.p2];
};

TCAD.TWO.Constraints.P2PDistance.prototype.SettableFields = {'d' : "Enter the distance"};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.P2PDistanceV = function(p1, p2, d) {
  this.p1 = p1;
  this.p2 = p2;
  this.d = d;
};

TCAD.TWO.Constraints.P2PDistanceV.prototype.aux = true;
TCAD.TWO.Constraints.P2PDistanceV.prototype.NAME = 'P2PDistanceV';
TCAD.TWO.Constraints.P2PDistanceV.prototype.UI_NAME = 'Distance Points';

TCAD.TWO.Constraints.P2PDistanceV.prototype.getSolveData = function() {
  var params = [];
  this.p1.collectParams(params);
  this.p2.collectParams(params);
  params.push(this.d);
  return [[this.NAME, params]];
};

// We don't serialize auxiliary constraints
//
//TCAD.TWO.Constraints.P2PDistanceV.prototype.serialize = function() {
//  return [this.NAME, [this.p1.id, this.p2.id, this.d.id]];
//};
//
//TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.P2PDistanceV.prototype.NAME] = function(refs, data) {
//  return new TCAD.TWO.Constraints.P2PDistanceV(refs(data[0]), refs(data[1]), refs(data[2]));
//};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.Radius = function(arc, d) {
  this.arc = arc;
  this.d = d;
};

TCAD.TWO.Constraints.Radius.prototype.NAME = 'Radius';
TCAD.TWO.Constraints.Radius.prototype.UI_NAME = 'Radius Value';


TCAD.TWO.Constraints.Radius.prototype.getSolveData = function(resolver) {
  return [['equalsTo', [this.arc.r], [resolver(this.d)]]];
};

TCAD.TWO.Constraints.Radius.prototype.serialize = function() {
  return [this.NAME, [this.arc.id, this.d]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.Radius.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.Radius(refs(data[0]), data[1]);
};

TCAD.TWO.Constraints.Radius.prototype.getObjects = function() {
  return [this.arc];
};

TCAD.TWO.Constraints.Radius.prototype.SettableFields = {'d' : "Enter the radius value"};

// ------------------------------------------------------------------------------------------------------------------ // 

/** @constructor */
TCAD.TWO.Constraints.RR = function(arc1, arc2) {
  this.arc1 = arc1;
  this.arc2 = arc2;
};

TCAD.TWO.Constraints.RR.prototype.NAME = 'RR';
TCAD.TWO.Constraints.RR.prototype.UI_NAME = 'Radius Equality';
//TCAD.TWO.Constraints.RR.prototype.reducible = true;


TCAD.TWO.Constraints.RR.prototype.getSolveData = function() {
  return [['equal', [this.arc1.r, this.arc2.r], []]];
};

TCAD.TWO.Constraints.RR.prototype.serialize = function() {
  return [this.NAME, [this.arc1.id, this.arc2.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.RR.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.RR(refs(data[0]), refs(data[1]));
};

TCAD.TWO.Constraints.RR.prototype.getObjects = function() {
  return [this.arc1, this.arc2];
};

// ------------------------------------------------------------------------------------------------------------------ // 

/** @constructor */
TCAD.TWO.Constraints.LL = function(line1, line2) {
  this.line1 = line1;
  this.line2 = line2;
  this.length = new TCAD.TWO.Ref(TCAD.math.distanceAB(line1.a, line1.b));
};

TCAD.TWO.Constraints.LL.prototype.NAME = 'LL';
TCAD.TWO.Constraints.LL.prototype.UI_NAME = 'Lines Equality';

TCAD.TWO.Constraints.LL.prototype.getSolveData = function() {
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

TCAD.TWO.Constraints.LL.prototype.serialize = function() {
  return [this.NAME, [this.line1.id, this.line2.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.LL.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.LL(refs(data[0]), refs(data[1]));
};

TCAD.TWO.Constraints.LL.prototype.getObjects = function() {
  return [this.line1, this.line2];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.Vertical = function(line) {
  this.line = line;
};

TCAD.TWO.Constraints.Vertical.prototype.NAME = 'Vertical';
TCAD.TWO.Constraints.Vertical.prototype.UI_NAME = 'Vertical';
//TCAD.TWO.Constraints.Vertical.prototype.reducible = true;

TCAD.TWO.Constraints.Vertical.prototype.getSolveData = function() {
  return [['equal', [this.line.a._x, this.line.b._x], []]];
};

TCAD.TWO.Constraints.Vertical.prototype.serialize = function() {
  return [this.NAME, [this.line.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.Vertical.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.Vertical(refs(data[0]));
};

TCAD.TWO.Constraints.Vertical.prototype.getObjects = function() {
  return [this.line];
};

// ------------------------------------------------------------------------------------------------------------------ // 

/** @constructor */
TCAD.TWO.Constraints.Horizontal = function(line) {
  this.line = line;
};

TCAD.TWO.Constraints.Horizontal.prototype.NAME = 'Horizontal';
TCAD.TWO.Constraints.Horizontal.prototype.UI_NAME = 'Horizontal';
//TCAD.TWO.Constraints.Horizontal.prototype.reducible = true;

TCAD.TWO.Constraints.Horizontal.prototype.getSolveData = function() {
  return [['equal', [this.line.a._y, this.line.b._y], []]];
};

TCAD.TWO.Constraints.Horizontal.prototype.serialize = function() {
  return [this.NAME, [this.line.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.Horizontal.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.Horizontal(refs(data[0]));
};

TCAD.TWO.Constraints.Horizontal.prototype.getObjects = function() {
  return [this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.Tangent = function(arc, line) {
  this.arc = arc;
  this.line = line;
};

TCAD.TWO.Constraints.Tangent.prototype.NAME = 'Tangent';
TCAD.TWO.Constraints.Tangent.prototype.UI_NAME = 'Tangent';

TCAD.TWO.Constraints.Tangent.prototype.getSolveData = function() {
  var params = [];
  this.arc.c.collectParams(params);
  this.line.collectParams(params);
  params.push(this.arc.r);
  return [['P2LDistanceV', params, []]];
};

TCAD.TWO.Constraints.Tangent.prototype.serialize = function() {
  return [this.NAME, [this.arc.id, this.line.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.Tangent.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.Tangent(refs(data[0]), refs(data[1]));
};

TCAD.TWO.Constraints.Tangent.prototype.getObjects = function() {
  return [this.arc, this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.PointOnLine = function(point, line) {
  this.point = point;
  this.line = line;
};

TCAD.TWO.Constraints.PointOnLine.prototype.NAME = 'PointOnLine';
TCAD.TWO.Constraints.PointOnLine.prototype.UI_NAME = 'Point On Line';

TCAD.TWO.Constraints.PointOnLine.prototype.getSolveData = function() {
  var params = [];
  this.point.collectParams(params);
  this.line.collectParams(params);
  return [['P2LDistance', params, [0]]];
};

TCAD.TWO.Constraints.PointOnLine.prototype.serialize = function() {
  return [this.NAME, [this.point.id, this.line.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.PointOnLine.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.PointOnLine(refs(data[0]), refs(data[1]));
};

TCAD.TWO.Constraints.PointOnLine.prototype.getObjects = function() {
  return [this.point, this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.PointOnArc = function(point, arc) {
  this.point = point;
  this.arc = arc;
};

TCAD.TWO.Constraints.PointOnArc.prototype.NAME = 'PointOnArc';
TCAD.TWO.Constraints.PointOnArc.prototype.UI_NAME = 'Point On Arc';

TCAD.TWO.Constraints.PointOnArc.prototype.getSolveData = function() {
  var params = [];
  this.point.collectParams(params);
  this.arc.c.collectParams(params);
  params.push(this.arc.r);
  return [['P2PDistanceV', params, []]];
};

TCAD.TWO.Constraints.PointOnArc.prototype.serialize = function() {
  return [this.NAME, [this.point.id, this.arc.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.PointOnArc.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.PointOnArc(refs(data[0]), refs(data[1]));
};

TCAD.TWO.Constraints.PointOnArc.prototype.getObjects = function() {
  return [this.point, this.arc];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.PointInMiddle = function(point, line) {
  this.point = point;
  this.line = line;
  this.length = new TCAD.TWO.Ref(TCAD.math.distanceAB(line.a, line.b) / 2);
};

TCAD.TWO.Constraints.PointInMiddle.prototype.NAME = 'PointInMiddle';
TCAD.TWO.Constraints.PointInMiddle.prototype.UI_NAME = 'Point In the Middle';

TCAD.TWO.Constraints.PointInMiddle.prototype.getSolveData = function() {
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

TCAD.TWO.Constraints.PointInMiddle.prototype.serialize = function() {
  return [this.NAME, [this.point.id, this.line.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.PointInMiddle.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.PointInMiddle(refs(data[0]), refs(data[1]));
};

TCAD.TWO.Constraints.PointInMiddle.prototype.getObjects = function() {
  return [this.point, this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.Symmetry = function(point, line) {
  this.point = point;
  this.line = line;
  this.length = new TCAD.TWO.Ref(TCAD.math.distanceAB(line.a, line.b) / 2);
};

TCAD.TWO.Constraints.Symmetry.prototype.NAME = 'Symmetry';
TCAD.TWO.Constraints.Symmetry.prototype.UI_NAME = 'Symmetry';

TCAD.TWO.Constraints.Symmetry.prototype.getSolveData = function(resolver) {
  var pointInMiddleData = TCAD.TWO.Constraints.PointInMiddle.prototype.getSolveData.call(this, [resolver]);
  var pointOnLineData = TCAD.TWO.Constraints.PointOnLine.prototype.getSolveData.call(this, [resolver]);
  return pointInMiddleData.concat(pointOnLineData);
};

TCAD.TWO.Constraints.Symmetry.prototype.serialize = function() {
  return [this.NAME, [this.point.id, this.line.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.Symmetry.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.Symmetry(refs(data[0]), refs(data[1]));
};

TCAD.TWO.Constraints.Symmetry.prototype.getObjects = function() {
  return [this.point, this.line];
};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.Angle = function(p1, p2, p3, p4, angle) {
  this.p1 = p1;
  this.p2 = p2;
  this.p3 = p3;
  this.p4 = p4;
  this._angle = new TCAD.TWO.Ref(0);
  this.angle = angle;
};

TCAD.TWO.Constraints.Angle.prototype.NAME = 'Angle';
TCAD.TWO.Constraints.Angle.prototype.UI_NAME = 'Lines Angle';

TCAD.TWO.Constraints.Angle.prototype.getSolveData = function(resolver) {
  this._angle.set(resolver(this.angle) / 180 * Math.PI);
  var params = [];
  this.p1.collectParams(params);
  this.p2.collectParams(params);
  this.p3.collectParams(params);
  this.p4.collectParams(params);
  params.push(this._angle);
  return [['angleConst', params, []]];
};

TCAD.TWO.Constraints.Angle.prototype.serialize = function() {
  return [this.NAME, [this.p1.id, this.p2.id, this.p3.id, this.p4.id, this.angle]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.Angle.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.Angle( refs(data[0]), refs(data[1]), refs(data[2]), refs(data[3]), data[4] );
};

TCAD.TWO.Constraints.Angle.prototype.getObjects = function() {
  var collector = new TCAD.TWO.Constraints.ParentsCollector();
  collector.check(this.p1);
  collector.check(this.p2);
  collector.check(this.p3);
  collector.check(this.p4);
  return collector.parents;
};

TCAD.TWO.Constraints.Angle.prototype.SettableFields = {'angle' : "Enter the angle value"};

// ------------------------------------------------------------------------------------------------------------------ //

/** @constructor */
TCAD.TWO.Constraints.LockConvex = function(c, a, t) {
  this.c = c;
  this.a = a;
  this.t = t;
};

TCAD.TWO.Constraints.LockConvex.prototype.NAME = 'LockConvex';
TCAD.TWO.Constraints.LockConvex.prototype.UI_NAME = 'Lock Convexity';

TCAD.TWO.Constraints.LockConvex.prototype.getSolveData = function() {
  var params = [];
  this.c.collectParams(params);
  this.a.collectParams(params);
  this.t.collectParams(params);
  return [['LockConvex', params, []]];
};

TCAD.TWO.Constraints.LockConvex.prototype.serialize = function() {
  return [this.NAME, [this.c.id, this.a.id, this.t.id]];
};

TCAD.TWO.Constraints.Factory[TCAD.TWO.Constraints.LockConvex.prototype.NAME] = function(refs, data) {
  return new TCAD.TWO.Constraints.LockConvex(refs(data[0]), refs(data[1]), refs(data[2]));
};

TCAD.TWO.Constraints.LockConvex.prototype.getObjects = function() {
  var collector = new TCAD.TWO.Constraints.ParentsCollector();
  collector.check(this.c);
  collector.check(this.a);
  collector.check(this.t);
  return collector.parents;
};

// ------------------------------------------------------------------------------------------------------------------ //