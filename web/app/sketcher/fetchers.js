TCAD.TWO.ParametricManager.prototype._fetchTwoPoints = function(objs) {
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

TCAD.TWO.ParametricManager.prototype._fetchPoints = function(objs) {
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

TCAD.TWO.ParametricManager.prototype._fetchArkCirc = function(objs, min) {
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

TCAD.TWO.ParametricManager.prototype._fetch = function(objs, types, min) {
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

TCAD.TWO.ParametricManager.prototype._fetchPointAndLine = function(objs) {

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

TCAD.TWO.ParametricManager.prototype._fetchLine = function(objs) {
  for (var i = 0; i < objs.length; ++i) {
    if (objs[i]._class == 'TCAD.TWO.Segment') {
      return objs[i];
    }
  }
  throw "Illegal Argument. Constraint requires a line."
};

TCAD.TWO.ParametricManager.prototype._fetchArcCircAndLine = function(objs) {

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

TCAD.TWO.ParametricManager.prototype._fetchTwoLines = function(objs) {
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

TCAD.TWO.utils._fetchSketchObjects = function(objs, silent, matching) {
  var fetched = [];
  for (var i = 0; i < objs.length; ++i) {
    for (var j = 0; j < matching.length; j++) {
      if (objs[i]._class ==  matching[j]) {
        fetched[j] = objs[i]; 
        matching[j] = null;
      }
    }
  }
  if (fetched.length != matching.length) {
    if (silent) {
      return null;
    } else {
      throw "Illegal Argument. " + matching + " required";
    }
  }
  return fetched;
};

