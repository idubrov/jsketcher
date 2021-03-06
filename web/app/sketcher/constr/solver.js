TCAD.parametric = {};

/** @constructor */
TCAD.parametric.Param = function(id, value, readOnly) {
  this.reset(value);
};

TCAD.parametric.Param.prototype.reset = function(value) {
  this.set(value);
  this.j = -1;
};

TCAD.parametric.Param.prototype.set = function(value) {
  this.value = value;
};

TCAD.parametric.Param.prototype.get = function() {
  return this.value;
};

TCAD.parametric.Param.prototype.nop = function() {};

/** @constructor */
TCAD.parametric.System = function(constraints) {
  this.constraints = constraints;
  this.params = [];
  for (var ci = 0; ci < constraints.length; ++ci) {
    var c = constraints[ci];
    for (var pi = 0; pi < c.params.length; ++pi) {
      var p = c.params[pi];
      if (p.j == -1) {
        p.j = this.params.length;
        this.params.push(p);
      }
    }
  }
};


TCAD.parametric.System.prototype.makeJacobian = function() {
  var jacobi = [];
  var i;
  var j;
  for (i=0; i < this.constraints.length; i++) {
    jacobi[i] = [];
    for (j=0; j < this.params.length; j++) {
      jacobi[i][j] = 0;
    }
  }
  for (i=0; i < this.constraints.length; i++) {
    var c = this.constraints[i];

    var cParams = c.params;
    var grad = [];
    c.gradient(grad);

    for (var p = 0; p < cParams.length; p++) {
      var param = cParams[p];
      j = param.j;
      jacobi[i][j] = grad[p];
    }
  }
  return jacobi;
};

TCAD.parametric.System.prototype.fillJacobian = function(jacobi) {
  for (i=0; i < this.constraints.length; i++) {
    var c = this.constraints[i];

    var cParams = c.params;
    var grad = [];
    c.gradient(grad);

    for (var p = 0; p < cParams.length; p++) {
      var param = cParams[p];
      j = param.j;
      jacobi[i][j] = grad[p];
    }
  }
  return jacobi;
};

TCAD.parametric.System.prototype.calcResidual = function(r) {

  var i=0;
  var err = 0.;

  for (i=0; i < this.constraints.length; i++) {
    var c = this.constraints[i];
    r[i] = c.error();
    err += r[i]*r[i];
  }

  err *= 0.5;
  return err;
};

TCAD.parametric.System.prototype.calcGrad_ = function(out) {
  var i;
  for (i = 0; i < out.length || i < this.params.length; ++i) {
    out[i][0] = 0;
  }

  for (i=0; i < this.constraints.length; i++) {
    var c = this.constraints[i];

    var cParams = c.params;
    var grad = [];
    c.gradient(grad);

    for (var p = 0; p < cParams.length; p++) {
      var param = cParams[p];
      var j = param.j;
      out[j][0] += this.constraints[i].error() * grad[p]; // (10.4)
    }
  }
};

TCAD.parametric.System.prototype.calcGrad = function(out) {
  var i;
  for (i = 0; i < out.length || i < this.params.length;  ++i) {
    out[i] = 0;
  }

  for (i=0; i < this.constraints.length; i++) {
    var c = this.constraints[i];

    var cParams = c.params;
    var grad = [];
    c.gradient(grad);

    for (var p = 0; p < cParams.length; p++) {
      var param = cParams[p];
      var j = param.j;
      out[j] += this.constraints[i].error() * grad[p]; // (10.4) 
    }
  }
};

TCAD.parametric.System.prototype.fillParams = function(out) {
  for (var p = 0; p < this.params.length; p++) {
    out[p] = this.params[p].get();
  }
};

TCAD.parametric.System.prototype.getParams = function() {
  var out = [];
  this.fillParams(out);
  return out;
};

TCAD.parametric.System.prototype.setParams = function(point) {
  for (var p = 0; p < this.params.length; p++) {
    this.params[p].set(point[p]);
  }
};

TCAD.parametric.System.prototype.error = function() {
  var error = 0;
  for (var i=0; i < this.constraints.length; i++) {
    error += Math.abs(this.constraints[i].error());
  }
  return error;
};

TCAD.parametric.System.prototype.errorSquare = function() {
  var error = 0;
  for (var i=0; i < this.constraints.length; i++) {
    var t = this.constraints[i].error();
    error += t * t;
  }
  return error * 0.5;
};

TCAD.parametric.System.prototype.getValues = function() {
  var values = [];
  for (var i=0; i < this.constraints.length; i++) {
    values[i] = this.constraints[i].error();
  }
  return values;
};

TCAD.parametric.wrapAux = function(constrs, locked) {

  var lockedSet = {};
  for (var i = 0; i < locked.length; i++) {
    lockedSet[locked[i].j] = true;
  }

  for (var i = 0; i < constrs.length; i++) {
    var c = constrs[i];
    var mask = [];
    var needWrap = false;
    for (var j = 0; j < c.params.length; j++) {
      var param = c.params[j];
      mask[j] = lockedSet[param.j] === true;
      needWrap = needWrap || mask[j];
    }
    if (needWrap) {
      var wrapper = new TCAD.constraints.ConstantWrapper(c, mask);
      constrs[i] = wrapper;
    }
  }
};

TCAD.parametric.lock2Equals2 = function(constrs, locked) {
  var _locked = [];
  for (var i = 0; i < locked.length; ++i) {
    _locked.push(new TCAD.constraints.EqualsTo([locked[i]], locked[i].get()));
  }
  return _locked;
};

TCAD.parametric.diagnose = function(sys) {
  if (sys.constraints.length == 0 || sys.params.length == 0) {
    return {
      conflict : false,
      dof : 0
    }
  }
  var jacobian = sys.makeJacobian();
  var qr = new TCAD.math.QR(jacobian);
  return {
    conflict : sys.constraints.length > qr.rank,
    dof : sys.params.length - qr.rank
  }
};

TCAD.parametric.prepare = function(constrs, locked, aux, alg) {

  var simpleMode = true;
  if (!simpleMode) {
    var lockingConstrs = TCAD.parametric.lock2Equals2(constrs, locked);
    Array.prototype.push.apply( constrs, lockingConstrs );
  }
            
  var sys = new TCAD.parametric.System(constrs);
  
  TCAD.parametric.wrapAux(constrs, aux);

  var model = function(point) {
    sys.setParams(point);
    return sys.getValues();
  };

  var jacobian = function(point) {
    sys.setParams(point);
    return sys.makeJacobian();
  };
  var nullResult = {
    evalCount : 0,
    error : 0,
    returnCode : 1
  };

  function solve(rough, alg) {
    //if (simpleMode) return nullResult;
    if (constrs.length == 0) return nullResult;
    if (sys.params.length == 0) return nullResult;
    switch (alg) {
      case 2:
        return TCAD.parametric.solve_lm(sys, model, jacobian, rough);
      case 1:
      default:    
        return optim.dog_leg(sys, rough);
    }
  }
  var systemSolver = {
    diagnose : function() {return TCAD.parametric.diagnose(sys)},
    error : function() {return sys.error()},
    solveSystem : solve,
    system : sys,
    updateLock : function(values) {
      for (var i = 0; i < values.length; ++i) {
        if (simpleMode) {
          locked[i].set(values[i]);
        } else {
          lockingConstrs[i].value = values[i];
        }
      }
    }
  };
  return systemSolver;
};

TCAD.parametric.solve_lm = function(sys, model, jacobian, rough) {
  var opt = new LMOptimizer(sys.getParams(), TCAD.math.vec(sys.constraints.length), model, jacobian);
  opt.evalMaximalCount = 100 * sys.params.length;
  var eps = rough ? 0.001 : 0.00000001;
  opt.init0(eps, eps, eps);
  var returnCode = 1;
  try {
    var res = opt.doOptimize();
  } catch (e) {
    returnCode = 2;
  }
  sys.setParams(res[0]);
  return {
    evalCount : opt.evalCount,
    error : sys.error(),
    returnCode : returnCode
  };
};

