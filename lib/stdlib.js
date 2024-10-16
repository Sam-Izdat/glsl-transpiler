/**
 * OpenGL/WebGL environment methods.
 *
 * @module  glsl-transpiler/lib/stdlib
 */

import { operators } from './operators.js';

// Returns the type of argument at index
// If argIndex is not specified, assume the first argument
function genType(argIndex) {
    if (Number.isInteger(argIndex)) {
        return function (node) {
            return this.process(node.children[argIndex + 1]).type;
        }
    }
    var node = argIndex;
    // node.children[0] is the function itself
    return this.process(node.children[1]).type;
}

/**
 * Types stubs
 */

function bool    (val) { return !!val; }
function int     (val) { return ti.i32(val); }
function uint    (val) { return ti.i32(val); }
function float   (val) { return ti.f32(val); }
function vec2    (x, y) { return [ti.f32(x), ti.f32(y)]; }
function vec3    (x, y, z) { return [ti.f32(x), ti.f32(y), ti.f32(z)]; }
function vec4    (x, y, z, w) { return [ti.f32(x), ti.f32(y), ti.f32(z), ti.f32(w)]; }
function mat2    (x) {
    if (x == null) x = 1;
    if (x.length === 4) return x;
    if (x.length === 2) return [x[0], 0, 0, x[1]];
    return ti.f32([x, 0, 0, x]);
}
function mat3    (x) {
    if (x == null) x = 1;
    if (x.length === 9) return x;
    if (x.length === 3) return [x[0], 0, 0, 0, x[1], 0, 0, 0, x[2]];
    return ti.f32([x, 0, 0, 0, x, 0, 0, 0, x])
}

function mat4    (x) {
    if (x == null) x = 1;
    if (x.length === 16) return x;
    if (x.length === 4) return [x[0], 0, 0, 0, 0, x[1], 0, 0, 0, 0, x[2], 0, 0, 0, 0, x[3]];
    return ti.f32([x, 0, 0, 0, 0, x, 0, 0, 0, 0, x, 0, 0, 0, 0, x]);
}


/**
 * Types operations.
 */
createOperations(vec2, 2);
createOperations(vec3, 3);
createOperations(vec4, 4);
createOperations(mat2, 4);

function createOperations(obj, len) {
    for (var operator in operators) {
        var comps = [];
        for (var i = 0; i < len; i++) {
            comps.push(`out[${i}] = a[${i}] ${operator} b[${i}]`);
        }

        obj[operators[operator]] = new Function('out', 'a', 'b',
            `${comps.join('\n')}\nreturn out`
        );

        obj[operators[operator]].source = new String(`function ${operators[operator]} (out, a, b) {
            ${comps.join(';\n')};\n
            return out;
        }`)

        obj[operators[operator]].toString =
            obj[operators[operator]].toSource = function () {
                return this.source
            }
    }
}


/**
 * Math
 */


function any           (x)           { return ti.dot(x, tm.ol(x)) > 0; }
function all           (x)           { return ti.dot(x, tm.ol(x)) / ti.len(x) >= 1; }
function log2          (x)           { return ti.log(x) / g.LN2; }
function exp2          (x)           { return ti.pow(2, x);  }
function log10         (x)           { return ti.log(x) / g.LN10; }
function normalize     (v)           { return ti.normalized(v); }
function next_pot      (x)           { return ti.pow(2, ti.ceil(tm.log2(x))); }
function mod           (n, m)        { return n - m * ti.floor(n / m); }
function lshift        (n, b)        { return n * ti.pow(2, b); }
function rshift        (n, b)        { return ti.floor(n / ti.pow(2, b)); }
// function reflect       (v, n)        { return 2. * ti.dot(n, v) * n - v; }
function clamp         (a, b, c)     { return ti.min(ti.max(a, b), c); }
function sign          (x)           { return ti.sgn(x); }
function mix           (a, b, w)     { return a * (1. - w) + b * w; }
function lerp          (a, b, w)     { return a * (1. - w) + b * w; }
function ilerp         (a, b, v)     { return (v - a) / (b - a); }
function remap         (i, o, v)     { return tm.lerp(o[0], o[1], tm.ilerp(i[0], i[1], v)); }   
function smoothstep    (e0, e1, x)   { return tm.cubic_polynomial(tm.saturate(tm.ilerp(e0, e1, x))); }
function fract         (x)           { return x - ti.floor(x); }
function radians       (d)           { return d * 0.017453292519943295; }
function degrees       (r)           { return r * 57.29577951308232; }
function inversesqrt   (x)           { return 1. / ti.sqrt(x); }


radians.type = genType;
degrees.type = genType;
// sin.type = genType;
// cos.type = genType;
// tan.type = genType;
// asin.type = genType;
// acos.type = genType;
// atan.type = genType;
// pow.type = genType;
// exp.type = genType;
// log.type = genType;
log2.type = genType;
exp2.type = genType;
// sqrt.type = genType;
inversesqrt.type = genType;
// abs.type = genType;
// floor.type = genType;
// ceil.type = genType;
// sign.type = genType;
fract.type = genType;
mod.type = genType;
// min.type = genType;
// max.type = genType;
clamp.type = genType;
mix.type = genType;
// step.type = genType(1);
smoothstep.type = genType(2);
// length.type = 'float';
// distance.type = 'float';
// dot.type = 'float';
// cross.type = 'vec3';
normalize.type = genType;

function faceforward(N, I, Nref) {
    if (Nref == null) Nref = N;

    var dot = 0;
    for (var i = 0; i < N.length; i++) {
        dot += Nref[i] * I[i];
    }

    return dot > 0 ? N.map(function (x) { return -x; }) : N;
}
faceforward.type = genType;

function reflect(I, N) {
    var dot = 0;
    for (var i = 0; i < N.length; i++) {
        dot += N[i] * I[i];
    }

    var out = Array(N.length);
    for (var i = 0; i < N.length; i++) {
        out[i] = I[i] - 2 * dot * N[i];
    }

    return out;
}
reflect.type = genType;

function refract(I, N, eta) {
    var dot = 0;
    for (var i = 0; i < N.length; i++) {
        dot += N[i] * I[i];
    }

    var k = 1 - eta * eta * (1 - dot * dot);

    var out = Array(N.length).fill(0);

    if (k > 0) {
        for (var i = 0; i < N.length; i++) {
            out[i] = eta * I[i] - (eta * dot + Math.sqrt(k)) * N[i];
        }
    }

    return out;
}
refract.type = genType;

/**
 * Vector relational functions
 */
function lessThan(x, y) {
    if (x.length) {
        var out = Array(x.length);
        for (var i = 0; i < x.length; i++) {
            out[i] = x[i] < y[i];
        }
        return out;
    }
    return x < y;
}

function lessThanEqual(x, y) {
    if (x.length) {
        var out = Array(x.length);
        for (var i = 0; i < x.length; i++) {
            out[i] = x[i] <= y[i];
        }
        return out;
    }
    return x <= y;
}

function greaterThan(x, y) {
    if (x.length) {
        var out = Array(x.length);
        for (var i = 0; i < x.length; i++) {
            out[i] = x[i] > y[i];
        }
        return out;
    }
    return x > y;
}

function greaterThanEqual(x, y) {
    if (x.length) {
        var out = Array(x.length);
        for (var i = 0; i < x.length; i++) {
            out[i] = x[i] >= y[i];
        }
        return out;
    }
    return x >= y;
}

function equal(x, y) {
    if (x.length) {
        var out = Array(x.length);
        for (var i = 0; i < x.length; i++) {
            out[i] = x[i] == y[i];
        }
        return out;
    }
    return x == y;
}

function notEqual(x, y) {
    if (x.length) {
        var out = Array(x.length);
        for (var i = 0; i < x.length; i++) {
            out[i] = x[i] != y[i];
        }
        return out;
    }
    return x != y;
}

// function any(x) {
//     return x.some(function (x) { return x; });
// }

// function all(x) {
//     return x.every(function (x) { return x; });
// }

function not(x) {
    if (x.length) {
        var out = Array(x.length);
        for (var i = 0; i < x.length; i++) {
            out[i] = !x[i];
        }
        return out;
    }
    return !x
}


/**
 * Matrices
 */
function matrixCompMult(x, y) {
    var out = Array(x.length);
    for (var i = 0; i < x.length; i++) {
        out[i] = x[i] * y[i];
    }
    return out;
}

function outerProduct(c, r) {
    var out = [];
    var l = c.length;
    for (var i = 0; i < c.length; i++) {
        for (var j = 0; j < r.length; j++) {
            out[i * l + j] = c[i] * r[j];
        }
    }
    return out;
}
outerProduct.type = function (node) {
    var child1Type = this.process(node.children[1]).type;
    var child2Type = this.process(node.children[2]).type;
    var dim1 = child1Type.slice(-1);
    var dim2 = child2Type.slice(-1);
    return `mat${dim1}x${dim2}`;
};

function transpose(m) {
    var l = m.length === 16 ? 4 : m.length === 9 ? 3 : 2;
    var out = Array(m.length);
    for (var i = 0; i < l; i++) {
        for (var j = 0; j < l; j++) {
            out[j * l + i] = m[i * l + j];
        }
    }
    return out;
}

function determinant(m) {
    if (m.length === 4) {
        return m[0] * m[3] - m[1] * m[2];
    }

    if (m.length === 9) {
        var a00 = m[0], a01 = m[1], a02 = m[2], a10 = m[3], a11 = m[4], a12 = m[5], a20 = m[6], a21 = m[7], a22 = m[8];

        return a00 * a11 * a22 + a01 * a12 * a20 + a02 * a10 * a21 - a02 * a11 * a20 - a01 * a10 * a22 - a00 * a12 * a21;
    }

    var a00 = m[0], a01 = m[1], a02 = m[2], a03 = m[3],
        a10 = m[4], a11 = m[5], a12 = m[6], a13 = m[7],
        a20 = m[8], a21 = m[9], a22 = m[10], a23 = m[11],
        a30 = m[12], a31 = m[13], a32 = m[14], a33 = m[15],

        b00 = a00 * a11 - a01 * a10,
        b01 = a00 * a12 - a02 * a10,
        b02 = a00 * a13 - a03 * a10,
        b03 = a01 * a12 - a02 * a11,
        b04 = a01 * a13 - a03 * a11,
        b05 = a02 * a13 - a03 * a12,
        b06 = a20 * a31 - a21 * a30,
        b07 = a20 * a32 - a22 * a30,
        b08 = a20 * a33 - a23 * a30,
        b09 = a21 * a32 - a22 * a31,
        b10 = a21 * a33 - a23 * a31,
        b11 = a22 * a33 - a23 * a32;

    return b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;
}
determinant.type = 'float';

//FIXME: optimize the method inclusion, per-matrix
//FIXME: inverse the dimensions of the input matrix: mat2x3 → mat3x2
function inverse(a) {
    var l = a.length;
    var out = Array(l);

    if (l === 4) {
        var a0 = a[0], a1 = a[1], a2 = a[2], a3 = a[3],

            det = a0 * a3 - a2 * a1;

        if (!det) {
            return out;
        }
        det = 1.0 / det;

        out[0] = a3 * det;
        out[1] = -a1 * det;
        out[2] = -a2 * det;
        out[3] = a0 * det;

        return out;
    }

    if (l === 9) {
        var a00 = a[0], a01 = a[1], a02 = a[2],
            a10 = a[3], a11 = a[4], a12 = a[5],
            a20 = a[6], a21 = a[7], a22 = a[8],

            b01 = a22 * a11 - a12 * a21,
            b11 = -a22 * a10 + a12 * a20,
            b21 = a21 * a10 - a11 * a20,

            det = a00 * b01 + a01 * b11 + a02 * b21;

        if (!det) {
            return out;
        }
        det = 1.0 / det;

        out[0] = b01 * det;
        out[1] = (-a22 * a01 + a02 * a21) * det;
        out[2] = (a12 * a01 - a02 * a11) * det;
        out[3] = b11 * det;
        out[4] = (a22 * a00 - a02 * a20) * det;
        out[5] = (-a12 * a00 + a02 * a10) * det;
        out[6] = b21 * det;
        out[7] = (-a21 * a00 + a01 * a20) * det;
        out[8] = (a11 * a00 - a01 * a10) * det;
        return out;
    }

    var a00 = a[0], a01 = a[1], a02 = a[2], a03 = a[3],
        a10 = a[4], a11 = a[5], a12 = a[6], a13 = a[7],
        a20 = a[8], a21 = a[9], a22 = a[10], a23 = a[11],
        a30 = a[12], a31 = a[13], a32 = a[14], a33 = a[15],

        b00 = a00 * a11 - a01 * a10,
        b01 = a00 * a12 - a02 * a10,
        b02 = a00 * a13 - a03 * a10,
        b03 = a01 * a12 - a02 * a11,
        b04 = a01 * a13 - a03 * a11,
        b05 = a02 * a13 - a03 * a12,
        b06 = a20 * a31 - a21 * a30,
        b07 = a20 * a32 - a22 * a30,
        b08 = a20 * a33 - a23 * a30,
        b09 = a21 * a32 - a22 * a31,
        b10 = a21 * a33 - a23 * a31,
        b11 = a22 * a33 - a23 * a32,

        det = b00 * b11 - b01 * b10 + b02 * b09 + b03 * b08 - b04 * b07 + b05 * b06;

    if (!det) {
        return out;
    }
    det = 1.0 / det;

    out[0] = (a11 * b11 - a12 * b10 + a13 * b09) * det;
    out[1] = (a02 * b10 - a01 * b11 - a03 * b09) * det;
    out[2] = (a31 * b05 - a32 * b04 + a33 * b03) * det;
    out[3] = (a22 * b04 - a21 * b05 - a23 * b03) * det;
    out[4] = (a12 * b08 - a10 * b11 - a13 * b07) * det;
    out[5] = (a00 * b11 - a02 * b08 + a03 * b07) * det;
    out[6] = (a32 * b02 - a30 * b05 - a33 * b01) * det;
    out[7] = (a20 * b05 - a22 * b02 + a23 * b01) * det;
    out[8] = (a10 * b10 - a11 * b08 + a13 * b06) * det;
    out[9] = (a01 * b08 - a00 * b10 - a03 * b06) * det;
    out[10] = (a30 * b04 - a31 * b02 + a33 * b00) * det;
    out[11] = (a21 * b02 - a20 * b04 - a23 * b00) * det;
    out[12] = (a11 * b07 - a10 * b09 - a12 * b06) * det;
    out[13] = (a00 * b09 - a01 * b07 + a02 * b06) * det;
    out[14] = (a31 * b01 - a30 * b03 - a32 * b00) * det;
    out[15] = (a20 * b03 - a21 * b01 + a22 * b00) * det;

    return out;
}

/**
 * mat * mat
 */
function matrixMult(m, n) {
    var l = m.length === 16 ? 4 : m.length === 9 ? 3 : 2;
    var out = Array(m.length);
    for (var i = 0; i < l; i++) {
        for (var j = 0; j < l; j++) {
            var sum = 0;
            for (var o = 0; o < l; o++) {
                sum += m[l * o + i] * n[j * l + o];
            }
            out[j * l + i] = sum;
        }
    }
    return out;
}


/**
 * Get texture value.
 * It has the output type of first arg.
 */
function texture(sampler, coord, bias) {
    var size = textureSize(sampler);
    var eps = 1e-6;
    var x = (Math.min(1 - eps, Math.max(0, coord[0])) * size[0]) | 0;
    var y = (Math.min(1 - eps, Math.max(0, coord[1])) * size[1]) | 0;
    var idx = y * 4 * size[0] + x * 4;
    if (sampler.data) {
        return sampler.data.slice(idx, idx + 4);
    }
    return sampler.slice(idx, idx + 4);
}
texture.include = ['textureSize'];
texture.type = function (node) {
    var samplerType = this.process(node.children[1]).type;
    return this.types[samplerType].type;
};

function textureSize(sampler, lod) {
    if (sampler.shape) return [sampler.shape[0], sampler.shape[1]];
    return [sampler.width, sampler.height];
};
textureSize.type = function (node) {
    var samplerType = this.process(node.children[1]).type;
    if (/1D/.test(samplerType)) return 'int';
    if (/2D|Cube/.test(samplerType)) return 'ivec2';
    return 'ivec3';
};

const stdlib = {}

stdlib.bool = bool;
stdlib.int = int;
stdlib.uint = uint;
stdlib.float = float;
stdlib.double = float;
stdlib.vec2 = vec2;
stdlib.vec3 = vec3;
stdlib.vec4 = vec4;
stdlib.dvec2 = vec2;
stdlib.dvec3 = vec3;
stdlib.dvec4 = vec4;
stdlib.ivec2 = vec2;
stdlib.ivec3 = vec3;
stdlib.ivec4 = vec4;
stdlib.uvec2 = vec2;
stdlib.uvec3 = vec3;
stdlib.uvec4 = vec4;
stdlib.mat2 = mat2;
stdlib.mat3 = mat3;
stdlib.mat4 = mat4;
stdlib.mat3x3 = mat3;
stdlib.mat4x4 = mat4;

stdlib.radians = radians;
stdlib.degrees = degrees;
// stdlib.sin = sin;
// stdlib.cos = cos;
// stdlib.tan = tan;
// stdlib.asin = asin;
// stdlib.acos = acos;
// stdlib.atan = atan;
// stdlib.pow = pow;
// stdlib.exp = exp;
// stdlib.log = log;
stdlib.log2 = log2;
stdlib.exp2 = exp2;
// stdlib.sqrt = sqrt;
stdlib.inversesqrt = inversesqrt;
// stdlib.abs = abs;
// stdlib.sign = sign;
// stdlib.floor = floor;
// stdlib.ceil = ceil;
stdlib.fract = fract;
stdlib.mod = mod;
// stdlib.min = min;
// stdlib.max = max;
stdlib.clamp = clamp;
stdlib.mix = mix;
// stdlib.step = step;
stdlib.smoothstep = smoothstep;
// stdlib.length = length;
// stdlib.distance = distance;
// stdlib.dot = dot;
// stdlib.cross = cross;
stdlib.faceforward = faceforward;
stdlib.normalize = normalize;
// stdlib.reflect = reflect;
// stdlib.refract = refract;
stdlib.lessThan = lessThan;
stdlib.lessThanEqual = lessThanEqual;
stdlib.greaterThan = greaterThan;
stdlib.greaterThanEqual = greaterThanEqual;
stdlib.equal = equal;
stdlib.notEqual = notEqual;
stdlib.any = any;
stdlib.all = all;
stdlib.not = not;
stdlib.matrixCompMult = matrixCompMult;
stdlib.matrixMult = matrixMult;
stdlib.outerProduct = outerProduct;
stdlib.transpose = transpose;
stdlib.determinant = determinant;
stdlib.inverse = inverse;

stdlib.texture1D =
    stdlib.texture2D =
    stdlib.texture3D =
    stdlib.textureCube =
    stdlib.shadow1D =
    stdlib.shadow2D =
    stdlib.shadow3D =
    stdlib.texture = texture;
stdlib.textureSize = textureSize;

export default stdlib;
