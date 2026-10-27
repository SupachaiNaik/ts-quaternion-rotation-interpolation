/**
 * ts-quaternion-rotation-interpolation
 *
 * Spherical linear interpolation between unit quaternions and
 * composition of rotation sequences for 3D orientations.
 *
 * Storage convention: [x, y, z, w] where (x, y, z) is the vector
 * part and w is the scalar part. Quaternion multiplication follows
 * the Hamilton product so that for unit quaternions q1 * q2
 * represents applying q2 first, then q1 (right-to-left, same as
 * column-major rotation matrices).
 */

/**
 * Plain quaternion value: [x, y, z, w].
 * @typedef {[number, number, number, number]} Quaternion
 */

/**
 * Returns the Euclidean norm (magnitude) of a quaternion.
 * Used internally to decide whether normalization is safe.
 *
 * @param {Quaternion} q
 * @returns {number}
 */
export function norm(q) {
  return Math.sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
}

/**
 * Returns a unit-length copy of q. If the input is the zero
 * quaternion (all components zero) we return the identity
 * quaternion [0,0,0,1] rather than NaN — the identity is the
 * only rotation that is a safe neutral element, and producing
 * NaN here would propagate silently through every later call.
 *
 * @param {Quaternion} q
 * @returns {Quaternion}
 */
export function normalize(q) {
  const n = norm(q);
  if (n < 1e-12) {
    return [0, 0, 0, 1];
  }
  const inv = 1 / n;
  return [q[0] * inv, q[1] * inv, q[2] * inv, q[3] * inv];
}

/**
 * Hamilton product of two quaternions.
 *
 * The sign convention is the standard one: (a + bi + cj + dk)
 * with i*i = j*j = k*k = i*j*k = -1. This makes the product
 * of two unit quaternions correspond to the composition of
 * their rotations.
 *
 * @param {Quaternion} a
 * @param {Quaternion} b
 * @returns {Quaternion}
 */
export function multiply(a, b) {
  const ax = a[0], ay = a[1], az = a[2], aw = a[3];
  const bx = b[0], by = b[1], bz = b[2], bw = b[3];
  return [
    aw * bx + ax * bw + ay * bz - az * by,
    aw * by - ax * bz + ay * bw + az * bx,
    aw * bz + ax * by - ay * bx + az * bw,
    aw * bw - ax * bx - ay * by - az * bz,
  ];
}

/**
 * Conjugate of a quaternion (negated vector part, same scalar).
 * For a unit quaternion the conjugate equals the inverse and
 * represents the opposite rotation.
 *
 * @param {Quaternion} q
 * @returns {Quaternion}
 */
export function conjugate(q) {
  return [-q[0], -q[1], -q[2], q[3]];
}

/**
 * Dot product of two quaternions, treating each as a 4-vector.
 * Used by slerp to decide whether to negate the second
 * quaternion to take the shorter arc.
 *
 * @param {Quaternion} a
 * @param {Quaternion} b
 * @returns {number}
 */
export function dot(a, b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + a[3] * b[3];
}

/**
 * Builds a unit quaternion representing a rotation of `angle`
 * radians about the unit axis `axis`. If the axis is not unit
 * length it is normalized first; if it is near zero we return
 * the identity, because a rotation about no axis is no rotation.
 *
 * @param {[number, number, number]} axis - rotation axis, need not be unit length
 * @param {number} angle - rotation angle in radians
 * @returns {Quaternion}
 */
export function fromAxisAngle(axis, angle) {
  const len = Math.sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
  if (len < 1e-12) {
    return [0, 0, 0, 1];
  }
  const half = angle / 2;
  const s = Math.sin(half);
  const inv = 1 / len;
  return [axis[0] * inv * s, axis[1] * inv * s, axis[2] * inv * s, Math.cos(half)];
}

/**
 * Converts a unit quaternion to its axis-angle representation.
 *
 * Edge case: when the rotation angle is effectively zero the
 * vector part vanishes and there is no well-defined axis. We
 * return the x-axis [1, 0, 0] in that case, because returning
 * a zero vector would be indistinguishable from "no input" and
 * would break downstream code that normalizes the axis.
 *
 * @param {Quaternion} q
 * @returns {{axis: [number, number, number], angle: number}}
 */
export function toAxisAngle(q) {
  // Work from a normalized copy so the result is independent of
  // the caller having passed a non-unit quaternion.
  const n = normalize(q);
  const x = n[0], y = n[1], z = n[2], w = n[3];

  const vLen = Math.sqrt(x * x + y * y + z * z);
  if (vLen < 1e-12) {
    // Near-identity rotation: angle is ~0, axis undefined.
    return { axis: [1, 0, 0], angle: 0 };
  }

  // Clamp w to the legal range for acos to avoid NaN from
  // floating-point drift slightly outside [-1, 1].
  const clamped = w < -1 ? -1 : w > 1 ? 1 : w;
  const angle = 2 * Math.acos(clamped);
  const inv = 1 / vLen;
  return { axis: [x * inv, y * inv, z * inv], angle };
}

/**
 * Spherical linear interpolation between two unit quaternions.
 *
 * t=0 returns a; t=1 returns b. The result is normalized before
 * being returned, so small floating-point drift in the
 * interpolation does not accumulate when slerp is called
 * repeatedly in a loop.
 *
 * The shorter-arc decision: if the dot product of a and b is
 * negative, we interpolate from a to -b. Both q and -q
 * represent the same orientation, so this picks the path that
 * travels the smaller angular distance and avoids the visually
 * jarring "long way around" spin.
 *
 * When the two quaternions are nearly identical the standard
 * slerp formula divides by sin(theta) ~ 0 and becomes unstable.
 * In that region we fall back to normalized linear interpolation
 * (nlerp), which is first-order accurate and perfectly stable.
 *
 * @param {Quaternion} a - start quaternion (need not be unit; normalized internally)
 * @param {Quaternion} b - end quaternion (need not be unit; normalized internally)
 * @param {number} t - interpolation parameter, clamped to [0, 1]
 * @returns {Quaternion}
 */
export function slerp(a, b, t) {
  // Clamp t so callers cannot accidentally extrapolate.
  const tc = t < 0 ? 0 : t > 1 ? 1 : t;

  const qa = normalize(a);
  let qb = normalize(b);

  let d = dot(qa, qb);

  // Take the shorter arc.
  if (d < 0) {
    qb = [-qb[0], -qb[1], -qb[2], -qb[3]];
    d = -d;
  }

  // If the quaternions are extremely close, nlerp is stable
  // where slerp is not. The threshold 1e-6 (about 0.00005
  // degrees of arc) is small enough that the visual difference
  // between nlerp and slerp is well below any perceptual or
  // numerical threshold.
  const THRESHOLD = 0.999999;
  if (d > THRESHOLD) {
    const r = [
      qa[0] + tc * (qb[0] - qa[0]),
      qa[1] + tc * (qb[1] - qa[1]),
      qa[2] + tc * (qb[2] - qa[2]),
      qa[3] + tc * (qb[3] - qa[3]),
    ];
    return normalize(r);
  }

  // Clamp d into the domain of acos for safety.
  const dc = d < -1 ? -1 : d > 1 ? 1 : d;
  const theta = Math.acos(dc);
  const sinTheta = Math.sin(theta);
  const s0 = Math.sin((1 - tc) * theta) / sinTheta;
  const s1 = Math.sin(tc * theta) / sinTheta;

  const r = [
    s0 * qa[0] + s1 * qb[0],
    s0 * qa[1] + s1 * qb[1],
    s0 * qa[2] + s1 * qb[2],
    s0 * qa[3] + s1 * qb[3],
  ];
  return normalize(r);
}

/**
 * Composes a sequence of rotations into a single unit quaternion.
 *
 * The rotations are applied in the order they appear: for
 * `compose([r1, r2, r3])` the combined rotation is r1 * r2 * r3,
 * meaning r3 is applied first to a vector and r1 last. This
 * matches the convention used by `multiply` and by most
 * column-major matrix libraries.
 *
 * An empty list returns the identity quaternion, since composing
 * no rotations leaves orientation unchanged. Each input is
 * normalized so that a non-unit input does not make the product
 * drift away from unit length.
 *
 * @param {Quaternion[]} quats
 * @returns {Quaternion}
 */
export function compose(quats) {
  if (!Array.isArray(quats) || quats.length === 0) {
    return [0, 0, 0, 1];
  }
  let acc = normalize(quats[0]);
  for (let i = 1; i < quats.length; i++) {
    acc = multiply(acc, normalize(quats[i]));
  }
  return normalize(acc);
}
