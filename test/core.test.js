import { test } from 'node:test';
import assert from 'node:assert/strict';

import {
  slerp,
  compose,
  normalize,
  multiply,
  conjugate,
  fromAxisAngle,
  toAxisAngle,
  norm,
  dot,
} from '../src/core.js';

const EPS = 1e-9;

function approxEqual(a, b, eps = EPS) {
  return Math.abs(a - b) < eps;
}

function quatApproxEqual(a, b, eps = EPS) {
  for (let i = 0; i < 4; i++) {
    if (!approxEqual(a[i], b[i], eps)) return false;
  }
  return true;
}

// Identity and two basis rotations reused across tests.
const IDENTITY = [0, 0, 0, 1];
const ROT_X_90 = fromAxisAngle([1, 0, 0], Math.PI / 2);
const ROT_Y_90 = fromAxisAngle([0, 1, 0], Math.PI / 2);

test('normalize returns identity for the zero quaternion', () => {
  assert.deepEqual(normalize([0, 0, 0, 0]), IDENTITY);
});

test('normalize scales a non-unit quaternion to unit length', () => {
  const q = normalize([0, 0, 3, 4]);
  assert.ok(approxEqual(norm(q), 1));
  assert.ok(quatApproxEqual(q, [0, 0, 0.6, 0.8]));
});

test('fromAxisAngle normalizes a non-unit axis', () => {
  const a = fromAxisAngle([0, 0, 5], Math.PI / 2);
  const b = fromAxisAngle([0, 0, 1], Math.PI / 2);
  assert.ok(quatApproxEqual(a, b));
});

test('fromAxisAngle with zero-length axis returns identity', () => {
  assert.deepEqual(fromAxisAngle([0, 0, 0], Math.PI / 2), IDENTITY);
});

test('fromAxisAngle at angle 0 returns identity', () => {
  assert.ok(quatApproxEqual(fromAxisAngle([1, 0, 0], 0), IDENTITY));
});

test('toAxisAngle recovers axis and angle for a 180-degree rotation about y', () => {
  const q = fromAxisAngle([0, 1, 0], Math.PI);
  const { axis, angle } = toAxisAngle(q);
  assert.ok(approxEqual(angle, Math.PI));
  assert.ok(approxEqual(Math.abs(axis[0]), 0, 1e-7));
  assert.ok(approxEqual(Math.abs(axis[1]), 1, 1e-7));
  assert.ok(approxEqual(Math.abs(axis[2]), 0, 1e-7));
});

test('toAxisAngle returns x-axis and zero angle for identity', () => {
  const { axis, angle } = toAxisAngle(IDENTITY);
  assert.deepEqual(axis, [1, 0, 0]);
  assert.ok(approxEqual(angle, 0));
});

test('toAxisAngle handles a non-unit input quaternion', () => {
  const q = fromAxisAngle([0, 0, 1], Math.PI / 3);
  const scaled = [q[0] * 5, q[1] * 5, q[2] * 5, q[3] * 5];
  const { axis, angle } = toAxisAngle(scaled);
  assert.ok(approxEqual(angle, Math.PI / 3));
  assert.ok(approxEqual(Math.abs(axis[0]), 0, 1e-7));
  assert.ok(approxEqual(Math.abs(axis[1]), 0, 1e-7));
  assert.ok(approxEqual(Math.abs(axis[2]), 1, 1e-7));
});

test('slerp endpoints: t=0 returns a, t=1 returns b', () => {
  const a = ROT_X_90;
  const b = ROT_Y_90;
  assert.ok(quatApproxEqual(slerp(a, b, 0), a));
  assert.ok(quatApproxEqual(slerp(a, b, 1), b));
});

test('slerp midpoint between identity and a 90-degree z rotation is a 45-degree z rotation', () => {
  const end = fromAxisAngle([0, 0, 1], Math.PI / 2);
  const mid = slerp(IDENTITY, end, 0.5);
  const expected = fromAxisAngle([0, 0, 1], Math.PI / 4);
  assert.ok(quatApproxEqual(mid, expected, 1e-7));
});

test('slerp takes the shorter arc when dot product is negative', () => {
  // end is -q, same orientation but the long way around.
  const end = fromAxisAngle([0, 0, 1], Math.PI / 4);
  const negEnd = [-end[0], -end[1], -end[2], -end[3]];
  const result = slerp(IDENTITY, negEnd, 0.5);
  const expected = fromAxisAngle([0, 0, 1], Math.PI / 8);
  assert.ok(quatApproxEqual(result, expected, 1e-7));
});

test('slerp between identical quaternions returns that quaternion', () => {
  const q = ROT_X_90;
  assert.ok(quatApproxEqual(slerp(q, q, 0.5), q));
});

test('slerp clamps t outside [0, 1]', () => {
  const a = ROT_X_90;
  const b = ROT_Y_90;
  assert.ok(quatApproxEqual(slerp(a, b, -0.5), a));
  assert.ok(quatApproxEqual(slerp(a, b, 1.5), b));
});

test('slerp returns a unit quaternion', () => {
  const a = ROT_X_90;
  const b = ROT_Y_90;
  const r = slerp(a, b, 0.37);
  assert.ok(approxEqual(norm(r), 1));
});

test('slerp handles non-unit inputs', () => {
  const end = fromAxisAngle([0, 0, 1], Math.PI / 2);
  const bigEnd = [end[0] * 3, end[1] * 3, end[2] * 3, end[3] * 3];
  const r = slerp(IDENTITY, bigEnd, 0.5);
  const expected = fromAxisAngle([0, 0, 1], Math.PI / 4);
  assert.ok(quatApproxEqual(r, expected, 1e-7));
});

test('multiply is associative for three quaternions', () => {
  const r1 = ROT_X_90;
  const r2 = ROT_Y_90;
  const r3 = fromAxisAngle([0, 0, 1], Math.PI / 3);
  const left = multiply(multiply(r1, r2), r3);
  const right = multiply(r1, multiply(r2, r3));
  assert.ok(quatApproxEqual(left, right, 1e-7));
});

test('conjugate of identity is identity', () => {
  assert.ok(quatApproxEqual(conjugate(IDENTITY), IDENTITY));
});

test('conjugate of a rotation negates the vector part', () => {
  const q = [1, 2, 3, 4];
  assert.deepEqual(conjugate(q), [-1, -2, -3, 4]);
});

test('compose of empty list returns identity', () => {
  assert.deepEqual(compose([]), IDENTITY);
});

test('compose of single quaternion normalizes it', () => {
  const q = fromAxisAngle([1, 0, 0], Math.PI / 3);
  const scaled = [q[0] * 2, q[1] * 2, q[2] * 2, q[3] * 2];
  const r = compose([scaled]);
  assert.ok(approxEqual(norm(r), 1));
  assert.ok(quatApproxEqual(r, q, 1e-7));
});

test('compose applies rotations right-to-left', () => {
  // compose([r1, r2]) == r1 * r2  (r2 applied first to a vector)
  const r1 = ROT_X_90;
  const r2 = ROT_Y_90;
  assert.ok(quatApproxEqual(compose([r1, r2]), multiply(r1, r2), 1e-7));
});

test('compose normalizes non-unit inputs so the result stays unit', () => {
  const a = [0, 0, 3, 4];
  const b = [0, 5, 0, 0];
  const r = compose([a, b]);
  assert.ok(approxEqual(norm(r), 1));
});

test('compose matches manual product of three rotations', () => {
  const r1 = ROT_X_90;
  const r2 = ROT_Y_90;
  const r3 = fromAxisAngle([0, 0, 1], Math.PI / 6);
  const expected = multiply(multiply(r1, r2), r3);
  assert.ok(quatApproxEqual(compose([r1, r2, r3]), expected, 1e-7));
});

test('round trip: quaternion to axis-angle and back is stable', () => {
  const original = fromAxisAngle([1, 1, 1], Math.PI / 5);
  const { axis, angle } = toAxisAngle(original);
  const rebuilt = fromAxisAngle(axis, angle);
  assert.ok(quatApproxEqual(rebuilt, original, 1e-7));
});
