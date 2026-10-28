# ts-quaternion-rotation-interpolation

Performs spherical linear interpolation between unit quaternions and composes rotation sequences for 3D orientations.

## Usage

```js
import {
  slerp,
  compose,
  fromAxisAngle,
  toAxisAngle,
  normalize,
  conjugate,
  multiply,
} from 'ts-quaternion-rotation-interpolation';

const a = fromAxisAngle([1, 0, 0], Math.PI / 2);   // 90° about X
const b = fromAxisAngle([0, 1, 0], Math.PI / 2);   // 90° about Y

const halfway = slerp(a, b, 0.5);                   // halfway orientation
const combined = compose([a, b]);                   // a * b: apply b first, then a

const { axis, angle } = toAxisAngle(halfway);      // recover axis + angle
```

A quaternion is a plain `[x, y, z, w]` array where `(x, y, z)` is the vector part and `w` is the scalar. The library assumes Hamilton-product multiplication, so `compose([r1, r2])` equals `r1 * r2` and applies `r2` first to a vector — the same right-to-left order as column-major rotation matrices.

## Why this exists

Interpolating 3D orientations with linear interpolation (lerp) produces rotation sequences that speed up and slow down as the quaternions drift away from and toward each other on the unit sphere. Spherical linear interpolation (slerp) travels along the great arc at constant angular velocity, which is what you want for camera paths, skeletal animation, and any case where the rate of rotation should be visually uniform.

The trade-off: slerp is more expensive than nlerp (normalized lerp) because it calls `acos` and `sin` per evaluation. For very small angular differences the two are visually indistinguishable, so this library falls back to nlerp below about 0.00005 degrees of arc to avoid dividing by `sin(theta) ≈ 0`. That threshold is far below any perceptual or numerical concern and keeps the implementation stable where naive slerp blows up.

## The awkward edge

Two quaternions `q` and `-q` represent the same orientation, so the dot product between two inputs can be negative even when the orientations are close. `slerp` checks the dot product and negates the second quaternion when needed, so it always takes the shorter arc. Without this, an animation can spin the long way around (nearly 720°) between two visually identical poses. This is the single most common bug in hand-written quaternion interpolation, and it is silent — the result is still a unit quaternion, just the wrong one.

## Exports

| Name | Signature | Notes |
|---|---|---|
| `slerp(a, b, t)` | `(Quaternion, Quaternion, number) => Quaternion` | `t` clamped to `[0, 1]`; shorter-arc aware; returns a unit quaternion. |
| `compose(quats)` | `(Quaternion[]) => Quaternion` | Right-to-left composition; empty list returns identity. |
| `fromAxisAngle(axis, angle)` | `([number, number, number], number) => Quaternion` | Axis normalized; zero-length axis returns identity. |
| `toAxisAngle(q)` | `(Quaternion) => { axis, angle }` | Identity returns `{ axis: [1,0,0], angle: 0 }` because the axis is otherwise undefined. |
| `normalize(q)` | `(Quaternion) => Quaternion` | Zero quaternion returns identity to avoid NaN propagation. |
| `conjugate(q)` | `(Quaternion) => Quaternion` | Negated vector part, same scalar. |
| `multiply(a, b)` | `(Quaternion, Quaternion) => Quaternion` | Hamilton product. |
