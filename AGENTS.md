# Notes for agents: fork changes

This is a fork of erincatto/box3d. Everything that differs from upstream `main` is listed here, with how it works and
where to look. Keep this file and the README's "Changes since upstream main" list in step with any further fork change.
Describe how the code works now, not its history.

## Spherical joint spring: swing and twist (`src/spherical_joint.c`)

**Problem.** The spring error was `b3DeltaQuatToRotation( relQ, targetRotation )`, about `2 * vec( q )` of the combined
rotation. With tilt and twist together the tilt part is scaled by `cos( twist / 2 )`, so a body tilted under gravity and
twisted 180 degrees about its z-axis gets no tilt torque. It settled around 75 to 78 degrees tilt (and could orbit
with low damping).

**Fix.** `b3ComputeSphericalSpringError( relQ, target )` returns the error in frame A as two parts:
- Tilt: the shortest rotation from the target's z-axis to the current z-axis, built from `1 + dot( a, b )` and
  `cross( a, b )` (half angle form, no `acos`). Exactly 180 degrees apart there is no unique axis, so an axis perpendicular
  to `a` is used.
- Twist: what remains of `relQ * conj( target )` after removing the swing is a rotation about `a`; its angle is
  taken the short way round.

Each part has magnitude `2 * sin( angle / 2 )`, the same as the old error for a pure tilt or a pure twist, so existing
spring behaviour without combined tilt and twist is unchanged. `b3SolveSphericalJoint` rotates the sum to the world
with `quatA` as before; the rest of the spring (softness, `rotationMass`, impulse) is untouched.

**Not changed.** `b3DeltaQuatToRotation` stays for the motor, weld and prismatic joints, whose angular springs hold
the whole orientation, twist included. Separate tilt and twist stiffness is not exposed.

**Tests.** The Houdini repo's `tests/core_joints.cpp` ("spherical spring"): a capsule on a sprung spherical joint,
tilted 70 degrees, twist 0 and 180, must be within 5 degrees of upright after 8 seconds.

## Thin point sets in `b3CreateHull` (`src/hull.c`)

**Problem.** The builder's tolerance scales with float epsilon and the size of the point set. A slab only a few times
`minRadius` thick (about 1.5e-5 of its size; failures seen up to about 1.5e-4 for a 10 unit slab) makes `MergeCoplanar`
judge the top and bottom faces coplanar. They share their whole boundary, so the first loop in
`b3HullBuilder_ConnectFaces` (`edgePrev->twin->face == twin->face`) never reaches a third face and never ends.

**Detection.** `ConnectFaces` returns `bool`. Before changing anything it counts the shared run of edges; reaching the
smaller face's edge count means the faces share everything, so it sets `b3HullBuilder::degenerate` and returns false.
`MergeConcave` and `MergeCoplanar` treat false as "no merge", `MergeFaces` returns, and `b3HullBuilder_Construct` stops
after the vertex that set the flag and returns false. A legitimate merge is never refused: the shared run is always
shorter than both faces' edge counts.

**Recovery.** `b3CreateHull` wraps `b3CreateHullInternal(..., bool* degenerate)`. If the build failed with `degenerate`
set, it retries once on `b3ThickenPoints`: every point is copied to both sides of the thinnest principal axis (Jacobi
eigen-solve of the covariance) by 2e-4 of the size, far above the failing band. The result is a thin but valid hull; if
the retry fails too, `NULL`. Nothing is thickened when the first build works, so existing hulls are unchanged.

**Not done.** Slabs that fail some other way (for example a needle thin on two axes) are not thickened further. No test
in this repo covers it; the Houdini repo's `tests/physics_shapes.cpp` is the place for one.

## Plane shape (`b3_planeShape`)

An infinite static half-space for ground collision. Used by the Houdini plugin (Ground Plane DOP), which can have many.

**Geometry.** The shape holds a `b3Plane` (`normal`, `offset`, normal must be a unit vector). The solid side is
`dot( normal, x ) <= offset` in the body frame; the normal points into free space and is the contact normal for
anything resting on it. Create with `b3CreatePlaneShape( bodyId, def, &plane )` (asserts `b3IsValidPlane`, returns
`b3_nullShapeId` if invalid or the body is not static). `b3Shape_GetPlane` reads it back. There is no setter: move the
body (`b3Body_SetTransform`) to move the plane. Per-plane friction and restitution come from the shape's material.

**Finite box.** The broad-phase needs a finite AABB, so `b3ComputePlaneAABB` returns a box around the plane's nearest
point to the body origin: `B3_PLANE_EXTENT` (`constants.h`, a quarter of `B3_HUGE`) along the plane, and that deep
below it. The dynamic tree rejects boxes wider than `B3_HUGE`, which is why the extent is a quarter. A tilted plane's
box is large in every axis, so every body keeps a (cheap) contact with it; a flat plane's box is thin. Nothing collides
outside the box.

**Code map** (all in `src/`):
- `plane.c`: everything plane specific: AABB, overlap, ray cast, shape cast, mover planes (`b3CollideMoverAndPlane`),
  the three manifold functions and `b3PlaneTimeOfImpact`.
- `shape.c` / `shape.h`: the union member `plane`, creation, and the switch cases (AABB, swept AABB, ray, shape cast,
  overlap, mover, time of impact, margin). The shape is static-only like compounds and height fields.
- `contact.c`: `b3AddType( b3_planeShape, sphere / capsule / hull )` makes the plane shape A of the pair;
  `b3ComputeConvexManifold` dispatches to the plane functions. Plane contacts are ordinary convex contacts (one manifold).
- `body.c`: `b3Body_SetType` refuses to make a body with a plane non-static.
- `sensor.c`: planes can be sensors (`b3OverlapPlane`).
- `recording_ops.inl`, `recording.[ch]`, `replay.[ch]`, `world_snapshot.c`: op `CreatePlaneShape` (0x47), the `PLANE`
  argument type and the snapshot geometry case.
- `physics_world.c`: the debug draw skips planes (nothing finite to draw).
- `types.h`: `b3_planeShape` is last before `b3_shapeTypeCount`, so values stored in recordings and snapshots stay valid.

**Contacts.** The manifold is in the plane's frame with normal = plane normal. Separation is the signed distance of
the other shape's surface; points are kept within `B3_SPECULATIVE_DISTANCE`, at the midpoint between surface and plane.
- Sphere: one point.
- Capsule: up to two points, one per end sphere, feature ids 0 and 1.
- Hull: every vertex within the speculative distance is a candidate. If there are more than the capacity, four are
  picked: the deepest, the one farthest from it, then the largest triangle on each side of that line. Feature id is the
  vertex index, so impulses warm start from step to step.

**Continuous collision.** `b3PlaneTimeOfImpact` is conservative advancement, not GJK: time advances by
(separation - `B3_LINEAR_SLOP`) divided by an upper bound on the speed along the normal (linear plus 1.25 x angular
x farthest point from the center of mass). It returns a hit at the first time within a quarter slop of the target, and
overlapped at fraction 0. At most 32 iterations; if it runs out it reports a hit at the last (still safe) time.

**Ray and shape casts.** Ray cast: a start inside the solid reports a hit at the origin with fraction 0 and a zero
normal (as spheres do). Shape cast: initial overlap is a miss; `canEncroach` is ignored.

**Not done.** Dynamic or kinematic planes, plane vs mesh / height field / compound children, a plane inside a compound,
a debug draw shape, a setter for the plane. Compound children cannot be planes.

**Tests.** The Houdini repo's `tests/core_plane.cpp` covers resting shapes, a raised and a tilted plane, a trough of two
planes, a 400 m/s sphere (no tunnelling) and refused input. The fork has no plane test of its own.
