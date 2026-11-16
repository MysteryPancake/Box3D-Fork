# Notes for agents: fork changes

This is a fork of erincatto/box3d. Everything that differs from upstream `main` is listed here, with how it works and
where to look. Keep this file and the README's "Changes since upstream main" list in step with any further fork change.
Describe how the code works now, not its history.

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

## World snapshots

A public, in-memory save / restore of a world: `b3World_CreateSnapshot`, `b3World_RestoreSnapshot`,
`b3WorldSnapshot_GetByteCount`, `b3DestroyWorldSnapshot` (`box3d.h`, implemented at the end of `world_snapshot.c`). User
docs are in `docs/recording.md`.

**How it works.** It reuses the recording player's keyframe serializer: `b3SerializeWorldImpl` and
`b3DeserializeIntoShell`. The player already restores into a populated world, so no new teardown code was needed. The
difference is how geometry is stored. The header flag `B3_SNAP_FLAG_HOST_GEOMETRY` makes mesh, height field and
compound shapes write the host's pointer instead of interning their bytes (a recording needs the bytes because it
leaves the process; a snapshot does not, and hashing a large mesh per snapshot would be slow). Hulls are still interned
into a private `b3Recording` registry held by the `b3WorldSnapshot`, because `b3AddHullToDatabase` clones from bytes. On
restore a temporary `b3RecReader` whose slots point at that registry supplies the hull bytes.

**Limits.** The pointers make a snapshot valid only inside the writing process, and only while the host keeps its mesh
data alive. Body, shape and joint `userData` are written as NULL, as in recordings. Restore refuses a locked or recording
world. Every snapshot copies each hull it uses: sharing hulls between snapshots would save memory (not done).

**Test.** `WorldSnapshotRestore` in `test/test_recording.c`: restore into the same world, then the next 60 steps match
the first run exactly. It uses hulls and spheres, not meshes.
