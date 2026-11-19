# Notes for agents: fork changes

This is a fork of erincatto/box3d. Everything that differs from upstream `main` is listed here, with how it works and
where to look. Keep this file and the README's "Changes since upstream main" list in step with any further fork change.
Describe how the code works now, not its history.

## Angular spring error: swing and twist (`src/math_internal.h`)

**Problem.** The spring error was `b3DeltaQuatToRotation( relQ, targetRotation )`, about `2 * vec( q )` of the combined
rotation. With tilt and twist together the tilt part is scaled by `cos( twist / 2 )`, so a body tilted under gravity and
twisted 180 degrees about its z-axis gets no tilt torque. It settled around 75 to 78 degrees tilt (and could orbit
with low damping).

**Fix.** `b3AngularSpringError( relQ, target )` returns the error in frame A as two parts:
- Tilt: the shortest rotation from the target's z-axis to the current z-axis, built from `1 + dot( a, b )` and
  `cross( a, b )` (half angle form, no `acos`). Exactly 180 degrees apart there is no unique axis, so an axis perpendicular
  to `a` is used.
- Twist: what remains of `relQ * conj( target )` after removing the swing is a rotation about `a`; its angle is
  taken the short way round.

Each part has magnitude `2 * sin( angle / 2 )`, the same as the old error for a pure tilt or a pure twist, so existing
spring behaviour without combined tilt and twist is unchanged. The spherical joint's spring (`spherical_joint.c`, rotated
to the world with `quatA`), the motor joint's angular spring (`motor_joint.c`) and the weld joint's soft angular constraint
(`weld_joint.c`, and its wide kernel `b3WeldAngularErrorW`) all use it. The rest of each spring (softness, `rotationMass`,
impulse) is untouched.

**Not changed.** `b3DeltaQuatToRotation` stays for the prismatic joint's rotation lock, a hard lock and not a spring.
Separate tilt and twist stiffness is not exposed.

**Tests.** `TestAngularSpringTwist` in `test/test_joint.c`: a capsule tilted 70 degrees must be within 5 degrees of upright
after 12 s on a sprung spherical joint with twist 90 and 180 free, and on a sprung motor and a weld joint (each over a
spherical joint that pins the twist) at 90 and 170 degrees. The Houdini repo's `tests/core_joints.cpp` ("spherical spring") covers the spherical case through the wrapper.

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

## Wide weld joint solver (`src/weld_joint_wide.inl`)

**Problem.** Joints were solved one at a time in scalar code. A scene with tens of thousands of weld joints (Glue,
Hard) spent most of its step in joint warm start, solve and relax.

**Kernel.** `weld_joint_wide.inl` is included by `wide4.c` and `wide8.c` after `contact_solver_wide.inl`, so it reuses
that file's wide types and `b3GatherBodies` / `b3ScatterBodies`, and the W4 / W8 choice follows `world->simdWidth` like
the contact solver (`B3_WIDE( name )` gives `nameW4` / `nameW8`). One joint per lane. `b3WeldLanes` holds one group of
`B3_SIMD_WIDTH` joints as plain per-lane float arrays. `b3SolveWeldLanes` / `b3WarmStartWeldLanes` are
`b3SolveWeldJoint` / `b3WarmStartWeldJoint` with the branches turned into lane masks: the quaternion polarity flip,
`fixedRotation`, the spring terms (only joints with a spring use them when `useBias` is false), the swing / twist
angular error (`b3WeldAngularErrorW`, identity target) and the singular check in the point constraint
(`b3SolvePointConstraintW`). Lanes past the end of a group have body index 0 (a static body), so their results are
discarded. Not built when `B3_SIMD_NONE`; `b3GetWeldJointBatchWidth` and `b3GetWeldGroupByteCount` then return 0 and the
scalar path runs.

**Two ways to feed it.** Both need the joints of a group to share no dynamic body, which holds within one graph color.
The overflow color stays scalar.
- Batches (`b3SolveWeldJoints`, `b3WarmStartWeldJoints`): `b3WarmStartJointsTask` / `b3SolveJointsTask` in `solver.c`
  collect the weld joints of a block, gather their lanes from `b3JointSim` on every call and write the impulses back.
  Used for colors that mix joint types. A leftover of one joint goes scalar.
- Groups (`b3PrepareWeldGroup`, `b3WarmStartWeldGroups`, `b3SolveWeldGroups`, `b3StoreWeldGroups`): used for a color whose
  joints are all welds (checked in the solver setup). Joint `i` of the color is lane `i % width` of group `i / width`.
  Each group is prepared straight into the color's `weldConstraints` array in the prepare stage
  (`b3PrepareJointsTask` calls `b3PrepareWeldGroup_Wide`), the substeps read only that array, and the new stage `b3_stageStoreWeldImpulses` writes the impulses back to the joint
  sims at the end (the next step, `b3GetJointReaction` and the feedback read them there). Memory is one stack
  allocation (`weld constraints`) holding the groups and the store stage's blocks. The color's graph blocks become
  `b3_graphWeldBlock` blocks counted in groups. 
  Force and torque thresholds are checked from the group solve (`b3CheckJointThresholds`, shared with the scalar path).

**Why groups.** The wide maths alone gave about 3x on solve, but warm start was bound by memory traffic through the large
`b3JointSim` structs. Groups cut the per-joint data from the union-sized sim to one contiguous array.

**Wide prepare.** `b3GatherWeldPrepareLane` (`joint.c`, plain C) reads one joint and its two bodies into a
`b3WeldPrepareLane`. `b3PrepareWeldGroup` (the `.inl`) writes those into the group's lane arrays and does the maths wide:
the frame quaternions and anchors, the angular mass (inverse of the summed inverse inertia, zero and fixed rotation when
singular) and the softness (`b3MakeSoftW`, a joint without its own spring uses the constraint softness). The joint sims
are only read: their prepared fields (`frameA`, `angularMass`, `indexA`, ...) are not updated, and nothing reads them for
a weld in a pure weld color. A prepare block prepares the groups whose first joint it holds, so a group that straddles
two blocks is prepared once, by the earlier block. Lanes past the end of a color use body index 0.

**Not done.** Other joint types (spherical, revolute, prismatic, motor, distance, wheel, parallel, pogo, mover) are
scalar. Mixed colors use batches, which still gather, and their prepare is scalar.

**Tests.** None in this repo. The Houdini repo's `tests/core_joints.cpp`, `core_glue.cpp` and `core_sidefx.cpp` run weld
joints; compare their results with the scalar path.

## Joint pair set for pair creation (`src/broad_phase.c`, `src/body.c`, `src/joint.c`)

**Problem.** `b3ShouldBodiesCollide` walked a body's joint list for every candidate shape pair. A pair rejected by a
joint with `collideConnected` false never becomes a contact, so it is not in `pairSet` and is tested again whenever
either shape moves. In a scene where every body has about 14 glue joints this was a large part of the pairs time.

**Fix.** `b3BroadPhase::jointPairSet` holds the key (`b3ShapePairKey( bodyIdA, bodyIdB, 0 )`) of every body pair joined
by at least one joint with `collideConnected` false. `b3ShouldBodiesCollide` returns true at once if either body has no
joints, otherwise it is one `b3ContainsKey`. It is kept in sync in `joint.c`: `b3AddJointPair` when such a joint is
created or has the flag cleared, `b3RemoveJointPairIfUnused` when it is destroyed or the flag is set. The removal walks
the shorter joint list for another joint joining the same two bodies, so two joints between one pair keep the key. The set
is saved and restored with the world snapshot next to `pairSet`. Pair creation only reads it, and it is only written
outside the step.

## Velocities in move events (`include/box3d/types.h`, `src/solver.c`, `src/solver_set.c`)

`b3BodyMoveEvent` has `linearVelocity` and `angularVelocity`, the body's velocities at the end of the step, so an
application that copies every moved body's state needs no per-body getter calls. `b3FinalizeBodiesTask` fills them from
the body state with the transform; the continuous collision pass (`solver.c`) updates the linear velocity of a fast body
it corrects (the angular velocity is not changed there), and `b3TrySleepIsland` (`solver_set.c`) sets both to zero for a
body that fell asleep, as `b3Body_GetLinearVelocity` returns for a sleeping body. Nothing else changes.

## Contact impulse events (`include/box3d/types.h`, `src/solver.c`, `src/contact_solver*.{c,inl}`, `src/physics_world.c`)

`b3World_EnableContactImpulseEvents( worldId, true )` makes `b3World_GetContactEvents` fill `impulseEvents` /
`impulseCount`: one `b3ContactImpulseEvent` per manifold point with `totalNormalImpulse > 0` of every contact in the
awake set, valid until the next step. An event has both shape and body ids, both bodies' user data, the point (A's
center plus `anchorA`, world space), the normal (A to B), the impulse and `awakeA` / `awakeB` (the body was in the awake
set when the event was made, which is before islands fall asleep, so a body that sleeps this step is still awake
here; a static body never is).

It follows the hit event path. The two store-impulse functions (`b3StoreImpulses` for overflow contacts, `b3StoreImpulsesW`
for the wide ones) set a bit in a per-worker `impulseEventBitSet` for each contact that carries impulse (only while the
world flag is on), `b3SolverStep` clears the bit sets with the hit ones, ORs the workers' sets together after the hit
events and walks the bits to build the events. Off, the only cost is the flag test in the store loops.

## Graph color count (`include/box3d/constants.h`)

`B3_GRAPH_COLOR_COUNT` is 32 (was 24), so dense joint networks overflow less. The overflow color is still the last one,
the dynamic colors are `B3_GRAPH_COLOR_COUNT - 4`, and the color name tables in `constraint_graph.c`, `joint.c` and
`physics_world.c`, `b3Counters::colorCounts` and the overflow color test scene follow the constant. More colors means
more solver sync points per substep.
