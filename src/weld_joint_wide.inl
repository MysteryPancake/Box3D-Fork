// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

// Wide (SIMD) weld joint warm start and solve. One joint per lane.
//
// This file is included after contact_solver_wide.inl by wide4.c and wide8.c. It reuses the wide
// vector, quaternion and symmetric matrix types and the body gather/scatter from that file.
//
// The joints in one batch must not share a dynamic body. The constraint graph colors guarantee that
// for joints in the same color. Overflow joints must use the scalar path.
//
// There are two ways to feed the kernels, both built on b3WeldLanes:
// - Batches (b3WarmStartWeldJoints, b3SolveWeldJoints) gather lanes from b3JointSim on every call. This
//   works for any mix of joint types in a color but touches the large joint sims each time.
// - Groups (b3PrepareWeldGroup, b3WarmStartWeldGroups, b3SolveWeldGroups, b3StoreWeldGroups) are used for
//   colors that contain only weld joints. The lanes are prepared once per step, in the prepare stage, straight
//   into a compact contiguous array. Joint i of the color is lane i % width of group i / width. The solver
//   iterations then read only that array and the impulses are written back to the joint sims at the end.
//
// The math matches b3SolveWeldJoint, with branches turned into lane masks.

#include "bitset.h"
#include "joint.h"

#if !defined( B3_SIMD_NONE )

typedef struct b3Vec3Lanes
{
	float x[B3_SIMD_WIDTH];
	float y[B3_SIMD_WIDTH];
	float z[B3_SIMD_WIDTH];
} b3Vec3Lanes;

typedef struct b3QuatLanes
{
	b3Vec3Lanes v;
	float s[B3_SIMD_WIDTH];
} b3QuatLanes;

typedef struct b3SymMatrix3Lanes
{
	float cxx[B3_SIMD_WIDTH];
	float cxy[B3_SIMD_WIDTH];
	float cxz[B3_SIMD_WIDTH];
	float cyy[B3_SIMD_WIDTH];
	float cyz[B3_SIMD_WIDTH];
	float czz[B3_SIMD_WIDTH];
} b3SymMatrix3Lanes;

static inline void b3SetLaneV( b3Vec3Lanes* lanes, int k, b3Vec3 v )
{
	lanes->x[k] = v.x;
	lanes->y[k] = v.y;
	lanes->z[k] = v.z;
}

static inline void b3SetLaneQ( b3QuatLanes* lanes, int k, b3Quat q )
{
	b3SetLaneV( &lanes->v, k, q.v );
	lanes->s[k] = q.s;
}

// The inertia matrices are symmetric, so only the unique entries are gathered
static inline void b3SetLaneM( b3SymMatrix3Lanes* lanes, int k, b3Matrix3 m )
{
	lanes->cxx[k] = m.cx.x;
	lanes->cxy[k] = m.cx.y;
	lanes->cxz[k] = m.cx.z;
	lanes->cyy[k] = m.cy.y;
	lanes->cyz[k] = m.cy.z;
	lanes->czz[k] = m.cz.z;
}

static inline b3Vec3W b3LoadLanesV( const b3Vec3Lanes* lanes )
{
	return (b3Vec3W){ b3LoadW( lanes->x ), b3LoadW( lanes->y ), b3LoadW( lanes->z ) };
}

static inline b3QuatW b3LoadLanesQ( const b3QuatLanes* lanes )
{
	return (b3QuatW){ b3LoadLanesV( &lanes->v ), b3LoadW( lanes->s ) };
}

static inline b3SymMatrix3W b3LoadLanesM( const b3SymMatrix3Lanes* lanes )
{
	return (b3SymMatrix3W){
		b3LoadW( lanes->cxx ), b3LoadW( lanes->cxy ), b3LoadW( lanes->cxz ),
		b3LoadW( lanes->cyy ), b3LoadW( lanes->cyz ), b3LoadW( lanes->czz ),
	};
}

static inline b3Vec3 b3GetLaneV( const b3Vec3Lanes* lanes, int k )
{
	return B3_LITERAL( b3Vec3 ){ lanes->x[k], lanes->y[k], lanes->z[k] };
}

static inline b3Vec3W b3BlendVW( b3Vec3W a, b3Vec3W b, b3FloatW mask )
{
	return (b3Vec3W){ b3BlendW( a.X, b.X, mask ), b3BlendW( a.Y, b.Y, mask ), b3BlendW( a.Z, b.Z, mask ) };
}

// a * b
static inline b3QuatW b3MulQuatW( b3QuatW a, b3QuatW b )
{
	b3Vec3W t1 = b3CrossW( a.V, b.V );
	b3Vec3W t2 = b3MulAddSVW( t1, a.S, b.V );
	b3Vec3W t3 = b3MulAddSVW( t2, b.S, a.V );
	return (b3QuatW){ t3, b3SubW( b3MulW( a.S, b.S ), b3DotW( a.V, b.V ) ) };
}

// inv(a) * b
static inline b3QuatW b3InvMulQuatW( b3QuatW a, b3QuatW b )
{
	b3Vec3W t1 = b3CrossW( b.V, a.V );
	b3Vec3W t2 = b3MulAddSVW( t1, a.S, b.V );
	b3Vec3W t3 = b3MulSubSVW( t2, b.S, a.V );
	return (b3QuatW){ t3, b3AddW( b3MulW( a.S, b.S ), b3DotW( a.V, b.V ) ) };
}

static inline b3QuatW b3ConjugateW( b3QuatW q )
{
	return (b3QuatW){ b3NegVW( q.V ), q.S };
}

// v + 2 * cross(q.v, cross(q.v, v) + q.s * v)
static inline b3Vec3W b3RotateVectorW( b3QuatW q, b3Vec3W v )
{
	b3Vec3W t1 = b3CrossW( q.V, v );
	b3Vec3W t2 = b3MulAddSVW( t1, q.S, v );
	b3Vec3W t3 = b3CrossW( q.V, t2 );
	return b3MulAddSVW( v, b3SplatW( 2.0f ), t3 );
}

// Wide b3AngularSpringError with an identity target, the weld joint case.
// The target z-axis is (0, 0, 1), so the swing axis is cross(z, b) = (-b.y, b.x, 0).
static inline b3Vec3W b3WeldAngularErrorW( b3QuatW relQ )
{
	b3FloatW zero = b3ZeroW();
	b3FloatW one = b3SplatW( 1.0f );
	b3Vec3W axisZ = { zero, zero, one };

	// b = rotate(relQ, z)
	b3Vec3W b = b3RotateVectorW( relQ, axisZ );

	// Half angle form of the tilt. Singular only when b = -z.
	b3FloatW w = b3AddW( one, b.Z );
	b3Vec3W swingAxis = b3CrossW( axisZ, b );

	b3FloatW lengthSq = b3AddW( b3DotW( swingAxis, swingAxis ), b3MulW( w, w ) );
	b3FloatW invLength = b3DivW( one, b3SqrtW( b3MaxW( lengthSq, b3SplatW( 1000.0f * FLT_MIN ) ) ) );
	b3QuatW swing = { b3MulSVW( invLength, swingAxis ), b3MulW( invLength, w ) };

	// Tilt of 180 degrees: any axis perpendicular to z will do. The scalar code picks +y.
	b3FloatW regular = b3GreaterThanW( w, b3SplatW( 1e-6f ) );
	b3Vec3W flipped = { zero, one, zero };
	swing.V = b3BlendVW( flipped, swing.V, regular );
	swing.S = b3BlendW( zero, swing.S, regular );

	// What is left after the swing is the twist about z.
	b3QuatW t = b3MulQuatW( b3ConjugateW( swing ), relQ );
	b3FloatW sign = b3BlendW( one, b3SplatW( -1.0f ), b3LessThanW( t.S, zero ) );
	b3FloatW twist = b3MulW( b3SplatW( 2.0f ), b3MulW( sign, t.V.Z ) );

	b3FloatW two = b3SplatW( 2.0f );
	return (b3Vec3W){ b3MulW( two, swing.V.X ), b3MulW( two, swing.V.Y ), b3AddW( b3MulW( two, swing.V.Z ), twist ) };
}

// Wide b3SolvePointConstraint
static inline b3Vec3W b3SolvePointConstraintW( b3FloatW mA, b3FloatW mB, b3SymMatrix3W iA, b3SymMatrix3W iB, b3Vec3W rA, b3Vec3W rB,
											   b3Vec3W rhs )
{
	b3Vec3W iAx = { iA.cxx, iA.cxy, iA.cxz };
	b3Vec3W iAy = { iA.cxy, iA.cyy, iA.cyz };
	b3Vec3W iAz = { iA.cxz, iA.cyz, iA.czz };
	b3Vec3W iBx = { iB.cxx, iB.cxy, iB.cxz };
	b3Vec3W iBy = { iB.cxy, iB.cyy, iB.cyz };
	b3Vec3W iBz = { iB.cxz, iB.cyz, iB.czz };

	b3FloatW axy = b3SubW( b3MulW( rA.Z, iA.cyy ), b3MulW( rA.Y, iA.cyz ) );
	b3FloatW axz = b3SubW( b3MulW( rA.Z, iA.cyz ), b3MulW( rA.Y, iA.czz ) );
	b3Vec3W ay = b3SubVW( b3MulSVW( rA.X, iAz ), b3MulSVW( rA.Z, iAx ) );
	b3Vec3W az = b3SubVW( b3MulSVW( rA.Y, iAx ), b3MulSVW( rA.X, iAy ) );

	b3FloatW bxy = b3SubW( b3MulW( rB.Z, iB.cyy ), b3MulW( rB.Y, iB.cyz ) );
	b3FloatW bxz = b3SubW( b3MulW( rB.Z, iB.cyz ), b3MulW( rB.Y, iB.czz ) );
	b3Vec3W by = b3SubVW( b3MulSVW( rB.X, iBz ), b3MulSVW( rB.Z, iBx ) );
	b3Vec3W bz = b3SubVW( b3MulSVW( rB.Y, iBx ), b3MulSVW( rB.X, iBy ) );

	b3FloatW m = b3AddW( mA, mB );

	b3FloatW kxx = b3AddW( b3AddW( b3SubW( b3MulW( axy, rA.Z ), b3MulW( axz, rA.Y ) ), b3SubW( b3MulW( bxy, rB.Z ), b3MulW( bxz, rB.Y ) ) ), m );
	b3FloatW kyy = b3AddW( b3AddW( b3SubW( b3MulW( ay.Z, rA.X ), b3MulW( ay.X, rA.Z ) ), b3SubW( b3MulW( by.Z, rB.X ), b3MulW( by.X, rB.Z ) ) ), m );
	b3FloatW kzz = b3AddW( b3AddW( b3SubW( b3MulW( az.X, rA.Y ), b3MulW( az.Y, rA.X ) ), b3SubW( b3MulW( bz.X, rB.Y ), b3MulW( bz.Y, rB.X ) ) ), m );
	b3FloatW kxy = b3AddW( b3SubW( b3MulW( ay.Y, rA.Z ), b3MulW( ay.Z, rA.Y ) ), b3SubW( b3MulW( by.Y, rB.Z ), b3MulW( by.Z, rB.Y ) ) );
	b3FloatW kxz = b3AddW( b3SubW( b3MulW( az.Y, rA.Z ), b3MulW( az.Z, rA.Y ) ), b3SubW( b3MulW( bz.Y, rB.Z ), b3MulW( bz.Z, rB.Y ) ) );
	b3FloatW kyz = b3AddW( b3SubW( b3MulW( az.Z, rA.X ), b3MulW( az.X, rA.Z ) ), b3SubW( b3MulW( bz.Z, rB.X ), b3MulW( bz.X, rB.Z ) ) );

	b3FloatW cxx = b3SubW( b3MulW( kyy, kzz ), b3MulW( kyz, kyz ) );
	b3FloatW cyy = b3SubW( b3MulW( kxx, kzz ), b3MulW( kxz, kxz ) );
	b3FloatW czz = b3SubW( b3MulW( kxx, kyy ), b3MulW( kxy, kxy ) );
	b3FloatW cxy = b3SubW( b3MulW( kxz, kyz ), b3MulW( kxy, kzz ) );
	b3FloatW cxz = b3SubW( b3MulW( kxy, kyz ), b3MulW( kyy, kxz ) );
	b3FloatW cyz = b3SubW( b3MulW( kxy, kxz ), b3MulW( kxx, kyz ) );

	b3FloatW det = b3AddW( b3AddW( b3MulW( kxx, cxx ), b3MulW( kxy, cxy ) ), b3MulW( kxz, cxz ) );
	b3FloatW valid = b3GreaterThanW( b3AbsW( det ), b3SplatW( 1000.0f * FLT_MIN ) );

	// Lanes with a singular matrix are replaced by zero below, so a bad reciprocal there is harmless
	b3FloatW invDet = b3DivW( b3SplatW( 1.0f ), det );

	b3Vec3W b = {
		b3MulW( invDet, b3AddW( b3MulW( cxx, rhs.X ), b3AddW( b3MulW( cxy, rhs.Y ), b3MulW( cxz, rhs.Z ) ) ) ),
		b3MulW( invDet, b3AddW( b3MulW( cxy, rhs.X ), b3AddW( b3MulW( cyy, rhs.Y ), b3MulW( cyz, rhs.Z ) ) ) ),
		b3MulW( invDet, b3AddW( b3MulW( cxz, rhs.X ), b3AddW( b3MulW( cyz, rhs.Y ), b3MulW( czz, rhs.Z ) ) ) ),
	};

	b3FloatW zero = b3ZeroW();
	return (b3Vec3W){ b3BlendW( zero, b.X, valid ), b3BlendW( zero, b.Y, valid ), b3BlendW( zero, b.Z, valid ) };
}

// Joint data for one batch, gathered from the joint sims
typedef struct b3WeldLanes
{
	// base 1 body indices, 0 is a static body or an unused lane
	_Alignas( B3_WIDE_ALIGNMENT ) int indexA[B3_SIMD_WIDTH];
	int indexB[B3_SIMD_WIDTH];

	float invMassA[B3_SIMD_WIDTH];
	float invMassB[B3_SIMD_WIDTH];
	b3SymMatrix3Lanes invIA;
	b3SymMatrix3Lanes invIB;
	b3SymMatrix3Lanes angularMass;
	b3QuatLanes frameQA;
	b3QuatLanes frameQB;
	b3Vec3Lanes frameA;
	b3Vec3Lanes frameB;
	b3Vec3Lanes deltaCenter;
	b3Vec3Lanes linearImpulse;
	b3Vec3Lanes angularImpulse;

	float fixedRotation[B3_SIMD_WIDTH];
	float linearHertz[B3_SIMD_WIDTH];
	float angularHertz[B3_SIMD_WIDTH];
	float linearBiasRate[B3_SIMD_WIDTH];
	float linearMassScale[B3_SIMD_WIDTH];
	float linearImpulseScale[B3_SIMD_WIDTH];
	float angularBiasRate[B3_SIMD_WIDTH];
	float angularMassScale[B3_SIMD_WIDTH];
	float angularImpulseScale[B3_SIMD_WIDTH];

	// 1 for lanes with a force or torque threshold. Only used by groups.
	float hasThreshold[B3_SIMD_WIDTH];
} b3WeldLanes;

static void b3SetWeldLane( b3WeldLanes* data, int k, const b3JointSim* base )
{
	const b3WeldJoint* joint = &base->weldJoint;
	B3_ASSERT( base->type == b3_weldJoint );

	// B3_NULL_INDEX is -1, so base 1 gives 0 for static bodies
	data->indexA[k] = joint->indexA + 1;
	data->indexB[k] = joint->indexB + 1;

	data->invMassA[k] = base->invMassA;
	data->invMassB[k] = base->invMassB;
	b3SetLaneM( &data->invIA, k, base->invIA );
	b3SetLaneM( &data->invIB, k, base->invIB );
	b3SetLaneM( &data->angularMass, k, joint->angularMass );
	b3SetLaneQ( &data->frameQA, k, joint->frameA.q );
	b3SetLaneQ( &data->frameQB, k, joint->frameB.q );
	b3SetLaneV( &data->frameA, k, joint->frameA.p );
	b3SetLaneV( &data->frameB, k, joint->frameB.p );
	b3SetLaneV( &data->deltaCenter, k, joint->deltaCenter );
	b3SetLaneV( &data->linearImpulse, k, joint->linearImpulse );
	b3SetLaneV( &data->angularImpulse, k, joint->angularImpulse );

	data->fixedRotation[k] = base->fixedRotation ? 1.0f : 0.0f;
	data->linearHertz[k] = joint->linearHertz;
	data->angularHertz[k] = joint->angularHertz;
	data->linearBiasRate[k] = joint->linearSpring.biasRate;
	data->linearMassScale[k] = joint->linearSpring.massScale;
	data->linearImpulseScale[k] = joint->linearSpring.impulseScale;
	data->angularBiasRate[k] = joint->angularSpring.biasRate;
	data->angularMassScale[k] = joint->angularSpring.massScale;
	data->angularImpulseScale[k] = joint->angularSpring.impulseScale;
	data->hasThreshold[k] = ( base->forceThreshold < FLT_MAX || base->torqueThreshold < FLT_MAX ) ? 1.0f : 0.0f;
}

// Unused lanes copy joint 0 so every lane holds finite numbers, but get body index 0 so their results are discarded
static void b3GatherWeldJoints( b3WeldLanes* data, b3JointSim* const* sims, int count )
{
	B3_ASSERT( 0 < count && count <= B3_SIMD_WIDTH );

	for ( int k = 0; k < B3_SIMD_WIDTH; ++k )
	{
		bool used = k < count;
		b3SetWeldLane( data, k, used ? sims[k] : sims[0] );
		if ( used == false )
		{
			data->indexA[k] = 0;
			data->indexB[k] = 0;
		}
	}
}

static void b3CopyWeldImpulsesToSims( const b3WeldLanes* data, b3JointSim* const* sims, int count )
{
	for ( int k = 0; k < count; ++k )
	{
		b3WeldJoint* joint = &sims[k]->weldJoint;
		joint->linearImpulse = b3GetLaneV( &data->linearImpulse, k );
		joint->angularImpulse = b3GetLaneV( &data->angularImpulse, k );
	}
}

static void b3WarmStartWeldLanes( const b3WeldLanes* lanes, b3StepContext* context )
{
	const b3WeldLanes* data = lanes;
	b3BodyState* states = context->states;
	b3BodyStateW bA = b3GatherBodies( states, data->indexA );
	b3BodyStateW bB = b3GatherBodies( states, data->indexB );

	b3FloatW mA = b3LoadW( data->invMassA );
	b3FloatW mB = b3LoadW( data->invMassB );
	b3SymMatrix3W iA = b3LoadLanesM( &data->invIA );
	b3SymMatrix3W iB = b3LoadLanesM( &data->invIB );
	b3Vec3W linearImpulse = b3LoadLanesV( &data->linearImpulse );
	b3Vec3W angularImpulse = b3LoadLanesV( &data->angularImpulse );

	b3Vec3W rA = b3RotateVectorW( bA.dq, b3LoadLanesV( &data->frameA ) );
	b3Vec3W rB = b3RotateVectorW( bB.dq, b3LoadLanesV( &data->frameB ) );

	bA.v = b3MulSubSVW( bA.v, mA, linearImpulse );
	bA.w = b3MulSubMVW( bA.w, iA, b3AddVW( b3CrossW( rA, linearImpulse ), angularImpulse ) );

	bB.v = b3MulAddSVW( bB.v, mB, linearImpulse );
	bB.w = b3MulAddMVW( bB.w, iB, b3AddVW( b3CrossW( rB, linearImpulse ), angularImpulse ) );

	b3ScatterBodies( states, data->indexA, &bA );
	b3ScatterBodies( states, data->indexB, &bB );
}

// Updates the impulses in the lanes
static void b3SolveWeldLanes( b3WeldLanes* lanes, b3StepContext* context, bool useBias )
{
	const b3WeldLanes* data = lanes;
	b3BodyState* states = context->states;
	b3BodyStateW bA = b3GatherBodies( states, data->indexA );
	b3BodyStateW bB = b3GatherBodies( states, data->indexB );

	b3FloatW zero = b3ZeroW();
	b3FloatW one = b3SplatW( 1.0f );
	b3FloatW allMask = b3EqualsW( zero, zero );

	b3FloatW mA = b3LoadW( data->invMassA );
	b3FloatW mB = b3LoadW( data->invMassB );
	b3SymMatrix3W iA = b3LoadLanesM( &data->invIA );
	b3SymMatrix3W iB = b3LoadLanesM( &data->invIB );

	b3Vec3W vA = bA.v, wA = bA.w, vB = bB.v, wB = bB.w;

	b3Vec3W linearImpulse = b3LoadLanesV( &data->linearImpulse );
	b3Vec3W angularImpulse = b3LoadLanesV( &data->angularImpulse );

	// angular constraint
	{
		b3QuatW quatA = b3MulQuatW( bA.dq, b3LoadLanesQ( &data->frameQA ) );
		b3QuatW quatB = b3MulQuatW( bB.dq, b3LoadLanesQ( &data->frameQB ) );

		// Keep the rotation angle in the range [-pi, pi]
		b3FloatW flip = b3LessThanW( b3AddW( b3DotW( quatA.V, quatB.V ), b3MulW( quatA.S, quatB.S ) ), zero );
		quatB.V = b3BlendVW( quatB.V, b3NegVW( quatB.V ), flip );
		quatB.S = b3BlendW( quatB.S, b3NegW( quatB.S ), flip );

		b3FloatW fixedMask = b3GreaterThanW( b3LoadW( data->fixedRotation ), zero );

		// Lanes that use the spring terms. Without bias only joints with a spring do.
		b3FloatW springMask = useBias ? allMask : b3GreaterThanW( b3LoadW( data->angularHertz ), zero );

		b3Vec3W bias = { zero, zero, zero };
		b3FloatW massScale = one;
		b3FloatW impulseScale = zero;

		if ( b3AnyTrueW( springMask ) )
		{
			b3QuatW relQ = b3InvMulQuatW( quatA, quatB );
			b3Vec3W error = b3WeldAngularErrorW( relQ );
			b3Vec3W c = b3RotateVectorW( quatA, error );

			b3FloatW biasRate = b3LoadW( data->angularBiasRate );
			bias = b3BlendVW( bias, b3MulSVW( biasRate, c ), springMask );
			massScale = b3BlendW( one, b3LoadW( data->angularMassScale ), springMask );
			impulseScale = b3BlendW( zero, b3LoadW( data->angularImpulseScale ), springMask );
		}

		b3Vec3W cdot = b3SubVW( wB, wA );
		b3SymMatrix3W angularMass = b3LoadLanesM( &data->angularMass );
		b3Vec3W m = b3MulMVW( angularMass, b3AddVW( cdot, bias ) );
		b3Vec3W impulse = b3MulSubSVW( b3MulSVW( b3NegW( massScale ), m ), impulseScale, angularImpulse );

		// No angular constraint when the rotation is fixed
		impulse = b3BlendVW( impulse, (b3Vec3W){ zero, zero, zero }, fixedMask );

		angularImpulse = b3AddVW( angularImpulse, impulse );

		wA = b3MulSubMVW( wA, iA, impulse );
		wB = b3MulAddMVW( wB, iB, impulse );
	}

	// linear constraint
	{
		b3Vec3W rA = b3RotateVectorW( bA.dq, b3LoadLanesV( &data->frameA ) );
		b3Vec3W rB = b3RotateVectorW( bB.dq, b3LoadLanesV( &data->frameB ) );

		b3Vec3W cdot = b3SubVW( b3AddVW( vB, b3CrossW( wB, rB ) ), b3AddVW( vA, b3CrossW( wA, rA ) ) );

		b3FloatW springMask = useBias ? allMask : b3GreaterThanW( b3LoadW( data->linearHertz ), zero );

		b3Vec3W bias = { zero, zero, zero };
		b3FloatW massScale = one;
		b3FloatW impulseScale = zero;

		if ( b3AnyTrueW( springMask ) )
		{
			b3Vec3W separation = b3AddVW( b3AddVW( b3SubVW( bB.dp, bA.dp ), b3SubVW( rB, rA ) ), b3LoadLanesV( &data->deltaCenter ) );

			b3FloatW biasRate = b3LoadW( data->linearBiasRate );
			bias = b3BlendVW( bias, b3MulSVW( biasRate, separation ), springMask );
			massScale = b3BlendW( one, b3LoadW( data->linearMassScale ), springMask );
			impulseScale = b3BlendW( zero, b3LoadW( data->linearImpulseScale ), springMask );
		}

		b3Vec3W b = b3SolvePointConstraintW( mA, mB, iA, iB, rA, rB, b3AddVW( cdot, bias ) );

		b3Vec3W impulse = b3MulSubSVW( b3MulSVW( b3NegW( massScale ), b ), impulseScale, linearImpulse );
		linearImpulse = b3AddVW( linearImpulse, impulse );

		vA = b3MulSubSVW( vA, mA, impulse );
		wA = b3MulSubMVW( wA, iA, b3CrossW( rA, impulse ) );
		vB = b3MulAddSVW( vB, mB, impulse );
		wB = b3MulAddMVW( wB, iB, b3CrossW( rB, impulse ) );
	}

	bA.v = vA;
	bA.w = wA;
	bB.v = vB;
	bB.w = wB;

	b3ScatterBodies( states, data->indexA, &bA );
	b3ScatterBodies( states, data->indexB, &bB );

	b3StoreVW( lanes->linearImpulse.x, lanes->linearImpulse.y, lanes->linearImpulse.z, linearImpulse );
	b3StoreVW( lanes->angularImpulse.x, lanes->angularImpulse.y, lanes->angularImpulse.z, angularImpulse );
}

void B3_WIDE( b3WarmStartWeldJoints )( b3JointSim* const* sims, int count, b3StepContext* context )
{
	b3WeldLanes data;
	b3GatherWeldJoints( &data, sims, count );
	b3WarmStartWeldLanes( &data, context );
}

void B3_WIDE( b3SolveWeldJoints )( b3JointSim* const* sims, int count, b3StepContext* context, bool useBias )
{
	b3WeldLanes data;
	b3GatherWeldJoints( &data, sims, count );
	b3SolveWeldLanes( &data, context, useBias );
	b3CopyWeldImpulsesToSims( &data, sims, count );
}

// Wide b3MakeSoft. Lanes with a zero hertz give zero softness, like the scalar code.
static inline void b3MakeSoftW( b3FloatW hertz, b3FloatW zeta, b3FloatW h, b3FloatW* biasRate, b3FloatW* massScale,
								b3FloatW* impulseScale )
{
	b3FloatW zero = b3ZeroW();
	b3FloatW one = b3SplatW( 1.0f );
	b3FloatW omega = b3MulW( b3SplatW( 2.0f * B3_PI ), hertz );
	b3FloatW hw = b3MulW( h, omega );
	b3FloatW a1 = b3AddW( b3MulW( b3SplatW( 2.0f ), zeta ), hw );
	b3FloatW a2 = b3MulW( hw, a1 );
	b3FloatW a3 = b3DivW( one, b3AddW( one, a2 ) );

	// a1 is zero only when the hertz is zero, which is masked out below
	b3FloatW valid = b3GreaterThanW( hertz, zero );
	*biasRate = b3BlendW( zero, b3DivW( omega, b3MaxW( a1, b3SplatW( FLT_MIN ) ) ), valid );
	*massScale = b3BlendW( zero, b3MulW( a2, a3 ), valid );
	*impulseScale = b3BlendW( zero, a3, valid );
}

// Prepare group g from the gathered lanes. This is b3PrepareWeldJoint plus the packing, without writing the joint sim.
void B3_WIDE( b3PrepareWeldGroup )( void* groups, int group, const b3WeldPrepareLane* in, int laneCount, float h, float invH,
									bool enableWarmStarting )
{
	b3WeldLanes* g = (b3WeldLanes*)groups + group;

	b3QuatLanes bodyQA, bodyQB, localQA, localQB;
	b3Vec3Lanes localPA, localPB, localCenterA, localCenterB;
	float constraintHertz[B3_SIMD_WIDTH], constraintDamping[B3_SIMD_WIDTH];
	float linearDamping[B3_SIMD_WIDTH], angularDamping[B3_SIMD_WIDTH];

	b3Quat identity = b3Quat_identity;
	for ( int k = 0; k < B3_SIMD_WIDTH; ++k )
	{
		if ( k < laneCount )
		{
			const b3WeldPrepareLane* lane = in + k;
			g->indexA[k] = lane->indexA + 1;
			g->indexB[k] = lane->indexB + 1;
			g->invMassA[k] = lane->invMassA;
			g->invMassB[k] = lane->invMassB;
			b3SetLaneM( &g->invIA, k, lane->invIA );
			b3SetLaneM( &g->invIB, k, lane->invIB );
			b3SetLaneV( &g->deltaCenter, k, lane->deltaCenter );
			b3SetLaneQ( &bodyQA, k, lane->bodyQA );
			b3SetLaneQ( &bodyQB, k, lane->bodyQB );
			b3SetLaneQ( &localQA, k, lane->localFrameA.q );
			b3SetLaneQ( &localQB, k, lane->localFrameB.q );
			b3SetLaneV( &localPA, k, lane->localFrameA.p );
			b3SetLaneV( &localPB, k, lane->localFrameB.p );
			b3SetLaneV( &localCenterA, k, lane->localCenterA );
			b3SetLaneV( &localCenterB, k, lane->localCenterB );
			constraintHertz[k] = lane->constraintHertz;
			constraintDamping[k] = lane->constraintDampingRatio;
			g->linearHertz[k] = lane->linearHertz;
			linearDamping[k] = lane->linearDampingRatio;
			g->angularHertz[k] = lane->angularHertz;
			angularDamping[k] = lane->angularDampingRatio;
			b3SetLaneV( &g->linearImpulse, k, enableWarmStarting ? lane->linearImpulse : b3Vec3_zero );
			b3SetLaneV( &g->angularImpulse, k, enableWarmStarting ? lane->angularImpulse : b3Vec3_zero );
			g->hasThreshold[k] = lane->hasThreshold ? 1.0f : 0.0f;
		}
		else
		{
			// Unused lane: a static body and benign numbers. The result is discarded.
			g->indexA[k] = 0;
			g->indexB[k] = 0;
			g->invMassA[k] = 0.0f;
			g->invMassB[k] = 0.0f;
			b3SetLaneM( &g->invIA, k, b3Mat3_zero );
			b3SetLaneM( &g->invIB, k, b3Mat3_zero );
			b3SetLaneV( &g->deltaCenter, k, b3Vec3_zero );
			b3SetLaneQ( &bodyQA, k, identity );
			b3SetLaneQ( &bodyQB, k, identity );
			b3SetLaneQ( &localQA, k, identity );
			b3SetLaneQ( &localQB, k, identity );
			b3SetLaneV( &localPA, k, b3Vec3_zero );
			b3SetLaneV( &localPB, k, b3Vec3_zero );
			b3SetLaneV( &localCenterA, k, b3Vec3_zero );
			b3SetLaneV( &localCenterB, k, b3Vec3_zero );
			constraintHertz[k] = 0.0f;
			constraintDamping[k] = 0.0f;
			g->linearHertz[k] = 0.0f;
			linearDamping[k] = 0.0f;
			g->angularHertz[k] = 0.0f;
			angularDamping[k] = 0.0f;
			b3SetLaneV( &g->linearImpulse, k, b3Vec3_zero );
			b3SetLaneV( &g->angularImpulse, k, b3Vec3_zero );
			g->hasThreshold[k] = 0.0f;
		}
	}

	b3FloatW zero = b3ZeroW();
	b3FloatW one = b3SplatW( 1.0f );

	// Angular mass is the inverse of the summed inverse inertia, and a singular sum means fixed rotation
	{
		b3SymMatrix3W iA = b3LoadLanesM( &g->invIA );
		b3SymMatrix3W iB = b3LoadLanesM( &g->invIB );
		b3FloatW xx = b3AddW( iA.cxx, iB.cxx );
		b3FloatW xy = b3AddW( iA.cxy, iB.cxy );
		b3FloatW xz = b3AddW( iA.cxz, iB.cxz );
		b3FloatW yy = b3AddW( iA.cyy, iB.cyy );
		b3FloatW yz = b3AddW( iA.cyz, iB.cyz );
		b3FloatW zz = b3AddW( iA.czz, iB.czz );

		b3FloatW cxx = b3SubW( b3MulW( yy, zz ), b3MulW( yz, yz ) );
		b3FloatW cxy = b3SubW( b3MulW( xz, yz ), b3MulW( xy, zz ) );
		b3FloatW cxz = b3SubW( b3MulW( xy, yz ), b3MulW( yy, xz ) );
		b3FloatW cyy = b3SubW( b3MulW( xx, zz ), b3MulW( xz, xz ) );
		b3FloatW cyz = b3SubW( b3MulW( xy, xz ), b3MulW( xx, yz ) );
		b3FloatW czz = b3SubW( b3MulW( xx, yy ), b3MulW( xy, xy ) );

		b3FloatW det = b3AddW( b3AddW( b3MulW( xx, cxx ), b3MulW( xy, cxy ) ), b3MulW( xz, cxz ) );
		b3FloatW threshold = b3SplatW( 1000.0f * FLT_MIN );
		b3FloatW valid = b3GreaterThanW( b3AbsW( det ), threshold );
		b3FloatW invDet = b3DivW( one, det );

		b3StoreW( g->angularMass.cxx, b3BlendW( zero, b3MulW( invDet, cxx ), valid ) );
		b3StoreW( g->angularMass.cxy, b3BlendW( zero, b3MulW( invDet, cxy ), valid ) );
		b3StoreW( g->angularMass.cxz, b3BlendW( zero, b3MulW( invDet, cxz ), valid ) );
		b3StoreW( g->angularMass.cyy, b3BlendW( zero, b3MulW( invDet, cyy ), valid ) );
		b3StoreW( g->angularMass.cyz, b3BlendW( zero, b3MulW( invDet, cyz ), valid ) );
		b3StoreW( g->angularMass.czz, b3BlendW( zero, b3MulW( invDet, czz ), valid ) );

		b3StoreW( g->fixedRotation, b3BlendW( zero, one, b3LessThanW( det, threshold ) ) );
	}

	// Joint frames with world space rotation, relative to the center of mass
	{
		b3QuatW qA = b3LoadLanesQ( &bodyQA );
		b3QuatW qB = b3LoadLanesQ( &bodyQB );

		b3QuatW frameQA = b3MulQuatW( qA, b3LoadLanesQ( &localQA ) );
		b3QuatW frameQB = b3MulQuatW( qB, b3LoadLanesQ( &localQB ) );
		b3StoreVW( g->frameQA.v.x, g->frameQA.v.y, g->frameQA.v.z, frameQA.V );
		b3StoreW( g->frameQA.s, frameQA.S );
		b3StoreVW( g->frameQB.v.x, g->frameQB.v.y, g->frameQB.v.z, frameQB.V );
		b3StoreW( g->frameQB.s, frameQB.S );

		b3Vec3W pA = b3RotateVectorW( qA, b3SubVW( b3LoadLanesV( &localPA ), b3LoadLanesV( &localCenterA ) ) );
		b3Vec3W pB = b3RotateVectorW( qB, b3SubVW( b3LoadLanesV( &localPB ), b3LoadLanesV( &localCenterB ) ) );
		b3StoreVW( g->frameA.x, g->frameA.y, g->frameA.z, pA );
		b3StoreVW( g->frameB.x, g->frameB.y, g->frameB.z, pB );
	}

	// Softness. A joint without its own spring uses the constraint softness.
	{
		b3FloatW hw = b3SplatW( h );
		b3FloatW hertz = b3MinW( b3LoadW( constraintHertz ), b3SplatW( 0.25f * invH ) );
		b3FloatW cBias, cMass, cImpulse;
		b3MakeSoftW( hertz, b3LoadW( constraintDamping ), hw, &cBias, &cMass, &cImpulse );

		b3FloatW linearHertz = b3LoadW( g->linearHertz );
		b3FloatW lBias, lMass, lImpulse;
		b3MakeSoftW( linearHertz, b3LoadW( linearDamping ), hw, &lBias, &lMass, &lImpulse );
		b3FloatW linearOwn = b3GreaterThanW( linearHertz, zero );
		b3StoreW( g->linearBiasRate, b3BlendW( cBias, lBias, linearOwn ) );
		b3StoreW( g->linearMassScale, b3BlendW( cMass, lMass, linearOwn ) );
		b3StoreW( g->linearImpulseScale, b3BlendW( cImpulse, lImpulse, linearOwn ) );

		b3FloatW angularHertz = b3LoadW( g->angularHertz );
		b3FloatW aBias, aMass, aImpulse;
		b3MakeSoftW( angularHertz, b3LoadW( angularDamping ), hw, &aBias, &aMass, &aImpulse );
		b3FloatW angularOwn = b3GreaterThanW( angularHertz, zero );
		b3StoreW( g->angularBiasRate, b3BlendW( cBias, aBias, angularOwn ) );
		b3StoreW( g->angularMassScale, b3BlendW( cMass, aMass, angularOwn ) );
		b3StoreW( g->angularImpulseScale, b3BlendW( cImpulse, aImpulse, angularOwn ) );
	}
}

int B3_WIDE( b3GetWeldGroupByteCount )( void )
{
	return (int)sizeof( b3WeldLanes );
}

void B3_WIDE( b3WarmStartWeldGroups )( void* groups, int startGroup, int groupCount, b3StepContext* context )
{
	const b3WeldLanes* data = (const b3WeldLanes*)groups + startGroup;
	for ( int g = 0; g < groupCount; ++g )
	{
		b3WarmStartWeldLanes( data + g, context );
	}
}

void B3_WIDE( b3SolveWeldGroups )( void* groups, b3JointSim* joints, int jointCount, int startGroup, int groupCount,
								   b3StepContext* context, bool useBias, b3BitSet* jointStateBitSet )
{
	b3WeldLanes* data = (b3WeldLanes*)groups + startGroup;
	for ( int g = 0; g < groupCount; ++g )
	{
		b3WeldLanes* lanes = data + g;
		b3SolveWeldLanes( lanes, context, useBias );

		if ( useBias && b3AnyTrueW( b3GreaterThanW( b3LoadW( lanes->hasThreshold ), b3ZeroW() ) ) )
		{
			int base = ( startGroup + g ) * B3_SIMD_WIDTH;
			for ( int k = 0; k < B3_SIMD_WIDTH && base + k < jointCount; ++k )
			{
				b3JointSim* sim = joints + base + k;
				sim->weldJoint.linearImpulse = b3GetLaneV( &lanes->linearImpulse, k );
				sim->weldJoint.angularImpulse = b3GetLaneV( &lanes->angularImpulse, k );
				b3CheckJointThresholds( sim, context, jointStateBitSet );
			}
		}
	}
}

// Write the impulses back so that the next step, the joint reaction and the feedback see them
void B3_WIDE( b3StoreWeldGroups )( void* groups, b3JointSim* joints, int jointCount, int startGroup, int groupCount )
{
	const b3WeldLanes* data = (const b3WeldLanes*)groups + startGroup;
	for ( int g = 0; g < groupCount; ++g )
	{
		const b3WeldLanes* lanes = data + g;
		int base = ( startGroup + g ) * B3_SIMD_WIDTH;
		for ( int k = 0; k < B3_SIMD_WIDTH && base + k < jointCount; ++k )
		{
			b3WeldJoint* joint = &joints[base + k].weldJoint;
			joint->linearImpulse = b3GetLaneV( &lanes->linearImpulse, k );
			joint->angularImpulse = b3GetLaneV( &lanes->angularImpulse, k );
		}
	}
}

#endif
