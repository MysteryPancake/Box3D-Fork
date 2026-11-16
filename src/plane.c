// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

// An infinite static half-space. The solid side is dot( normal, x ) <= offset, so the normal
// points into free space. The shape can only live on a static body.
//
// The broad-phase needs a finite box, so the half-space is clipped to B3_PLANE_EXTENT around the
// plane origin. Nothing farther than that out along the plane, or deeper than that below it,
// collides.

#include "manifold.h"
#include "math_internal.h"
#include "shape.h"

#include "box3d/base.h"
#include "box3d/collision.h"
#include "box3d/constants.h"

#include <float.h>
#include <math.h>

b3AABB b3ComputePlaneAABB( const b3Plane* shape, b3Transform transform )
{
	b3Vec3 normal = b3RotateVector( transform.q, shape->normal );
	b3Vec3 foot = b3TransformPoint( transform, b3MulSV( shape->offset, shape->normal ) );

	// An oriented box of half extent H along the plane and H / 2 along the normal, centered below
	// the plane. The tangent frame is never built: a unit tangent u has |u_i| + |v_i| at most
	// sqrt( 2 * ( 1 - n_i^2 ) ) for the two tangents u and v.
	float h = B3_PLANE_EXTENT;
	b3Vec3 center = b3MulSub( foot, 0.5f * h, normal );

	b3Vec3 n2 = { normal.x * normal.x, normal.y * normal.y, normal.z * normal.z };
	b3Vec3 extent;
	extent.x = h * sqrtf( 2.0f * b3MaxFloat( 0.0f, 1.0f - n2.x ) ) + 0.5f * h * fabsf( normal.x );
	extent.y = h * sqrtf( 2.0f * b3MaxFloat( 0.0f, 1.0f - n2.y ) ) + 0.5f * h * fabsf( normal.y );
	extent.z = h * sqrtf( 2.0f * b3MaxFloat( 0.0f, 1.0f - n2.z ) ) + 0.5f * h * fabsf( normal.z );

	return ( b3AABB ){ b3Sub( center, extent ), b3Add( center, extent ) };
}

// The smallest signed distance from the proxy surface to the plane. Negative when overlapped.
static float b3PlaneProxySeparation( const b3Plane* plane, const b3ShapeProxy* proxy, b3Transform proxyToPlane, int* index )
{
	float minSeparation = FLT_MAX;
	int minIndex = 0;
	for ( int i = 0; i < proxy->count; ++i )
	{
		b3Vec3 p = b3TransformPoint( proxyToPlane, proxy->points[i] );
		float separation = b3Dot( plane->normal, p ) - plane->offset;
		if ( separation < minSeparation )
		{
			minSeparation = separation;
			minIndex = i;
		}
	}

	if ( index != NULL )
	{
		*index = minIndex;
	}

	return minSeparation - proxy->radius;
}

bool b3OverlapPlane( const b3Plane* shape, b3Transform shapeTransform, const b3ShapeProxy* proxy )
{
	// The proxy is in the frame of the plane's body, so only the shape transform is needed
	B3_UNUSED( shapeTransform );
	return b3PlaneProxySeparation( shape, proxy, b3Transform_identity, NULL ) < B3_OVERLAP_SLOP;
}

b3CastOutput b3RayCastPlane( const b3Plane* shape, const b3RayCastInput* input )
{
	B3_ASSERT( b3IsValidRay( input ) );
	b3CastOutput output = { 0 };

	float originSeparation = b3Dot( shape->normal, input->origin ) - shape->offset;
	if ( originSeparation <= 0.0f )
	{
		// initial overlap
		output.point = input->origin;
		output.hit = true;
		return output;
	}

	float speed = b3Dot( shape->normal, input->translation );
	if ( speed >= 0.0f )
	{
		// moving away or parallel
		return output;
	}

	float fraction = originSeparation / -speed;
	if ( fraction > input->maxFraction )
	{
		return output;
	}

	output.fraction = fraction;
	output.point = b3MulAdd( input->origin, fraction, input->translation );
	output.normal = shape->normal;
	output.hit = true;
	return output;
}

b3CastOutput b3ShapeCastPlane( const b3Plane* shape, const b3ShapeCastInput* input )
{
	b3CastOutput output = { 0 };

	int index;
	float separation = b3PlaneProxySeparation( shape, &input->proxy, b3Transform_identity, &index );
	if ( separation < 0.0f )
	{
		// Initial overlap is treated as a miss
		return output;
	}

	float speed = b3Dot( shape->normal, input->translation );
	if ( speed >= 0.0f )
	{
		return output;
	}

	float fraction = separation / -speed;
	if ( fraction > input->maxFraction )
	{
		return output;
	}

	// The contact point on the plane under the lowest proxy point at the time of impact
	b3Vec3 p = b3MulAdd( input->proxy.points[index], fraction, input->translation );
	float height = b3Dot( shape->normal, p ) - shape->offset;

	output.fraction = fraction;
	output.point = b3MulSub( p, height, shape->normal );
	output.normal = shape->normal;
	output.hit = true;
	return output;
}

int b3CollideMoverAndPlane( b3PlaneResult* result, const b3Plane* shape, const b3Capsule* mover )
{
	float s1 = b3Dot( shape->normal, mover->center1 ) - shape->offset;
	float s2 = b3Dot( shape->normal, mover->center2 ) - shape->offset;
	b3Vec3 closest = s1 < s2 ? mover->center1 : mover->center2;
	float separation = b3MinFloat( s1, s2 ) - mover->radius;
	if ( separation > 0.0f )
	{
		return 0;
	}

	// The plane faces the mover, offset by the mover radius so the separation is of the surface
	b3Plane plane = { shape->normal, shape->offset + mover->radius };
	b3Vec3 point = b3MulSub( closest, b3MinFloat( s1, s2 ), shape->normal );
	*result = ( b3PlaneResult ){ plane, point, 0, 0, 0 };
	return 1;
}

// Conservative advancement of a sweeping proxy toward the plane. The plane is static in its own
// frame. Time moves forward by the gap divided by an upper bound on the speed of the closest
// approach, so a step never passes the first time the separation reaches the target.
b3TOIOutput b3PlaneTimeOfImpact( const b3Plane* shape, const b3Sweep* sweepA, const b3ShapeProxy* proxyB, const b3Sweep* sweepB,
								 float maxFraction )
{
	b3TOIOutput output = { 0 };
	output.state = b3_toiStateSeparated;
	output.fraction = maxFraction;

	b3Transform xfA = {
		.p = b3Sub( sweepA->c1, b3RotateVector( sweepA->q1, sweepA->localCenter ) ),
		.q = sweepA->q1,
	};

	const float target = B3_LINEAR_SLOP;
	const float tolerance = 0.25f * B3_LINEAR_SLOP;

	// Bound the speed along the plane normal, per unit of fraction.
	b3Vec3 worldNormal = b3RotateVector( xfA.q, shape->normal );
	float linearSpeed = fabsf( b3Dot( worldNormal, b3Sub( sweepB->c2, sweepB->c1 ) ) );

	b3Quat dq = b3InvMulQuat( sweepB->q1, sweepB->q2 );
	float angle = 2.0f * atan2f( b3Length( dq.v ), fabsf( dq.s ) );

	float maxRadius = 0.0f;
	for ( int i = 0; i < proxyB->count; ++i )
	{
		maxRadius = b3MaxFloat( maxRadius, b3Distance( proxyB->points[i], sweepB->localCenter ) );
	}

	// Allow for the rotation being interpolated non-uniformly
	float speedBound = linearSpeed + 1.25f * angle * maxRadius;
	if ( speedBound < FLT_EPSILON )
	{
		return output;
	}

	float t = 0.0f;
	for ( int iteration = 0; iteration < 32; ++iteration )
	{
		b3Transform xfB = b3GetSweepTransform( sweepB, t );
		b3Transform relative = b3InvMulTransforms( xfA, xfB );

		int index;
		float separation = b3PlaneProxySeparation( shape, proxyB, relative, &index );
		output.distanceIterations = iteration + 1;

		if ( separation < target + tolerance )
		{
			if ( t == 0.0f )
			{
				output.state = b3_toiStateOverlapped;
				output.fraction = 0.0f;
				return output;
			}

			// Hit
			b3Vec3 p = b3TransformPoint( relative, proxyB->points[index] );
			float height = b3Dot( shape->normal, p ) - shape->offset;
			output.state = b3_toiStateHit;
			output.fraction = t;
			output.point = b3TransformPoint( xfA, b3MulSub( p, height, shape->normal ) );
			output.normal = worldNormal;
			output.distance = separation;
			return output;
		}

		t += ( separation - target ) / speedBound;
		if ( t >= maxFraction )
		{
			return output;
		}
	}

	// Out of iterations: the last time is still before contact
	if ( t > 0.0f )
	{
		b3Transform xfB = b3GetSweepTransform( sweepB, t );
		b3Transform relative = b3InvMulTransforms( xfA, xfB );
		int index;
		float separation = b3PlaneProxySeparation( shape, proxyB, relative, &index );
		b3Vec3 p = b3TransformPoint( relative, proxyB->points[index] );
		float height = b3Dot( shape->normal, p ) - shape->offset;
		output.state = b3_toiStateHit;
		output.fraction = t;
		output.point = b3TransformPoint( xfA, b3MulSub( p, height, shape->normal ) );
		output.normal = worldNormal;
		output.distance = separation;
	}

	return output;
}

// The contact point is the midpoint between the plane and the surface of the shape, as for the
// other convex pairs. The normal points from the plane to the shape.
static void b3SetPlanePoint( b3LocalManifoldPoint* point, const b3Plane* plane, b3Vec3 center, float radius, b3FeaturePair pair )
{
	float centerSeparation = b3Dot( plane->normal, center ) - plane->offset;
	float separation = centerSeparation - radius;
	point->point = b3MulSub( center, radius + 0.5f * separation, plane->normal );
	point->separation = separation;
	point->pair = pair;
}

void b3CollidePlaneAndSphere( b3LocalManifold* manifold, int capacity, const b3Plane* planeA, const b3Sphere* sphereB,
							  b3Transform transformBtoA )
{
	manifold->pointCount = 0;

	if ( capacity < 1 )
	{
		return;
	}

	b3Vec3 center = b3TransformPoint( transformBtoA, sphereB->center );
	float separation = b3Dot( planeA->normal, center ) - planeA->offset - sphereB->radius;
	if ( separation > B3_SPECULATIVE_DISTANCE )
	{
		return;
	}

	manifold->normal = planeA->normal;
	manifold->pointCount = 1;
	b3SetPlanePoint( manifold->points + 0, planeA, center, sphereB->radius, b3FeaturePair_single );
}

void b3CollidePlaneAndCapsule( b3LocalManifold* manifold, int capacity, const b3Plane* planeA, const b3Capsule* capsuleB,
							   b3Transform transformBtoA )
{
	manifold->pointCount = 0;

	if ( capacity < 1 )
	{
		return;
	}

	b3Vec3 centers[2] = {
		b3TransformPoint( transformBtoA, capsuleB->center1 ),
		b3TransformPoint( transformBtoA, capsuleB->center2 ),
	};

	float radius = capsuleB->radius;
	int count = 0;
	b3LocalManifoldPoint points[2];
	for ( int i = 0; i < 2; ++i )
	{
		float separation = b3Dot( planeA->normal, centers[i] ) - planeA->offset - radius;
		if ( separation > B3_SPECULATIVE_DISTANCE )
		{
			continue;
		}

		b3FeaturePair pair = b3MakeFeaturePair( b3_featureShapeB, i, b3_featureShapeB, i );
		b3SetPlanePoint( points + count, planeA, centers[i], radius, pair );
		count += 1;
	}

	if ( count == 0 )
	{
		return;
	}

	manifold->normal = planeA->normal;
	for ( int i = 0; i < count && i < capacity; ++i )
	{
		manifold->points[i] = points[i];
	}
	manifold->pointCount = b3MinInt( count, capacity );
}

void b3CollidePlaneAndHull( b3LocalManifold* manifold, int capacity, const b3Plane* planeA, const b3HullData* hullB,
							b3Transform transformBtoA )
{
	manifold->pointCount = 0;

	if ( capacity < 1 )
	{
		return;
	}

	const b3Vec3* vertices = b3GetHullPoints( hullB );
	int vertexCount = hullB->vertexCount;
	b3Vec3 normal = planeA->normal;

	// Hull vertices that are near or below the plane
	b3Vec3 points[B3_MAX_HULL_VERTICES];
	float separations[B3_MAX_HULL_VERTICES];
	int indices[B3_MAX_HULL_VERTICES];
	int count = 0;
	for ( int i = 0; i < vertexCount && count < B3_MAX_HULL_VERTICES; ++i )
	{
		b3Vec3 p = b3TransformPoint( transformBtoA, vertices[i] );
		float separation = b3Dot( normal, p ) - planeA->offset;
		if ( separation <= B3_SPECULATIVE_DISTANCE )
		{
			points[count] = p;
			separations[count] = separation;
			indices[count] = i;
			count += 1;
		}
	}

	if ( count == 0 )
	{
		return;
	}

	// Pick at most four points that span the most area. The choice depends only on the geometry,
	// and a vertex index is the feature, so points persist from step to step.
	int picks[B3_MAX_HULL_VERTICES];
	int pickCount = 0;

	if ( count <= capacity )
	{
		for ( int i = 0; i < count; ++i )
		{
			picks[pickCount++] = i;
		}
	}
	else
	{
		// Deepest point
		int first = 0;
		for ( int i = 1; i < count; ++i )
		{
			if ( separations[i] < separations[first] )
			{
				first = i;
			}
		}

		// Farthest from the first
		int second = first == 0 ? 1 : 0;
		float maxDistanceSqr = -1.0f;
		for ( int i = 0; i < count; ++i )
		{
			if ( i == first )
			{
				continue;
			}

			float distanceSqr = b3DistanceSquared( points[i], points[first] );
			if ( distanceSqr > maxDistanceSqr )
			{
				maxDistanceSqr = distanceSqr;
				second = i;
			}
		}

		picks[pickCount++] = first;
		picks[pickCount++] = second;

		// Largest triangle on each side of the line through the first two
		b3Vec3 edge = b3Sub( points[second], points[first] );
		int third = -1;
		int fourth = -1;
		float maxArea = 0.0f;
		float minArea = 0.0f;
		for ( int i = 0; i < count; ++i )
		{
			if ( i == first || i == second )
			{
				continue;
			}

			float area = b3Dot( normal, b3Cross( edge, b3Sub( points[i], points[first] ) ) );
			if ( area > maxArea )
			{
				maxArea = area;
				third = i;
			}
			else if ( area < minArea )
			{
				minArea = area;
				fourth = i;
			}
		}

		if ( third >= 0 )
		{
			picks[pickCount++] = third;
		}

		if ( fourth >= 0 )
		{
			picks[pickCount++] = fourth;
		}

		// Everything is collinear: the two end points of the edge are enough
		pickCount = b3MinInt( pickCount, capacity );
	}

	manifold->normal = normal;
	manifold->pointCount = pickCount;
	for ( int i = 0; i < pickCount; ++i )
	{
		int k = picks[i];
		b3LocalManifoldPoint* point = manifold->points + i;

		// The point is the midpoint of the vertex and its projection onto the plane
		point->point = b3MulSub( points[k], 0.5f * separations[k], normal );
		point->separation = separations[k];
		point->pair = b3MakeFeaturePair( b3_featureShapeB, indices[k], b3_featureShapeA, 0 );
	}
}
