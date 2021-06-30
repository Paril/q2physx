#pragma once

#include <PxPhysicsAPI.h>

using namespace physx;

void PhysX_Init(void);

void PhysX_Shutdown(void);

void PhysX_Simulate();

void PhysX_MakePlayer(edict_t *player);

void PhysX_UpdatePlayer(edict_t *player);

void PhysX_DestroyPlayer(edict_t *player);

PxConvexMesh *PhysX_CookMD2(const char *name, int frame = 0);

void PhysX_MakeBrushModel(edict_t *bmodel);

void PhysX_DestroyBrushModel(edict_t *bmodel);

PxRigidDynamic *PhysX_ConvertToPhysX(edict_t *self, const PxGeometry &geometry, const float &density);

extern PxScene *scene;

extern PxPhysics *physics;

extern PxMaterial *material;

extern PxMaterial *highFrictionMaterial;

inline PxQuat QuatFromEuler(vec3_t euler)
{
	// cache trig values of half radian values
	double cx = cos( euler[2] * M_PI/360 );
	double cy = cos( euler[0] * M_PI/360 );
	double cz = cos( euler[1] * M_PI/360 );
	double sx = sin( euler[2] * M_PI/360 );
	double sy = sin( euler[0] * M_PI/360 );
	double sz = sin( euler[1] * M_PI/360 );

	// cache products
	double cc = cx*cz;
	double cs = cx*sz;
	double sc = sx*cz;
	double ss = sx*sz;

	return PxQuat(
		(cy*sc) - (sy*cs),
		(cy*ss) + (sy*cc),
		(cy*cs) - (sy*sc),
		(cy*cc) + (sy*ss)
	).getNormalized();
}

inline PxTransform TransformFromQuake(edict_t *ent)
{
	return PxTransform(PxVec3(ent->s.origin[0], ent->s.origin[1], ent->s.origin[2]), QuatFromEuler(ent->s.angles));
}

inline void QuakeFromTransform(const PxTransform &xf, vec3_t origin, vec3_t angles)
{
	VectorSet(origin, xf.p.x, xf.p.y, xf.p.z);
	PxVec3 axis;
	float angle;
	xf.q.toRadiansAndUnitAxis(angle, axis);
	auto &q = xf.q;

	// roll (x-axis rotation)
	double sinr_cosp = 2 * (q.w * q.x + q.y * q.z);
	double cosr_cosp = 1 - 2 * (q.x * q.x + q.y * q.y);
	angles[2] = atan2(sinr_cosp, cosr_cosp) * (180 / M_PI);

	// pitch (y-axis rotation)
	double sinp = 2 * (q.w * q.y - q.z * q.x);
	if (fabs(sinp) >= 1)
		angles[0] = copysign(M_PI / 2, sinp) * (180 / M_PI); // use 90 degrees if out of range
	else
		angles[0] = asin(sinp) * (180 / M_PI);

	// yaw (z-axis rotation)
	double siny_cosp = 2 * (q.w * q.z + q.x * q.y);
	double cosy_cosp = 1 - 2 * (q.y * q.y + q.z * q.z);
	angles[1] = atan2(siny_cosp, cosy_cosp) * (180 / M_PI);
}