#include "g_local.h"
#include "physx.h"

class PxQuake2ErrorCallback : public PxErrorCallback
{
public:
	void reportError(PxErrorCode::Enum code, const char* message, const char* file, int line)
	{
		gi.dprintf("PhysX (%i): %s (%s:%i)\n", code, message, file, line);
	}
};

static PxQuake2ErrorCallback physxErrorCallback;

static PxDefaultAllocator physxAllocatorCallback;

static PxFoundation *foundation;

static PxCooking *cooking;

static PxPvd *pvd;

PxPhysics *physics;

PxScene *scene;

static PxCpuDispatcher *dispatcher;

PxMaterial *material, *highFrictionMaterial;

const char *PVD_HOST = "localhost";

#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <string>
#include <memory>
#include <thread>
#include <variant>

static std::jthread physxThread;

struct PakFile
{
	std::shared_ptr<std::ifstream>	file;
	uint32_t		offset;
	uint32_t		size;
};

static std::unordered_map<std::string, std::variant<std::filesystem::path, PakFile>> files;

static std::string &replace(std::string &s, const std::string &from, const std::string &to)
{
    for(size_t pos = 0; (pos = s.find(from, pos)) != std::string::npos; pos += to.size())
        s.replace(pos, from.size(), to);
    return s;
}

static void AddGameDirectory(const char *dir)
{
	std::filesystem::path baseq2(dir);

	// add all the paks
	for (auto &f : std::filesystem::directory_iterator(baseq2))
	{
		if (!f.is_regular_file() || f.path().extension().compare(".pak") != 0)
			continue;

		auto pak = std::make_shared<std::ifstream>();
		pak->open(f.path(), std::ifstream::in | std::ifstream::binary);
		pak->seekg(4, std::ios_base::beg);
		uint32_t offset, length;
		pak->read((char *) &offset, sizeof(offset));
		pak->read((char *) &length, sizeof(length));

		pak->seekg(offset);

		uint32_t numFiles = length / 64;

		for (uint32_t i = 0; i < numFiles; i++)
		{
			char name[56];
			pak->read(name, sizeof(name));
			pak->read((char *) &offset, sizeof(offset));
			pak->read((char *) &length, sizeof(length));

			files[name] = PakFile { pak, offset, length };
		}
	}

	// add loose files
	auto recurseLooseFiles = [](auto recurseLooseFiles, std::filesystem::path path, std::filesystem::path base) -> void
	{
		for (auto &f : std::filesystem::directory_iterator(path))
		{
			if (f.is_directory())
			{
				recurseLooseFiles(recurseLooseFiles, f, base);
				continue;
			}
			else if (f.path().extension().compare(".pak") == 0)
				continue;

			auto quakePath = f.path().lexically_relative(base).string();
			replace(quakePath, "\\", "/");
			files[quakePath] = f.path();
		}
	};

	recurseLooseFiles(recurseLooseFiles, baseq2, baseq2);
}

static void Fs_Init()
{
	AddGameDirectory("baseq2");
	AddGameDirectory(gi.cvar("game", "", 0)->string);
}

const std::shared_ptr<std::ifstream> Fs_Find(const char *name)
{
	auto it = files.find(name);

	if (it == files.end())
		return NULL;

	auto &variant = (*it).second;

	if (std::holds_alternative<std::filesystem::path>(variant))
	{
		auto stream = std::make_shared<std::ifstream>();
		stream->open(std::get<std::filesystem::path>(variant), std::ifstream::in | std::ifstream::binary);
		return stream;
	}

	auto &pak = std::get<PakFile>(variant);
	pak.file->seekg(pak.offset);
	return pak.file;
}

#include <vector>

#include <tuple>
// function has to live in the std namespace 
// so that it is picked up by argument-dependent name lookup (ADL).
namespace std{
    namespace
    {

        // Code from boost
        // Reciprocal of the golden ratio helps spread entropy
        //     and handles duplicates.
        // See Mike Seymour in magic-numbers-in-boosthash-combine:
        //     https://stackoverflow.com/questions/4948780

        template <class T>
        inline void hash_combine(std::size_t& seed, T const& v)
        {
            seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
        }

        // Recursive template code derived from Matthieu M.
        template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
        struct HashValueImpl
        {
          static void apply(size_t& seed, Tuple const& tuple)
          {
            HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
            hash_combine(seed, get<Index>(tuple));
          }
        };

        template <class Tuple>
        struct HashValueImpl<Tuple,0>
        {
          static void apply(size_t& seed, Tuple const& tuple)
          {
            hash_combine(seed, get<0>(tuple));
          }
        };
    }

    template <typename ... TT>
    struct hash<std::tuple<TT...>> 
    {
        size_t
        operator()(std::tuple<TT...> const& tt) const
        {                                              
            size_t seed = 0;                             
            HashValueImpl<std::tuple<TT...> >::apply(seed, tt);    
            return seed;                                 
        }                                              

    };
}

static std::unordered_map<std::tuple<std::string, int>, PxConvexMesh *> cookedModels;

PxConvexMesh *PhysX_CookMD2(const char *name, int frame)
{
	auto exist = cookedModels.find(std::make_tuple(name, frame));

	if (exist != cookedModels.end())
		return (*exist).second;

	auto pf = Fs_Find(name);
	int32_t frameSize;
	std::streamoff start = pf->tellg();
	pf->seekg(start + 16);
	pf->read((char *) &frameSize, sizeof(frameSize));

	pf->seekg(start + 24);

	int32_t numVertices;
	pf->read((char *) &numVertices, sizeof(numVertices));

	int32_t ofsFrames;
	pf->seekg(start + 56);
	pf->read((char *) &ofsFrames, sizeof(ofsFrames));

	pf->seekg(start + ofsFrames + (frameSize * frame));
	vec3_t scale, translate;
	pf->read((char *) scale, sizeof(scale));
	pf->read((char *) translate, sizeof(translate));
	char fn[16];
	pf->read(fn, sizeof(fn));

	std::vector<PxVec3> vecs;
	vecs.reserve(numVertices);

	for (int i = 0; i < numVertices; i++)
	{
		uint8_t v[4];
		pf->read((char *) v, 4);
		vecs.push_back(PxVec3((v[0] * scale[0]) + translate[0], (v[1] * scale[1]) + translate[1], (v[2] * scale[2]) + translate[2]));
	}

	PxConvexMeshDesc convexDesc;
	convexDesc.points.count     = numVertices;
	convexDesc.points.stride    = sizeof(PxVec3);
	convexDesc.points.data      = vecs.data();
	convexDesc.flags            = PxConvexFlag::eCOMPUTE_CONVEX;
	convexDesc.vertexLimit		= 32;

	PxDefaultMemoryOutputStream buf;
	if(!cooking->cookConvexMesh(convexDesc, buf))
		return nullptr;

	PxDefaultMemoryInputData input(buf.getData(), buf.getSize());
	PxConvexMesh* convexMesh = physics->createConvexMesh(input);

	cookedModels.insert(std::make_pair(std::make_tuple(name, frame), convexMesh));

	return convexMesh;
}

struct BSPLump
{
	uint32_t offset, size;
};

struct BSPPlane
{
	vec3_t	normal;
	float	distance;
	int32_t	type;
};

struct BSPBrushSide
{
    uint16_t plane;
    int16_t texInfo;
};

struct BSPBrush
{
    uint32_t firstSide;
    uint32_t numSides;
    uint32_t contents;
};

struct BSPModel
{
	vec3_t mins, maxs, origin;
	uint32_t headNode;
	uint32_t firstFace, numFaces;
};

struct BSPNode
{
	uint32_t plane;
	int32_t children[2];
	uint16_t mins[3];
	uint16_t maxs[3];
	uint16_t firstFace;
	uint16_t numFaces;
};

struct BSPLeaf
{
	uint32_t contents;
	uint16_t cluster;
	uint16_t area;
	uint16_t mins[3];
	uint16_t maxs[3];
	uint16_t firstLeafFace;
	uint16_t numLeafFaces;
	uint16_t firstLeafBrush;
	uint16_t numLeafBrushes;
};

template<typename T>
std::vector<T> ReadBSPLump(std::shared_ptr<std::ifstream> stream, fpos_t start, const BSPLump &lump)
{
	std::vector<T> data;

	data.resize(lump.size / sizeof(T));
	stream->seekg(start + lump.offset);
	stream->read((char *) data.data(), lump.size);

	return data;
}

#include <unordered_set>

std::vector<std::vector<PxConvexMeshGeometry>> bspMeshes;

void PhysX_CookBSP(std::shared_ptr<std::ifstream> stream)
{
	std::streamoff start = stream->tellg();
	stream->seekg(start + 8);
	BSPLump lumps[16];

	stream->read((char *) lumps, sizeof(lumps));

	// planes
	auto planes = ReadBSPLump<BSPPlane>(stream, start, lumps[1]);
	
	// nodes
	auto nodes = ReadBSPLump<BSPNode>(stream, start, lumps[4]);
	
	// leaves
	auto leaves = ReadBSPLump<BSPLeaf>(stream, start, lumps[8]);
	
	// leaves
	auto leafbrushes = ReadBSPLump<uint16_t>(stream, start, lumps[10]);
	
	// models
	auto models = ReadBSPLump<BSPModel>(stream, start, lumps[13]);

	// brush sides
	auto brushSides = ReadBSPLump<BSPBrushSide>(stream, start, lumps[15]);

	// brushes
	auto brushes = ReadBSPLump<BSPBrush>(stream, start, lumps[14]);

	std::vector<PxVec4> equations;
	std::vector<PxVec3> verts;

	auto isPointInsidePlanes = [&equations, &verts] (const PxVec3 &point, float margin) {
		int numbrushes = equations.size();
		for (int i = 0; i < numbrushes; i++)
		{
			auto &N1 = equations[i];
			float dist = N1.getXYZ().dot(point) + N1.w - margin;
			if (dist > 0.f)
				return false;
		}
		return true;
	};
	
	auto getVerticesFromPlaneEquations = [&equations, &verts, isPointInsidePlanes] () {
		int numbrushes = equations.size();

		// brute force:
		for (int i = 0; i < numbrushes; i++) {
			auto &N1 = equations[i];

			for (int j = i + 1; j < numbrushes; j++) {
				auto &N2 = equations[j];

				for (int k = j + 1; k < numbrushes; k++) {
					auto &N3 = equations[k];
				
					auto n2n3 = N2.getXYZ().cross(N3.getXYZ());
					auto n3n1 = N3.getXYZ().cross(N1.getXYZ());
					auto n1n2 = N1.getXYZ().cross(N2.getXYZ());

					if (!((n2n3.magnitudeSquared() > 0.0001f) &&
							(n3n1.magnitudeSquared() > 0.0001f) &&
							(n1n2.magnitudeSquared() > 0.0001f)))
						continue;

					// point P out of 3 plane equations:

					// 	     d1 ( N2 * N3 ) + d2 ( N3 * N1 ) + d3 ( N1 * N2 )  
					// P =  -------------------------------------------------------------------------  
					//    N1 . ( N2 * N3 )  
					float quotient = N1.getXYZ().dot(n2n3);
					if (fabsf(quotient) <= 0.000001f)
						continue;

					quotient = -1.f / quotient;
					n2n3 *= N1.w;
					n3n1 *= N2.w;
					n1n2 *= N3.w;
					auto potentialVertex = n2n3;
					potentialVertex += n3n1;
					potentialVertex += n1n2;
					potentialVertex *= quotient;

					// check if inside, and replace supportingVertexOut if needed
					if (isPointInsidePlanes(potentialVertex, 0.01f))
						verts.push_back(potentialVertex);
				}
			}
		}
	};

	auto recurseModelNode = [&nodes, &leaves, &leafbrushes] (int32_t node, auto &brushIDs, auto &recurse) {
		if (node < 0)
		{
			int32_t leaf = -(node + 1);

			for (int32_t l = 0; l < leaves[leaf].numLeafBrushes; l++)
				brushIDs.insert(leafbrushes[leaves[leaf].firstLeafBrush + l]);

			return;
		}
		
		recurse(nodes[node].children[0], brushIDs, recurse);
		recurse(nodes[node].children[1], brushIDs, recurse);
	};

	for (auto &m : models)
	{
		auto n = m.headNode;
		std::unordered_set<int32_t> brushIDs;
		std::vector<PxConvexMeshGeometry> meshes;

		recurseModelNode(n, brushIDs, recurseModelNode);

		for (auto &brushID : brushIDs)
		{
			auto &b = brushes[brushID];

			if (!(b.contents & MASK_SOLID))
				continue;

			for (int s = 0; s < b.numSides; s++)
			{
				auto side = brushSides[b.firstSide + s];
				auto plane = planes[side.plane];
				equations.push_back({ plane.normal[0], plane.normal[1], plane.normal[2], -plane.distance });
			}

			getVerticesFromPlaneEquations();

			if (verts.size())
			{
				PxConvexMeshDesc convexDesc;
				convexDesc.points.count     = verts.size();
				convexDesc.points.stride    = sizeof(PxVec3);
				convexDesc.points.data      = verts.data();
				convexDesc.flags            = PxConvexFlag::eCOMPUTE_CONVEX;
				convexDesc.vertexLimit		= 32;

				PxDefaultMemoryOutputStream buf;
				if (cooking->cookConvexMesh(convexDesc, buf))
				{
					PxDefaultMemoryInputData input(buf.getData(), buf.getSize());
					PxConvexMesh *convexMesh = physics->createConvexMesh(input);
					meshes.push_back(PxConvexMeshGeometry(convexMesh));
				}
			}

			equations.clear();
			verts.clear();
		}

		bspMeshes.push_back(std::move(meshes));
	}
}

#include <chrono>

using namespace std::chrono_literals;

#include <atomic>

const int32_t physicSubsteps = 10;
const float physicTime = FRAMETIME / physicSubsteps;

static std::atomic<int32_t> numRan(0);

static void PhysX_Thread(std::stop_token stoken)
{
	while (true)
	{
		if (stoken.stop_requested())
			return;
		else if (numRan < physicSubsteps)
		{
			scene->lockWrite();
			scene->simulate(physicTime);
			scene->fetchResults(true);
			numRan++;
			scene->unlockWrite();
		}
	}
}

physx::PxFilterFlags PxQuake2CollisionFilterShader(
	physx::PxFilterObjectAttributes attributes0, physx::PxFilterData filterData0,
	physx::PxFilterObjectAttributes attributes1, physx::PxFilterData filterData1,
	physx::PxPairFlags& retPairFlags, const void* constantBlock, physx::PxU32 constantBlockSize)
{
	return PxDefaultSimulationFilterShader(attributes0, filterData0, attributes1, filterData1, retPairFlags, constantBlock, constantBlockSize) | physx::PxFilterFlag::eCALLBACK;
}

static class PxQuake2FilterCallbackType : public PxSimulationFilterCallback
{
	PxFilterFlags	pairFound(	PxU32 pairID,
		PxFilterObjectAttributes attributes0, PxFilterData filterData0, const PxActor* a0, const PxShape* s0,
		PxFilterObjectAttributes attributes1, PxFilterData filterData1, const PxActor* a1, const PxShape* s1,
		PxPairFlags& pairFlags)
	{
		if (a0->userData && a1->userData)
		{
			edict_t *a0e = (edict_t *) a0->userData;
			edict_t *a1e = (edict_t *) a1->userData;

			if (a0e->owner == a1e || a1e->owner == a0e)
				return PxFilterFlag::eKILL;
		}

		return PxFilterFlags();
	}

	void			pairLost(	PxU32 pairID,
		PxFilterObjectAttributes attributes0,
		PxFilterData filterData0,
		PxFilterObjectAttributes attributes1,
		PxFilterData filterData1,
		bool objectRemoved)
	{
	}

	bool			statusChange(PxU32& pairID, PxPairFlags& pairFlags, PxFilterFlags& filterFlags)
	{
		return false;
	}
} PxQuake2FilterCallback;

void PhysX_Init()
{
	srand(time(NULL));

	Fs_Init();

	foundation = PxCreateFoundation(PX_PHYSICS_VERSION, physxAllocatorCallback, physxErrorCallback);
	if(!foundation)
		gi.error((char *)"PxCreateFoundation failed!");

	PxTolerancesScale tolerances;
	tolerances.length = 32.f;
	tolerances.speed = sv_gravity->value;

	cooking = PxCreateCooking(PX_PHYSICS_VERSION, *foundation, PxCookingParams(tolerances));
	if (!cooking)
		gi.error((char *)"PxCreateCooking failed!");

	bool recordMemoryAllocations = true;

	//pvd = PxCreatePvd(*foundation);
	//PxPvdTransport* transport = PxDefaultPvdSocketTransportCreate(PVD_HOST, 5425, 10000);
	//PxPvdTransport* transport = PxDefaultPvdFileTransportCreate("physx/transport.pxd2");
	//pvd->connect(*transport,PxPvdInstrumentationFlag::eDEBUG);

	physics = PxCreateBasePhysics(PX_PHYSICS_VERSION, *foundation,
		tolerances, recordMemoryAllocations, pvd);
	if(!physics)
		gi.error((char *)"PxCreatePhysics failed!");

	PxSceneDesc sceneDesc(physics->getTolerancesScale());
	sceneDesc.gravity = PxVec3(0.0f, 0.0f, -sv_gravity->value);
	dispatcher = PxDefaultCpuDispatcherCreate(4);
	sceneDesc.cpuDispatcher	= dispatcher;
	sceneDesc.filterShader	= PxQuake2CollisionFilterShader;
	sceneDesc.filterCallback = &PxQuake2FilterCallback;
	scene = physics->createScene(sceneDesc);
	if(!scene)
		gi.error((char *)"createScene failed!");

	material = physics->createMaterial(0.5f, 0.5f, 0.f);
	material->setFlag(PxMaterialFlag::eIMPROVED_PATCH_FRICTION, true);

	highFrictionMaterial = physics->createMaterial(0.7f, 0.7f, 0.8f);
	highFrictionMaterial->setFlag(PxMaterialFlag::eIMPROVED_PATCH_FRICTION, true);

	PhysX_CookBSP(Fs_Find(std::format("maps/{}.bsp", level.mapname).c_str()));

	physxThread = std::jthread(PhysX_Thread);
}

PxRigidDynamic *PhysX_ConvertToPhysX(edict_t *self, const PxGeometry &geometry, const float &density)
{
	self->movetype = MOVETYPE_PHYSX;

	auto box = PxCreateDynamic(*physics, TransformFromQuake(self), geometry, *highFrictionMaterial, density);
	box->setLinearVelocity(PxVec3(self->velocity[0], self->velocity[1], self->velocity[2]));
	box->setAngularVelocity(PxVec3(self->avelocity[0], self->avelocity[1], self->avelocity[2]));
	box->userData = self;
	self->physicsBody = box;
	scene->lockWrite();
	scene->addActor(*box);
	scene->unlockWrite();

	return box;
}

void PhysX_MakeBrushModel(edict_t *bmodel)
{
	if (bmodel->physicsBody)
		return;

	int id = bmodel->s.modelindex - 1;

	if (bmodel->s.modelindex == 1)
	{
		auto body = physics->createRigidStatic(PxTransform(PxVec3(0, 0, 0), PxIdentity));

		for (auto &brush : bspMeshes[id])
		{
			auto shape = physics->createShape(brush, *material, true);
			body->attachShape(*shape);
			shape->release();
		}

		bmodel->physicsBody = body;
		
		scene->lockWrite();
		scene->addActor(*body);
		scene->unlockWrite();
	}
	else
	{
		auto body = physics->createRigidDynamic(PxTransform(PxVec3(0, 0, 0), PxIdentity));
		body->setRigidBodyFlag(PxRigidBodyFlag::eKINEMATIC, true);

		for (auto &brush : bspMeshes[id])
		{
			auto shape = physics->createShape(brush, *material, true);
			body->attachShape(*shape);
			shape->release();
		}

		bmodel->physicsBody = body;
		
		scene->lockWrite();
		scene->addActor(*body);
		scene->unlockWrite();
	}
}

void PhysX_DestroyBrushModel(edict_t *bmodel)
{
	if (!bmodel->physicsBody)
		return;
	
	scene->lockWrite();
	scene->removeActor(*((PxActor *)bmodel->physicsBody));
	scene->unlockWrite();
}

#define lengthof(array) (sizeof(array) / sizeof(array[0]))

void PhysX_MakePlayer(edict_t *player)
{
	auto body = PxCreateKinematic(*physics, 
		PxTransform(PxVec3(player->s.origin[0], player->s.origin[1], player->s.origin[2]), PxQuat(PxHalfPi, PxVec3(0, 0, 1))),
		PxCapsuleGeometry(player->maxs[0] / 3, (player->maxs[2] - player->mins[2]) / 2),
		*material,
		0.5,
		PxTransform(PxVec3(0, 0, -(player->maxs[2] + player->mins[2]) - 4)));
	player->physicsBody = body;
	body->userData = player;
	scene->lockWrite();
	scene->addActor(*body);
	scene->unlockWrite();
}

void PhysX_UpdatePlayer(edict_t *player)
{
	assert(player->physicsBody);

	auto body = (PxRigidDynamic *)player->physicsBody;
	scene->lockWrite();
	body->setKinematicTarget(PxTransform(PxVec3(player->s.origin[0], player->s.origin[1], player->s.origin[2]), PxQuat(PxHalfPi, PxVec3(0, 0, 1))));
	scene->unlockWrite();
}

void PhysX_DestroyPlayer(edict_t *player)
{
	assert(player->physicsBody);
	
	scene->lockWrite();
	scene->removeActor(*(PxRigidBody *)player->physicsBody);
	((PxRigidBody *)player->physicsBody)->release();
	scene->unlockWrite();
	player->physicsBody = nullptr;
}

void PhysX_Simulate()
{
	scene->lockWrite();
	if (numRan != physicSubsteps)
	{
		int numLeft = physicSubsteps - numRan;
		gi.dprintf("Didn't run %i frames, catching up with %i\n", physicSubsteps, numLeft);
		for (int i = 0; i < numLeft; i++)
		{
			scene->simulate(physicTime);
			scene->fetchResults(true);
		}
	}

	numRan = 0;

	for (int i = game.maxclients + 1; i < globals.num_edicts; i++)
	{
		auto *barrel = &g_edicts[i];

		if (!barrel->physicsBody)
			continue;

		if (barrel->solid == SOLID_BSP)
		{
			auto box = ((PxRigidDynamic *)barrel->physicsBody);

			box->setKinematicTarget(PxTransform(PxVec3(barrel->s.origin[0], barrel->s.origin[1], barrel->s.origin[2]), QuatFromEuler(barrel->s.angles)));
		}
		else
		{
			auto box = ((PxRigidActor *)barrel->physicsBody);

			QuakeFromTransform(box->getGlobalPose(), barrel->s.origin, barrel->s.angles);

			gi.linkentity(barrel);
		}
	}
	scene->unlockWrite();
}

void PhysX_Shutdown()
{
	if (physxThread.joinable())
	{
		physxThread.request_stop();
		physxThread.join();
	}

	if (cooking)
		cooking->release();
	if (physics)
		physics->release();
	if (pvd)
		pvd->release();
	if (foundation)
		foundation->release();
}