// from https://github.com/wonkoRT
//  at gist https://gist.github.com/wonkoRT/4a2d07cc292126bddf10/
// cpp port by @saephoed
// ogre based (working; some minor bugs; no ui) port of

// see http://experilous.com/1/planet-generator/2014-09-28/planet-generator.js
// see http://experilous.com/1/blog/post/procedural-planet-generation
// see twitter @AndyGainey

/// usage somehow like ///
/*
void generateWorld()
{
  // get params
  
  world_ = new simu::World(mGUI, mPlatform, mRenderWindow, sceneMgr_);
	worldNode_ = sceneMgr_->getRootSceneNode()->createChildSceneNode();
	surfaceNode_ = worldNode_->createChildSceneNode();
	plateBoundariesNode_ = worldNode_->createChildSceneNode();
	plateMovementsNode_ = worldNode_->createChildSceneNode();
	airCurrentsNode_ = worldNode_->createChildSceneNode();
	planet_ = new simu::Planet();
	planet_->renderData.surface = sceneMgr_->createManualObject("surfaceObj");
	planet_->renderData.plateBoundaries = sceneMgr_->createManualObject("plateBoundariesObj");
	planet_->renderData.plateMovements = sceneMgr_->createManualObject("plateMovementsObj");
	planet_->renderData.airCurrents = sceneMgr_->createManualObject("airCurrentsObj");

	world_->generatePlanetAsynchronous(planet_, originalSeed, seed, subdivisions, distortionRate, plateCount, oceanicRate, heatLevel, moistureLevel);

	surfaceNode_->attachObject(planet_->renderData.surface);
	plateBoundariesNode_->attachObject(planet_->renderData.plateBoundaries);
	plateMovementsNode_->attachObject(planet_->renderData.plateMovements);
	airCurrentsNode_->attachObject(planet_->renderData.airCurrents);

	// update stats
}

/// ogre material ///
material experilousworld/planetmat1
{
   technique
   {
      pass
      {
         diffuse vertexcolour
         specular vertexcolour
         ambient vertexcolour
         lighting on
      }
   }
}
*/


/// hpp ///

#ifndef INCLUDED_World
#define INCLUDED_World

#include <list>
#include <vector>

#include <Ogre.h>
#include "OgreVector3.h"
#include "OgreQuaternion.h"


namespace simu {

	struct Border;
	struct Tile;
	struct Plate;


	struct XorShift128
	{
		unsigned long x_;
		unsigned long y_;
		unsigned long z_;
		unsigned long w_;

		XorShift128(unsigned long x = 0, unsigned long y = 0, unsigned long z = 0, unsigned long w = 0);

		unsigned long next();
		double unit();
		double unitInclusive();
		unsigned long integer(unsigned long min, unsigned long max);
		unsigned long integerExclusive(unsigned long min, unsigned long max);
		double real(double min, double max);
		double realInclusive(double min, double max);
		void reseed(unsigned long x, unsigned long y, unsigned long z, unsigned long w);
	};


	struct UidObj
	{
		size_t id;
		virtual bool operator ==(const UidObj& other) const
		{
			return id == other.id;
		}

		UidObj(size_t p_id)
			: id(p_id)
		{}
	};

	struct Corner : public UidObj
	{
		Ogre::Vector3 position;
		double area;
		std::vector<Corner*> corners;
		std::vector<Border*> borders;
		std::vector<Tile*> tiles;
		double distanceToPlateRoot;
		double distanceToPlateBoundary;
		bool betweenPlates;
		double pressure;
		double shear;
		double elevation;
		Ogre::Vector3 airCurrent;
		double airCurrentSpeed; //kilometers per hour
		std::vector<double> airCurrentOutflows;
		// AirHeat
		double temperature;
		// 2do: temporaries for temperature calc. ?!
		double airHeat;
		double newAirHeat;
		double heat;
		double heatAbsorption;
		double maxHeat;
		// AirMoisture
		double moisture;
		// 2do: temporaries for temperature calc. ?!
		double airMoisture;
		double newAirMoisture;
		double precipitation;
		double precipitationRate;
		double maxPrecipitation;

		Corner(size_t p_id, const Ogre::Vector3& p_position, size_t cornerCount, size_t borderCount, size_t tileCount);
		Ogre::Vector3 vectorTo(const Corner& corner);
		std::string toString();
	};
	struct Tile : public UidObj
	{
		Ogre::Vector3 position;
		Ogre::Sphere boundingSphere;
		Ogre::Vector3 normal;
		Ogre::Vector3 averagePosition;
		double area;
		double elevation;
		std::vector<Corner*> corners;
		std::vector<Border*> borders;
		std::vector<Tile*> tiles;
		Plate* plate;
		double temperature;
		double moisture;
		std::string biome;
		Ogre::Vector3 plateMovement;

		Tile(size_t p_id, const Ogre::Vector3& p_position, size_t cornerCount, size_t borderCount, size_t tileCount);
		bool intersectRay(const Ogre::Ray& ray);
		std::string toString();
	};
	struct Border : public UidObj
	{
		std::vector<Corner*> corners;
		std::vector<Border*> borders;
		std::vector<Tile*> tiles;
		Ogre::Vector3 midpoint;
		bool betweenPlates;

		Border(size_t p_id, size_t cornerCount, size_t borderCount, size_t tileCount);
		Corner& oppositeCorner(const Corner& corner);
		const Corner& oppositeCorner(const Corner& corner) const;
		Tile& oppositeTile(const Tile& tile);
		const Tile& oppositeTile(const Tile& tile) const;
		double length() const;
		bool isLandBoundary() const;
		std::string toString() const;
	};




	// 2do attention; order of members matters because of initializer values, see f.e. http://stackoverflow.com/questions/11516657/c-structure-initialization
	struct Face
	{
		std::vector<size_t> n;
		std::vector<size_t> e;

		Ogre::Vector3 centroid;
		Ogre::Sphere boundingSphere;
		std::vector<Tile*> children;
	};
	struct Edge
	{
		std::vector<size_t> n;
		std::vector<size_t> f;

		std::vector<size_t> subdivided_n;
		std::vector<size_t> subdivided_e;
	};
	struct Node
	{
		Ogre::Vector3 p;

		std::vector<size_t> e;
		std::vector<size_t> f;
	};
	struct Mesh
	{
		std::vector<Node*> nodes;
		std::vector<Edge*> edges;
		std::vector<Face*> faces;
		std::vector<Corner*> corners;
		std::vector<Border*> borders;
	};

	struct Whorl
	{
		Whorl()
			: center()
			, strength(0)
			, radius(0)
		{}
		Whorl(Ogre::Vector3&& p_center, double p_strength, double p_radius)
			: center(p_center)
			, strength(p_strength)
			, radius(p_radius)
		{}
		Ogre::Vector3 center;
		double strength;
		double radius;
	};

	struct ElevationBorderOrigin
	{
		typedef std::function<double(double, double, double, double, double, double)> CalculateElevationFunc;
		ElevationBorderOrigin()
			: corner(nullptr)
			, pressure(0)
			, shear(0)
			, plate(nullptr)
		{}
		ElevationBorderOrigin(Corner* p_corner, double p_pressure, double p_shear, Plate* p_plate, CalculateElevationFunc& p_calculateElevation)
			: corner(p_corner)
			, pressure(p_pressure)
			, shear(p_shear)
			, plate(p_plate)
			, calculateElevation(p_calculateElevation)
		{}
		Corner* corner;
		double pressure;
		double shear;
		Plate* plate;
		CalculateElevationFunc calculateElevation;
	};
	struct ElevationBorder
	{
		ElevationBorder()
			: border(nullptr)
			, corner(nullptr)
			, nextCorner(nullptr)
			, distanceToPlateBoundary(0)
		{}
		ElevationBorder(ElevationBorderOrigin& p_origin, Border* p_border, Corner* p_corner, Corner* p_nextCorner, double p_distanceToPlateBoundary)
			: origin(p_origin)
			, border(p_border)
			, corner(p_corner)
			, nextCorner(p_nextCorner)
			, distanceToPlateBoundary(p_distanceToPlateBoundary)
		{}
		ElevationBorderOrigin origin; // 2do: maybe wrong storage because there are 2 distinctive creation types, so is the ref one wrong?
		Border* border;
		Corner* corner;
		Corner* nextCorner;
		double distanceToPlateBoundary;
	};
	struct Stress
	{
		Stress()
			: pressure(0)
			, shear(0)
		{}
		Stress(double p_pressure, double p_shear)
			: pressure(p_pressure)
			, shear(p_shear)
		{}
		double pressure;
		double shear;
	};
	// 2do attention end





	struct Plate
	{
		Ogre::ColourValue color;
		Ogre::Vector3 driftAxis;
		double driftRate;
		double spinRate;
		double elevation;
		bool oceanic;
		Corner* root;
		std::vector<Tile*> tiles;
		std::vector<Corner*> boundaryCorners;
		std::vector<Border*> boundaryBorders;
		double area;
		double circumference;

		Plate(const Ogre::ColourValue& p_color, const Ogre::Vector3& p_driftAxis, double p_driftRate, double p_spinRate, double p_elevation, bool p_oceanic, Corner& p_root)
			: color(p_color)
			, driftAxis(p_driftAxis)
			, driftRate(p_driftRate)
			, spinRate(p_spinRate)
			, elevation(p_elevation)
			, oceanic(p_oceanic)
			, root(&p_root)
			, tiles()
			, boundaryCorners()
			, boundaryBorders()
			, area(0)
			, circumference(0)
		{}
		Ogre::Vector3 calculateMovement(const Ogre::Vector3& position);
	};

	struct SpatialPartition
	{
		Ogre::Sphere boundingSphere;
		std::vector<SpatialPartition> partitions;
		std::vector<Tile*> tiles;

		SpatialPartition()
			: boundingSphere()
			, partitions()
			, tiles()
		{}
		SpatialPartition(const Ogre::Sphere& p_boundingSphere, const std::vector<SpatialPartition>& p_partitions, const std::vector<Tile*>& p_tiles)
			: boundingSphere(p_boundingSphere)
			, partitions(p_partitions)
			, tiles(p_tiles)
		{}
		bool intersectRay(const Ogre::Ray& ray);
	};

	struct Topology
	{
		std::vector<Corner> corners;
		std::vector<Border> borders;
		std::vector<Tile> tiles;
	};


	template<typename T>
	struct ValueStats
	{
		ValueStats()
		{
			reset();
		}
		T min, avg, max;
		void reset()
		{
			min = std::numeric_limits<T>::max();
			max = std::numeric_limits<T>::lowest();
			avg = 0;
		};
	};
	struct PlanetStatistics
	{
		struct Corners
		{
			size_t count;
			ValueStats<double> airCurrent;
			ValueStats<double> elevation;
			ValueStats<double> temperature;
			ValueStats<double> moisture;
			ValueStats<double> distanceToPlateBoundary;
			ValueStats<double> distanceToPlateRoot;
			ValueStats<double> pressure;
			ValueStats<double> shear;
			size_t doublePlateBoundaryCount;
			size_t triplePlateBoundaryCount;
			size_t innerLandBoundaryCount;
			size_t outerLandBoundaryCount;
		};
		Corners corners;

		struct Borders
		{
			size_t count;
			ValueStats<double> length;
			size_t plateBoundaryCount;
			double plateBoundaryPercentage;
			size_t landBoundaryCount;
			double landBoundaryPercentage;
		};
		Borders borders;

		struct Tiles
		{
			size_t count;
			double totalArea;
			ValueStats<double> area;
			ValueStats<double> elevation;
			ValueStats<double> temperature;
			ValueStats<double> moisture;
			ValueStats<double> plateMovement;
			std::map<std::string, size_t> biomeCounts;
			std::map<std::string, double> biomeAreas;
			size_t pentagonCount;
			size_t hexagonCount;
			size_t heptagonCount;
		};
		Tiles tiles;

		struct Plates
		{
			size_t count;
			ValueStats<size_t> tileCount;
			ValueStats<double> area;
			ValueStats<double> boundaryElevation;
			ValueStats<size_t> boundaryBorders;
			ValueStats<double> circumference;
		};
		Plates plates;
	};

	struct SurfaceRenderObj
	{
		std::vector<Ogre::ColourValue> terrainColors;
		std::vector<Ogre::ColourValue> plateColors;
		std::vector<Ogre::ColourValue> elevationColors;
		std::vector<Ogre::ColourValue> temperatureColors;
		std::vector<Ogre::ColourValue> moistureColors;
	};
	struct RenderData
	{
		Ogre::ManualObject* surface;
		Ogre::ManualObject* plateBoundaries;
		Ogre::ManualObject* plateMovements;
		Ogre::ManualObject* airCurrents;
	};

	struct Planet
	{
		Planet() {}

		unsigned long seed;
		unsigned long originalSeed;
		Topology topology;

		std::vector<Plate> plates;
		SpatialPartition partition;
		RenderData renderData;
		PlanetStatistics statistics;
	};


	namespace tools
	{
		Ogre::Vector3 slerp(const Ogre::Vector3& p0, const Ogre::Vector3& p1, double t);
		Ogre::Vector3 randomUnitVector(XorShift128& random);
		Ogre::Quaternion randomQuaternion(XorShift128& random);
		Ogre::Vector3 projectOnVector(const Ogre::Vector3& v1, const Ogre::Vector3& v2);
		bool intersectRayWithSphere(const Ogre::Ray& ray, const Ogre::Sphere& sphere);
		double calculateTriangleArea(const Ogre::Vector3& pa, const Ogre::Vector3& pb, const Ogre::Vector3& pc);
		double adjustRange(double value, double oldMin, double oldMax, double newMin, double newMax);
		void setLength(Ogre::Vector3& v, double len);
		Ogre::Vector3 setLength(Ogre::Vector3&& v, double len);
		Ogre::ColourValue ocv(unsigned int rgb);
	}


	class World
	{
	public:
		World(MyGUI::Gui* gui, MyGUI::OgrePlatform* platform, Ogre::RenderWindow* renderWindow, Ogre::SceneManager* sceneMgr);
		~World();

		void generatePlanetAsynchronous(Planet*& planet, unsigned long originalSeed, unsigned long seed, size_t icosahedronSubdivision, double topologyDistortionRate, size_t plateCount, double oceanicRate, double heatLevel, double moistureLevel);

	private:
		MyGUI::Gui* mGUI;
		MyGUI::OgrePlatform* mPlatform;
		Ogre::RenderWindow* mRenderWindow;
		Ogre::SceneManager* sceneMgr_;

	private:
		typedef std::function<bool(const Node&, const Node&, const Node&, const Node&)> RotationPredicateType;

		void generatePlanet(size_t icosahedronSubdivision, double topologyDistortionRate, size_t plateCount, double oceanicRate, double heatLevel, double moistureLevel, XorShift128& random, Planet& planet);
		Mesh generatePlanetMesh(size_t icosahedronSubdivision, double topologyDistortionRate, XorShift128& random);
		void generateIcosahedron(Mesh& ret);
		void generateSubdividedIcosahedron(size_t degree, Mesh& ret);
		size_t getEdgeOppositeFaceIndex(const Edge& edge, size_t faceIndex);
		size_t getFaceOppositeNodeIndex(const Face& face, const Edge& edge);
		size_t findNextFaceIndex(const Mesh& mesh, size_t nodeIndex, size_t faceIndex);
		bool conditionalRotateEdge(const Mesh& mesh, size_t edgeIndex, RotationPredicateType& predicate);
		Ogre::Vector3 calculateFaceCentroid(const Ogre::Vector3& pa, const Ogre::Vector3& pb, const Ogre::Vector3& pc);
		bool distortMesh(const Mesh& mesh, size_t degree, XorShift128& random);
		double relaxMesh(const Mesh& mesh, double multiplier);
		bool generatePlanetTopology(const Mesh& mesh, Topology& ret);
		void generatePlanetPartition(std::vector<Tile>& tiles, SpatialPartition& rootPartition);
		void generatePlanetTerrain(Planet& planet, size_t plateCount, double oceanicRate, double heatLevel, double moistureLevel, XorShift128& random);
		void generatePlanetTectonicPlates(Topology& topology, size_t plateCount, double oceanicRate, XorShift128& random, std::vector<Plate>& plates);
		void calculateCornerDistancesToPlateRoot(std::vector<Plate>& plates);
		void generatePlanetElevation(Topology& topology, std::vector<Plate>& plates);
		void identifyBoundaryBorders(std::vector<Border>& borders);
		void collectBoundaryCorners(std::vector<Corner>& corners, std::vector<Corner*>& boundaryCorners);
		void calculatePlateBoundaryStress(const std::vector<Corner*>& boundaryCorners, std::vector<size_t>& boundaryCornerInnerBorderIndexes);
		void calculateStress(const Ogre::Vector3& movement0, const Ogre::Vector3& movement1, const Ogre::Vector3& boundaryVector, const Ogre::Vector3& boundaryNormal, Stress& stress);
		void blurPlateBoundaryStress(const std::vector<Corner*>& boundaryCorners, size_t stressBlurIterations, double stressBlurCenterWeighting);
		void populateElevationBorderQueue(const std::vector<Corner*>& boundaryCorners, const std::vector<size_t>& boundaryCornerInnerBorderIndexes, std::list<ElevationBorder>& elevationBorderQueue);
		void processElevationBorderQueue(std::list<ElevationBorder>& elevationBorderQueue, std::function<bool(const ElevationBorder&, const ElevationBorder&)> elevationBorderQueueSorter);
		void calculateTileAverageElevations(std::vector<Tile>& tiles);
		void generatePlanetWeather(Topology& topology, SpatialPartition& partition, double heatLevel, double moistureLevel, XorShift128& random);
		void generateAirCurrentWhorls(double planetRadius, XorShift128& random, std::vector<Whorl>& whorls);
		void calculateAirCurrents(std::vector<Corner>& corners, const std::vector<Whorl>& whorls, double planetRadius);
		void initializeAirHeat(std::vector<Corner>& corners, double heatLevel, std::list<Corner*>& activeCorners, double& airHeat);
		double processAirHeat(std::list<Corner*>& activeCorners);
		void calculateTemperature(std::vector<Corner>& corners, std::vector<Tile>& tiles, double planetRadius);
		void initializeAirMoisture(std::vector<Corner>& corners, double moistureLevel, std::list<Corner*>& activeCorners, double& airMoisture);
		double processAirMoisture(std::list<Corner*>& activeCorners);
		void calculateMoisture(std::vector<Corner>& corners, std::vector<Tile>& tiles);
		void generatePlanetBiomes(std::vector<Tile>& tiles, double planetRadius);
		void generatePlanetRenderData(Topology& topology, XorShift128& random, RenderData& renderData);
		void buildSurfaceRenderObject(std::vector<Tile>& tiles, XorShift128& random, Ogre::ManualObject& mo);
		void buildPlateBoundariesRenderObject(std::vector<Border>& borders, Ogre::ManualObject& mo);
		void buildPlateMovementsRenderObject(std::vector<Tile>& tiles, Ogre::ManualObject& mo);
		void buildAirCurrentsRenderObject(std::vector<Corner>& corners, Ogre::ManualObject& mo);
		void buildArrow(Ogre::ManualObject& mo, const Ogre::Vector3& position, const Ogre::Vector3& direction, const Ogre::Vector3& normal, double baseWidth, const Ogre::ColourValue& color);
		void buildTileWedge(std::vector<Face>& f, size_t b, size_t s, size_t t, size_t n);
		void buildTileWedgeColors(std::vector<Ogre::ColourValue>& f, const Ogre::ColourValue& c, const Ogre::ColourValue& bc);
		void generatePlanetStatistics(Topology& topology, std::vector<Plate>& plates, PlanetStatistics& planetStatistics);
	};


}

#endif // #ifndef INCLUDED_World



/// cpp ///

#include "World.hpp"
#include "OgreRay.h"
#include "OgreSphere.h"

using namespace simu::tools;


namespace simu {


	template <typename T>
	typename T::iterator getAdvancedIt(typename T& value, size_t off)
	{
		T::iterator it = value.begin();
		std::advance(it, off);
		return it;
	}
	template <typename T>
	typename T::value_type& getAdvancedItVal(typename T& value, size_t off)
	{
		auto it = value.begin();
		std::advance(it, off);
		return *it;
	}

	XorShift128::XorShift128(unsigned long x, unsigned long y, unsigned long z, unsigned long w)
			//: x_(x ? x : 123456789)
			//, y_(y ? y : 362436069)
			//, z_(z ? z : 521288629)
			//, w_(w ? w : 88675123)
		{
			reseed(x, y, z, w);
		}

	unsigned long XorShift128::next()
		{
			unsigned long t = x_ ^ (x_ << 11) & 0x7FFFFFFF;
			x_ = y_;
			y_ = z_;
			z_ = w_;
			w_ = (w_ ^ (w_ >> 19)) ^ (t ^ (t >> 8));
			return w_;
		};

	double XorShift128::unit()
		{
			return (double)next() / (double)0x80000000;
		};

	double XorShift128::unitInclusive()
		{
			return (double)next() / (double)0x7FFFFFFF;
		};

	unsigned long XorShift128::integer(unsigned long min, unsigned long max)
		{
			return integerExclusive(min, max + 1);
		};

	unsigned long XorShift128::integerExclusive(unsigned long min, unsigned long max)
		{
			return std::floor(unit() * (max - min)) + min;
		};

	double XorShift128::real(double min, double max)
		{
			return unit() * (max - min) + min;
		};

	double XorShift128::realInclusive(double min, double max)
		{
			return unitInclusive() * (max - min) + min;
		};

	void XorShift128::reseed(unsigned long x, unsigned long y, unsigned long z, unsigned long w)
		{
			x_ = (x ? x : 123456789);
			y_ = (y ? y : 362436069);
			z_ = (z ? z : 521288629);
			w_ = (w ? w : 88675123);
		};



	template <typename T, typename F>
	T lerp(const T& from, const T& to, const F f)
	{
		return from + (to - from) * f;
	}

	namespace tools {

		Ogre::Vector3 slerp(const Ogre::Vector3& p0, const Ogre::Vector3& p1, double t)
		{
			double omega = std::acos(p0.dotProduct(p1));
			return p0 * std::sin((1 - t) * omega) + p1 * std::sin(t * omega) / std::sin(omega);
		}

		Ogre::Vector3 randomUnitVector(XorShift128& random)
		{
			double theta = random.real(0, M_PI * 2);
			double phi = std::acos(random.realInclusive(-1, 1));
			double sinPhi = std::sin(phi);
			return Ogre::Vector3(
				std::cos(theta) * sinPhi,
				std::sin(theta) * sinPhi,
				std::cos(phi));
		}

		Ogre::Quaternion randomQuaternion(XorShift128& random)
		{
			double theta = random.real(0, M_PI * 2);
			double phi = std::acos(random.realInclusive(-1, 1));
			double sinPhi = std::sin(phi);
			double gamma = random.real(0, M_PI * 2);
			double sinGamma = std::sin(gamma);
			return Ogre::Quaternion(
				std::cos(theta) * sinPhi * sinGamma,
				std::sin(theta) * sinPhi * sinGamma,
				std::cos(phi) * sinGamma,
				std::cos(gamma));
		}

		Ogre::Vector3 projectOnVector(const Ogre::Vector3& v1, const Ogre::Vector3& v2)
		{
			auto ang = v1.angleBetween(v2);
			auto ret = v1 / v1.length() * std::cos(ang.valueRadians());
			return ret;
		}

		bool intersectRayWithSphere(const Ogre::Ray& ray, const Ogre::Sphere& sphere)
		{
			auto v1 = sphere.getCenter() - ray.getOrigin();

			auto v2 = projectOnVector(v1, ray.getDirection());

			double d = v1.distance(v2);
			return d <= sphere.getRadius();
		}

		double calculateTriangleArea(const Ogre::Vector3& pa, const Ogre::Vector3& pb, const Ogre::Vector3& pc)
		{
			auto vab = pb - pa;
			auto vac = pc - pa;
			auto faceNormal = vab.crossProduct(vac);
			auto vabNormal = faceNormal.crossProduct(vab);
			vabNormal.normalise();
			auto plane = Ogre::Plane(vabNormal, pa);
			auto height = plane.getDistance(pc);
			auto width = vab.length();
			double area = width * height * 0.5;
			return area;
		}

		double adjustRange(double value, double oldMin, double oldMax, double newMin, double newMax)
		{
			return (value - oldMin) / (oldMax - oldMin) * (newMax - newMin) + newMin;
		}

		void setLength(Ogre::Vector3& v, double len)
		{
			v.normalise();
			v *= len;
		}
		Ogre::Vector3 setLength(Ogre::Vector3&& v, double len)
		{
			v.normalise();
			v *= len;
			return v;
		}

		Ogre::ColourValue ocv(unsigned int rgb)
		{
			Ogre::ColourValue ret;
			ret.setAsARGB(0xFF000000 | rgb);
			return ret;
		}
	}


	template <typename T, typename F>
	bool findIndex(const T& container, const F& value, size_t& idx)
	{
		auto it = std::find(container.begin(), container.end(), value);
		if (it == container.end())
			return false;
		idx = std::distance(container.begin(), it);
		return true;
	}

	template <typename T, typename F>
	bool removeIfFindIndex(T& container, const F& value)
	{
		auto it = std::find(container.begin(), container.end(), value);
		if (it == container.end())
			return false;
		//container.remove(it);
		container.erase(it);
		return true;
	}









	Corner::Corner(size_t p_id, const Ogre::Vector3& p_position, size_t cornerCount, size_t borderCount, size_t tileCount)
		: UidObj(p_id)
		, position(p_position)
		, area(0)
		, corners()
		, borders()
		, tiles()
		, distanceToPlateRoot(0)
		, distanceToPlateBoundary(0)
		, betweenPlates(false)
		, pressure(0)
		, shear(0)
		, elevation(0)
		, airCurrent()
		, airCurrentSpeed(0)
		, airCurrentOutflows()
		, temperature(0)
		, airHeat(0)
		, newAirHeat(0)
		, heat(0)
		, heatAbsorption(0)
		, maxHeat(0)
		, moisture(0)
		, airMoisture(0)
		, newAirMoisture(0)
		, precipitation(0)
		, precipitationRate(0)
		, maxPrecipitation(0)
	{
		corners.resize(cornerCount);
		borders.resize(borderCount);
		tiles.resize(tileCount);
	}
	Ogre::Vector3 Corner::vectorTo(const Corner& corner)
	{
		return corner.position - position;
	};
	std::string Corner::toString()
	{
		std::stringstream ss;
		ss << "Corner " << id << " < " << position.x << ", " << position.y << ", " << position.z << " >";
		return ss.str();
	};






	Tile::Tile(size_t p_id, const Ogre::Vector3& p_position, size_t cornerCount, size_t borderCount, size_t tileCount)
		: UidObj(p_id)
		, position(p_position)
		, boundingSphere()
		, normal()
		, averagePosition()
		, area(0)
		, elevation(0)
		, corners()
		, borders()
		, tiles()
		, plate(nullptr)
		, temperature(0)
		, moisture(0)
		, biome()
		, plateMovement()
	{
		corners.resize(cornerCount);
		borders.resize(borderCount);
		tiles.resize(tileCount);
	}
	bool Tile::intersectRay(const Ogre::Ray& ray)
	{
		if (!intersectRayWithSphere(ray, boundingSphere)) return false;

		Ogre::Plane surface(normal, averagePosition);
		if (surface.getDistance(ray.getOrigin()) <= 0) return false;

		auto denominator = surface.normal.dotProduct(ray.getDirection());
		if (denominator == 0) return false;

		auto t = -(ray.getOrigin().dotProduct(surface.normal) - std::abs(surface.getDistance(Ogre::Vector3::ZERO))) / denominator;
		auto point = ray.getDirection() * t + ray.getOrigin();

		auto origin = Ogre::Vector3::ZERO;
		for (auto i = 0; i < corners.size(); ++i)
		{
			auto j = (i + 1) % corners.size();
			Ogre::Plane side(corners[j]->position, corners[i]->position, origin);

			if (side.getDistance(point) < 0) return false;
		}

		return true;
	};
	std::string Tile::toString()
	{
		std::stringstream ss;
		ss << "Tile " << id << " (" << tiles.size() << " Neighbors) < " << position.x << ", " << position.y << ", " << position.z << " >";
		return ss.str();
	};






	Border::Border(size_t p_id, size_t cornerCount, size_t borderCount, size_t tileCount)
		: UidObj(p_id)
		, corners()
		, borders()
		, tiles()
		, midpoint()
		, betweenPlates(false)
	{
		corners.resize(cornerCount);
		borders.resize(borderCount);
		tiles.resize(tileCount);
	}
	Corner& Border::oppositeCorner(const Corner& corner)
	{
		return (*corners[0] == corner) ? *corners[1] : *corners[0];
	};
	const Corner& Border::oppositeCorner(const Corner& corner) const
	{
		return (*corners[0] == corner) ? *corners[1] : *corners[0];
	};
	Tile& Border::oppositeTile(const Tile& tile)
	{
		return (*tiles[0] == tile) ? *tiles[1] : *tiles[0];
	};
	const Tile& Border::oppositeTile(const Tile& tile) const
	{
		return (*tiles[0] == tile) ? *tiles[1] : *tiles[0];
	};
	double Border::length() const
	{
		return corners[0]->position.distance(corners[1]->position);
	};
	bool Border::isLandBoundary() const
	{
		return (tiles[0]->elevation > 0) != (tiles[1]->elevation > 0);
	};
	std::string Border::toString() const
	{
		std::stringstream ss;
		ss << "Border " << id;
		return ss.str();
	};







	Ogre::Vector3 Plate::calculateMovement(const Ogre::Vector3& position)
	{
		auto movement = setLength(driftAxis.crossProduct(position), driftRate * projectOnVector(position, driftAxis).distance(position));
		auto tmp = setLength(root->position.crossProduct(position), spinRate * projectOnVector(position, root->position).distance(position));
		movement += tmp;
		return movement;
	};




	bool SpatialPartition::intersectRay(const Ogre::Ray& ray)
	{
		if (intersectRayWithSphere(ray, boundingSphere))
		{
			for (auto& partition : partitions)
			{
				auto intersection = partition.intersectRay(ray);
				if (intersection != false)
				{
					return true;
				}
			}

			for (auto& tile : tiles)
			{
				if (tile->intersectRay(ray))
				{
					return true;
				}
			}
		}

		return false;
	};













	enum eValues
	{
		subdivisions = 0,
		distortionLevel,
		plateCount,
		oceanicRate,
		heatLevel,
		moistureLevel,
		seed,
	};



World::World(MyGUI::Gui* gui, MyGUI::OgrePlatform* platform, Ogre::RenderWindow* renderWindow, Ogre::SceneManager* sceneMgr)
	: mGUI(gui)
	, mPlatform(platform)
	, mRenderWindow(renderWindow)
	, sceneMgr_(sceneMgr)
{
}

World::~World()
{
}



void World::generatePlanetAsynchronous(Planet*& planet, unsigned long originalSeed, unsigned long seed, size_t icosahedronSubdivision, double topologyDistortionRate, size_t plateCount, double oceanicRate, double heatLevel, double moistureLevel)
{
	XorShift128 random(seed);

	generatePlanet(icosahedronSubdivision, topologyDistortionRate, plateCount, oceanicRate, heatLevel, moistureLevel, random, *planet);
	planet->seed = seed;
	planet->originalSeed = originalSeed;
}

void World::generatePlanet(size_t icosahedronSubdivision, double topologyDistortionRate, size_t plateCount, double oceanicRate, double heatLevel, double moistureLevel, XorShift128& random, Planet& planet)
{
	auto mesh = generatePlanetMesh(icosahedronSubdivision, topologyDistortionRate, random);

	generatePlanetTopology(mesh, planet.topology);

	generatePlanetPartition(planet.topology.tiles, planet.partition);
	generatePlanetTerrain(planet, plateCount, oceanicRate, heatLevel, moistureLevel, random);
	generatePlanetRenderData(planet.topology, random, planet.renderData);
	generatePlanetStatistics(planet.topology, planet.plates, planet.statistics);
}

Mesh World::generatePlanetMesh(size_t icosahedronSubdivision, double topologyDistortionRate, XorShift128& random)
{
	Mesh mesh;
	generateSubdividedIcosahedron(icosahedronSubdivision, mesh);


	// Distorting Triangle Mesh
	double totalDistortion = std::ceil(mesh.edges.size() * topologyDistortionRate);
	for (size_t remainingIterations = 6; remainingIterations > 0;)
	{
		double iterationDistortion = std::floor(totalDistortion / remainingIterations);
		totalDistortion -= iterationDistortion;
		distortMesh(mesh, iterationDistortion, random);
		relaxMesh(mesh, 0.5);
		--remainingIterations;
	}

	
	// Relaxing Triangle Mesh
	double averageNodeRadius = std::sqrt(4 * M_PI / mesh.nodes.size());
	double minShiftDelta = averageNodeRadius / 50000 * mesh.nodes.size();
	double maxShiftDelta = averageNodeRadius / 50 * mesh.nodes.size();

	double priorShift = relaxMesh(mesh, 0.5);
	for (int i = 0; i < 300; i++)
	{
		double currentShift = relaxMesh(mesh, 0.5);
		double shiftDelta = std::abs(currentShift - priorShift);
		if (shiftDelta < minShiftDelta)
			break;
		priorShift = currentShift;
	}


	// Calculating Triangle Centroids
	for (auto pface : mesh.faces)
	{
		auto& face = *pface;
		auto& p0 = mesh.nodes[face.n[0]]->p;
		auto& p1 = mesh.nodes[face.n[1]]->p;
		auto& p2 = mesh.nodes[face.n[2]]->p;
		face.centroid = calculateFaceCentroid(p0, p1, p2);
		face.centroid.normalise();
	}


	// Reordering Triangle Nodes
	int i = 0;
	for (auto pnode : mesh.nodes)
	{
		auto& node = *pnode;
		size_t faceIndex = node.f[0];
		for (int j = 1; j < node.f.size() - 1; ++j)
		{
			faceIndex = findNextFaceIndex(mesh, i, faceIndex);
			size_t k;
			bool ok = findIndex(node.f, faceIndex, k);
			assert(ok);
			node.f[k] = node.f[j];
			node.f[j] = faceIndex;
		}
		++i;
	}

	return mesh;
}

void World::generateIcosahedron(Mesh& ret)
{
	double phi = (1 + std::sqrt(5.)) / 2;
	double du = 1 / std::sqrt(phi * phi + 1.0);
	double dv = phi * du;

	ret.nodes =
	{
		new Node{ Ogre::Vector3(0, +dv, +du) },
		new Node{ Ogre::Vector3(0, +dv, -du) },
		new Node{ Ogre::Vector3(0, -dv, +du) },
		new Node{ Ogre::Vector3(0, -dv, -du) },
		new Node{ Ogre::Vector3(+du, 0, +dv) },
		new Node{ Ogre::Vector3(-du, 0, +dv) },
		new Node{ Ogre::Vector3(+du, 0, -dv) },
		new Node{ Ogre::Vector3(-du, 0, -dv) },
		new Node{ Ogre::Vector3(+dv, +du, 0) },
		new Node{ Ogre::Vector3(+dv, -du, 0) },
		new Node{ Ogre::Vector3(-dv, +du, 0) },
		new Node{ Ogre::Vector3(-dv, -du, 0) },
	};

	ret.edges =
	{
		new Edge{ { 0, 1 } },
		new Edge{ { 0, 4 } },
		new Edge{ { 0, 5 } },
		new Edge{ { 0, 8 } },
		new Edge{ { 0, 10 } },
		new Edge{ { 1, 6 } },
		new Edge{ { 1, 7 } },
		new Edge{ { 1, 8 } },
		new Edge{ { 1, 10 } },
		new Edge{ { 2, 3 } },
		new Edge{ { 2, 4 } },
		new Edge{ { 2, 5 } },
		new Edge{ { 2, 9 } },
		new Edge{ { 2, 11 } },
		new Edge{ { 3, 6 } },
		new Edge{ { 3, 7 } },
		new Edge{ { 3, 9 } },
		new Edge{ { 3, 11 } },
		new Edge{ { 4, 5 } },
		new Edge{ { 4, 8 } },
		new Edge{ { 4, 9 } },
		new Edge{ { 5, 10 } },
		new Edge{ { 5, 11 } },
		new Edge{ { 6, 7 } },
		new Edge{ { 6, 8 } },
		new Edge{ { 6, 9 } },
		new Edge{ { 7, 10 } },
		new Edge{ { 7, 11 } },
		new Edge{ { 8, 9 } },
		new Edge{ { 10, 11 } },
	};

	ret.faces =
	{
		new Face{ { 0, 1, 8 }, { 0, 7, 3 }, },
		new Face{ { 0, 4, 5 }, { 1, 18, 2 }, },
		new Face{ { 0, 5, 10 }, { 2, 21, 4 }, },
		new Face{ { 0, 8, 4 }, { 3, 19, 1 }, },
		new Face{ { 0, 10, 1 }, { 4, 8, 0 }, },
		new Face{ { 1, 6, 8 }, { 5, 24, 7 }, },
		new Face{ { 1, 7, 6 }, { 6, 23, 5 }, },
		new Face{ { 1, 10, 7 }, { 8, 26, 6 }, },
		new Face{ { 2, 3, 11 }, { 9, 17, 13 }, },
		new Face{ { 2, 4, 9 }, { 10, 20, 12 }, },
		new Face{ { 2, 5, 4 }, { 11, 18, 10 }, },
		new Face{ { 2, 9, 3 }, { 12, 16, 9 }, },
		new Face{ { 2, 11, 5 }, { 13, 22, 11 }, },
		new Face{ { 3, 6, 7 }, { 14, 23, 15 }, },
		new Face{ { 3, 7, 11 }, { 15, 27, 17 }, },
		new Face{ { 3, 9, 6 }, { 16, 25, 14 }, },
		new Face{ { 4, 8, 9 }, { 19, 28, 20 }, },
		new Face{ { 5, 11, 10 }, { 22, 29, 21 }, },
		new Face{ { 6, 9, 8 }, { 25, 28, 24 }, },
		new Face{ { 7, 10, 11 }, { 26, 29, 27 }, },
	};

	ret.nodes.reserve(ret.edges.size() * /* 2do: guess */ 3);
	for (size_t i = 0; i < ret.edges.size(); ++i)
		for (size_t j = 0; j < ret.edges[i]->n.size(); ++j)
			ret.nodes[j]->e.push_back(i);

	for (size_t i = 0; i < ret.faces.size(); ++i)
		for (size_t j = 0; j < ret.faces[i]->n.size(); ++j)
			ret.nodes[j]->f.push_back(i);

	ret.edges.reserve(ret.faces.size() * /* 2do: guess */ 3);
	for (size_t i = 0; i < ret.faces.size(); ++i)
		for (size_t j = 0; j < ret.faces[i]->e.size(); ++j)
			ret.edges[j]->f.push_back(i);
}


void World::generateSubdividedIcosahedron(size_t degree, Mesh& ret)
{
	auto& nodes = ret.nodes;
	auto& edges = ret.edges;
	auto& faces = ret.faces;

	Mesh icosahedron;
	generateIcosahedron(icosahedron);

	nodes.reserve(icosahedron.nodes.size() + icosahedron.edges.size() * (degree - 1) + icosahedron.faces.size() * (degree - 1)*(degree - 1) / 2);
	for (auto& node : icosahedron.nodes)
	{
		nodes.push_back(new Node{ node->p });
	}


	edges.reserve(icosahedron.edges.size() + icosahedron.edges.size()*(degree - 1) + 3 * icosahedron.faces.size() * (degree - 1)*(degree - 1)/2);

	for (auto& pedge : icosahedron.edges)
	{
		auto& edge = *pedge;

		edge.subdivided_n.clear();
		edge.subdivided_e.clear();
		auto& n0 = *icosahedron.nodes[edge.n[0]];
		auto& n1 = *icosahedron.nodes[edge.n[1]];
		auto& p0 = n0.p;
		auto& p1 = n1.p;
		auto delta = p1 - p0;
		nodes[edge.n[0]]->e.push_back(edges.size());
		size_t priorNodeIndex = edge.n[0];
		for (auto s = 1; s < degree; ++s)
		{
			size_t edgeIndex = edges.size();
			size_t nodeIndex = nodes.size();
			edge.subdivided_e.push_back(edgeIndex);
			edge.subdivided_n.push_back(nodeIndex);
			edges.push_back(new Edge{ { priorNodeIndex, nodeIndex } });
			priorNodeIndex = nodeIndex;
			nodes.push_back(new Node{ slerp(p0, p1, (double)s / degree), { edgeIndex, edgeIndex + 1 } });
		}
		edge.subdivided_e.push_back(edges.size());
		nodes[edge.n[1]]->e.push_back(edges.size());
		edges.push_back(new Edge{ { priorNodeIndex, edge.n[1] } });
	}

	faces.reserve(2 * icosahedron.faces.size() * (degree - 1)*(degree - 1) / 2);

	for (auto& pface : icosahedron.faces)
	{
		auto& face = *pface;

		auto& edge0 = *icosahedron.edges[face.e[0]];
		auto& edge1 = *icosahedron.edges[face.e[1]];
		auto& edge2 = *icosahedron.edges[face.e[2]];
		auto& point0 = icosahedron.nodes[face.n[0]]->p;
		auto& point1 = icosahedron.nodes[face.n[1]]->p;
		auto& point2 = icosahedron.nodes[face.n[2]]->p;
		auto delta = point1 - point0;

		typedef std::function<size_t(size_t)> GetIdxFun;
		auto getEdgeNode0 = (face.n[0] == edge0.n[0]
			? GetIdxFun([&](size_t k) -> size_t { return edge0.subdivided_n[k]; })
			: GetIdxFun([&](size_t k) -> size_t { return edge0.subdivided_n[degree - 2 - k]; }));
		auto getEdgeNode1 = (face.n[1] == edge1.n[0]
			? GetIdxFun([&](size_t k) -> size_t { return edge1.subdivided_n[k]; })
			: GetIdxFun([&](size_t k) -> size_t { return edge1.subdivided_n[degree - 2 - k]; }));
		auto getEdgeNode2 = (face.n[0] == edge2.n[0]
			? GetIdxFun([&](size_t k) -> size_t { return edge2.subdivided_n[k]; })
			: GetIdxFun([&](size_t k) -> size_t { return edge2.subdivided_n[degree - 2 - k]; }));

		std::vector<size_t> faceNodes; // (1 + edge0.subdivided_n.size() + 1 + (degree - 1)*(degree - 1) / 2 + 1);
		faceNodes.reserve(1 + edge0.subdivided_n.size() + 1 + (degree - 1)*(degree - 1) / 2 + 1);
		faceNodes.push_back(face.n[0]);
		for (auto j = 0; j < edge0.subdivided_n.size(); ++j)
			faceNodes.push_back(getEdgeNode0(j));
		faceNodes.push_back(face.n[1]);
		for (auto s = 1; s < degree; ++s)
		{
			faceNodes.push_back(getEdgeNode2(s - 1));
			auto& p0 = nodes[getEdgeNode2(s - 1)]->p;
			auto& p1 = nodes[getEdgeNode1(s - 1)]->p;
			for (auto t = 1; t < degree - s; ++t)
			{
				faceNodes.push_back(nodes.size());
				nodes.push_back(new Node{ slerp(p0, p1, (double)t / (degree - s)) });
			}
			faceNodes.push_back(getEdgeNode1(s - 1));
		}
		faceNodes.push_back(face.n[2]);

		auto getEdgeEdge0 = (face.n[0] == edge0.n[0])
			? GetIdxFun([&](size_t k) -> size_t { return edge0.subdivided_e[k]; })
			: GetIdxFun([&](size_t k) -> size_t { return edge0.subdivided_e[degree - 1 - k]; });
		auto getEdgeEdge1 = (face.n[1] == edge1.n[0])
			? GetIdxFun([&](size_t k) -> size_t { return edge1.subdivided_e[k]; })
			: GetIdxFun([&](size_t k) -> size_t { return edge1.subdivided_e[degree - 1 - k]; });
		auto getEdgeEdge2 = (face.n[0] == edge2.n[0])
			? GetIdxFun([&](size_t k) -> size_t { return edge2.subdivided_e[k]; })
			: GetIdxFun([&](size_t k) -> size_t { return edge2.subdivided_e[degree - 1 - k]; });

		std::vector<size_t> faceEdges0; // (degree + (degree - 1)*(degree - 1) / 2);
		faceEdges0.reserve(degree + (degree - 1)*(degree - 1) / 2);
		for (auto j = 0; j < degree; ++j)
			faceEdges0.push_back(getEdgeEdge0(j));
		auto nodeIndex = degree + 1;
		for (auto s = 1; s < degree; ++s)
		{
			for (auto t = 0; t < degree - s; ++t)
			{
				faceEdges0.push_back(edges.size());
				auto pedge = new Edge{ { faceNodes[nodeIndex], faceNodes[nodeIndex + 1] } };
				auto& edge = *pedge;
				nodes[edge.n[0]]->e.push_back(edges.size());
				nodes[edge.n[1]]->e.push_back(edges.size());
				edges.push_back(pedge);
				++nodeIndex;
			}
			++nodeIndex;
		}

		std::vector<size_t> faceEdges1; // (degree * (degree - 1) / 2);
		faceEdges1.reserve(degree * (degree - 1) / 2);
		nodeIndex = 1;
		for (auto s = 0; s < degree; ++s)
		{
			for (auto t = 1; t < degree - s; ++t)
			{
				faceEdges1.push_back(edges.size());
				auto pedge = new Edge{ { faceNodes[nodeIndex], faceNodes[nodeIndex + degree - s] } };
				auto& edge = *pedge;
				nodes[edge.n[0]]->e.push_back(edges.size());
				nodes[edge.n[1]]->e.push_back(edges.size());
				edges.push_back(pedge);
				++nodeIndex;
			}
			faceEdges1.push_back(getEdgeEdge1(s));
			nodeIndex += 2;
		}

		std::vector<size_t> faceEdges2; // (degree * (degree - 1) / 2);
		faceEdges2.reserve(degree * (degree - 1) / 2);
		nodeIndex = 1;
		for (auto s = 0; s < degree; ++s)
		{
			faceEdges2.push_back(getEdgeEdge2(s));
			for (auto t = 1; t < degree - s; ++t)
			{
				faceEdges2.push_back(edges.size());
				auto pedge = new Edge{ { faceNodes[nodeIndex], faceNodes[nodeIndex + degree - s + 1] } };
				auto& edge = *pedge;
				nodes[edge.n[0]]->e.push_back(edges.size());
				nodes[edge.n[1]]->e.push_back(edges.size());
				edges.push_back(pedge);
				++nodeIndex;
			}
			nodeIndex += 2;
		}

		nodeIndex = 0;
		auto edgeIndex = 0;
		for (auto s = 0; s < degree; ++s)
		{
			for (auto t = 1; t < degree - s + 1; ++t)
			{
				auto psubFace = new Face{
					{ faceNodes[nodeIndex], faceNodes[nodeIndex + 1], faceNodes[nodeIndex + degree - s + 1] },
					{ faceEdges0[edgeIndex], faceEdges1[edgeIndex], faceEdges2[edgeIndex] } };
				auto& subFace = *psubFace;
				nodes[subFace.n[0]]->f.push_back(faces.size());
				nodes[subFace.n[1]]->f.push_back(faces.size());
				nodes[subFace.n[2]]->f.push_back(faces.size());
				edges[subFace.e[0]]->f.push_back(faces.size());
				edges[subFace.e[1]]->f.push_back(faces.size());
				edges[subFace.e[2]]->f.push_back(faces.size());
				faces.push_back(psubFace);
				++nodeIndex;
				++edgeIndex;
			}
			++nodeIndex;
		}

		nodeIndex = 1;
		edgeIndex = 0;
		for (auto s = 1; s < degree; ++s)
		{
			for (auto t = 1; t < degree - s + 1; ++t)
			{
				auto psubFace = new Face{
						{ faceNodes[nodeIndex], faceNodes[nodeIndex + degree - s + 2], faceNodes[nodeIndex + degree - s + 1], },
						{ faceEdges2[edgeIndex + 1], faceEdges0[edgeIndex + degree - s + 1], faceEdges1[edgeIndex], }, };
				auto& subFace = *psubFace;
				nodes[subFace.n[0]]->f.push_back(faces.size());
				nodes[subFace.n[1]]->f.push_back(faces.size());
				nodes[subFace.n[2]]->f.push_back(faces.size());
				edges[subFace.e[0]]->f.push_back(faces.size());
				edges[subFace.e[1]]->f.push_back(faces.size());
				edges[subFace.e[2]]->f.push_back(faces.size());
				faces.push_back(psubFace);
				++nodeIndex;
				++edgeIndex;
			}
			nodeIndex += 2;
			edgeIndex += 1;
		}
	}
}






size_t World::getEdgeOppositeFaceIndex(const Edge& edge, size_t faceIndex)
{
	if (edge.f[0] == faceIndex) return edge.f[1];
	if (edge.f[1] == faceIndex) return edge.f[0];
	throw "Given face is not part of given edge.";
}

size_t World::getFaceOppositeNodeIndex(const Face& face, const Edge& edge)
{
	if (face.n[0] != edge.n[0] && face.n[0] != edge.n[1]) return 0;
	if (face.n[1] != edge.n[0] && face.n[1] != edge.n[1]) return 1;
	if (face.n[2] != edge.n[0] && face.n[2] != edge.n[1]) return 2;
	throw std::exception("Cannot find node of given face that is not also a node of given edge.");
}

size_t World::findNextFaceIndex(const Mesh& mesh, size_t nodeIndex, size_t faceIndex)
{
	auto node = mesh.nodes[nodeIndex];
	auto& face = *mesh.faces[faceIndex];
	size_t nodeFaceIndex;
	bool ok = findIndex(face.n, nodeIndex, nodeFaceIndex);
	assert(ok);
	auto& edge = *mesh.edges[face.e[(nodeFaceIndex + 2) % 3]];
	return getEdgeOppositeFaceIndex(edge, faceIndex);
}

bool World::conditionalRotateEdge(const Mesh& mesh, size_t edgeIndex, RotationPredicateType& predicate)
{
	auto& edge = *mesh.edges[edgeIndex];
	auto& face0 = *mesh.faces[edge.f[0]];
	auto& face1 = *mesh.faces[edge.f[1]];
	size_t farNodeFaceIndex0;
	size_t farNodeFaceIndex1;
	farNodeFaceIndex0 = getFaceOppositeNodeIndex(face0, edge);
	farNodeFaceIndex1 = getFaceOppositeNodeIndex(face1, edge);
	auto newNodeIndex0 = face0.n[farNodeFaceIndex0];
	auto oldNodeIndex0 = face0.n[(farNodeFaceIndex0 + 1) % 3];
	auto newNodeIndex1 = face1.n[farNodeFaceIndex1];
	auto oldNodeIndex1 = face1.n[(farNodeFaceIndex1 + 1) % 3];
	auto& oldNode0 = *mesh.nodes[oldNodeIndex0];
	auto& oldNode1 = *mesh.nodes[oldNodeIndex1];
	auto& newNode0 = *mesh.nodes[newNodeIndex0];
	auto& newNode1 = *mesh.nodes[newNodeIndex1];
	auto newEdgeIndex0 = face1.e[(farNodeFaceIndex1 + 2) % 3];
	auto newEdgeIndex1 = face0.e[(farNodeFaceIndex0 + 2) % 3];
	auto& newEdge0 = *mesh.edges[newEdgeIndex0];
	auto& newEdge1 = *mesh.edges[newEdgeIndex1];

	if (!predicate(oldNode0, oldNode1, newNode0, newNode1)) return false;

	removeIfFindIndex(oldNode0.e, edgeIndex);
	removeIfFindIndex(oldNode1.e, edgeIndex);
	newNode0.e.push_back(edgeIndex);
	newNode1.e.push_back(edgeIndex);

	edge.n[0] = newNodeIndex0;
	edge.n[1] = newNodeIndex1;

	removeIfFindIndex(newEdge0.f, edge.f[1]);
	removeIfFindIndex(newEdge1.f, edge.f[0]);
	newEdge0.f.push_back(edge.f[0]);
	newEdge1.f.push_back(edge.f[1]);

	removeIfFindIndex(oldNode0.f, edge.f[1]);
	removeIfFindIndex(oldNode1.f, edge.f[0]);
	newNode0.f.push_back(edge.f[1]);
	newNode1.f.push_back(edge.f[0]);

	face0.n[(farNodeFaceIndex0 + 2) % 3] = newNodeIndex1;
	face1.n[(farNodeFaceIndex1 + 2) % 3] = newNodeIndex0;

	face0.e[(farNodeFaceIndex0 + 1) % 3] = newEdgeIndex0;
	face1.e[(farNodeFaceIndex1 + 1) % 3] = newEdgeIndex1;
	face0.e[(farNodeFaceIndex0 + 2) % 3] = edgeIndex;
	face1.e[(farNodeFaceIndex1 + 2) % 3] = edgeIndex;

	return true;
}

Ogre::Vector3 World::calculateFaceCentroid(const Ogre::Vector3& pa, const Ogre::Vector3& pb, const Ogre::Vector3& pc)
{
	auto centroid = pa + pb + pc;
	centroid /= 3;
	return centroid;
}

bool World::distortMesh(const Mesh& mesh, size_t degree, XorShift128& random)
{
	auto totalSurfaceArea = 4. * M_PI;
	auto idealFaceArea = totalSurfaceArea / (double)mesh.faces.size();
	auto idealEdgeLength = std::sqrt(idealFaceArea * 4. / std::sqrt(3.));
	auto idealFaceHeight = idealEdgeLength * std::sqrt(3.) / 2.;

	auto rotationPredicate = RotationPredicateType([&](const Node& oldNode0, const Node& oldNode1, const Node& newNode0, const Node& newNode1) -> bool
	{
		if (newNode0.f.size() >= 7 ||
			newNode1.f.size() >= 7 ||
			oldNode0.f.size() <= 5 ||
			oldNode1.f.size() <= 5) return false;
		double oldEdgeLength = oldNode0.p.distance(oldNode1.p);
		if (oldEdgeLength == 0)
			return false;
		double newEdgeLength = newNode0.p.distance(newNode1.p);
		auto ratio = oldEdgeLength / newEdgeLength;
		if (ratio >= 2 || ratio <= 0.5) return false;
		auto v0 = (oldNode1.p - oldNode0.p) / oldEdgeLength;
		auto v1 = (newNode0.p - oldNode0.p); v1.normalise();
		auto v2 = (newNode1.p - oldNode0.p); v2.normalise();
		if (v0.dotProduct(v1) < 0.2 || v0.dotProduct(v2) < 0.2) return false;
		v0 *= -1.;
		auto v3 = (newNode0.p - oldNode1.p); v3.normalise();
		auto v4 = (newNode1.p - oldNode1.p); v4.normalise();
		if (v0.dotProduct(v3) < 0.2 || v0.dotProduct(v4) < 0.2) return false;
		return true;
	});

	for (auto i = 0; i < degree; i++)
	{
		auto consecutiveFailedAttempts = 0;
		auto edgeIndex = random.integerExclusive(0, mesh.edges.size());
		while (!conditionalRotateEdge(mesh, edgeIndex, rotationPredicate))
		{
			if (++consecutiveFailedAttempts >= mesh.edges.size()) return false;
			edgeIndex = (edgeIndex + 1) % mesh.edges.size();
		}
	};

	return true;
}

double World::relaxMesh(const Mesh& mesh, double multiplier)
{
	double  totalSurfaceArea = 4. * M_PI;
	double  idealFaceArea = totalSurfaceArea / mesh.faces.size();
	double  idealEdgeLength = std::sqrt(idealFaceArea * 4. / std::sqrt(3.));
	double  idealDistanceToCentroid = idealEdgeLength * std::sqrt(3) / 3. * 0.9;

	std::vector<Ogre::Vector3> pointShifts(mesh.nodes.size());
	for (auto i = 0; i < mesh.nodes.size(); ++i)
		pointShifts[i] = Ogre::Vector3(0, 0, 0);

	auto i = 0;
	for (auto& pface : mesh.faces)
	{
		auto& face = *pface;

		auto& n0 = *mesh.nodes[face.n[0]];
		auto& n1 = *mesh.nodes[face.n[1]];
		auto& n2 = *mesh.nodes[face.n[2]];
		auto& p0 = n0.p;
		auto& p1 = n1.p;
		auto& p2 = n2.p;
		auto e0 = p1.distance(p0) / idealEdgeLength;
		auto e1 = p2.distance(p1) / idealEdgeLength;
		auto e2 = p0.distance(p2) / idealEdgeLength;
		auto centroid = calculateFaceCentroid(p0, p1, p2); centroid.normalise();
		auto v0 = centroid - p0;
		auto v1 = centroid - p1;
		auto v2 = centroid - p2;
		auto length0 = v0.length();
		auto length1 = v1.length();
		auto length2 = v2.length();
		v0 *= (multiplier * (length0 - idealDistanceToCentroid) / length0);
		v1 *= (multiplier * (length1 - idealDistanceToCentroid) / length1);
		v2 *= (multiplier * (length2 - idealDistanceToCentroid) / length2);
		pointShifts[face.n[0]] += v0; // does this change the tmp ret or the ref?
		pointShifts[face.n[1]] += v1;
		pointShifts[face.n[2]] += v2;

		++i;
	}

	auto origin = Ogre::Vector3(0, 0, 0);
	Ogre::Plane plane;
	for (auto i = 0; i < mesh.nodes.size(); ++i)
	{
		plane.redefine(mesh.nodes[i]->p, origin); //.setFromNormalAndCoplanarPoint(mesh.nodes[i].p, origin);
		pointShifts[i] = mesh.nodes[i]->p + plane.projectVector(pointShifts[i]); // .projectPoint(pointShifts[i]); 
		pointShifts[i].normalise();
	}

	std::vector<double> rotationSupressions(mesh.nodes.size());
	for (auto i = 0; i < mesh.nodes.size(); ++i)
		rotationSupressions[i] = 0;

	for (auto pedge : mesh.edges)
	{
		auto& edge = *pedge;
		auto& oldPoint0 = mesh.nodes[edge.n[0]]->p;
		auto& oldPoint1 = mesh.nodes[edge.n[1]]->p;
		auto& newPoint0 = pointShifts[edge.n[0]];
		auto& newPoint1 = pointShifts[edge.n[1]];
		auto oldVector = oldPoint1 - oldPoint0; oldVector.normalise();
		auto newVector = newPoint1 - newPoint0; newVector.normalise();
		auto suppression = (1 - oldVector.dotProduct(newVector)) * 0.5;
		rotationSupressions[edge.n[0]] = std::max(rotationSupressions[edge.n[0]], suppression);
		rotationSupressions[edge.n[1]] = std::max(rotationSupressions[edge.n[1]], suppression);
	}

	double totalShift = 0;
	i = 0;
	for (auto pnode : mesh.nodes)
	{
		auto& node = *pnode;

		auto& point = node.p;
		auto delta = point;
		point = lerp(point, pointShifts[i], 1. - std::sqrt(rotationSupressions[i]));
		point.normalise();
		delta -= point;
		totalShift += delta.length();

		++i;
	}

	return totalShift;
}

bool World::generatePlanetTopology(const Mesh& mesh, Topology& ret)
{
	auto& corners = ret.corners;
	auto& borders = ret.borders;
	auto& tiles = ret.tiles;
	size_t i;

	i = 0;
	corners.reserve(mesh.faces.size());
	for (auto& pface : mesh.faces)
	{
		auto& face = *pface;
		corners.push_back(Corner(i, face.centroid * 1000, face.e.size(), face.e.size(), face.n.size()));
		++i;
	}

	i = 0;
	borders.reserve(mesh.edges.size());
	for (auto pedge : mesh.edges)
	{
		auto& edge = *pedge;
		borders.push_back(Border(i, 2, 4, 2)); //edge.f.size(), mesh.faces[edge.f[0]].e.length + mesh.faces[edge.f[1]].e.length - 2, edge.n.length
		++i;
	}

	i = 0;
	tiles.reserve(mesh.nodes.size());
	for (auto pnode : mesh.nodes)
	{
		auto& node = *pnode;
		tiles.push_back(Tile(i, node.p * 1000, node.f.size(), node.e.size(), node.e.size()));
		++i;
	}

	i = 0;
	for (auto& corner : corners)
	{
		auto& face = *mesh.faces[i];
		for (auto j = 0; j < face.e.size(); ++j)
		{
			corner.borders[j] = &borders[face.e[j]]; // 2do: assign ptr or val? and if ptr, no realloc may happen in any std::vector !!! how to assure that?
		}
		for (auto j = 0; j < face.n.size(); ++j)
		{
			corner.tiles[j] = &tiles[face.n[j]];
		}
		++i;
	}

	i = 0;
	for (auto& border : borders)
	{
		auto& edge = *mesh.edges[i];
		auto averageCorner = Ogre::Vector3(0, 0, 0);
		auto n = 0;
		//border.corners.reserve(edge.f.size());
		//border.corners.clear();
		// 2do: how2 reserve border.borders
		for (auto j = 0; j < edge.f.size(); ++j)
		{
			auto& corner = corners[edge.f[j]];
			averageCorner += corner.position;
			border.corners[j] = &corner;
			for (auto pcornerBorder : corner.borders)
			{
				auto& cornerBorder = *pcornerBorder;
				if (pcornerBorder != &border) border.borders[n++] = pcornerBorder;
			}
		}
		border.midpoint = averageCorner * (1. / border.corners.size());
		for (auto j = 0; j < edge.n.size(); ++j)
		{
			border.tiles[j] = &tiles[edge.n[j]];
		}
		++i;
	}

	for (auto& corner : corners)
	{
		for (auto j = 0; j < corner.borders.size(); ++j)
		{
			corner.corners[j] = &corner.borders[j]->oppositeCorner(corner);
		}
	}

	i = 0;
	for (auto& tile : tiles)
	{
		auto& node = *mesh.nodes[i];
		for (auto j = 0; j < node.f.size(); ++j)
		{
			tile.corners[j] = &corners[node.f[j]];
		}
		for (auto j = 0; j < node.e.size(); ++j)
		{
			auto& border = borders[node.e[j]];
			if (border.tiles[0] == &tile)
			{
				for (auto k = 0; k < tile.corners.size(); ++k)
				{
					auto& corner0 = *tile.corners[k];
					auto& corner1 = *tile.corners[(k + 1) % tile.corners.size()];
					if (border.corners[1] == &corner0 && border.corners[0] == &corner1)
					{
						border.corners[0] = &corner0;
						border.corners[1] = &corner1;
					}
					else if (border.corners[0] != &corner0 || border.corners[1] != &corner1)
					{
						continue;
					}
					tile.borders[k] = &border;
					tile.tiles[k] = &border.oppositeTile(tile);
					break;
				}
			}
			else
			{
				for (auto k = 0; k < tile.corners.size(); ++k)
				{
					auto& corner0 = *tile.corners[k];
					auto& corner1 = *tile.corners[(k + 1) % tile.corners.size()];
					if (border.corners[0] == &corner0 && border.corners[1] == &corner1)
					{
						border.corners[1] = &corner0;
						border.corners[0] = &corner1;
					}
					else if (border.corners[1] != &corner0 || border.corners[0] != &corner1)
					{
						continue;
					}
					tile.borders[k] = &border;
					tile.tiles[k] = &border.oppositeTile(tile);
					break;
				}
			}
		}

		tile.averagePosition = Ogre::Vector3(0, 0, 0);
		for (auto j = 0; j < tile.corners.size(); ++j)
		{
			tile.averagePosition += tile.corners[j]->position;
		}
		tile.averagePosition *= (1. / tile.corners.size());

		double maxDistanceToCorner = 0;
		for (auto j = 0; j < tile.corners.size(); ++j)
		{
			maxDistanceToCorner = std::max(maxDistanceToCorner, (double)tile.corners[j]->position.distance(tile.averagePosition));
		}

		double area = 0;
		for (auto j = 0; j < tile.borders.size(); ++j)
		{
			area += calculateTriangleArea(tile.position, tile.borders[j]->corners[0]->position, tile.borders[j]->corners[1]->position);
		}
		tile.area = area;

		tile.normal = tile.position;
		tile.normal.normalise();

		tile.boundingSphere = Ogre::Sphere(tile.averagePosition, maxDistanceToCorner);
		++i;
	}

	for (auto& corner : corners)
	{
		corner.area = 0;
		for (auto j = 0; j < corner.tiles.size(); ++j)
		{
			corner.area += corner.tiles[j]->area / corner.tiles[j]->corners.size();
		}
	}

	return true;
}

void World::generatePlanetPartition(std::vector<Tile>& tiles, SpatialPartition& rootPartition)
{
	Mesh icosahedron;
	generateIcosahedron(icosahedron);


	for (auto pface : icosahedron.faces)
	{
		auto& face = *pface;

		auto p0 = icosahedron.nodes[face.n[0]]->p * 1000;
		auto p1 = icosahedron.nodes[face.n[1]]->p * 1000;
		auto p2 = icosahedron.nodes[face.n[2]]->p * 1000;
		auto center = (p0 + p1 + p2) / 3.;
		auto radius = std::max(center.distance(p0), std::max(center.distance(p2), center.distance(p2)));
		face.boundingSphere = Ogre::Sphere(center, radius);
		//face.children = [];
	}


	std::vector<Tile*> unparentedTiles;
	double maxDistanceFromOrigin = 0;
	for (auto& tile : tiles)
	{
		maxDistanceFromOrigin = std::max(maxDistanceFromOrigin, (double)(tile.boundingSphere.getCenter().length() + tile.boundingSphere.getRadius()));

		bool parentFound = false;
		for (auto pface : icosahedron.faces)
		{
			auto& face = *pface;
			double distance = tile.boundingSphere.getCenter().distance(face.boundingSphere.getCenter()) + tile.boundingSphere.getRadius();
			if (distance < face.boundingSphere.getRadius())
			{
				face.children.push_back(&tile);
				parentFound = true;
				break;
			}
		}
		if (!parentFound)
		{
			unparentedTiles.push_back(&tile);
		}
	}

	rootPartition.boundingSphere = Ogre::Sphere(Ogre::Vector3(0, 0, 0), maxDistanceFromOrigin);
	rootPartition.tiles = unparentedTiles;
	for (auto pface : icosahedron.faces)
	{
		auto& face = *pface;
		rootPartition.partitions.push_back(SpatialPartition(face.boundingSphere, {}, face.children));
		// 2do
		//delete face.boundingSphere;
		//delete face.children;
	}
}

void World::generatePlanetTerrain(Planet& planet, size_t plateCount, double oceanicRate, double heatLevel, double moistureLevel, XorShift128& random)
{
	generatePlanetTectonicPlates(planet.topology, plateCount, oceanicRate, random, planet.plates);
	generatePlanetElevation(planet.topology, planet.plates);
	generatePlanetWeather(planet.topology, planet.partition, heatLevel, moistureLevel, random);
	generatePlanetBiomes(planet.topology.tiles, 1000);
}

void World::generatePlanetTectonicPlates(Topology& topology, size_t plateCount, double oceanicRate, XorShift128& random, std::vector<Plate>& plates)
{
	struct Plateless
	{
		Tile* tile;
		Plate* plate;
	};
	std::list<Plateless> plateless;

	plates.reserve(plateCount); // 2do: hack: vector may not grow anymore or all pointers r invalid; faster though :/

	auto failedCount = 0;
	while (plates.size() < plateCount && failedCount < 10000)
	{
		auto& corner = topology.corners[random.integerExclusive(0, topology.corners.size())];
		bool adjacentToExistingPlate = false;
		for (auto i = 0; i < corner.tiles.size(); ++i)
		{
			if (corner.tiles[i]->plate)
			{
				adjacentToExistingPlate = true;
				failedCount += 1;
				break;
			}
		}
		if (adjacentToExistingPlate) continue;

		failedCount = 0;

		bool oceanic = (random.unit() < oceanicRate);
		plates.push_back(Plate(
			Ogre::ColourValue(random.realInclusive(0, 1), random.realInclusive(0, 1), random.realInclusive(0, 1)),
			randomUnitVector(random),
			random.realInclusive(-M_PI / 30, M_PI / 30),
			random.realInclusive(-M_PI / 30, M_PI / 30),
			oceanic ? random.realInclusive(-0.8, -0.3) : random.realInclusive(0.1, 0.5),
			oceanic,
			corner));
		auto pplate = &plates.back();
		auto& plate = plates.back();

		for (auto i = 0; i < corner.tiles.size(); ++i)
		{
			corner.tiles[i]->plate = &plate;
			plate.tiles.push_back(corner.tiles[i]);
		}

		for (auto i = 0; i < corner.tiles.size(); ++i)
		{
			auto& tile = *corner.tiles[i];
			for (auto j = 0; j < tile.tiles.size(); ++j)
			{
				auto& adjacentTile = *tile.tiles[j];
				if (!adjacentTile.plate)
				{
					plateless.push_back(Plateless{ &adjacentTile, pplate });
				}
			}
		}
	}



	while (plateless.size() > 0)
	{
		size_t index = (size_t)std::floor(std::pow(random.unit(), 2) * plateless.size());
		auto indexit = getAdvancedIt(plateless, index);
		auto ptile = indexit->tile;
		auto pplate = indexit->plate;
		auto& tile = *ptile;
		auto& plate = *pplate;
		plateless.erase(indexit);
		if (!tile.plate)
		{
			tile.plate = &plate;
			plate.tiles.push_back(ptile);
			for (auto j = 0; j < tile.tiles.size(); ++j)
			{
				if (!tile.tiles[j]->plate)
				{
					plateless.push_back(Plateless{ tile.tiles[j], &plate });
				}
			}
		}
	}


	calculateCornerDistancesToPlateRoot(plates);
}

void World::calculateCornerDistancesToPlateRoot(std::vector<Plate>& plates)
{
	struct DistanceCorner //2do: order matters
	{
		Corner* corner;
		double distanceToPlateRoot;
	};
	std::list<DistanceCorner> distanceCornerQueue;
	for (auto i = 0; i < plates.size(); ++i)
	{
		auto& corner = *plates[i].root;
		corner.distanceToPlateRoot = 0;
		for (auto j = 0; j < corner.corners.size(); ++j)
		{
			distanceCornerQueue.push_back(DistanceCorner{ corner.corners[j], corner.borders[j]->length() });
		}
	}

	while (distanceCornerQueue.size() > 0)
	{
		size_t iEnd = distanceCornerQueue.size();
		for (auto i = 0; i < iEnd; ++i)
		{
			auto& front = getAdvancedItVal(distanceCornerQueue, i);
			auto& corner = *front.corner;
			auto distanceToPlateRoot = front.distanceToPlateRoot;
			if (!corner.distanceToPlateRoot || corner.distanceToPlateRoot > distanceToPlateRoot)
			{
				corner.distanceToPlateRoot = distanceToPlateRoot;
				for (auto j = 0; j < corner.corners.size(); ++j)
				{
					distanceCornerQueue.push_back(DistanceCorner{ corner.corners[j], distanceToPlateRoot + corner.borders[j]->length() });
				}
			}
		}
		auto itend = distanceCornerQueue.begin();
		std::advance(itend, iEnd);
		distanceCornerQueue.erase(distanceCornerQueue.begin(), itend);

		distanceCornerQueue.sort([](DistanceCorner const& left, DistanceCorner const& right) -> bool 
			{ return left.distanceToPlateRoot </*-*/ right.distanceToPlateRoot; });
	}
}

void World::generatePlanetElevation(Topology& topology, std::vector<Plate>& plates)
{
	std::vector<Corner*> boundaryCorners;
	std::vector<size_t> boundaryCornerInnerBorderIndexes;
	std::list<ElevationBorder> elevationBorderQueue;
	std::function<bool(const ElevationBorder&, const ElevationBorder&)> elevationBorderQueueSorter = 
		[](const ElevationBorder& left, const ElevationBorder& right) -> bool { return left.distanceToPlateBoundary < /*-*/ right.distanceToPlateBoundary; };

	identifyBoundaryBorders(topology.borders);
	collectBoundaryCorners(topology.corners, boundaryCorners);
	calculatePlateBoundaryStress(boundaryCorners, boundaryCornerInnerBorderIndexes);
	blurPlateBoundaryStress(boundaryCorners, 3, 0.4);
	populateElevationBorderQueue(boundaryCorners, boundaryCornerInnerBorderIndexes, elevationBorderQueue);
	processElevationBorderQueue(elevationBorderQueue, elevationBorderQueueSorter);
	calculateTileAverageElevations(topology.tiles);
}

void World::identifyBoundaryBorders(std::vector<Border>& borders)
{
	for (auto& border : borders)
	{
		if (border.tiles[0]->plate != border.tiles[1]->plate)
		{
			border.betweenPlates = true;
			border.corners[0]->betweenPlates = true;
			border.corners[1]->betweenPlates = true;
			border.tiles[0]->plate->boundaryBorders.push_back(&border);
			border.tiles[1]->plate->boundaryBorders.push_back(&border);
		}
	}
}

void World::collectBoundaryCorners(std::vector<Corner>& corners, std::vector<Corner*>& boundaryCorners)
{
	for (auto& corner : corners)
	{
		if (corner.betweenPlates)
		{
			boundaryCorners.push_back(&corner);
			corner.tiles[0]->plate->boundaryCorners.push_back(&corner);
			if (corner.tiles[1]->plate != corner.tiles[0]->plate) corner.tiles[1]->plate->boundaryCorners.push_back(&corner);
			if (corner.tiles[2]->plate != corner.tiles[0]->plate && corner.tiles[2]->plate != corner.tiles[1]->plate) corner.tiles[2]->plate->boundaryCorners.push_back(&corner);
		}
	}
}

const size_t UdefIdx = std::numeric_limits<size_t>::max();

void World::calculatePlateBoundaryStress(const std::vector<Corner*>& boundaryCorners, std::vector<size_t>& boundaryCornerInnerBorderIndexes)
{
	boundaryCornerInnerBorderIndexes.resize(boundaryCorners.size());
	size_t i = 0;
	for (auto pcorner : boundaryCorners)
	{
		auto& corner = *pcorner;
		corner.distanceToPlateBoundary = 0;

		Border* innerBorder = nullptr;
		size_t innerBorderIndex;
		for (size_t j = 0; j < corner.borders.size(); ++j)
		{
			auto pborder = corner.borders[j];
			auto& border = *pborder;
			if (!border.betweenPlates)
			{
				innerBorder = &border;
				innerBorderIndex = j;
				break;
			}
		}

		if (innerBorder)
		{
			boundaryCornerInnerBorderIndexes[i] = innerBorderIndex;
			auto& outerBorder0 = *corner.borders[(innerBorderIndex + 1) % corner.borders.size()];
			auto& outerBorder1 = *corner.borders[(innerBorderIndex + 2) % corner.borders.size()];
			auto& farCorner0 = outerBorder0.oppositeCorner(corner);
			auto& farCorner1 = outerBorder1.oppositeCorner(corner);
			auto& plate0 = *innerBorder->tiles[0]->plate;
			auto& plate1 = *(outerBorder0.tiles[0]->plate != &plate0 ? outerBorder0.tiles[0]->plate : outerBorder0.tiles[1]->plate);
			auto boundaryVector = farCorner0.vectorTo(farCorner1);
			auto boundaryNormal = boundaryVector.crossProduct(corner.position);
			Stress stress;
			calculateStress(plate0.calculateMovement(corner.position), plate1.calculateMovement(corner.position), boundaryVector, boundaryNormal, stress);
			corner.pressure = stress.pressure;
			corner.shear = stress.shear;
		}
		else
		{
			boundaryCornerInnerBorderIndexes[i] = UdefIdx;
			auto& plate0 = *corner.tiles[0]->plate;
			auto& plate1 = *corner.tiles[1]->plate;
			auto& plate2 = *corner.tiles[2]->plate;
			auto boundaryVector0 = corner.corners[0]->vectorTo(corner);
			auto boundaryVector1 = corner.corners[1]->vectorTo(corner);
			auto boundaryVector2 = corner.corners[2]->vectorTo(corner);
			auto boundaryNormal0 = boundaryVector0.crossProduct(corner.position);
			auto boundaryNormal1 = boundaryVector1.crossProduct(corner.position);
			auto boundaryNormal2 = boundaryVector2.crossProduct(corner.position);
			Stress stress0, stress1, stress2;
			calculateStress(plate0.calculateMovement(corner.position), plate1.calculateMovement(corner.position), boundaryVector0, boundaryNormal0, stress0);
			calculateStress(plate1.calculateMovement(corner.position), plate2.calculateMovement(corner.position), boundaryVector1, boundaryNormal1, stress1);
			calculateStress(plate2.calculateMovement(corner.position), plate0.calculateMovement(corner.position), boundaryVector2, boundaryNormal2, stress2);

			corner.pressure = (stress0.pressure + stress1.pressure + stress2.pressure) / 3;
			corner.shear = (stress0.shear + stress1.shear + stress2.shear) / 3;
		}

		++i;
	}
}

void World::calculateStress(const Ogre::Vector3& movement0, const Ogre::Vector3& movement1, const Ogre::Vector3& boundaryVector, const Ogre::Vector3& boundaryNormal, Stress& stress)
{
	auto relativeMovement = movement0 - movement1;
	auto pressureVector = projectOnVector(relativeMovement, boundaryNormal);
	double pressure = pressureVector.length();
	if (pressureVector.dotProduct(boundaryNormal) > 0) pressure = -pressure;
	double shear = projectOnVector(relativeMovement, boundaryVector).length();

	stress.pressure = 2 / (1 + std::exp(-pressure / 30)) - 1;
	stress.shear = 2 / (1 + std::exp(-shear / 30)) - 1;
}

void World::blurPlateBoundaryStress(const std::vector<Corner*>& boundaryCorners, size_t stressBlurIterations, double stressBlurCenterWeighting)
{
	std::vector<double> newCornerPressure(boundaryCorners.size());
	std::vector<double> newCornerShear(boundaryCorners.size());
	for (size_t i = 0; i < stressBlurIterations; ++i)
	{
		for (size_t j = 0; j < boundaryCorners.size(); ++j)
		{
			auto& corner = *boundaryCorners[j];
			double averagePressure = 0;
			double averageShear = 0;
			size_t neighborCount = 0;
			for (size_t k = 0; k < corner.corners.size(); ++k)
			{
				auto& neighbor = *corner.corners[k];
				if (neighbor.betweenPlates)
				{
					averagePressure += neighbor.pressure;
					averageShear += neighbor.shear;
					++neighborCount;
				}
			}
			newCornerPressure[j] = corner.pressure * stressBlurCenterWeighting + (averagePressure / neighborCount) * (1 - stressBlurCenterWeighting);
			newCornerShear[j] = corner.shear * stressBlurCenterWeighting + (averageShear / neighborCount) * (1 - stressBlurCenterWeighting);
		}

		for (size_t j = 0; j < boundaryCorners.size(); ++j)
		{
			auto& corner = *boundaryCorners[j];
			if (corner.betweenPlates)
			{
				corner.pressure = newCornerPressure[j];
				corner.shear = newCornerShear[j];
			}
		}
	}
}

double calculateCollidingElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear);
double calculateSuperductingElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear);
double calculateSubductingElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear);
double calculateDivergingElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear);
double calculateShearingElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear);
double calculateDormantElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear);

void World::populateElevationBorderQueue(const std::vector<Corner*>& boundaryCorners, const std::vector<size_t>& boundaryCornerInnerBorderIndexes, std::list<ElevationBorder>& elevationBorderQueue)
{
	size_t i = 0;
	for (auto pcorner : boundaryCorners)
	{
		auto& corner = *pcorner;

		size_t innerBorderIndex = boundaryCornerInnerBorderIndexes[i];
		if (innerBorderIndex != UdefIdx)
		{
			auto& innerBorder = *corner.borders[innerBorderIndex];
			auto& outerBorder0 = *corner.borders[(innerBorderIndex + 1) % corner.borders.size()];
			auto& plate0 = *innerBorder.tiles[0]->plate;
			auto& plate1 = *(outerBorder0.tiles[0]->plate != &plate0 ? outerBorder0.tiles[0]->plate : outerBorder0.tiles[1]->plate);

			ElevationBorderOrigin::CalculateElevationFunc calculateElevation;

			if (corner.pressure > 0.3)
			{
				corner.elevation = std::max(plate0.elevation, plate1.elevation) + corner.pressure;
				if (plate0.oceanic == plate1.oceanic)
					calculateElevation = &calculateCollidingElevation;
				else if (plate0.oceanic)
					calculateElevation = &calculateSubductingElevation;
				else
					calculateElevation = &calculateSuperductingElevation;
			}
			else if (corner.pressure < -0.3)
			{
				corner.elevation = std::max(plate0.elevation, plate1.elevation) - corner.pressure / 4;
				calculateElevation = &calculateDivergingElevation;
			}
			else if (corner.shear > 0.3)
			{
				corner.elevation = std::max(plate0.elevation, plate1.elevation) + corner.shear / 8;
				calculateElevation = &calculateShearingElevation;
			}
			else
			{
				corner.elevation = (plate0.elevation + plate1.elevation) / 2;
				calculateElevation = &calculateDormantElevation;
			}

			auto& nextCorner = innerBorder.oppositeCorner(corner);
			if (!nextCorner.betweenPlates)
			{
				elevationBorderQueue.push_back(ElevationBorder(
					ElevationBorderOrigin( &corner, corner.pressure, corner.shear, &plate0, calculateElevation ),
					&innerBorder,
					&corner,
					&nextCorner,
					innerBorder.length()
				));
			}
		}
		else
		{
			auto& plate0 = *corner.tiles[0]->plate;
			auto& plate1 = *corner.tiles[1]->plate;
			auto& plate2 = *corner.tiles[2]->plate;

			//2do: this was here, but was maybe just an artefact; or is it some javascript global??   elevation = 0;

			if (corner.pressure > 0.3)
			{
				corner.elevation = std::max(plate0.elevation, std::max(plate1.elevation, plate2.elevation)) + corner.pressure;
			}
			else if (corner.pressure < -0.3)
			{
				corner.elevation = std::max(plate0.elevation, std::max(plate1.elevation, plate2.elevation)) + corner.pressure / 4;
			}
			else if (corner.shear > 0.3)
			{
				corner.elevation = std::max(plate0.elevation, std::max(plate1.elevation, plate2.elevation)) + corner.shear / 8;
			}
			else
			{
				corner.elevation = (plate0.elevation + plate1.elevation + plate2.elevation) / 3;
			}
		}

		++i;
	}
}

double calculateCollidingElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear)
{
	double t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
	if (t < 0.5)
	{
		t = t / 0.5;
		return plateElevation + std::pow(t - 1, 2) * (boundaryElevation - plateElevation);
	}
	else
	{
		return plateElevation;
	}
}

double calculateSuperductingElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear)
{
	double t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
	if (t < 0.2)
	{
		t = t / 0.2;
		return boundaryElevation + t * (plateElevation - boundaryElevation + pressure / 2);
	}
	else if (t < 0.5)
	{
		t = (t - 0.2) / 0.3;
		return plateElevation + std::pow(t - 1, 2) * pressure / 2;
	}
	else
	{
		return plateElevation;
	}
}

double calculateSubductingElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear)
{
	double t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
	return plateElevation + std::pow(t - 1, 2) * (boundaryElevation - plateElevation);
}

double calculateDivergingElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear)
{
	double t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
	if (t < 0.3)
	{
		t = t / 0.3;
		return plateElevation + std::pow(t - 1, 2) * (boundaryElevation - plateElevation);
	}
	else
	{
		return plateElevation;
	}
}

double calculateShearingElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear)
{
	double t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
	if (t < 0.2)
	{
		t = t / 0.2;
		return plateElevation + std::pow(t - 1, 2) * (boundaryElevation - plateElevation);
	}
	else
	{
		return plateElevation;
	}
}

double calculateDormantElevation(double distanceToPlateBoundary, double distanceToPlateRoot, double boundaryElevation, double plateElevation, double pressure, double shear)
{
	double t = distanceToPlateBoundary / (distanceToPlateBoundary + distanceToPlateRoot);
	double elevationDifference = boundaryElevation - plateElevation;
	double a = 2 * elevationDifference;
	double b = -3 * elevationDifference;
	return t * t * elevationDifference * (2 * t - 3) + boundaryElevation;
}

void World::processElevationBorderQueue(std::list<ElevationBorder>& elevationBorderQueue, std::function<bool(const ElevationBorder&, const ElevationBorder&)> elevationBorderQueueSorter)
{
	while (elevationBorderQueue.size() > 0)
	{
		size_t iEnd = elevationBorderQueue.size();
		auto it = elevationBorderQueue.begin();
		for (auto i = 0; i < iEnd; ++i)
		{
			auto& front = *it;
			++it;
			auto& corner = *front.nextCorner;
			if (!corner.elevation)
			{
				corner.distanceToPlateBoundary = front.distanceToPlateBoundary;
				corner.elevation = front.origin.calculateElevation(
					corner.distanceToPlateBoundary,
					corner.distanceToPlateRoot,
					front.origin.corner->elevation,
					front.origin.plate->elevation,
					front.origin.pressure,
					front.origin.shear);

				for (auto j = 0; j < corner.borders.size(); ++j)
				{
					auto& border = *corner.borders[j];
					if (!border.betweenPlates)
					{
						auto& nextCorner = *corner.corners[j];
						double distanceToPlateBoundary = corner.distanceToPlateBoundary + border.length();
						if (!nextCorner.distanceToPlateBoundary || nextCorner.distanceToPlateBoundary > distanceToPlateBoundary)
						{
							elevationBorderQueue.push_back(ElevationBorder(
								front.origin, // 2do: maybe wrong storage because there are 2 distinctive creation types, so is the ref one wrong?
								&border,
								&corner,
								&nextCorner,
								distanceToPlateBoundary
							));
						}
					}
				}
			}
		}
		//elevationBorderQueue.splice(0, iEnd);
		auto itend = elevationBorderQueue.begin();
		std::advance(itend, iEnd);
		elevationBorderQueue.erase(elevationBorderQueue.begin(), itend);

		elevationBorderQueue.sort(elevationBorderQueueSorter);
	}
}

void World::calculateTileAverageElevations(std::vector<Tile>& tiles)
{
	for (auto& tile: tiles)
	{
		double elevation = 0;
		for (auto j = 0; j < tile.corners.size(); ++j)
		{
			elevation += tile.corners[j]->elevation;
		}
		tile.elevation = elevation / tile.corners.size();
	}
}

void World::generatePlanetWeather(Topology& topology, SpatialPartition& partition, double heatLevel, double moistureLevel, XorShift128& random)
{
	double planetRadius = 1000;
	std::vector<Whorl> whorls;
	std::list<Corner*> activeCorners;
	double totalHeat;
	double remainingHeat;
	double totalMoisture;
	double remainingMoisture;

	
	generateAirCurrentWhorls(planetRadius, random, whorls);
	calculateAirCurrents(topology.corners, whorls, planetRadius);


	initializeAirHeat(topology.corners, heatLevel, activeCorners, totalHeat);
	remainingHeat = totalHeat;
	double consumedHeat;
	do {
		consumedHeat = processAirHeat(activeCorners);
		remainingHeat -= consumedHeat;
	}
	while (remainingHeat > 0 && consumedHeat >= 0.0001); // action.loop(1 - remainingHeat / totalHeat);
	
	calculateTemperature(topology.corners, topology.tiles, planetRadius);


	initializeAirMoisture(topology.corners, moistureLevel, activeCorners, totalMoisture);
	remainingMoisture = totalMoisture;
	double consumedMoisture;
	do
	{
		consumedMoisture = processAirMoisture(activeCorners);
		remainingMoisture -= consumedMoisture;
	}
	while (remainingMoisture > 0 && consumedMoisture >= 0.0001); // action.loop(1 - remainingMoisture / totalMoisture);

	calculateMoisture(topology.corners, topology.tiles);
}

void World::generateAirCurrentWhorls(double planetRadius, XorShift128& random, std::vector<Whorl>& whorls)
{
	double direction = random.integer(0, 1) ? 1 : -1;
	double layerCount = random.integer(4, 7);
	double circumference = M_PI * 2 * planetRadius;
	double fullRevolution = M_PI * 2;
	double baseWhorlRadius = circumference / (2 * (layerCount - 1));

	whorls.reserve(1 + layerCount * (size_t)(std::ceil(planetRadius * fullRevolution / baseWhorlRadius)) / 2 + 1);

	Ogre::Matrix3 m1, m2;
	m1.FromAngleAxis(Ogre::Vector3(1, 0, 0), Ogre::Radian(random.realInclusive(0, fullRevolution / (2 * (layerCount + 4)))));
	m2.FromAngleAxis(Ogre::Vector3(0, 1, 0), Ogre::Radian(random.real(0, fullRevolution)));
	whorls.push_back(Whorl{
		Ogre::Vector3(0, planetRadius, 0) * m1 * m2,
		random.realInclusive(fullRevolution / 36, fullRevolution / 24) * direction,
		random.realInclusive(baseWhorlRadius * 0.8, baseWhorlRadius * 1.2)
	});

	for (size_t i = 1; i < layerCount - 1; ++i)
	{
		direction = -direction;
		double baseTilt = i / (layerCount - 1) * fullRevolution / 2;
		size_t layerWhorlCount = (size_t)(std::ceil((std::sin(baseTilt) * planetRadius * fullRevolution) / baseWhorlRadius));
		Ogre::Matrix3 m3, m4;
		m3.FromAngleAxis(Ogre::Vector3(1, 0, 0), Ogre::Radian(baseTilt));
		for (size_t j = 0; j < layerWhorlCount; ++j)
		{
			m4.FromAngleAxis(Ogre::Vector3(0, 1, 0), Ogre::Radian(fullRevolution * (j + (i % 2) / 2) / layerWhorlCount));
			whorls.push_back(Whorl{
				Ogre::Vector3(0, planetRadius, 0) * m1 * m2 * m3 * m4,
				random.realInclusive(fullRevolution / 48, fullRevolution / 32) * direction,
				random.realInclusive(baseWhorlRadius * 0.8, baseWhorlRadius * 1.2)
			});
		}
	}

	direction = -direction;
	Ogre::Matrix3 m3;
	m3.FromAngleAxis(Ogre::Vector3(1, 0, 0), Ogre::Radian(fullRevolution / 2));
	whorls.push_back(Whorl{
		Ogre::Vector3(0, planetRadius, 0) * m1 * m2 * m3,
		random.realInclusive(fullRevolution / 36, fullRevolution / 24) * direction,
		random.realInclusive(baseWhorlRadius * 0.8, baseWhorlRadius * 1.2)
	});
}

void World::calculateAirCurrents(std::vector<Corner>& corners, const std::vector<Whorl>& whorls, double planetRadius)
{
	for (auto& corner : corners)
	{
		Ogre::Vector3 airCurrent(0, 0, 0);
		double weight = 0;
		for (auto& whorl : whorls)
		{
			auto angle = whorl.center.angleBetween(corner.position);
			double distance = angle.valueRadians() * planetRadius;
			if (distance < whorl.radius)
			{
				double normalizedDistance = distance / whorl.radius;
				double whorlWeight = 1 - normalizedDistance;
				double whorlStrength = planetRadius * whorl.strength * whorlWeight * normalizedDistance;
				Ogre::Vector3 whorlCurrent = setLength(whorl.center.crossProduct(corner.position), whorlStrength);
				airCurrent += whorlCurrent;
				weight += whorlWeight;
			}
		}
		airCurrent /= weight;
		corner.airCurrent = airCurrent;
		corner.airCurrentSpeed = airCurrent.length(); //kilometers per hour

		corner.airCurrentOutflows.reserve(corner.borders.size());
		auto airCurrentDirection = airCurrent.normalisedCopy();
		double outflowSum = 0;
		for (auto pcornerCorner : corner.corners)
		{
			auto vector = corner.vectorTo(*pcornerCorner).normalisedCopy();
			auto dot = vector.dotProduct(airCurrentDirection);
			if (dot > 0)
			{
				corner.airCurrentOutflows.push_back(dot);
				outflowSum += dot;
			}
			else
			{
				corner.airCurrentOutflows.push_back(0);
			}
		}

		if (outflowSum > 0)
		{
			for (auto j = 0; j < corner.borders.size(); ++j)
			{
				corner.airCurrentOutflows[j] /= outflowSum;
			}
		}
	}
}

void World::initializeAirHeat(std::vector<Corner>& corners, double heatLevel, std::list<Corner*>& activeCorners, double& airHeat)
{
	activeCorners.clear();
	airHeat = 0;
	for (auto& corner : corners)
	{
		corner.airHeat = corner.area * heatLevel;
		corner.newAirHeat = 0;
		corner.heat = 0;

		corner.heatAbsorption = 0.1 * corner.area / std::max(0.1, std::min(corner.airCurrentSpeed, 1.0));
		if (corner.elevation <= 0)
		{
			corner.maxHeat = corner.area;
		}
		else
		{
			corner.maxHeat = corner.area;
			corner.heatAbsorption *= 2;
		}

		activeCorners.push_back(&corner);
		airHeat += corner.airHeat;
	}
}

double World::processAirHeat(std::list<Corner*>& activeCorners)
{
	double consumedHeat = 0;
	auto activeCornerCount = activeCorners.size();

	auto it = activeCorners.begin();
	for (auto i = 0; i < activeCornerCount; ++i, ++it)
	{
		auto& corner = *(*it);
		if (corner.airHeat == 0) continue;

		double heatChange = std::max(0.0, std::min(corner.airHeat, corner.heatAbsorption * (1 - corner.heat / corner.maxHeat)));
		corner.heat += heatChange;
		consumedHeat += heatChange;
		double heatLoss = corner.area * (corner.heat / corner.maxHeat) * 0.02;
		heatChange = std::min(corner.airHeat, heatChange + heatLoss);

		double remainingCornerAirHeat = corner.airHeat - heatChange;
		corner.airHeat = 0;

		for (auto j = 0; j < corner.corners.size(); ++j)
		{
			auto outflow = corner.airCurrentOutflows[j];
			if (outflow > 0)
			{
				corner.corners[j]->newAirHeat += remainingCornerAirHeat * outflow;
				activeCorners.push_back(corner.corners[j]);
			}
		}
	}

	activeCorners.erase(activeCorners.begin(), it);

	for (auto pcorner : activeCorners)
	{
		pcorner->airHeat = pcorner->newAirHeat;
		pcorner->newAirHeat = 0;
	}

	return consumedHeat;
}

void World::calculateTemperature(std::vector<Corner>& corners, std::vector<Tile>& tiles, double planetRadius)
{
	for (auto& corner : corners)
	{
		double latitudeEffect = std::sqrt(1 - std::abs(corner.position.y) / planetRadius);
		double elevationEffect = 1 - std::pow(std::max(0.0, std::min(corner.elevation * 0.8, 1.0)), 2.0);
		double normalizedHeat = corner.heat / corner.area;
		corner.temperature = (latitudeEffect * elevationEffect * 0.7 + normalizedHeat * 0.3) * 5. / 3 - 2. / 3;
		// 2do
		//delete corner.airHeat;
		//delete corner.newAirHeat;
		//delete corner.heat;
		//delete corner.maxHeat;
		//delete corner.heatAbsorption;
	}

	for (auto& tile : tiles)
	{
		tile.temperature = 0;
		for (auto pcorner : tile.corners)
		{
			tile.temperature += pcorner->temperature;
		}
		tile.temperature /= tile.corners.size();
	}
}

void World::initializeAirMoisture(std::vector<Corner>& corners, double moistureLevel, std::list<Corner*>& activeCorners, double& airMoisture)
{
	airMoisture = 0;
	for (auto& corner : corners)
	{
		corner.airMoisture = (corner.elevation > 0) ? 0 : corner.area * moistureLevel * std::max(0.0, std::min(0.5 + corner.temperature * 0.5, 1.0));
		corner.newAirMoisture = 0;
		corner.precipitation = 0;

		corner.precipitationRate = 0.0075 * corner.area / std::max(0.1, std::min(corner.airCurrentSpeed, 1.0));
		corner.precipitationRate *= 1 + (1 - std::max(0.0, std::min(corner.temperature, 1.0))) * 0.1;
		if (corner.elevation > 0)
		{
			corner.precipitationRate *= 1 + corner.elevation * 0.5;
			corner.maxPrecipitation = corner.area * (0.25 + std::max(0.0, std::min(corner.elevation, 1.0)) * 0.25);
		}
		else
		{
			corner.maxPrecipitation = corner.area * 0.25;
		}

		activeCorners.push_back(&corner);
		airMoisture += corner.airMoisture;
	}
}

double World::processAirMoisture(std::list<Corner*>& activeCorners)
{
	double consumedMoisture = 0;
	auto activeCornerCount = activeCorners.size();
	auto it = activeCorners.begin();
	for (auto i = 0; i < activeCornerCount; ++i, ++it)
	{
		auto& corner = *(*it);
		if (corner.airMoisture == 0) continue;

		double moistureChange = std::max(0.0, std::min(corner.airMoisture, corner.precipitationRate * (1.0 - corner.precipitation / corner.maxPrecipitation)));
		corner.precipitation += moistureChange;
		consumedMoisture += moistureChange;
		double moistureLoss = corner.area * (corner.precipitation / corner.maxPrecipitation) * 0.02;
		moistureChange = std::min(corner.airMoisture, moistureChange + moistureLoss);

		double remainingCornerAirMoisture = corner.airMoisture - moistureChange;
		corner.airMoisture = 0;

		for (auto j = 0; j < corner.corners.size(); ++j)
		{
			double outflow = corner.airCurrentOutflows[j];
			if (outflow > 0)
			{
				corner.corners[j]->newAirMoisture += remainingCornerAirMoisture * outflow;
				activeCorners.push_back(corner.corners[j]);
			}
		}
	}

	activeCorners.erase(activeCorners.begin(), it);

	for (auto pcorner : activeCorners)
	{
		pcorner->airMoisture = pcorner->newAirMoisture;
		pcorner->newAirMoisture = 0;
	}

	return consumedMoisture;
}

void World::calculateMoisture(std::vector<Corner>& corners, std::vector<Tile>& tiles)
{
	for (auto& corner : corners)
	{
		corner.moisture = corner.precipitation / corner.area / 0.5;
		// 2do
		//delete corner.airMoisture;
		//delete corner.newAirMoisture;
		//delete corner.precipitation;
		//delete corner.maxPrecipitation;
		//delete corner.precipitationRate;
	}

	for (auto& tile : tiles)
	{
		tile.moisture = 0;
		for (auto pcorner : tile.corners)
		{
			tile.moisture += pcorner->temperature;
		}
		tile.moisture /= tile.corners.size();
	}
}

void World::generatePlanetBiomes(std::vector<Tile>& tiles, double planetRadius)
{
	for (auto& tile : tiles)
	{
		auto elevation = std::max(0.0, tile.elevation);
		auto latitude = std::abs(tile.position.y / planetRadius);
		auto temperature = tile.temperature;
		auto moisture = tile.moisture;

		if (elevation <= 0)
		{
			if (temperature > 0)
			{
				tile.biome = "ocean";
			}
			else
			{
				tile.biome = "oceanGlacier";
			}
		}
		else if (elevation < 0.6)
		{
			if (temperature > 0.75)
			{
				if (moisture < 0.25)
				{
					tile.biome = "desert";
				}
				else
				{
					tile.biome = "rainForest";
				}
			}
			else if (temperature > 0.5)
			{
				if (moisture < 0.25)
				{
					tile.biome = "rocky";
				}
				else if (moisture < 0.50)
				{
					tile.biome = "plains";
				}
				else
				{
					tile.biome = "swamp";
				}
			}
			else if (temperature > 0)
			{
				if (moisture < 0.25)
				{
					tile.biome = "plains";
				}
				else if (moisture < 0.50)
				{
					tile.biome = "grassland";
				}
				else
				{
					tile.biome = "deciduousForest";
				}
			}
			else
			{
				if (moisture < 0.25)
				{
					tile.biome = "tundra";
				}
				else
				{
					tile.biome = "landGlacier";
				}
			}
		}
		else if (elevation < 0.8)
		{
			if (temperature > 0)
			{
				if (moisture < 0.25)
				{
					tile.biome = "tundra";
				}
				else
				{
					tile.biome = "coniferForest";
				}
			}
			else
			{
				tile.biome = "tundra";
			}
		}
		else
		{
			if (temperature > 0 || moisture < 0.25)
			{
				tile.biome = "mountain";
			}
			else
			{
				tile.biome = "snowyMountain";
			}
		}
	}
}

void World::generatePlanetRenderData(Topology& topology, XorShift128& random, RenderData& renderData)
{
	buildSurfaceRenderObject(topology.tiles, random, *renderData.surface);
	buildPlateBoundariesRenderObject(topology.borders, *renderData.plateBoundaries);
	buildPlateMovementsRenderObject(topology.tiles, *renderData.plateMovements);
	buildAirCurrentsRenderObject(topology.corners, *renderData.airCurrents);
}


void doBuildTileWedge(Ogre::ManualObject& mo, size_t b, size_t s, size_t t)
{
	mo.triangle(b + s + 2, b, b + t + 2);
	mo.triangle(b + s + 1, b + t + 2, b + t + 1);
	mo.triangle(b + s + 1, b + s + 2, b + t + 2);
}

void World::buildSurfaceRenderObject(std::vector<Tile>& tiles, XorShift128& random, Ogre::ManualObject& mo) // SurfaceRenderObj& surfaceRenderObj)
{
	mo.begin("experilousworld/planetmat1", Ogre::RenderOperation::OT_TRIANGLE_LIST);

	auto baseIndex = 0;
	for (auto& tile : tiles)
	{
		auto colorDeviance = Ogre::ColourValue(random.unit(), random.unit(), random.unit());
		Ogre::ColourValue terrainColor;
		if (tile.elevation <= 0)
		{
			auto normalizedElevation = std::min(-tile.elevation, 1.0);
			if (tile.biome == "ocean") terrainColor = lerp(lerp(ocv(0x0066FF), ocv(0x0044BB), std::min(-tile.elevation, 1.0)), colorDeviance, 0.10);
			else if (tile.biome == "oceanGlacier") terrainColor = lerp(ocv(0xDDEEFF), colorDeviance, 0.10);
			else terrainColor = ocv(0xFF00FF);
		}
		else if (tile.elevation < 0.6)
		{
			auto normalizedElevation = tile.elevation / 0.6;
			if (tile.biome == "desert") terrainColor = lerp(lerp(ocv(0xDDDD77), ocv(0xBBBB55), normalizedElevation), colorDeviance, 0.10);
			else if (tile.biome == "rainForest") terrainColor = lerp(lerp(ocv(0x44DD00), ocv(0x229900), normalizedElevation), colorDeviance, 0.20);
			else if (tile.biome == "rocky") terrainColor = lerp(lerp(ocv(0xAA9977), ocv(0x887755), normalizedElevation), colorDeviance, 0.15);
			else if (tile.biome == "plains") terrainColor = lerp(lerp(ocv(0x99BB44), ocv(0x667722), normalizedElevation), colorDeviance, 0.10);
			else if (tile.biome == "grassland") terrainColor = lerp(lerp(ocv(0x77CC44), ocv(0x448822), normalizedElevation), colorDeviance, 0.15);
			else if (tile.biome == "swamp") terrainColor = lerp(lerp(ocv(0x77AA44), ocv(0x446622), normalizedElevation), colorDeviance, 0.25);
			else if (tile.biome == "deciduousForest") terrainColor = lerp(lerp(ocv(0x33AA22), ocv(0x116600), normalizedElevation), colorDeviance, 0.10);
			else if (tile.biome == "tundra") terrainColor = lerp(lerp(ocv(0x9999AA), ocv(0x777788), normalizedElevation), colorDeviance, 0.15);
			else if (tile.biome == "landGlacier") terrainColor = lerp(ocv(0xDDEEFF), colorDeviance, 0.10);
			else terrainColor = ocv(0xFF00FF);
		}
		else if (tile.elevation < 0.8)
		{
			auto normalizedElevation = (tile.elevation - 0.6) / 0.2;
			if (tile.biome == "tundra") terrainColor = lerp(lerp(ocv(0x777788), ocv(0x666677), normalizedElevation), colorDeviance, 0.10);
			else if (tile.biome == "coniferForest") terrainColor = lerp(lerp(ocv(0x338822), ocv(0x116600), normalizedElevation), colorDeviance, 0.10);
			else if (tile.biome == "snow") terrainColor = lerp(lerp(ocv(0xEEEEEE), ocv(0xDDDDDD), normalizedElevation), colorDeviance, 0.10);
			else if (tile.biome == "mountain") terrainColor = lerp(lerp(ocv(0x555544), ocv(0x444433), normalizedElevation), colorDeviance, 0.05);
			else terrainColor = ocv(0xFF00FF);
		}
		else
		{
			auto normalizedElevation = std::min((tile.elevation - 0.8) / 0.5, 1.0);
			if (tile.biome == "mountain") terrainColor = lerp(lerp(ocv(0x444433), ocv(0x333322), normalizedElevation), colorDeviance, 0.05);
			else if (tile.biome == "snowyMountain") terrainColor = lerp(lerp(ocv(0xDDDDDD), ocv(0xFFFFFF), normalizedElevation), colorDeviance, 0.10);
			else terrainColor = ocv(0xFF00FF);
		}

		auto plateColor = tile.plate->color;

		Ogre::ColourValue elevationColor;
		if (tile.elevation <= 0) elevationColor = lerp(ocv(0x224488), ocv(0xAADDFF), std::max(0.0, std::min((tile.elevation + 3. / 4) / (3. / 4), 1.)));
		else if (tile.elevation < 0.75) elevationColor = lerp(ocv(0x997755), ocv(0x553311), std::max(0.0, std::min((tile.elevation) / (3. / 4), 1.)));
		else elevationColor = lerp(ocv(0x553311), ocv(0x222222), std::max(0.0, std::min((tile.elevation - 3. / 4) / (1. / 2), 1.0)));

		Ogre::ColourValue temperatureColor;
		if (tile.temperature <= 0) temperatureColor = lerp(ocv(0x0000FF), ocv(0xBBDDFF), std::max(0.0, std::min((tile.temperature + 2. / 3) / (2. / 3), 1.)));
		else temperatureColor = lerp(ocv(0xFFFF00), ocv(0xFF0000), std::max(0.0, std::min((tile.temperature) / (3.0 / 3), 1.)));

		auto moistureColor = lerp(ocv(0xFFCC00), ocv(0x0066FF), std::max(0.0, std::min(tile.moisture, 1.)));




		mo.position(tile.averagePosition);
		mo.normal(tile.normal);
		mo.colour(terrainColor); // this essentially replaces setSurfaceRenderMode on the fly or at least tries to. 2do..
		for (auto j = 0; j < tile.corners.size(); ++j)
		{
			auto& cornerPosition = tile.corners[j]->position;

			mo.position(cornerPosition);
			mo.normal(tile.normal);
			mo.colour(terrainColor);

			mo.position((tile.averagePosition - cornerPosition) * 0.1 + cornerPosition); // planetGeometry.vertices.push(tile.averagePosition.clone().sub(cornerPosition).multiplyScalar(0.1).add(cornerPosition));
			mo.normal(tile.normal); // buildTileWedge(mo, baseIndex, i0, i1, tile.normal);
			mo.colour(terrainColor * 0.5);

			size_t i0 = j * 2;
			size_t i1 = ((j + 1) % tile.corners.size()) * 2;
			doBuildTileWedge(mo, (size_t)baseIndex, i0, i1);
		}

		baseIndex += 1 + tile.corners.size() * 2;
	}

	mo.end();
}

void World::buildPlateBoundariesRenderObject(std::vector<Border>& borders, Ogre::ManualObject& mo)
{
	mo.begin("experilousworld/planetmat1", Ogre::RenderOperation::OT_TRIANGLE_LIST);

	auto baseIndex = mo.getCurrentVertexCount();
	for (auto& border : borders)
	{
		if (border.betweenPlates)
		{
			auto normal = border.midpoint.normalisedCopy();
			auto offset = normal * 1;

			auto& borderPoint0 = border.corners[0]->position;
			auto& borderPoint1 = border.corners[1]->position;
			auto& tilePoint0 = border.tiles[0]->averagePosition;
			auto& tilePoint1 = border.tiles[1]->averagePosition;
			auto nl = (border.tiles[0]->normal + border.tiles[1]->normal) / 2;

			auto pressure = std::max(-1.0, std::min((border.corners[0]->pressure + border.corners[1]->pressure) / 2, 1.0));
			auto shear = std::max(0.0, std::min((border.corners[0]->shear + border.corners[1]->shear) / 2, 1.0));
			auto innerColor = (pressure <= 0) ? Ogre::ColourValue(1 + pressure, 1.0, 0.0) : Ogre::ColourValue(1.0, 1 - pressure, 0.0);
			auto outerColor = Ogre::ColourValue(0.0, shear / 2, shear);

			mo.position(borderPoint0 + offset);
			mo.normal(nl);
			mo.colour(innerColor);

			mo.position(borderPoint1 + offset);
			mo.normal(nl);
			mo.colour(innerColor);

			mo.position((tilePoint0 - borderPoint0) * 0.2 + borderPoint0 + offset);
			mo.normal(nl);
			mo.colour(outerColor);

			mo.position((tilePoint0 - borderPoint1) * 0.2 + borderPoint1 + offset);
			mo.normal(nl);
			mo.colour(outerColor);

			mo.position((tilePoint1 - borderPoint0) * 0.2 + borderPoint0 + offset);
			mo.normal(nl);
			mo.colour(outerColor);

			mo.position((tilePoint1 - borderPoint1) * 0.2 + borderPoint1 + offset);
			mo.normal(nl);
			mo.colour(outerColor);


			mo.triangle(baseIndex + 0, baseIndex + 2, baseIndex + 1);
			mo.triangle(baseIndex + 1, baseIndex + 2, baseIndex + 3);
			mo.triangle(baseIndex + 1, baseIndex + 5, baseIndex + 0);
			mo.triangle(baseIndex + 0, baseIndex + 5, baseIndex + 4);

			baseIndex += 6;
		}
	}

	mo.end();
}

void World::buildPlateMovementsRenderObject(std::vector<Tile>& tiles, Ogre::ManualObject& mo)
{
	mo.begin("experilousworld/planetmat1", Ogre::RenderOperation::OT_TRIANGLE_LIST);

	for (auto& tile : tiles)
	{
		auto& plate = *tile.plate;
		auto movement = plate.calculateMovement(tile.position);
		auto plateMovementColor = Ogre::ColourValue(1 - plate.color.r, 1 - plate.color.g, 1 - plate.color.b);

		buildArrow(mo, tile.position * 1.002, movement * 0.5, tile.position.normalisedCopy(), std::min((double)movement.length(), 4.0), plateMovementColor);

		tile.plateMovement = movement;
	}

	mo.end();
}

void World::buildAirCurrentsRenderObject(std::vector<Corner>& corners, Ogre::ManualObject& mo)
{
	mo.begin("experilousworld/planetmat1", Ogre::RenderOperation::OT_TRIANGLE_LIST);

	for (auto& corner : corners)
	{
		buildArrow(mo, corner.position * 1.002, corner.airCurrent * 0.5, corner.position.normalisedCopy(), std::min((double)corner.airCurrent.length(), 4.0), Ogre::ColourValue(1,1,1,1));
	}

	mo.end();
}

void World::buildArrow(Ogre::ManualObject& mo, const Ogre::Vector3& position, const Ogre::Vector3& direction, const Ogre::Vector3& normal, double baseWidth, const Ogre::ColourValue& color)
{
	if (direction.squaredLength() == 0) return;

	auto sideOffset = setLength(direction.crossProduct(normal), baseWidth / 2);
	auto baseIndex = mo.getCurrentVertexCount();

	mo.position(position + sideOffset);
	mo.normal(normal);
	mo.colour(color);

	mo.position(position + direction);
	mo.normal(normal);
	mo.colour(color);

	mo.position(position - sideOffset);
	mo.normal(normal);
	mo.colour(color);

	mo.triangle(baseIndex, baseIndex + 1, baseIndex + 2);
}


template <typename T>
void updateMinMaxAvg(ValueStats<T>& stats, T value)
{
	stats.min = std::min(stats.min, value);
	stats.max = std::max(stats.max, value);
	stats.avg += value;
}

void World::generatePlanetStatistics(Topology& topology, std::vector<Plate>& plates, PlanetStatistics& planetStatistics)
{
	auto& statistics = planetStatistics;

	statistics.corners.count = topology.corners.size();
	statistics.corners.airCurrent.reset();
	statistics.corners.elevation.reset();
	statistics.corners.temperature.reset();
	statistics.corners.moisture.reset();
	statistics.corners.distanceToPlateBoundary.reset();
	statistics.corners.distanceToPlateRoot.reset();
	statistics.corners.pressure.reset();
	statistics.corners.shear.reset();
	statistics.corners.doublePlateBoundaryCount = 0;
	statistics.corners.triplePlateBoundaryCount = 0;
	statistics.corners.innerLandBoundaryCount = 0;
	statistics.corners.outerLandBoundaryCount = 0;

	for (auto& corner : topology.corners)
	{
		updateMinMaxAvg(statistics.corners.airCurrent, (double)corner.airCurrent.length());
		updateMinMaxAvg(statistics.corners.elevation, corner.elevation);
		updateMinMaxAvg(statistics.corners.temperature, corner.temperature);
		updateMinMaxAvg(statistics.corners.moisture, corner.moisture);
		updateMinMaxAvg(statistics.corners.distanceToPlateBoundary, corner.distanceToPlateBoundary);
		updateMinMaxAvg(statistics.corners.distanceToPlateRoot, corner.distanceToPlateRoot);
		if (corner.betweenPlates)
		{
			updateMinMaxAvg(statistics.corners.pressure, corner.pressure);
			updateMinMaxAvg(statistics.corners.shear, corner.shear);
			if (!corner.borders[0]->betweenPlates || !corner.borders[1]->betweenPlates || !corner.borders[2]->betweenPlates)
			{
				statistics.corners.doublePlateBoundaryCount += 1;
			}
			else
			{
				statistics.corners.triplePlateBoundaryCount += 1;
			}
		}
		auto landCount = ((corner.tiles[0]->elevation > 0) ? 1 : 0) + ((corner.tiles[1]->elevation > 0) ? 1 : 0) + ((corner.tiles[2]->elevation > 0) ? 1 : 0);
		if (landCount == 2)
		{
			statistics.corners.innerLandBoundaryCount += 1;
		}
		else if (landCount == 1)
		{
			statistics.corners.outerLandBoundaryCount += 1;
		}
		if (corner.corners.size() != 3) throw "Corner has as invalid number of neighboring corners.";
		if (corner.borders.size() != 3) throw "Corner has as invalid number of borders.";
		if (corner.tiles.size() != 3) throw "Corner has as invalid number of tiles.";
	}

	statistics.corners.airCurrent.avg /= statistics.corners.count;
	statistics.corners.elevation.avg /= statistics.corners.count;
	statistics.corners.temperature.avg /= statistics.corners.count;
	statistics.corners.moisture.avg /= statistics.corners.count;
	statistics.corners.distanceToPlateBoundary.avg /= statistics.corners.count;
	statistics.corners.distanceToPlateRoot.avg /= statistics.corners.count;
	statistics.corners.pressure.avg /= (statistics.corners.doublePlateBoundaryCount + statistics.corners.triplePlateBoundaryCount);
	statistics.corners.shear.avg /= (statistics.corners.doublePlateBoundaryCount + statistics.corners.triplePlateBoundaryCount);

	statistics.borders.count = topology.borders.size();
	statistics.borders.length.reset();
	statistics.borders.plateBoundaryCount = 0,
	statistics.borders.plateBoundaryPercentage = 0,
	statistics.borders.landBoundaryCount = 0;
	statistics.borders.landBoundaryPercentage = 0;

	for (auto& border : topology.borders)
	{
		auto length = border.length();
		updateMinMaxAvg(statistics.borders.length, length);
		if (border.betweenPlates)
		{
			statistics.borders.plateBoundaryCount += 1;
			statistics.borders.plateBoundaryPercentage += length;
		}
		if (border.isLandBoundary())
		{
			statistics.borders.landBoundaryCount += 1;
			statistics.borders.landBoundaryPercentage += length;
		}
		if (border.corners.size() != 2) throw "Border has as invalid number of corners.";
		if (border.borders.size() != 4) throw "Border has as invalid number of neighboring borders.";
		if (border.tiles.size() != 2) throw "Border has as invalid number of tiles.";
	}

	statistics.borders.plateBoundaryPercentage /= statistics.borders.length.avg;
	statistics.borders.landBoundaryPercentage /= statistics.borders.length.avg;
	statistics.borders.length.avg /= statistics.borders.count;

	statistics.tiles.count = topology.tiles.size();
	statistics.tiles.totalArea = 0;
	statistics.tiles.area.reset();
	statistics.tiles.elevation.reset();
	statistics.tiles.temperature.reset();
	statistics.tiles.moisture.reset();
	statistics.tiles.plateMovement.reset();
	statistics.tiles.biomeCounts.clear();
	statistics.tiles.biomeAreas.clear();
	statistics.tiles.pentagonCount = 0;
	statistics.tiles.hexagonCount = 0;
	statistics.tiles.heptagonCount = 0;

	for (auto& tile : topology.tiles)
	{
		updateMinMaxAvg(statistics.tiles.area, tile.area);
		updateMinMaxAvg(statistics.tiles.elevation, tile.elevation);
		updateMinMaxAvg(statistics.tiles.temperature, tile.temperature);
		updateMinMaxAvg(statistics.tiles.moisture, tile.moisture);
		updateMinMaxAvg(statistics.tiles.plateMovement, (double)tile.plateMovement.length());
		if (!statistics.tiles.biomeCounts[tile.biome]) statistics.tiles.biomeCounts[tile.biome] = 0;
		statistics.tiles.biomeCounts[tile.biome] += 1;
		if (!statistics.tiles.biomeAreas[tile.biome]) statistics.tiles.biomeAreas[tile.biome] = 0;
		statistics.tiles.biomeAreas[tile.biome] += tile.area;
		if (tile.tiles.size() == 5) statistics.tiles.pentagonCount += 1;
		else if (tile.tiles.size() == 6) statistics.tiles.hexagonCount += 1;
		else if (tile.tiles.size() == 7) statistics.tiles.heptagonCount += 1;
		else throw "Tile has an invalid number of neighboring tiles.";
		if (tile.tiles.size() != tile.borders.size()) throw "Tile has a neighbor and border count that do not match.";
		if (tile.tiles.size() != tile.corners.size()) throw "Tile has a neighbor and corner count that do not match.";
	}

	statistics.tiles.totalArea = statistics.tiles.area.avg;
	statistics.tiles.area.avg /= statistics.tiles.count;
	statistics.tiles.elevation.avg /= statistics.tiles.count;
	statistics.tiles.temperature.avg /= statistics.tiles.count;
	statistics.tiles.moisture.avg /= statistics.tiles.count;
	statistics.tiles.plateMovement.avg /= statistics.tiles.count;

	statistics.plates.count = plates.size();
	statistics.plates.tileCount.reset();
	statistics.plates.area.reset();
	statistics.plates.boundaryElevation.reset();
	statistics.plates.boundaryBorders.reset();
	statistics.plates.circumference.reset();

	for (auto& plate : plates)
	{
		updateMinMaxAvg(statistics.plates.tileCount, plate.tiles.size());
		plate.area = 0;
		for (auto ptile : plate.tiles)
		{
			auto& tile = *ptile;
			plate.area += tile.area;
		}
		updateMinMaxAvg(statistics.plates.area, plate.area);
		double elevation = 0;
		for (auto pcorner : plate.boundaryCorners)
		{
			auto& corner = *pcorner;
			elevation += corner.elevation;
		}
		updateMinMaxAvg(statistics.plates.boundaryElevation, elevation / plate.boundaryCorners.size());
		updateMinMaxAvg(statistics.plates.boundaryBorders, plate.boundaryBorders.size());
		plate.circumference = 0;
		for (auto pborder : plate.boundaryBorders)
		{
			auto& border = *pborder;
			plate.circumference += border.length();
		}
		updateMinMaxAvg(statistics.plates.circumference, plate.circumference);
	}

	statistics.plates.tileCount.avg /= statistics.plates.count;
	statistics.plates.area.avg /= statistics.plates.count;
	statistics.plates.boundaryElevation.avg /= statistics.plates.count;
	statistics.plates.boundaryBorders.avg /= statistics.plates.count;
	statistics.plates.circumference.avg /= statistics.plates.count;
}

}



