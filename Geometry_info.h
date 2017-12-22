#ifndef GEOMETRY_INFO_H
#define GEOMETRY_INFO_H

#include<math.h>
#include<valarray>
#include<vector>
#include<map>
#include<iostream>
#include<limits>
#include<assert.h>

using namespace std;
typedef std::valarray<double> Vector;

// creat a package for the geometry
namespace rui
{
	typedef enum GeometryType
	{
		NULL_GEOMETRY,
		RECTANGLE,
		SEGMENTED_LINE,
	};

	class Segment;
	class ConvexPolygon;

	class Point
	{
	public:
		double x;
		double y;
		double z;
		int id;
		Point();
		Point(double x, double y);
		Point(double x, double y, double z);

		void setX(double W);
		void setY(double W);
		void setZ(double W);
		void set(double W, double WW);				// Set two demension
		void set(double W, double WW, double WWW); // Set three demension
		void set(const Point p);
		void set(const Point *P);

		//////////////////////////////define bool identification for points////////////
		bool operator ==(Point P) const;
		bool operator !=(Point P) const;

		bool operator < (Point P) const;
		bool operator > (Point P) const;

		//////////////////////////// define operator for point////////////////////////
		Point operator -(Point P) const;
		Point operator-(Vector P) const;
		Point operator +(Point P) const;
		Point operator+(Vector P) const;
		Point operator/(double P) const;
		Point operator*(double P) const;
		double operator*(Point P) const;
		double operator*(Vector P) const;
		Point operator ^(Point P) const;
		Point operator^(Vector P) const;

		void operator+=(Point P);
		double norm()  const;    //Set MATLAB norm() function
		double sqNorm() const;
		void print() const;
		double angle() const;
		double & operator[](size_t i);
		double operator[](size_t i) const;
	};

	////////////////////define line in the geometry///////////

	class Line
	{
	protected:
		Point p;
		Point v;
	public:
		Line(Point origin, Point vector);

		bool intersects(const Line *l) const;
		bool intersects(const Segment *s) const;
		bool intersects(const Geometry *g) const;
		bool on(const Point *P) const;

		vector<Point> intersection(const Geometry * g) const;
		Point intersection(const Line *l) const;
		Point intersection(const Segment *l) const;

		const Point * vector() const;
		const Point * origin() const;

		Point projection(const Point *P) const;
	};


	//////////////////////Segment of the domain ///////////////////

	class Segment
	{
	protected:
		Point f;
		Point s;
		Point mid;
		Point vec;

	public:
		Segment(const Point P0, const Point P1);
		Segment();
		virtual ~Segment();

		bool intersections(const Line *l) const;
		bool intersections(const Segment *s) const;
		bool intersections(const Geometry *g) const;
		bool on(const Point *P) const;

		void setFirst(Point P);
		void setFirst(double x, double y);
		void setSecond(Point P);
		void setSecond(double x, double y);
		void set(Point P0, Point P1);
		void set(double x0, double y0, double x1, double y1);

		void print() const;

		const Point * first() const;
		const Point * second() const;

		Point * first();
		Point * second();

		Point intersection(const Line *l) const;
		Point intersection(const Segment *l) const;
		vector<Point> intersection(const Geometry * g) const;
	};

	/////////////////////// Geometry of the domain//////////////////
	class Geometry
	{
	protected:
		std::valarray<Point *> *inPoints;
		bool sampled;

		Point center;

		virtual void computeCemter() = 0;

		GeometryType gType;

	public:

		Geometry();
		Geometry(size_t numPoints);
		virtual ~Geometry() { delete inPoints; }

		virtual std::valarray<Point *> * getBoundingPoints() const = 0;
		virtual Point* getBoundingPoint(size_t i) const = 0;
		virtual std::valarray<Point *> * getInPoints() const;
		virtual Point* getInPoint(size_t i) const;
		virtual Point * getPoint(size_t i) const = 0;
		virtual const Point * getCenter() const;
		virtual Point * getCenter();
		virtual void project(Point *) const = 0;

		GeometryType getGeometryType() const;

		virtual double getRadius() const = 0;

		virtual void sampleBoundingSurface(size_t num_points) = 0;
		virtual void sampleSurface(size_t num_points) = 0;
		virtual bool in(const Point p) = 0;

		virtual bool in(const Point *p) const{ return false; }

		virtual size_t size() const = 0;
		virtual double area() const = 0;
		virtual size_t sides() const { return 3; }
		virtual bool is1D() const = 0;
	};


	class PointSet
	{
	protected:
		valarray<Point *>* boundingPoints;
		size_t chullEndPos;

	public:
		PointSet();
		PointSet(size_t npoints);

		virtual ~PointSet() { delete boundingPoints; }
		double x(size_t i);
		double y(size_t i);
		double z(size_t i);

		void setX(size_t i, double v);
		void setY(size_t i, double v);
		void setZ(size_t i, double v);

		void set(size_t i, Point *p);
		void set(size_t i, double x, double y);
		void set(size_t i, double x, double y, double z);

		Point * operator [](size_t i);
		Point * operator [](size_t i) const;
		Point * getPoint(size_t i) const;
		Point * getPoint(size_t i);

		virtual bool in(const Point p);
		virtual size_t size() const;

		void removePoint(size_t index);
		Point * computeCenter();

		ConvexPolygon * convexHull();

		typedef Point** iterator;

		iterator begin() const;
		iterator end() const;

	};

	class NonConvexGeometry : public PointSet, public Geometry
	{
	protected:
		valarray<Point *> orderedSet;
		vector<size_t> stopPos;

	public:
		NonConvexGeometry();				 // constructor
		NonConvexGeometry(size_t numPoints);
		NonConvexGeometry(PointSet * p);
		virtual ~NonConvexGeometry(){ };	//destructor

		virtual valarray<Point *> * getBoundingPoints() const;
		virtual Point* getBoundingPoint(size_t i) const;
		virtual Point* getPoint(size_t i) const;
		virtual size_t size() const;
		virtual double area() const = 0;

		virtual void project(Point *) const = 0;

	};
}
#endif
