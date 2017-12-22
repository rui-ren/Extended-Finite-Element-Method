// define the geometry: Triangle. Circle, Segmentline
// the geometry of this file is used to define the discontinuity of the crack or void
// Rui Dec-3-2017
// reference:  stephane Bordas, phu vinh Nguyen, etc

#ifndef __GEOMETRY_2D_H_
#define __GEOMETRY_2D_H_

#include "geometry_info.h"
namespace rui
{

	class Rectangle    // define rectangle.
	{
	protected:
		double size_y;
		double size_x;
		size_t numberOfPointsAlongX;
		size_t numberOfPointsAlongY;

		virtual void computeCenter();

	public:
		Rectangle(double x, double y, double originX, double originY);
		Rectangle(double x, double y, Point &center);
		Rectangle();
		virtual ~Rectangle() { };

		virtual void sampleBoundingSurface(size_t num_points);
		virtual void sampleSurface(size_t num_points);

		virtual double width() const;
		virtual double height() const;
		virtual double area() const;
		virtual double getRadius() const;

		virtual size_t sides() const { return 4; }

		virtual bool is1D() const;
		virtual bool in(const Point p);
		virtual bool in(const Point *p) const;

		virtual void project(Point *) const;

	};

	class SegmentedLine : public NonConvexGeometry       
	{
	protected:
		virtual void computeCenter();
	public:
		SegmentedLine(PointSet * points);
		virtual ~SegmentedLine() { };

		virtual void sampleBoundingSurface(size_t num_points);

		virtual void sampleSurface(size_t num_points);

		virtual bool in(const Point v);

		virtual double area() const { return 0; }

		virtual Point * getHead() const;

		virtual Point * getTail() const;

		virtual bool is1D() const;

		virtual void project(Point *) const;

		virtual double getRadius() const;
	};
};

#endif
