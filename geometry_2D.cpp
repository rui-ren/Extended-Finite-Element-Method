#include "geometry_2D.h"
//reference:  stephane Bordas, phu vinh Nguyen, etc

// Define pi in C++
#ifndef M_PI    
#define M_PI     3.14159265358
using namespace rui;


Rectangle::Rectangle(double x, double y, double originX, double originY) : ConvexGeometry(4), size_y(y), size_x(x)
{
	gType = RECTANGLE;
	this->center = Point(originX, originY);
	(*boundingPoints)[0] = new Point(originX - 0.5*x, originY + 0.5*y);
	(*boundingPoints)[1] = new Point(originX - 0.5*x, originY - 0.5*y);
	(*boundingPoints)[2] = new Point(originX + 0.5*x, originY - 0.5*y);
	(*boundingPoints)[3] = new Point(originX + 0.5*x, originY + 0.5*y);
}

Rectangle::Rectangle(double x, double y, Point &center) : ConvexGeometry(4), size_y(y), size_x(x)
{
	gType = RECTANGLE;
	this->center = center;
	(*boundingPoints)[0] = new Point(center.x - 0.5*x, center.y + 0.5*y);
	(*boundingPoints)[1] = new Point(center.x - 0.5*x, center.y - 0.5*y);
	(*boundingPoints)[2] = new Point(center.x + 0.5*x, center.y - 0.5*y);
	(*boundingPoints)[3] = new Point(center.x + 0.5*x, center.y + 0.5*y);
}

Rectangle::Rectangle() : ConvexGeometry(4), size_y(2), size_x(2)
{
	gType = RECTANGLE;
	this->center = Point(0, 0);
	(*boundingPoints)[0] = new Point(-1, 1);
	(*boundingPoints)[1] = new Point(-1, -1);
	(*boundingPoints)[2] = new Point(1, -1);
	(*boundingPoints)[3] = new Point(1, 1);
}


void Rectangle::computeCenter()
{
	for (size_t i = 0; i < this->size(); i++)
		this->center += *this->getPoint(i);

	this->center = this->center / this->size();
}

double  Rectangle::getRadius() const
{
	return sqrt(width()*width()*0.25 + height()*height()*0.25);
}

bool Rectangle::in(const Point p)
{
	if (p.x < getCenter()->x - 0.5*size_x)
		return false;
	if (p.x > getCenter()->x + 0.5*size_x)
		return false;
	if (p.y > getCenter()->y + 0.5*size_y)
		return false;
	if (p.y < getCenter()->y - 0.5*size_y)
		return false;

	return true;
}
bool Rectangle::in(const Point * p) const
{
	if (p->x < getCenter()->x - 0.5*size_x)
		return false;
	if (p->x > getCenter()->x + 0.5*size_x)
		return false;
	if (p->y > getCenter()->y + 0.5*size_y)
		return false;
	if (p->y < getCenter()->y - 0.5*size_y)
		return false;

	return true;
}

double Rectangle::area() const
{
	return this->size_x*this->size_y;
}

bool Rectangle::is1D() const
{
	return false;
}

double Rectangle::width() const
{
	return size_x;
}

double Rectangle::height() const
{
	return size_y;
}

void Rectangle::project(Point * p) const
{

	Line l(*getCenter(), *p);

	for (size_t i = 0; i < boundingPoints->size(); i++)
	{
		Segment test(*getBoundingPoint(i), *getBoundingPoint((i + 1) % boundingPoints->size()));

		if (l.intersects(&test))
		{
			Point inter = l.intersection(&test);
			p->x = inter.x;
			p->y = inter.y;
			return;
		}
	}
}

void Rectangle::sampleBoundingSurface(size_t num_points)
{
	assert(num_points % 4 == 0);
	double perimeter = 2 * (size_x + size_y);

	double distanceBetweenPoints = perimeter / num_points;

	this->numberOfPointsAlongX = static_cast<size_t>(std::ceil(size_x / distanceBetweenPoints) + 1);
	double distanceBetweenPointsAlongX = size_x / (this->numberOfPointsAlongX - 1);

	this->numberOfPointsAlongY = static_cast<size_t>(std::ceil(size_y / distanceBetweenPoints) + 1);
	double distanceBetweenPointsAlongY = size_y / (this->numberOfPointsAlongY - 1);

	num_points = ((numberOfPointsAlongX)* 2 + (numberOfPointsAlongY)* 2 - 4);

	boundingPoints->resize(num_points);

	for (size_t i = 0; i < numberOfPointsAlongY; i++)
	{
		(*boundingPoints)[i] = new Point(center.x - 0.5*size_x, center.y + 0.5*size_y - i*distanceBetweenPointsAlongY);
	}
	for (size_t i = 1; i < numberOfPointsAlongX; i++)
	{
		(*boundingPoints)[numberOfPointsAlongY + i - 1] = new Point(center.x - 0.5*size_x + i*distanceBetweenPointsAlongX,
			getCenter()->y - 0.5*size_y);
	}
	for (size_t i = 1; i < numberOfPointsAlongY; i++)
	{
		(*boundingPoints)[numberOfPointsAlongX + numberOfPointsAlongY + i - 2] = new Point(center.x + 0.5*size_x,
			center.y - 0.5*size_y + i*distanceBetweenPointsAlongY);
	}
	for (size_t i = 1; i < numberOfPointsAlongX - 1; i++)
	{
		assert(2 * numberOfPointsAlongY + numberOfPointsAlongX + i - 3< num_points);
		(*boundingPoints)[2 * numberOfPointsAlongY + numberOfPointsAlongX + i - 3] = new Point(
			center.x + 0.5*size_x - i*distanceBetweenPointsAlongX,
			center.y + 0.5*size_y);

	}
}

void Rectangle::sampleSurface(size_t num_points)
{
	if (this->size() == 4)
	{
		this->Rectangle::sampleBoundingSurface(num_points);
	}

	size_t nip = static_cast<size_t>((numberOfPointsAlongX - 2)*(numberOfPointsAlongY - 2));

	inPoints->resize(nip);

	if (nip > 0)
	{
		for (size_t i = 0; i < this->numberOfPointsAlongX - 2; i++)
		{
			for (size_t j = 0; j < this->numberOfPointsAlongY - 2; j++)
			{
				(*inPoints)[i*(numberOfPointsAlongY - 2) + j] = new Point((*this->boundingPoints)[numberOfPointsAlongY + i]->x, (*this->boundingPoints)[j + 1]->y);
			}
		}
	}
}

SegmentedLine::SegmentedLine(PointSet * points) : NonConvexGeometry(1)
{
	gType = SEGMENTED_LINE;
	this->center = Point();
	boundingPoints->resize(points->size());
	std::copy(points->begin(), points->end(), &(*boundingPoints)[0]);
}

void SegmentedLine::computeCenter()
{
}

double SegmentedLine::getRadius() const
{
	//! \todo make it do something
	return 0;
}

Point * SegmentedLine::getHead() const
{
	return (*this->begin());
}

Point * SegmentedLine::getTail() const
{
	return (*this->end() - 1);
}

void SegmentedLine::sampleBoundingSurface(size_t num_points)
{
}

void SegmentedLine::project(Point *p) const
{
	std::map<double, Segment> m;
	for (size_t i = 0; i < boundingPoints->size() - 1; i++)
	{
		Segment s(*getBoundingPoint(i), *getBoundingPoint(i + 1));
		m[squareDist(p, s.midPoint())] = s;
	}
	Line l(*m.rbegin()->second.first(), *m.rbegin()->second.vector());
	Point proj = l.projection(p);
	p->x = proj.x;
	p->y = proj.y;
}

void SegmentedLine::sampleSurface(size_t num_points)
{
}

bool SegmentedLine::in(const Point v)
{
	for (size_t i = 0; i < boundingPoints->size(); i++)
	if ((*(*boundingPoints)[i]) == v)
		return true;
	return false;
}

bool SegmentedLine::is1D() const
{
	return true;
}
#endif