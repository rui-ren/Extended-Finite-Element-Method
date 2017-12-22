#include "Geometry_info.h"
#include "Matrix.h"
#include <vector>

using namespace rui;

typedef std::valarray<double> Vector;

Point::Point()
{
	this->x = 0;
	this->y = 0;
	this->z = 0;
	this->id = -1;
}

double Point::angle() const
{
	assert(z == 0);
	return atan2(y, x);
}

Point::Point(double x, double y)
{
	//(*(dynamic_cast< valarray<double> *>(this))) ;
	this->x = x;
	this->y = y;
	this->z = 0;
	this->id = -1;
}

Point::Point(double x, double y, double z)
{
	//(*(dynamic_cast< valarray<double> *>(this))) ;
	this->x = x;
	this->y = y;
	this->z = z;
	this->id = -1;
}

void Point::print() const
{
	std::cout << " (" << x << ", " << y << ") ";
}

double Point::norm() const
{
	return sqrt(x*x + y*y + z*z);
}

double Point::sqNorm() const
{
	return x*x + y*y + z*z;
}
void Point::setX(double v)
{
	x = v;
}

void Point::setY(double v)
{
	y = v;
}

void Point::setZ(double v)
{
	z = v;
}

void Point::set(const Point p)
{
	x = p.x;
	y = p.y;
	z = p.z;
}

void Point::set(const Point * p)
{
	x = p->x;
	y = p->y;
	z = p->z;
}

void Point::set(double v, double vv)
{
	x = v;
	y = vv;
}

void Point::set(double v, double vv, double vvv)
{
	x = v;
	y = vv;
	z = vvv;
}

bool Point::operator==(Point p) const
{
	return  std::abs(x - p.x) < 2 * std::numeric_limits<double>::epsilon() &&
		std::abs(y - p.y) < 2 * std::numeric_limits<double>::epsilon() && std::abs(z - p.z) < 2 * std::numeric_limits<double>::epsilon();
}

bool Point::operator!=(Point p) const
{
	return   std::abs(x - p.x) >  2 * std::numeric_limits<double>::epsilon() ||
		std::abs(y - p.y) >  2 * std::numeric_limits<double>::epsilon() || std::fabs(z - p.z) >  2 * std::numeric_limits<double>::epsilon();
}

Point Point::operator-(Point p) const
{
	Point ret((*this));
	ret.x -= p.x;
	ret.y -= p.y;
	ret.z -= p.z;
	return ret;
}

Point Point::operator-(Vector p) const
{
	Point ret((*this));
	ret.x -= p[0];
	ret.y -= p[1];
	ret.z -= p[2];
	return ret;
}

Point Point::operator+(Point p) const
{
	Point ret((*this));
	ret.x += p.x;
	ret.y += p.y;
	ret.z += p.z;
	return ret;
}

void Point::operator+=(Point p)
{
	x += p.x;
	y += p.y;
	z += p.z;
}

Point Point::operator+(Vector p) const
{
	Point ret((*this));
	ret.x += p[0];
	ret.y += p[1];
	ret.z += p[2];
	return ret;
}

Point Point::operator/(double p) const
{
	Point ret((*this));
	ret.x /= p;
	ret.y /= p;
	ret.z /= p;
	return ret;
}


bool Point::operator <(Point p) const
{
	return (y < p.y) || ((y == p.y) && (x < p.x) || (y == p.y) && (x == p.x) && (z < p.z));
}

bool Point::operator >(Point p) const
{
	return (y > p.y) || ((y == p.y) && (x > p.x) || (y == p.y) && (x == p.x) && (z > p.z));
}

Point Point::operator*(double p)  const
{
	Point ret((*this));
	ret.x *= p;
	ret.y *= p;
	ret.z *= p;
	return ret;
}

double Point::operator*(Point p) const
{
	return x*p.x + y*p.y + z*p.z;
}

double Point::operator*(Vector p) const
{
	return x*p[0] + y*p[1] + z*p[2];
}

Point Point::operator^(Point p) const
{
	Point ret;
	ret.x = y*p.z - z*p.y;
	ret.y = z*p.x - x*p.z;
	ret.z = x*p.y - y*p.x;

	return ret;
}

Point Point::operator^(Vector p) const
{
	Point ret;
	ret.x = y*p[2] - z*p[1];
	ret.y = z*p[0] - x*p[2];
	ret.z = x*p[1] - y*p[0];

	return ret;
}

PointSet::PointSet()
{
	this->boundingPoints = new std::valarray<Point *>(0);
	this->chullEndPos = 0;
}


std::valarray<Point *> * NonConvexGeometry::getBoundingPoints() const
{
	return this->boundingPoints;
}

Point* NonConvexGeometry::getBoundingPoint(size_t i) const
{
	return (*this->boundingPoints)[i];
}

Point* Geometry::getInPoint(size_t i) const
{
	return (*this->inPoints)[i];
}

double & Point::operator[](size_t i)
{
	return (*(&x + i));
}
double Point::operator[](size_t i) const
{
	return (*(&x + i));
}

std::valarray<Point* > * Geometry::getInPoints() const
{
	return this->inPoints;
}


const Point* Geometry::getCenter() const
{
	return &center;
}

Point * Geometry::getCenter()
{
	return &center;
}

PointSet::PointSet(size_t npoints)
{
	this->boundingPoints = new std::valarray<Point *>(npoints);
	for (size_t i = 0; i < npoints; i++)
		(*this->boundingPoints)[i] = NULL;
	this->chullEndPos = 0;
};


size_t NonConvexGeometry::size() const
{
	return 	boundingPoints->size() + inPoints->size();
}

double PointSet::x(size_t i)
{
	return (*boundingPoints)[i]->x;
}

double PointSet::y(size_t i)
{
	return (*boundingPoints)[i]->y;
}

double PointSet::z(size_t i)
{
	return (*boundingPoints)[i]->z;
}

void PointSet::setX(size_t i, double vv)
{
	(*boundingPoints)[i]->Point::setX(vv);
}

void PointSet::setY(size_t i, double vv)
{
	(*boundingPoints)[i]->Point::setY(vv);
}

void PointSet::setZ(size_t i, double vv)
{
	(*boundingPoints)[i]->Point::setZ(vv);
}

void PointSet::set(size_t i, Point * p)
{
	delete (*boundingPoints)[i];
	(*boundingPoints)[i] = p;
}

void PointSet::set(size_t i, double x, double y)
{
	int id = -1;
	if ((*boundingPoints)[i] != NULL)
		id = (*boundingPoints)[i]->id;

	delete (*boundingPoints)[i];

	(*boundingPoints)[i] = new Point(x, y);
	(*boundingPoints)[i]->id = id;
}

void PointSet::set(size_t i, double x, double y, double z)
{
	int id = -1;
	if ((*boundingPoints)[i] != NULL)
		id = (*boundingPoints)[i]->id;

	delete (*boundingPoints)[i];

	(*boundingPoints)[i] = new Point(x, y, z);
	(*boundingPoints)[i]->id = id;
}

Point * PointSet::operator[](size_t i)
{
	return (*boundingPoints)[i];
}

Point * PointSet::operator[](size_t i) const
{
	return (*boundingPoints)[i];
}

Point * PointSet::getPoint(size_t i) const
{
	return (*boundingPoints)[i];
}

Point * PointSet::getPoint(size_t i)
{
	return (*boundingPoints)[i];
}

PointSet::iterator PointSet::begin() const
{
	return &(*boundingPoints)[0];
}

PointSet::iterator PointSet::end() const
{
	return  &(*boundingPoints)[boundingPoints->size()];
}



size_t PointSet::size() const
{
	return  boundingPoints->size();
}

Point * PointSet::computeCenter()
{
	Point * ret = new Point();
	for (size_t i = 0; i < boundingPoints->size(); i++)
	{
		ret->x += (*boundingPoints)[i]->x / boundingPoints->size();
		ret->y += (*boundingPoints)[i]->y / boundingPoints->size();
	}
	return ret;
}

void PointSet::removePoint(size_t index)
{
	std::valarray<Point *> n(size() - 1);
	//std::copy((*this)[0], (*this)[index-1], (*n)[0]) ;
	for (size_t i = 0; i < index; i++)
	{
		n[i] = (*boundingPoints)[i];
	}
	for (size_t i = index + 1; i < size(); i++)
	{
		n[i] = (*boundingPoints)[i];
	}

	boundingPoints->resize(n.size());
	(*boundingPoints) = n;
}


Geometry::Geometry() : gType(NULL_GEOMETRY)
{
	sampled = false;
	inPoints = new std::valarray<Point * >(0);
}

Geometry::Geometry(size_t numPoints) : gType(NULL_GEOMETRY)
{
	inPoints = new std::valarray<Point * >(numPoints);
}

GeometryType Geometry::getGeometryType() const
{
	return gType;
}

// ConvexPolygon::ConvexPolygon()
// {
// 	std::cout << "calling default ConvexPolygon constructor" << std::endl ;
// }


NonConvexGeometry::NonConvexGeometry() : PointSet(1)
{
	orderedSet.resize(1);
	orderedSet[0] = (*this->inPoints)[0];
}

NonConvexGeometry::NonConvexGeometry(size_t numPoints) : PointSet(numPoints)
{
	orderedSet.resize(numPoints);
	for (size_t i = 0; i < numPoints; i++)
		orderedSet[i] = (*this->inPoints)[i];
}

NonConvexGeometry::NonConvexGeometry(PointSet * p)
{
	this->boundingPoints = new std::valarray<Point *>(p->size());
	std::copy(p->begin(), p->end(), begin());//&boundingPoints[0]) ;
}

Point * NonConvexGeometry::getPoint(size_t i) const
{
	if (i < inPoints->size())
		return (*inPoints)[i];
	return (*boundingPoints)[i - inPoints->size()];
}

Line::Line(Point origin, Point vector)
{
	p = origin;
	v = vector;
}

bool Line::intersects(const Line *l) const
{
	return v.x * l->vector()->y - v.y * l->vector()->x != 0;
}

const Point * Line::vector() const
{
	return dynamic_cast<const Point *>(&v);
}

Point Line::intersection(const Line *l) const
{
	double t = 0;
	if (v.x != 0 && v.y != 0)
		t = ((p.y - l->origin()->y) + v.y / v.x * (p.x - l->origin()->x)) / (v.x / v.y - l->vector()->y);
	else if (v.x == 0)
		t = (p.x - l->origin()->x) / l->vector()->x;
	else if (v.y == 0)
		t = (p.y - l->origin()->y) / l->vector()->y;

	return (*l->origin()) + (*l->vector())*t;
}


Point Line::projection(const Point *m) const
{
	Line l((*m), Point(-v.y, v.x));
	return l.intersection(this);
}


Segment::Segment(const Point p0, const Point p1)
{
	f = p0;
	s = p1;
	mid = p0*0.5 + p1*0.5;
	vec = p1 - p0;
}


Segment::Segment()
{
	f = Point(0, 0);
	s = Point(0, 0);
	mid = f*0.5 + s*0.5;
	vec = f - s;
}
Segment::~Segment()
{
}


void Segment::print() const
{
	std::cout << "[ (" << f.x << ", " << f.y << ") ; (" << s.x << ", " << s.y << ") ]" << std::endl;
}


void Segment::setFirst(Point p)
{
	f = p;
	mid = f*0.5 + s*0.5;
	vec = s - f;
}

void Segment::setFirst(double x, double y)
{
	f.set(x, y);
	mid = f*0.5 + s*0.5;
	vec = s - f;
}

void Segment::setSecond(Point p)
{
	s = p;
	mid = f*0.5 + s*0.5;
	vec = s - f;
}

void Segment::setSecond(double x, double y)
{
	s.set(x, y);
	mid = f*0.5 + s*0.5;
	vec = s - f;
}

void Segment::set(Point p0, Point p1)
{
	f.set(p0);
	s.set(p1);
	mid = f*0.5 + s*0.5;
	vec = s - f;
}

void Segment::set(double x0, double y0, double x1, double y1)
{
	f.set(x0, y0);
	s.set(x1, y1);
	mid = f*0.5 + s*0.5;
	vec = s - f;
}

const Point * Segment::first() const
{
	return &f;
}

const Point * Segment::second() const
{
	return &s;
}

Point * Segment::first()
{
	return &f;
}

Point * Segment::second()
{
	return &s;
}

bool isInTriangle(const Point test, const Point p0, const Point p1, const Point p2)
{
	return isOnTheSameSide(test, p0, p1, p2) && isOnTheSameSide(test, p1, p0, p2) && isOnTheSameSide(test, p2, p1, p2);
}

bool isOnTheSameSide(const Point test, const Point witness, const Point f0, const Point f1)
{
	return (((f1.x - f0.x)*(test.y - f0.y) - (f1.y - f0.y)*(test.x - f0.x)) *
		((f1.x - f0.x)*(witness.y - f0.y) - (f1.y - f0.y)*(witness.x - f0.x)) > -2 * std::numeric_limits<double>::epsilon());
}

double dist(const Point v1, const Point v2)
{
	return sqrt((v2.x - v1.x)*(v2.x - v1.x) + (v2.y - v1.y)*(v2.y - v1.y));
}

double squareDist(const Point v1, const Point v2)
{
	return (v2.x - v1.x)*(v2.x - v1.x) + (v2.y - v1.y)*(v2.y - v1.y);
}

double squareDist(const Point *v1, const Point *v2)
{
	return (v2->x - v1->x)*(v2->x - v1->x) + (v2->y - v1->y)*(v2->y - v1->y);
}

