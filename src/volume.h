
#pragma once

#include "initial3d.h"

namespace vol {

	//namespace initial3d2 = initial3d; //for asthetic purposes

	class Sphere;

	class Volume {
	private:
	public:
		virtual bool contains(const initial3d::vec3d &) = 0;
		virtual initial3d::vec3d getMin() = 0;
		virtual initial3d::vec3d getMax() = 0;
		virtual Sphere * getSphere() = 0;
		virtual ~Volume() { }
	};

	class Sphere : public Volume {
	private:
		initial3d::vec3d position;
		double radius;
	public:
		virtual bool contains(const initial3d::vec3d &p) {
			using namespace initial3d;
			return +(p - position) < radius;
		}
		virtual initial3d::vec3d getMin() {
			using namespace initial3d;
			return position - vec3d::one() * radius;
		}
		virtual initial3d::vec3d getMax() {
			using namespace initial3d;
			return position + vec3d::one() * radius;
		}
		virtual Sphere * getSphere() {
			return new Sphere(position, radius);
		}

		double getRadius() {return radius;}
		initial3d::vec3d getPosition() {return position;}

		Sphere(initial3d::vec3d p, double r) : position(p), radius(r) {}
	};

	class Plane : public Volume {
	private:
		initial3d::vec3d normal;
		double cutoff;
	public:
		virtual bool contains(const initial3d::vec3d &p) {
			using namespace initial3d;
			return p * normal > cutoff;//ben's black magic... DO NOT TOUCH!!!
		}
		virtual initial3d::vec3d getMin() {
			using namespace initial3d;
			return vec3d::one() * -math::inf<double>();
		}
		virtual initial3d::vec3d getMax() {
			using namespace initial3d;
			return vec3d::one() * math::inf<double>();
		}
		virtual Sphere * getSphere() {
			using namespace initial3d;
			return new Sphere(vec3d::zero(), math::inf<double>());
		}

		Plane(initial3d::vec3d p, initial3d::vec3d n) {
			using namespace initial3d;
			normal = ~n;
			cutoff = normal * p;
		}
	};

	class Box : public Volume {
	private:	
		initial3d::vec3d minimum;
		initial3d::vec3d maximum;

	public:
		virtual bool contains(const initial3d::vec3d &p) {
			using namespace initial3d;
			const vec3d &minimum = this->minimum;
			const vec3d &maximum = this->maximum;
			return p.x() >= minimum.x()
				&& p.y() >= minimum.y()
				&& p.z() >= minimum.z()
				&& p.x() <= maximum.x()
				&& p.y() <= maximum.y()
				&& p.z() <= maximum.z();
		}
		virtual initial3d::vec3d getMin() {
			return minimum;
		}
		virtual initial3d::vec3d getMax() {
			return maximum;
		}
		virtual Sphere * getSphere() {
			using namespace initial3d;
			double radius = +(maximum - minimum) * 0.5;
			return new Sphere(minimum + vec3d::one() * radius, radius);
		}

		Box(initial3d::vec3d min, initial3d::vec3d max) : minimum(min), maximum(max) {}
	};

	class Cylinder : public Volume {
	private:
		initial3d::vec3d position;
		double height;
		double bottom;
		double top;
	public:
		virtual bool contains(const initial3d::vec3d &p) {
			using namespace initial3d;
			vec3d pos = p - position;
			if (pos.z() < 0.0 || pos.z()> height )
				return false;
			double radius = bottom + (top - bottom) * (pos.z() / height);
			return radius >= sqrt(pos.x()*pos.x() + pos.y()*pos.y());
		}
		virtual initial3d::vec3d getMin() {
			using namespace initial3d;
			return vec3d(-std::max(bottom, top), -std::max(bottom, top), 0) + position;
		}
		virtual initial3d::vec3d getMax() {
			using namespace initial3d;
			return vec3d(std::max(bottom, top), std::max(bottom, top), height) + position;
		}
		virtual Sphere * getSphere() {
			using namespace initial3d;
			return new Sphere(vec3d(0, 0, height/2), + (getMax() - getMin()) * 0.5);
		}
		//specified where the top of the cone is the position
		Cylinder(initial3d::vec3d p, double h, double r) : position(p), height(h), bottom(r), top(r) {}
		Cylinder(initial3d::vec3d p, double h, double b, double t) : position(p), height(h), bottom(b), top(t) {}
	};

	class Cone : public Volume {
	private:
		initial3d::vec3d position;
		double height;
		double base;
	public:
		virtual bool contains(const initial3d::vec3d &p) {
			using namespace initial3d;
			vec3d pos = p - position;
			if (pos.z() < 0.0 || pos.z() > height)
				return false;
			double radius = base * (1-(pos.z() / height));
			return radius >= sqrt(pos.x()*pos.x() + pos.y()*pos.y());
		}
		virtual initial3d::vec3d getMin() {
			using namespace initial3d;
			return vec3d(-base, -base, 0) + position;
		}
		virtual initial3d::vec3d getMax() {
			using namespace initial3d;
			return vec3d(base, base, height) + position;
		}
		virtual Sphere * getSphere() {
			using namespace initial3d;
			return new Sphere(vec3d(0, 0, height/2), + (getMax() - getMin()) * 0.5);
		}
		Cone(initial3d::vec3d pos, double h, double r) : position(pos), height(h), base(r) {}
	};


	//Pair classes

	class Union : public Volume {
	private:
		Volume *left;
		Volume *right;
	public:
		virtual bool contains(const initial3d::vec3d &p) {
			using namespace initial3d;
			return left->contains(p) || right->contains(p);
		}
		virtual initial3d::vec3d getMin() {
			using namespace initial3d;
			return vec3d::negative_extremes(left->getMin(), right->getMin());
		}
		virtual initial3d::vec3d getMax() {
			using namespace initial3d;
			return vec3d::positive_extremes(left->getMax(), right->getMax());
		}
		virtual Sphere * getSphere() {
			using namespace initial3d;
			//TODO
			return NULL;
		}
		Union(Volume *l, Volume *r) : left(l), right(r) {}
	};

	class Intersection : public Volume {
	private:
		Volume *left;
		Volume *right;
	public:
		virtual bool contains(const initial3d::vec3d &p) {
			using namespace initial3d;
			return left->contains(p) && right->contains(p);
		}
		virtual initial3d::vec3d getMin() {
			using namespace initial3d;
			return vec3d::positive_extremes(left->getMin(), right->getMin());
		}
		virtual initial3d::vec3d getMax() {
			using namespace initial3d;
			return vec3d::negative_extremes(left->getMax(), right->getMax());
		}
		virtual Sphere * getSphere() {
			using namespace initial3d;
			double radius = +(getMax() - getMin()) * 0.5;
			return new Sphere(getMin() + vec3d::one() * radius, radius);
		}
		Intersection(Volume *l, Volume *r) : left(l), right(r) {}
	};

	class Subtraction : public Volume {
	private:
		Volume *left;
		Volume *right;
	public:
		virtual bool contains(const initial3d::vec3d &p) {
			using namespace initial3d;
			return left->contains(p) && !(right->contains(p));
		}
		virtual initial3d::vec3d getMin() {
			return left->getMin();
		}
		virtual initial3d::vec3d getMax() {
			return left->getMax();
		}
		virtual Sphere * getSphere() {
			using namespace initial3d;
			double radius = +(getMax() - getMin()) * 0.5;
			return new Sphere(getMin() + vec3d::one() * radius, radius);
		}
		Subtraction(Volume *l, Volume *r) : left(l), right(r) {}
	};


	//space manipulation classes

	class Translation : public Volume {
	private:
		Volume *volume;
		initial3d::vec3d translation;
	public:
		virtual bool contains(const initial3d::vec3d &p) {
			using namespace initial3d;
			return volume->contains(p - translation);
		}
		virtual initial3d::vec3d getMin() {
			using namespace initial3d;
			return volume->getMin() + translation;
		}
		virtual initial3d::vec3d getMax() {
			using namespace initial3d;
			return volume->getMax() + translation;
		}
		virtual Sphere * getSphere() {
			using namespace initial3d;
			double radius = +(getMax() - getMin()) * 0.5;
			return new Sphere(getMin() + vec3d::one() * radius, radius);
		}
		Translation(Volume *v, initial3d::vec3d trans) : translation(trans) {}
	};


	class Rotation : public Volume {
	private:
		Volume *volume;
		initial3d::quatd rotation;
		initial3d::vec3d center;
	public:
		virtual bool contains(const initial3d::vec3d &p) {
			using namespace initial3d;
			return volume->contains(((!rotation) * (p - center)) + center);
		}
		virtual initial3d::vec3d getMin() {
			using namespace initial3d;
			vec3d min = volume->getMin();
			vec3d max = volume->getMax();

			vec3d xyz = (rotation * (vec3d(min.x(), min.y(), min.z()) - center)) + center;
			vec3d Xyz = (rotation * (vec3d(max.x(), min.y(), min.z()) - center)) + center;
			vec3d xYz = (rotation * (vec3d(min.x(), max.y(), min.z()) - center)) + center;
			vec3d xyZ = (rotation * (vec3d(min.x(), min.y(), max.z()) - center)) + center;
			vec3d XYz = (rotation * (vec3d(max.x(), max.y(), min.z()) - center)) + center;
			vec3d XyZ = (rotation * (vec3d(max.x(), min.y(), max.z()) - center)) + center;
			vec3d xZY = (rotation * (vec3d(min.x(), max.y(), max.z()) - center)) + center;
			vec3d XZY = (rotation * (vec3d(max.x(), max.y(), max.z()) - center)) + center;

			return vec3d::negative_extremes(xyz,
					vec3d::negative_extremes(Xyz,
					vec3d::negative_extremes(xYz,
					vec3d::negative_extremes(xyZ,
					vec3d::negative_extremes(XYz,
					vec3d::negative_extremes(XyZ,
					vec3d::negative_extremes(xZY,XZY
				)))))));
		}
		virtual initial3d::vec3d getMax() {
			using namespace initial3d;
			vec3d min = volume->getMin();
			vec3d max = volume->getMax();

			vec3d xyz = (rotation * (vec3d(min.x(), min.y(), min.z()) - center)) + center;
			vec3d Xyz = (rotation * (vec3d(max.x(), min.y(), min.z()) - center)) + center;
			vec3d xYz = (rotation * (vec3d(min.x(), max.y(), min.z()) - center)) + center;
			vec3d xyZ = (rotation * (vec3d(min.x(), min.y(), max.z()) - center)) + center;
			vec3d XYz = (rotation * (vec3d(max.x(), max.y(), min.z()) - center)) + center;
			vec3d XyZ = (rotation * (vec3d(max.x(), min.y(), max.z()) - center)) + center;
			vec3d xZY = (rotation * (vec3d(min.x(), max.y(), max.z()) - center)) + center;
			vec3d XZY = (rotation * (vec3d(max.x(), max.y(), max.z()) - center)) + center;

			return vec3d::positive_extremes(xyz,
					vec3d::positive_extremes(Xyz,
					vec3d::positive_extremes(xYz,
					vec3d::positive_extremes(xyZ,
					vec3d::positive_extremes(XYz,
					vec3d::positive_extremes(XyZ,
					vec3d::positive_extremes(xZY,XZY
				)))))));
		}
		virtual Sphere * getSphere() {
			return new Sphere(rotation * volume->getSphere()->getPosition(), volume->getSphere()->getRadius());
		}
		Rotation(Volume *v, initial3d::quatd r) : volume(v), rotation(r), center(initial3d::vec3d::zero()) {}
		Rotation(Volume *v, initial3d::quatd r, initial3d::vec3d c) : volume(v), rotation(r), center(c) {}
	};
}