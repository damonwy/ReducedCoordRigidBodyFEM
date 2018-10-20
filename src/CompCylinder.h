#pragma once
#include "Comp.h"

class CompCylinder : public Comp
{
public:
	CompCylinder();
	CompCylinder(std::shared_ptr<Body> parent, double r, double h);
	virtual ~CompCylinder();

	void load(const std::string &RESOURCE_DIR, std::string shape);
	void init();
	void update();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P)const;
	void setTransform(Eigen::Matrix4d E);

protected:
	double m_r;
	double m_h;

	Eigen::Matrix4d E_wi;	// Where the component is wrt world
	Eigen::Matrix4d E_ji;	// Where the component is wrt body

	std::shared_ptr<Shape> m_shape;
	std::shared_ptr<Body> m_parent;
};