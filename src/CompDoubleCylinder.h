#pragma once
#include "Comp.h"

class CompDoubleCylinder : public Comp
{
public:
	CompDoubleCylinder();
	CompDoubleCylinder(std::shared_ptr<Body> parentA, double rA, std::shared_ptr<Body> parentB, double rB);
	virtual ~CompDoubleCylinder();

	void init();
	void update();
	void load(const std::string &RESOURCE_DIR, std::string shapeA, std::string shapeB);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P)const;
	void setTransformA(Eigen::Matrix4d E);
	void setTransformB(Eigen::Matrix4d E);

protected:
	double m_rA;
	double m_hA;

	double m_rB;
	double m_hB;

	std::shared_ptr<Body> m_parentA;
	std::shared_ptr<Body> m_parentB;

	Eigen::Matrix4d E_wiA;	// Where the component is wrt world
	Eigen::Matrix4d E_jiA;	// Where the component is wrt body

	Eigen::Matrix4d E_wiB;	// Where the component is wrt world
	Eigen::Matrix4d E_jiB;	// Where the component is wrt body

	std::shared_ptr<Shape> m_shapeA;
	std::shared_ptr<Shape> m_shapeB;
};