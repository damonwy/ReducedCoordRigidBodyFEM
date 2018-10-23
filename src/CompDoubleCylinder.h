#pragma once
#include "Comp.h"
#include "Node.h"
#include "Vector.h"

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
	double getRadiusA() { return m_rA; }
	double getRadiusB() { return m_rB; }
	std::shared_ptr<Node> getOriginA() { return m_OA; }
	std::shared_ptr<Node> getOriginB() { return m_OB; }
	std::shared_ptr<Vector> getZAxisA() { return m_ZA; }
	std::shared_ptr<Vector> getZAxisB() { return m_ZB; }
	void setZAxisA(std::shared_ptr<Vector> ZA) { m_ZA = ZA; }
	void setOriginA(std::shared_ptr<Node> OA) { m_OA = OA; }
	void setZAxisB(std::shared_ptr<Vector> ZB) { m_ZB = ZB; }
	void setOriginB(std::shared_ptr<Node> OB) { m_OB = OB; }

protected:
	double m_rA;
	double m_hA;

	double m_rB;
	double m_hB;

	std::shared_ptr<Body> m_parentA;
	std::shared_ptr<Body> m_parentB;

	std::shared_ptr<Vector> m_ZA;	//Z axis
	std::shared_ptr<Node> m_OA;      // Origin
	std::shared_ptr<Vector> m_ZB;	//Z axis
	std::shared_ptr<Node> m_OB;      // Origin

	Eigen::Matrix4d E_wiA;	// Where the component is wrt world
	Eigen::Matrix4d E_jiA;	// Where the component is wrt body

	Eigen::Matrix4d E_wiB;	// Where the component is wrt world
	Eigen::Matrix4d E_jiB;	// Where the component is wrt body

	std::shared_ptr<Shape> m_shapeA;
	std::shared_ptr<Shape> m_shapeB;
};