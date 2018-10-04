#pragma once
#ifndef MUSCLEMASS_SRC_SOFTBODY_H_
#define MUSCLEMASS_SRC_SOFTBODY_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class MatrixStack;
class Program;
class Node;
class Body;
class FaceTriangle;
class Tetrahedron;

class SoftBody {

public:
	SoftBody();
	SoftBody(double density, double young, double poisson);
	virtual ~SoftBody();

	virtual void load(const std::string &RESOURCE_DIR, const std::string &MESH_NAME);
	virtual void init();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;

	void updatePosNor();
	void computeEnergies(Eigen::Vector3d grav, double &T, double &V);

private:
	std::vector<std::shared_ptr<Node> > m_nodes;
	std::vector<std::shared_ptr<FaceTriangle> > m_trifaces;
	std::vector<std::shared_ptr<Tetrahedron> > m_tets;

	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;

	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;

	double m_young;
	double m_poisson;
	double m_density;

};


#endif // MUSCLEMASS_SRC_SOFTBODY_H_