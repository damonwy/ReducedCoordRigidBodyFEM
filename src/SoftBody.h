#pragma once
#ifndef MUSCLEMASS_SRC_SOFTBODY_H_
#define MUSCLEMASS_SRC_SOFTBODY_H_

#define EIGEN_USE_MKL_ALL

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "MLCommon.h"

class MatrixStack;
class Program;
class Node;
class Body;
class FaceTriangle;
class Tetrahedron;
class Vector;
class Line;

typedef Eigen::Triplet<double> T;
class SoftBody {

public:
	SoftBody();
	SoftBody(double density, double young, double poisson, Material material);
	virtual ~SoftBody() {}

	void insertAdditionalPoints() {}
	virtual void load(const std::string &RESOURCE_DIR, const std::string &MESH_NAME, const std::string &TETGEN_FLAGS);
	virtual void init();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	void updatePosNor();

	virtual void countDofs(int &nm, int &nr);
	virtual void computeJacobian(Eigen::MatrixXd &J);
	virtual void computeJacobianSparse(std::vector<T> &J_);
	virtual void computeMass(Eigen::MatrixXd &M);
	virtual void computeMassSparse(std::vector<T> &M_);
	virtual Energy computeEnergies(Vector3d grav, Energy ener);

	void computeForce(Vector3d grav, Eigen::VectorXd &f);
	void computeStiffness(Eigen::MatrixXd &K);
	void computeStiffnessSparse(std::vector<T> &K_);
	void computeForceDamping(Eigen::VectorXd &f, Eigen::MatrixXd &D);
	void computeForceDampingSparse(Eigen::VectorXd &f, std::vector<T> &D_);

	virtual void gatherDofs(Eigen::VectorXd &y, int nr);
	virtual Eigen::VectorXd gatherDDofs(Eigen::VectorXd &ydot, int nr);
	virtual void scatterDofs(Eigen::VectorXd &y, int nr);
	virtual void scatterDDofs(Eigen::VectorXd &ydot, int nr);

	inline void setColor(Vector3f color) { m_color = color; }
	void setAttachments(int id, std::shared_ptr<Body> body);
	void setAttachmentsByLine(Vector3d dir, Vector3d orig, std::shared_ptr<Body> body);
	void setAttachmentsByLine(std::shared_ptr<Line> l, std::shared_ptr<Body> body);
	void setAttachmentsByLine(std::shared_ptr<Line> l);
	void setAttachmentsByXYSurface(double z, double range, Vector2d xrange, Vector2d yrange, std::shared_ptr<Body> body);
	void setAttachmentsByYZSurface(double x, double range, Vector2d yrange, Vector2d zrange, std::shared_ptr<Body> body);
	void setAttachmentsByXZSurface(double y, double range, Vector2d xrange, Vector2d zrange, std::shared_ptr<Body> body);

	void setAttachmentsByXZCircle(double y, double range, Vector2d O, double r, std::shared_ptr<Body> body);
	void setAttachmentsByYZCircle(double x, double range, Vector2d O, double r, std::shared_ptr<Body> body);

	void setSlidingNodes(int id, std::shared_ptr<Body> body, Vector3d init_dir);
	void setSlidingNodesByXYSurface(double z, Vector2d xrange, Vector2d yrange, double dir, std::shared_ptr<Body> body);
	void setSlidingNodesByYZSurface(double x, Vector2d yrange, Vector2d zrange, double dir, std::shared_ptr<Body> body);
	void setSlidingNodesByXZSurface(double y, Vector2d xrange, Vector2d zrange, double dir, std::shared_ptr<Body> body);
	void setSlidingNodesByYZCircle(double x, double range_x, Vector2d O, double r, std::shared_ptr<Body> body);

	inline void setInvertiblity(bool isInvertible) { m_isInvertible = isInvertible; }
	inline void setFloor(double floor_y) { m_floor_y = floor_y; m_isCollisionWithFloor = true; }
	inline void setYthreshold(double collision_y) { m_collision_y = collision_y; }
	void setDamping(double damping) { m_damping = damping; }
	inline bool getInvertiblity() { return m_isInvertible; }

	inline const std::vector<std::shared_ptr<Node> > & getNodes() const { return m_nodes; }
	inline const std::vector<std::shared_ptr<Tetrahedron> > & getTets() const { return m_tets; }
	inline const std::vector<std::shared_ptr<FaceTriangle> > & getFaces() const { return m_trifaces; }

	void transform(Vector3d dx);
	void transform(Matrix4d E);
	// attached 
	std::vector<std::shared_ptr<Node> > m_attach_nodes;
	std::vector<std::shared_ptr<Body> > m_attach_bodies;
	std::vector<Vector3d> m_r;

	// sliding
	std::vector<std::shared_ptr<Node> > m_sliding_nodes;
	std::vector<std::shared_ptr<Body> > m_sliding_bodies;
	std::vector<Vector3d> m_r_sliding;
	std::vector<std::shared_ptr<Vector>> m_normals_sliding;

	// nodes for comparing to sliding
	std::vector<std::shared_ptr<Node> > m_compared_nodes;

	std::shared_ptr<SoftBody> next;
	std::vector<std::shared_ptr<FaceTriangle> > m_trifaces;
	bool m_isInverted;
	
	int m_npotentialcols;		// how many nodes need to check


protected:
	int m_type;
	bool m_isInvertible;
	bool m_isGravity;
	bool m_isElasticForce;
	bool m_isCollisionWithFloor;
	Material m_material;
	Vector3f m_color;

	double m_floor_y;
	double m_collision_y;	// nodes on the mesh surface and whose y(1) below this, need to do collision detection

	std::vector<std::shared_ptr<Node> > m_nodes;	
	std::vector<std::shared_ptr<Tetrahedron> > m_tets;

	std::vector<unsigned int> eleBuf;
	std::vector<float> posBuf;
	std::vector<float> norBuf;
	std::vector<float> texBuf;

	unsigned eleBufID;
	unsigned posBufID;
	unsigned norBufID;
	unsigned texBufID;

	double m_damping;
	double m_young;
	double m_poisson;
	double m_density;
	double m_mass;

	virtual void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	virtual void computeStiffnessSparse_(std::vector<T> &K_);
	virtual void computeStiffness_(Eigen::MatrixXd &K);
	virtual void computeForce_(Vector3d grav, Eigen::VectorXd &f);
	
};


#endif // MUSCLEMASS_SRC_SOFTBODY_H_