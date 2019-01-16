#include "Scene.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Node.h"
#include "Shape.h"
#include "Program.h"
#include "Joint.h"
#include "MatlabDebug.h"
#include "Vector.h"
#include "JsonEigen.h"
#include "World.h"
#include "Solver.h"
#include "SolverDense.h"
#include "SolverSparse.h"
//#include "Spring.h"
#include "Deformable.h"
#include "DeformableSpring.h"

#include "SoftBody.h"
#include "MeshEmbedding.h"

#include "BrenderManager.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

#include <unsupported/Eigen/MatrixFunctions> // TODO: avoid using this later, write a func instead

BrenderManager *brender;

Scene::Scene() :
	t(0.0),
	h(1e-2),
	grav(0.0, 0.0, 0.0),
	time_step(0),
	drawHz(10)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{	
	//read a JSON file
	ifstream i(RESOURCE_DIR + "input.json");
	i >> js;
	i.close();

	// Units: meters, kilograms, seconds
	h = js["h"];
	Eigen::from_json(js["grav"], grav);
	drawHz = js["drawHz"];

	m_world = make_shared<World>(FINGERS);//_INVERTIBLE
	m_world->load(RESOURCE_DIR);

	//m_solver = make_shared<SolverDense>(m_world, REDMAX_EULER);
	m_solver = make_shared<SolverSparse>(m_world, REDMAX_EULER, PARDISO_LU);

	brender = BrenderManager::getInstance();
	brender->add(m_world);	
	brender->setExportDir("D:/Research/Muscles/Projects/ReducedCoordRigidBodyFEM/resources/hand/");
	
	m_world->export_part = 0;

	brender->exportBrender(t);
	m_world->export_part = 1;
}


void Scene::init()
{
	count = 0;
	m_world->init();
	VectorXd y0, y1;
	y0.resize(2 * m_world->nr);
	y0.setZero();
	
	//y1 = m_solver->dynamics(y0);
	y.resize(2 * m_world->nr);
	y.setZero();
	m_world->getJoint0()->reparam();
	m_world->getJoint0()->gatherDofs(y, m_world->nr);
	m_world->getDeformable0()->gatherDofs(y, m_world->nr);
	m_world->getSoftBody0()->gatherDofs(y, m_world->nr);
	m_world->getMeshEmbedding0()->gatherDofs(y, m_world->nr);
	//m_solution = m_solver->solve();
	//vec_to_file(m_solution->t, "t");
	//mat_to_file(m_solution->y, "y");

	//tk = m_solution->t(0);
	drawH = 1.0 / drawHz;
	search_idx = 0;

}

void Scene::reset()
{
	
}

void Scene::solve() {
	

}

int torend = 0;
void Scene::step()
{	
	//int n_steps = m_solution->getNsteps();
	t += 0.01;
	//int output_idx;
	//double s;
	VectorXd ys;

	y = m_solver->dynamics(y);
	//m_world->getJoint0()->reparam();
	//m_world->getJoint0()->gatherDofs(y, m_world->nr);
	m_world->update();
	m_world->incrementTime();
	count++;
	if (count == 99) {
		cout << count << endl;
	}
	
	//if(tk < m_solution->t(n_steps-1)) {
	//	m_solution->searchTime(tk, search_idx, output_idx, s);
	//	search_idx = output_idx;
	//	ys = (1 - s)* m_solution->y.row(output_idx) + s * m_solution->y.row(output_idx + 1);

	//	m_world->getJoint0()->scatterDofs(ys, m_world->nr);
	//	m_world->getSpring0()->scatterDofs(ys, m_world->nr);
	//	m_world->getSoftBody0()->scatterDofs(ys, m_world->nr);
	//	tk = tk + drawH;
	//}
	//else {
	//	// reset
	//	tk = m_solution->t(0);
	//}	
	//if (t > 130.0 && t < 150.0) 
	//{
	//	if (torend % 1 == 0) {
	//		//brender->exportBrender(t);
	//	}

	// }
	if (torend % 5 == 0) {
		brender->exportBrender(t);
	}
	torend++;
	if (t > 50.0) {
		m_world->export_part = 2;
		brender->exportBrender(t);
		exit(1);
	}
}

void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, const shared_ptr<Program> progSoft, shared_ptr<MatrixStack> P) const
{
	m_world->draw(MV, prog, progSimple, progSoft, P);
	
}

void Scene::toggleCoarseMesh() {
	bool isOn = m_world->getMeshEmbedding0()->m_isCoarseMesh;
	
	m_world->getMeshEmbedding0()->toggleDrawingCoarseMesh(!isOn);

}