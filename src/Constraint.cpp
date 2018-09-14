#include "Constraint.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Constraint::Constraint() {

}

Constraint::Constraint(int _nconEM, int _nconER, int _nconIM, int _nconIR) :
nconEM(_nconEM),
nconER(_nconER),
nconIM(_nconIM),
nconIR(_nconIR),
activeM(false),
activeR(false)
{

}

void Constraint::update() {

}

void Constraint::draw() {

}

void Constraint::countDofs(int &nem, int &ner, int &nim, int &nir) {
	// Counts DOFs
	idxEM = nem;
	idxER = ner;
	idxIM = nim;
	idxIR = nir;

	nem += nconEM;
	ner += nconER;
	nim += nconIM;
	nir += nconIR;
}

void Constraint::computeJacEqM(MatrixXd &Gm, MatrixXd &Gmdot, VectorXd &gm) {

	computeJacEqM_(Gm, Gmdot, gm);
	if (next != nullptr) {
		next->computeJacEqM(Gm, Gmdot, gm);
	}
}


void Constraint::computeJacEqM_(MatrixXd &Gm, MatrixXd &Gmdot, VectorXd &gm) {

}

void Constraint::computeJacEqR(MatrixXd &Gr, MatrixXd &Grdot, VectorXd &gr) {

	computeJacEqM_(Gr, Grdot, gr);
	if (next != nullptr) {
		next->computeJacEqM(Gr, Grdot, gr);
	}
}


void Constraint::computeJacEqR_(MatrixXd &Gr, MatrixXd &Grdot, VectorXd &gr) {


}

void Constraint::computeJacIneqM(Eigen::MatrixXd &Cm, Eigen::MatrixXd &Cmdot, Eigen::VectorXd &cm) {
	computeJacIneqM_(Cm, Cmdot, cm);
	if (next != nullptr) {
		next->computeJacIneqM(Cm, Cmdot, cm);
	}
}

void Constraint::computeJacIneqM_(Eigen::MatrixXd &Cm, Eigen::MatrixXd &Cmdot, Eigen::VectorXd &cm) {


}

void Constraint::computeJacIneqR(Eigen::MatrixXd &Cr, Eigen::MatrixXd &Crdot, Eigen::VectorXd &cr) {
	computeJacIneqR_(Cr, Crdot, cr);
	if (next != nullptr) {
		next->computeJacIneqR(Cr, Crdot, cr);
	}

}

void Constraint::computeJacIneqR_(Eigen::MatrixXd &Cr, Eigen::MatrixXd &Crdot, Eigen::VectorXd &cr) {


}

Constraint::~Constraint() {

}