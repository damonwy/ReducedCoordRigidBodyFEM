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

void Constraint::getActiveList(std::vector<int> &listM, std::vector<int> &listR) {
	// Gets list of active inequality indices
	if (activeM) {
		listM.push_back(idxIM);
	}
	if (activeR) {
		listR.push_back(idxIR);
	}
	if (next != nullptr) {
		next->getActiveList(listM, listR);
	}
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

void Constraint::scatterForceEqM(Eigen::MatrixXd Gmt, Eigen::VectorXd lm) {
	if (nconEM > 0) {
		fcon = -Gmt.block(idxQ, idxEM, nQ, nconEM) * lm.segment(idxEM, nconEM);
	}
	else {
		fcon.resize(nQ);
		fcon.setZero();
	}
	scatterForceEqM_();
	if (next != nullptr) {
		next->scatterForceEqM(Gmt, lm);
	}
}

void Constraint::scatterForceEqR(Eigen::MatrixXd Grt, Eigen::VectorXd lr) {
	if (nconER > 0) {
		fcon = -Grt.block(idxQ, idxER, nQ, nconER) * lr.segment(idxER, nconER);
	}
	else {
		fcon.resize(nQ);
		fcon.setZero();
	}
	scatterForceEqR_();
	if (next != nullptr) {
		next->scatterForceEqR(Grt, lr);
	}
}

void Constraint::scatterForceIneqR(Eigen::MatrixXd Crt, Eigen::VectorXd lr) {
	if (nconIR > 0) {
		fcon = -Crt.block(idxQ, idxIR, nQ, nconIR) * lr.segment(idxIR, nconIR);
	}
	else {
		fcon.resize(nQ);
		fcon.setZero();
	}
	scatterForceIneqR_();
	if (next != nullptr) {
		next->scatterForceIneqR(Crt, lr);
	}
}

void Constraint::scatterForceIneqM(Eigen::MatrixXd Cmt, Eigen::VectorXd lm) {
	if (nconIM > 0) {
		fcon = -Cmt.block(idxQ, idxIM, nQ, nconIM) * lm.segment(idxEM, nconIM);
	}
	else {
		fcon.resize(nQ);
		fcon.setZero();
	}
	scatterForceIneqM_();
	if (next != nullptr) {
		next->scatterForceIneqM(Cmt, lm);
	}
}

void Constraint::scatterForceEqM_() {


}

void Constraint::scatterForceEqR_() {


}

void Constraint::scatterForceIneqR_() {


}

void Constraint::scatterForceIneqM_() {


}

Constraint::~Constraint() {

}