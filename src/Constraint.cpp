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


Constraint::~Constraint() {

}