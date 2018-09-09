#include "Wrench.h"

#include <iostream>
#include <fstream>
#include <json.hpp>

#include "Body.h"
using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Wrench::Wrench(Body *body) {



}

Wrench::~Wrench() {

}


void Wrench::setWrenchLocal(Vector6d wrenchLocal) {


}

void Wrench::setWrenchWorld(Vector6d wrenchWorld) {


}