#include "Comp.h"

#include <fstream>
#include <json.hpp>

#include "Joint.h"
#include "SE3.h"
#include "Shape.h"
#include "MatrixStack.h"
#include "Program.h"
#include "Body.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Comp::Comp() {

}

Comp:: ~Comp() {

}

void Comp::load(const std::string &RESOURCE_DIR, std::string shape) {
}


void Comp::init() {
}

void Comp::update() {
	
}

void Comp::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P)const {

}



