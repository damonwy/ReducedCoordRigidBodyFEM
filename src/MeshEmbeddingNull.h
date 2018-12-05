#pragma once
// MeshEmbeddingNull 

#ifndef REDUCEDCOORD_SRC_MESHEMBEDDINGNULL_H_
#define REDUCEDCOORD_SRC_MESHEMBEDDINGNULL_H_

#include "MeshEmbedding.h"
class MeshEmbeddingNull : public MeshEmbedding {

public:
	MeshEmbeddingNull() : MeshEmbedding() {}
	virtual ~MeshEmbeddingNull() {}

};

#endif // REDUCEDCOORD_SRC_MESHEMBEDDINGNULL_H_