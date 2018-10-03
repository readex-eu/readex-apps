
#ifndef SRC_BASIS_MATRICES_MATRIXTYPE_H_
#define SRC_BASIS_MATRICES_MATRIXTYPE_H_

namespace espreso {

enum class MatrixType : int {
	REAL_SYMMETRIC_POSITIVE_DEFINITE = 0,
	REAL_SYMMETRIC_INDEFINITE = 1,
	REAL_UNSYMMETRIC = 2
};

}


#endif /* SRC_BASIS_MATRICES_MATRIXTYPE_H_ */
