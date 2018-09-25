/*!
 * @file    Mesh.cpp
 * @author  Michal Merta 
 * @date    July 4, 2013
 * 
 */

#ifdef MESH_H
namespace bem4i {

template<class LO,class SC>
Mesh<LO, SC>::Mesh() {
}

template<class LO, class SC>
Mesh<LO,SC>::Mesh(const Mesh& orig) {
}

template<class LO, class SC>
Mesh<LO,SC>::~Mesh() {
}


// specialization for LO=int, SC=double
//
//template<>
//Mesh<int,std::complex<double> >::Mesh() {
//}
//
//template<>
//Mesh<int,std::complex<double> >::Mesh(const Mesh& orig) {
//}
//
//template<>
//Mesh<int,std::complex<double> >::~Mesh() {
//}

}

#endif
