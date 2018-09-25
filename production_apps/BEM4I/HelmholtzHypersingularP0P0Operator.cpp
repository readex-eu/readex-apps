/*!
 * @file    HelmholtzHypersingularP0P0Operator.cpp
 * @author  Jan Zapletal
 * @author  Michal Merta
 * @date    July 23, 2015
 */

#ifdef HELMHOLTZHYPERSINGULARP0P0OPERATOR_H

namespace bem4i {

template<class LO, class SC>
HelmholtzHypersingularP0P0Operator<LO, SC>::
HelmholtzHypersingularP0P0Operator( ) {
}

template<class LO, class SC>
HelmholtzHypersingularP0P0Operator<LO, SC>::
HelmholtzHypersingularP0P0Operator(
    const HelmholtzHypersingularP0P0Operator & orig
    ) {
}

template<class LO, class SC>
HelmholtzHypersingularP0P0Operator<LO, SC>::HelmholtzHypersingularP0P0Operator(
    LinearOperator< LO, SC > * H1,
    LinearOperator< LO, SC > * H2,
    SurfaceMesh3D<LO, SC> &mesh
    ) {

  this->H2 = H2;
  
  Vector<LO, SC> ones(mesh.getNElements());
  ones.setAll(1.0);
  this->diagElems1 = new Vector<LO, SC>(mesh.getNElements());
  this->diagElems2 = new Vector<LO, SC>(mesh.getNElements());

  H1->apply(ones, *diagElems1);
  H2->apply(ones, *diagElems2);

}

template<class LO, class SC>
HelmholtzHypersingularP0P0Operator<LO, SC>::
~HelmholtzHypersingularP0P0Operator( ) {
  delete diagElems1;
  delete diagElems2;
}

template<class LO, class SC>
void HelmholtzHypersingularP0P0Operator<LO, SC>::apply(
    const Vector<LO, SC> & x,
    Vector<LO, SC> & y,
    bool transA,
    SC alpha,
    SC beta
    ) {
  
  //Vector<LO, SC> tmp1(y.getLength());
  Vector<LO, SC> tmp2(y.getLength());
  Vector<LO, SC> tmp3(y.getLength());
  this->H2->apply(x, y, false, - (SCVT) 1.0*alpha, beta);
  
  for (LO i = 0 ; i< y.getLength(); i++) {
    tmp2.set(i, -alpha*x.get(i) * diagElems1->get(i));
  }
  for (LO i = 0 ; i< y.getLength(); i++) {
    tmp3.set(i, alpha*x.get(i) * diagElems2->get(i));
  }
  
  y.add(tmp2);
  y.add(tmp3);
  
}

}

#endif
