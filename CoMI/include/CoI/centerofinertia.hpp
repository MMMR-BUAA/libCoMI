//
// Copyright (c) 2015-2019 CNRS INRIA
//

#ifndef __CoMI_CoI_centerofinertia_hpp__
#define __CoMI_CoI_centerofinertia_hpp__

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/centroidal.hpp"

namespace pinocchio
{

  ///
  /// \brief Computes the CoI dynamics.
  ///
  /// \tparam Scalar The scalar type.
  /// \tparam Options Eigen Alignment options.
  /// \tparam JointCollection Collection of Joint types.
  ///
  /// \param[in] model The model structure of the rigid body system.
  /// \param[in] data The data structure of the rigid body system.
  ///
  /// \returns The whole multibody force.
  ///
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  inline void computeCenterofInertiaDynamics(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_com);
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, typename ConfigVectorType, typename TangentVectorType>
  inline void computeCenterofInertiaDynamics(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_com,
                                    const Eigen::MatrixBase<ConfigVectorType> & q,
                                    const Eigen::MatrixBase<TangentVectorType> & v)
  {
    forwardKinematics(model,data,q,v);
    return computeCenterofInertiaDynamics(model,data,A_com,F_com);
  }
  
} // namespace pinocchio

/* --- Details -------------------------------------------------------------------- */
#include "CoMI/CoI/centerofinertia.hxx"

#endif // ifndef __CoMI_CoI_centerofinertia_hpp__
