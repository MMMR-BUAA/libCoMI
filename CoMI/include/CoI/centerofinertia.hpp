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
  /// \brief Computes the nonlinear force of the CoI dynamics wrt. to the world system or floating body.
  ///
  /// \tparam Scalar: The scalar type.
  /// \tparam Options: Eigen Alignment options.
  /// \tparam JointCollection: Collection of Joint types.
  /// \tparam ConfigVectorType: Config Vector of Joint  types.
  /// \tparam TangentVectorType: Tangent Vector of Joint Velocity types.
  /// \tparam gra_flag: 1 is considered gravity
  /// \tparam floating: 1 is wrt. to the floating body(the first rigid), and the first joint is floating!.
  ///
  /// \param[in] model:The model structure of the rigid body system.
  /// \param[in] data:The data structure of the rigid body system.
  /// \param[in] q:The data of joint.
  /// \param[in] v:The data of joint velocity.
  /// \param[out] nlf_com:The nonlinear force of the CoI dynamics wrt. to the world system.
  ///
  template<typename Scalar, int Options,template<typename,int> class JointCollectionTpl, 
           int gra_flag = 1>
  inline void computeCenterofInertiaDynamics_nl(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com);                                  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, 
           typename ConfigVectorType, typename TangentVectorType,
           int gra_flag = 1>
  inline void computeCenterofInertiaDynamics_nl(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const Eigen::MatrixBase<ConfigVectorType> & q,
                                    const Eigen::MatrixBase<TangentVectorType> & v,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com)
  {
    forwardKinematics(model,data,q,v);
    computeCenterofInertiaDynamics_nl<Scalar,Options,JointCollectionTpl,gra_flag>(model,data,nlf_com);
  }


  ///
  /// \brief Computes the forward CoI dynamics wrt. to the world system.
  ///
  /// \tparam Scalar: The scalar type.
  /// \tparam Options: Eigen Alignment options.
  /// \tparam JointCollection: Collection of Joint types.
  /// \tparam ConfigVectorType: Config Vector of Joint  types.
  /// \tparam TangentVectorType: Tangent Vector of Joint Velocity types.
  /// \tparam gra_flag: 1 is considered gravity 
  /// \tparam floating: 1 is wrt. to the floating body(the first rigid), and the first joint is floating!.
  ///
  /// \param[in] model:The model structure of the rigid body system.
  /// \param[in] data:The data structure of the rigid body system.
  /// \param[in] nlf_com:The nonlinear force of the CoI dynamics wrt. to the world system.
  /// \param[in] F_com:The whole body force of the CoI dynamics wrt. to the world system.
  /// \param[in] q:The data of joint.
  /// \param[in] v:The data of joint velocity.
  /// \param[out] A_com:The whole body accelation of the CoI dynamics wrt. to the world system.  
  ///
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl,
           int floating_flag = 1>
  inline void computeCenterofInertiaDynamics_forward(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_com,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com);

  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, 
           typename ConfigVectorType, typename TangentVectorType,
           int gra_flag = 1, int floating_flag = 1>
  inline void computeCenterofInertiaDynamics_forward(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_com,
                                    const Eigen::MatrixBase<ConfigVectorType> & q,
                                    const Eigen::MatrixBase<TangentVectorType> & v,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com)
  {
    assert((model.njoints<floating_flag) && "floating body is not exist.");
    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force nlf;
    computeCenterofInertiaDynamics_nl<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType,gra_flag>(model,data,q,v,nlf);
    computeCenterofInertiaDynamics_forward<Scalar,Options,JointCollectionTpl,floating_flag>(model,data,nlf,F_com,A_com);
  }


  ///
  /// \brief Computes the inverse CoI dynamics wrt. to the world system.
  ///
  /// \tparam Scalar: The scalar type.
  /// \tparam Options: Eigen Alignment options.
  /// \tparam JointCollection: Collection of Joint types.
  /// \tparam ConfigVectorType: Config Vector of Joint  types.
  /// \tparam TangentVectorType: Tangent Vector of Joint Velocity types.
  /// \tparam gra_flag: 1 is considered gravity
  /// \tparam floating: 1 is wrt. to the floating body(the first rigid), and the first joint is floating!.
  ///
  /// \param[in] model:The model structure of the rigid body system.
  /// \param[in] data:The data structure of the rigid body system.
  /// \param[in] nlf_com:The nonlinear force of the CoI dynamics wrt. to the world system.
  /// \param[in] A_com:The whole body accelation of the CoI dynamics wrt. to the world system.  
  /// \param[in] q:The data of joint.
  /// \param[in] v:The data of joint velocity.
  /// \param[out] F_com:The whole body force of the CoI dynamics wrt. to the world system.  
  ///
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl,
           int floating_flag = 1>
  inline void computeCenterofInertiaDynamics_inverse(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_com);

  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl,
           typename ConfigVectorType, typename TangentVectorType,
           int gra_flag = 1, int floating_flag = 1>
  inline void computeCenterofInertiaDynamics_inverse(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com,
                                    const Eigen::MatrixBase<ConfigVectorType> & q,
                                    const Eigen::MatrixBase<TangentVectorType> & v,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_com)
  {
    assert((model.njoints<floating_flag) && "floating body is not exist.");
    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force nlf;
    computeCenterofInertiaDynamics_nl<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType,gra_flag>(model,data,q,v,nlf);
    computeCenterofInertiaDynamics_inverse<Scalar,Options,JointCollectionTpl,floating_flag>(model,data,nlf,A_com,F_com);
  }
  
} // namespace pinocchio

/* --- Details -------------------------------------------------------------------- */
#include "CoMI/CoI/centerofinertia.hxx"

#endif // ifndef __CoMI_CoI_centerofinertia_hpp__
