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
  /// \brief Computes the nonlinear force of the CoI dynamics wrt. the world system.
  ///
  /// \tparam Scalar: The scalar type.
  /// \tparam Options: Eigen Alignment options.
  /// \tparam JointCollection: Collection of Joint types.
  /// \tparam ConfigVectorType: Config Vector of Joint  types.
  /// \tparam TangentVectorType: Tangent Vector of Joint Velocity types.
  /// \tparam gra_flag: 1 is considered gravity
  ///
  /// \param[in] model:The model structure of the rigid body system.
  /// \param[in] data:The data structure of the rigid body system.
  /// \param[in] q:The data of joint.
  /// \param[in] v:The data of joint velocity.
  /// \param[in] V_s:The whole velocity of CoI dynamic wrt. the world system.
  /// \param[out] nlf_com:The nonlinear force of the CoI dynamics wrt. the world system.
  ///
  /*template<int gra_flag = 1,
           typename Scalar, int Options,template<typename,int> class JointCollectionTpl>
  inline void computeCenterofInertiaDynamics_nl(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com);*/  
  template<typename Scalar, int Options,template<typename,int> class JointCollectionTpl>
  inline void computeCenterofInertiaDynamics_nl(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion &V_s,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com);                                                             
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, 
           typename ConfigVectorType, typename TangentVectorType>
  inline void computeCenterofInertiaDynamics_nl(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const Eigen::MatrixBase<ConfigVectorType> & q,
                                    const Eigen::MatrixBase<TangentVectorType> & v,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion &V_s,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com)
  {
    forwardKinematics(model,data,q,v);
    //data.oMi is the SE3 of all joint in the world;
    //data.liMi is the SE3 of all joint in self coordinate;
    //data.v is the spatial of all joint wrt the world in self coordinate;
    computeCenterofInertiaDynamics_nl<Scalar,Options,JointCollectionTpl>(model,data,V_s,nlf_com);
  }


  ///
  /// \brief Computes the forward CoI dynamics wrt. the floating system.
  ///
  /// \tparam Scalar: The scalar type.
  /// \tparam Options: Eigen Alignment options.
  /// \tparam JointCollection: Collection of Joint types.
  /// \tparam ConfigVectorType: Config Vector of Joint  types.
  /// \tparam TangentVectorType: Tangent Vector of Joint Velocity types.
  /// \tparam gra_flag: 1 is considered gravity 
  /// \tparam floating_i: 1 is wrt. the floating system(the first rigid), and the first joint is floating!.
  ///
  /// \param[in] model:The model structure of the rigid body system.
  /// \param[in] data:The data structure of the rigid body system.
  /// \param[in] V_s:The whole velocity of CoI dynamic wrt. the world system.
  /// \param[in] nlf_com:The nonlinear force of the CoI dynamics wrt. the body system.
  /// \param[in] F_s:The whole body force of the CoI dynamics wrt. the world system.
  /// \param[in] q:The data of joint.
  /// \param[in] v:The data of joint velocity.
  /// \param[out] A_com:The whole body accelation of the CoI dynamics wrt. the body system.  
  ///
  template<int gravity_flag = 1,int floating_i = 1,
           typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  inline void computeCenterofInertiaDynamics_forward(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & V_s,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_s,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com);

  template<int gravity_flag = 1, int floating_i = 1,
           typename Scalar, int Options, template<typename,int> class JointCollectionTpl, 
           typename ConfigVectorType, typename TangentVectorType>
  inline void computeCenterofInertiaDynamics_forward(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_s,
                                    const Eigen::MatrixBase<ConfigVectorType> & q,
                                    const Eigen::MatrixBase<TangentVectorType> & v,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com)
  {
    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force nlf;
    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion V_s;
    computeCenterofInertiaDynamics_nl<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType>(model,data,q,v,V_s,nlf);
    computeCenterofInertiaDynamics_forward<gravity_flag,floating_i,Scalar,Options,JointCollectionTpl>(model,data,V_s,nlf,F_s,A_com);
  }


  ///
  /// \brief Computes the inverse CoI dynamics wrt. the floating body.
  ///
  /// \tparam Scalar: The scalar type.
  /// \tparam Options: Eigen Alignment options.
  /// \tparam JointCollection: Collection of Joint types.
  /// \tparam ConfigVectorType: Config Vector of Joint  types.
  /// \tparam TangentVectorType: Tangent Vector of Joint Velocity types.
  /// \tparam gra_flag: 1 is considered gravity
  /// \tparam floating_i: 1 is wrt. the floating system(the first rigid), and the first joint is floating!.
  ///
  /// \param[in] model:The model structure of the rigid body system.
  /// \param[in] data:The data structure of the rigid body system.
  /// \param[in] V_s:The whole velocity of CoI dynamic wrt. the world system.
  /// \param[in] nlf_com:The nonlinear force of the CoI dynamics wrt. the body system.
  /// \param[in] A_com:The whole body accelation of the CoI dynamics wrt. the body system.
  /// \param[in] q:The data of joint.
  /// \param[in] v:The data of joint velocity.
  /// \param[out] F_s:The whole body force of the CoI dynamics wrt. the world system.  
  ///
  template<int gravity_flag = 1,int floating_i = 0,
           typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  inline void computeCenterofInertiaDynamics_inverse(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & V_s,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_s);

  template<int gravity_flag = 1,int floating_i = 0,
           typename Scalar, int Options, template<typename,int> class JointCollectionTpl,
           typename ConfigVectorType, typename TangentVectorType>
  inline void computeCenterofInertiaDynamics_inverse(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com,
                                    const Eigen::MatrixBase<ConfigVectorType> & q,
                                    const Eigen::MatrixBase<TangentVectorType> & v,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_s)
  {
    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force nlf;
    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion V_s;
    computeCenterofInertiaDynamics_nl<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType>(model,data,q,v,V_s,nlf);
    computeCenterofInertiaDynamics_inverse<gravity_flag,floating_i,Scalar,Options,JointCollectionTpl>(model,data,V_s,nlf,A_com,F_s);
  }
  
} // namespace pinocchio

/* --- Details -------------------------------------------------------------------- */
#include "CoMI/CoI/centerofinertia.hxx"

#endif // ifndef __CoMI_CoI_centerofinertia_hpp__
