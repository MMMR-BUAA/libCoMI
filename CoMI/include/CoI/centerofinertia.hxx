//
// Copyright (c) 2015-2023 CNRS INRIA
//

#ifndef __CoMI_CoI_centerofinertia_hxx__
#define __CoMI_CoI_centerofinertia_hxx__

#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/spatial/act-on-set.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/check.hpp"
#include "pinocchio/spatial/symmetric3.hpp"
#include "pinocchio/spatial/skew.hpp"

/// @cond DEV

namespace pinocchio
{

  namespace details
  {
    /// \brief get the generated velocity of the momentum V=M-1*H
    ///
    /// \tparam Scalar: The scalar type.
    /// \tparam Options: Eigen Alignment options.
    ///
    /// \param[in] I:The generated inertia.
    /// \param[in] f:The momentum.
    /// \param[out] v:The generated velocity.
    ///
    template <typename Scalar, int Options>
    inline void CoM_Minvh(const InertiaTpl<Scalar, Options> &I, const ForceTpl<Scalar, Options> &f, MotionTpl<Scalar, Options> &v)
    {
      typedef Eigen::Matrix<Scalar,6,1,Options> Vector6;
      typedef Symmetric3Tpl<Scalar,Options> Symmetric3;
      typedef Eigen::Matrix<Scalar,3,3,Options> Matrix3;
      //get the rotate inertia
      Vector6 rI=I.inertia().data();
      //get the determinant of rotate inertia
      Scalar detrI=rI(0)*(rI(2)*rI(5)-rI(4)*rI(4))-rI(1)*rI(1)*rI(5)+rI(3)*(rI(1)*rI(4)*2-rI(2)*rI(3));
      assert(!Eigen::internal::isApprox(detrI,0.0) && "the determinant of rotate inertia is zero!");
      //get the inverse matrix of the generated inertia using the element.
      Matrix3 rIInv;
      rIInv(0,0)=(rI(2)*rI(5)-rI(4)*rI(4))/detrI;
      rIInv(1,0)=(rI(3)*rI(4)-rI(1)*rI(5))/detrI;rIInv(1,1)=(rI(0)*rI(5)-rI(3)*rI(3))/detrI;
      rIInv(2,0)=(rI(1)*rI(4)-rI(2)*rI(3))/detrI;rIInv(2,1)=(rI(1)*rI(3)-rI(0)*rI(4))/detrI;rIInv(2,2)=(rI(0)*rI(2)-rI(1)*rI(1))/detrI;
      rIInv(0,1)=rIInv(1,0);rIInv(0,2)=rIInv(2,0);rIInv(1,2)=rIInv(2,1);
      //get the generated velocity
      Matrix3 sp=skew(I.lever());
      v.linear()=(Matrix3::Identity()/I.mass()-sp*rIInv*sp)*f.linear()+sp*rIInv*f.angular();
      v.angular()=-sp*rIInv*f.linear()+rIInv*f.angular();
    }

    template <typename Scalar, int Options>
    inline MotionTpl<Scalar,Options> CoM_Minvh(const InertiaTpl<Scalar,Options> &I, const ForceTpl<Scalar,Options> &f)
    {
      MotionTpl<Scalar,Options> v;
      CoM_Minvh<Scalar,Options>(I,f,v);
      return v;
    }

  }
  /// Computes the Iterative algorithm of CoI dynamics.
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  struct DCcrbaCoIStep
  : public fusion::JointUnaryVisitorBase< DCcrbaCoIStep<Scalar,Options,JointCollectionTpl> >
  {
    typedef ModelTpl<Scalar,Options,JointCollectionTpl> Model;
    typedef DataTpl<Scalar,Options,JointCollectionTpl> Data;
    
    typedef boost::fusion::vector<const Model &,
                                  Data &
                                  > ArgsType;
    
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data)
    {
      typedef typename Data::Inertia Inertia;
      typedef typename Data::Force   Force;
      typedef typename SizeDepType<JointModel::NV>::template ColsReturn<typename Data::Matrix6x>::Type ColsBlock;
      
      const JointIndex & i = jmodel.id();
      const JointIndex & parent = model.parents[i];
      const Inertia & Y = data.oYcrb[i];
      const Force & oh = data.oh[i];
      const typename Inertia::Matrix6 & doYcrb = data.doYcrb[i];
      
      data.oYcrb[parent] += Y;
      if(parent > 0)
        data.oh[parent]+=oh;
        data.doYcrb[parent] += doYcrb;
    }
    
  }; // struct DCcrbaCoIStep

  //Computes the nonlinear force of the CoI dynamics wrt. to the world system.
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl, int gra_flag>
  inline void computeCenterofInertiaDynamics_nl(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com)
  {
    assert(model.check(data) && "data is not consistent with model.");
    
    typedef ModelTpl<Scalar,Options,JointCollectionTpl> Model;
    typedef typename Model::JointIndex JointIndex;
    typedef DataTpl<Scalar,Options,JointCollectionTpl> Data;
    typedef typename Data::Motion Motion;
    typedef typename Data::Inertia Inertia;
    typedef typename Data::Force   Force;

    //compute every body and its time derivate on the world.
    data.oYcrb[0].setZero();
    data.doYcrb[0].setZero();
    for(JointIndex i=1; i<(JointIndex)(model.njoints); ++i)
    {
      data.oYcrb[i] = data.oinertias[i] = data.oMi[i].act(model.inertias[i]);
      data.ov[i] = data.oMi[i].act(data.v[i]); // v_i expressed in the world frame
      data.oh[i]=data.oYcrb[i]*data.ov[i];
      data.doYcrb[i] = data.oYcrb[i].variation(data.ov[i]);
    }
    //compute whole body and its time derivate on the world.
    typedef DCcrbaCoIStep<Scalar,Options,JointCollectionTpl> PassCoI;
    for(JointIndex i=(JointIndex)(model.njoints-1); i>0; --i)
    {
      PassCoI::run(model.joints[i],data.joints[i],
                 typename PassCoI::ArgsType(model,data));
    }
    //return to the joint h
    for(JointIndex i=1; i<(JointIndex)(model.njoints); ++i)
    {
      data.h[i]=data.oMi[i].actInv(data.oh[i]);
      data.com[i] = data.oYcrb[i].lever();
      data.mass[i] = data.oYcrb[i].mass();
    }
    //obtain the CoI parameters.
    data.mass[0] = data.oYcrb[0].mass();
    data.com[0] = data.oYcrb[0].lever();
    data.hg=data.oh[0];
    data.Ig=data.oYcrb[0];
    data.Itmp=data.doYcrb[0];
    Motion V_com=details::CoM_Minvh<Scalar,Options>(data.Ig,data.hg);
    //compute forward centerofinertia M*dot_V+dot_M*V+ad(V)*M*V=Fout on the whole.
    if(gra_flag)
    {
      nlf_com.linear()=model.gravity981*data.mass[0];
      nlf_com.angular()=data.com[0].cross(nlf_com.linear());
    }
    nlf_com+=data.hg.motionAction(V_com);
    nlf_com.linear()-=(data.Itmp.template block<3,3>(Inertia::LINEAR,Inertia::LINEAR) * V_com.linear());
    nlf_com.linear()-=(data.Itmp.template block<3,3>(Inertia::LINEAR,Inertia::ANGULAR) * V_com.angular());
    nlf_com.angular()-=(data.Itmp.template block<3,3>(Inertia::ANGULAR,Inertia::LINEAR) * V_com.linear());
    nlf_com.angular()-=(data.Itmp.template block<3,3>(Inertia::ANGULAR,Inertia::ANGULAR) * V_com.angular());
  }
  //Computes the forward CoI dynamics wrt. to the world system.
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  inline void computeCenterofInertiaDynamics_forward(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_com,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com)
  {
    A_com=details::CoM_Minvh<Scalar,Options>(data.Ig,F_com+nlf_com);
  }
  //Computes the inverse CoI dynamics wrt. to the world system.
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  inline void computeCenterofInertiaDynamics_inverse(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                    DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & nlf_com,
                                    const typename DataTpl<Scalar,Options,JointCollectionTpl>::Motion & A_com,
                                    typename DataTpl<Scalar,Options,JointCollectionTpl>::Force & F_com)
  {
    F_com=data.Ig*A_com-nlf_com;
  }
} // namespace pinocchio

/// @endcond

#endif // ifndef __CoMI_CoI_centerofinertia_hxx__
