#include <iostream>
#include <sstream>
#include <fcntl.h>
#include "pinocchio/parsers/urdf.hpp"
#include "pinocchio/parsers/sample-models.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/algorithm/centroidal.hpp"
#include "pinocchio/algorithm/centroidal-derivatives.hpp"
#include "CoMI/CoI/centerofinertia.hpp"

using namespace std;

/******************************************************define*****************************************************/




/********************************************main***********************************************************************/
int main(int argc, char ** argv)
{
    std::cout<<"/**********test start*****************/"<<std::endl;
    using namespace pinocchio;
    
    //urdf机器人模型文件位置
    const std::string urdf_filename=std::string("/home/cjw/MMMR_code/libCoMI/test/model/urdf/robot_series_two.urdf");
    std::cout<<urdf_filename<<std::endl;

    //加载机器人模型
    pinocchio::Model model;
    pinocchio::JointModelFreeFlyer root_joint;
    pinocchio::urdf::buildModel(urdf_filename,root_joint,model);
    //pinocchio::urdf::buildModel(urdf_filename,model);
    std::cout << "model name: " << model.name << std::endl;

    //模型数据
    Data data(model);

    //产生关节数据
    Eigen::VectorXd q = pinocchio::neutral(model);
    std::cout << "q: " << q.transpose() << std::endl;
    Eigen::VectorXd v = Eigen::VectorXd::Zero(model.nv);
    std::cout << "v: " << v.transpose() << std::endl;
    Eigen::VectorXd a = Eigen::VectorXd::Zero(model.nv);
    std::cout << "a: " << a.transpose() << std::endl;
 

    //计算正运动学
    forwardKinematics(model, data, q,v);
    //打印正运动学结果
    for (JointIndex joint_id = 0; joint_id < (JointIndex)model.njoints; ++joint_id)
    {
        std::cout << std::setw(24) << std::left
                  << model.names[joint_id] << ": "
                  << std::fixed << std::setprecision(2)
                  << data.oMi[joint_id].translation().transpose()<< ": "
                  << std::fixed << std::setprecision(2)
                  << data.ov[joint_id].linear().transpose()<< "+"
                  << data.ov[joint_id].angular().transpose()
                  << std::endl;
    }
    //计算动力学
    //const Eigen::VectorXd & tau = pinocchio::rnea(model,data,q,v,a);
    //std::cout << "tau = " << tau.transpose() << std::endl;
    //pinocchio::computeCentroidalMomentum(model,data,q,v);
    pinocchio::Force  Hg=pinocchio::computeCentroidalMomentumTimeVariation(model,data,q,v,a);
    std::cout << "data.dhg = " << data.dhg.linear().transpose()<< data.dhg.angular().transpose() << std::endl;
    std::cout << "data.hg = " << data.hg.linear().transpose()<< data.hg.angular().transpose() << std::endl;
    std::cout << "data.com = " << data.com[0].transpose() << std::endl;
    std::cout << "data.vcom = " << data.vcom[0].transpose() << std::endl;
    Eigen::Matrix<double,6,Eigen::Dynamic> dh_dq=Eigen::MatrixXd::Zero(6,model.nv);
    Eigen::Matrix<double,6,Eigen::Dynamic> dhdot_dq=Eigen::MatrixXd::Zero(6,model.nv);
    Eigen::Matrix<double,6,Eigen::Dynamic> dhdot_dv=Eigen::MatrixXd::Zero(6,model.nv);
    Eigen::Matrix<double,6,Eigen::Dynamic> dhdot_da=Eigen::MatrixXd::Zero(6,model.nv);
    pinocchio::computeCentroidalDynamicsDerivatives(model,data,q,v,a,dh_dq,dhdot_dq,dhdot_dv,dhdot_da);
    std::cout << "data.dh_dq = " << dh_dq<< std::endl;
    std::cout << "data.dhdot_dq = " << dhdot_dq<< std::endl;
    std::cout << "data.dhdot_dv = " << dhdot_dv<< std::endl;
    std::cout << "data.dhdot_da = " << dhdot_da<< std::endl;
    pinocchio::Motion A_com=data.oa[0];
    pinocchio::Force  F_com;
    pinocchio::computeCenterofInertiaDynamics(model,data,A_com,F_com);
    std::cout << "F_com = " << F_com.linear().transpose()<<"+"<<F_com.angular().transpose() << std::endl;



    std::cout<<"/**********test end*****************/"<<std::endl;
    return 0;
}