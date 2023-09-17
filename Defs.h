#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

typedef double FP;
typedef vector<FP> fVector;

const FP r_f0 = 0.99999,

f_base = 0.7,
f_max = 0.94,

We_base = 0.,
We_max = 2.5 * pow(10., 15.),

Wa_base = 0.,
Wa_max = 1. * pow(10., 16.),

Ws_base = 0.,
Ws_max = 1. * pow(10., 16.),

Wn_base = 1. * pow(10., 11.),
Wn_max = 1. * pow(10., 12.),

J = 0.1,
W = 1.5 * pow(10., 2.),
M_mh = 1.406 * pow(10., -7.),
k_mh = 1.,
del_mh = 1.,
M_mr = 1.24 * pow(10., -6.),
k_mr = 1.,
del_mr = 1.,
k_c = 4.714 * pow(10., -1.),
x_c = 4.,
x_h = 3.,
x_r = 1.,
M_h = 7.4011 * pow(10., 2.),
del_h = 1.433,
alf_L = 48.46,
alf__L = 13.17,
h_v = 1.94 * pow(10., -1.),
b_L = 60.27,
gam_L = 74.4,
alf_v = 6.782 * pow(10., 2.),
alf__v = 6.585,
w = 4.54 * pow(10., -1.),
//w = 0., 
b_v = 60.27,
gam_v = 74.4,
gam_r = 2.232 * pow(10., 2.),
P = 180.,
b_0 = 0.,
m_l = 1.,
m_v = 2.,
M_r = 2.876 * pow(10., 7.),
R_chol_L = 3400.,
R_chol_v = 3100.,
M_c = 2.099 * pow(10., 5.),
del_c = 2.679,
ro_v = 2.521 * pow(10., -1.),
sig_l = 6.19 * pow(10., -7.),
sig_v = 1.561 * pow(10., -7.),
v_v = 11.05,
v_l = 2.786,

b_p = 60.27,
alf_p = 7.299,
alf__p = 2.232,
k_mp = 1.,
x_p = 1.,
del_mp = 1.,
m_p_const = 1.,

pE0_ = 1.,		
//pE0_ = 1.5 * pow(10., 14.),		
v_p = 3.27 * pow(10., 13.) / pE0_,

M_mp = 0.272 * pow(10., -6.),
M_p = 183.083 * pow(10., 18.) / pE0_,
gam_p = 2.232 * pow(10., 2.),

//anti-PCSK9 1
e_p = 0.223 * pow(10., -9.) * pE0_,
e__p = 22.3,

//anti-PCSK9 2
//e_p = 0.223 * pow(10., -8.) * pE0_,
//e__p = 11.16,

//anti-PCSK9 3
//e_p = 0.223 * pow(10., -10.) * pE0_,
//e__p = 22.3,

//statins
S0_ = 8.21 * pow(10., 16.),
SE0_ = 1. * pow(10., 12.),
e_S = 0.223 * pow(10., -9.),
e__S = 223.,
CLs = 223. * pow(10., -3.),

//siRNA
N_E0_ = 1.,		
//N_E0_ = pow(10., 11.),	
CLn = 44.64,
gam_N = 22.32,
del_N = 4.464,
alf_N = 0.223 * pow(10., -9.),
alf__N = 22.32,
del_mN = 1.,
del_Np = 100.;

void stepRungeKutt(FP& m_h_0, FP& m_h, FP& m_r_0, FP& m_r, FP& h_0, FP& h, FP& l_E_0, FP& l_E, FP& l_RB_0, FP& l_RB, FP& l_I_0, FP& l_I, FP& v_E_0, FP& v_E, FP& v_RB_0, FP& v_RB, FP& v_I_0, FP& v_I, FP& r_f_0, FP& r_f, FP& r_I_0, FP& r_I, FP& c_0, FP& c, FP& m_p_0, FP& m_p, FP& p_I_0, FP& p_I, FP& p_E_0, FP& p_E, FP& p_RB_0, FP& p_RB, FP& A_E_0, FP& A_E, FP& p_AB_0, FP& p_AB, FP& S_E_0, FP& S_E, FP& S_i_0, FP& S_i, FP& S_ih_0, FP& S_ih, FP& N_E_0, FP& N_E, FP& N_I_0, FP& N_I, FP& m_N_0, FP& m_N, FP& m_Np_0, FP& m_Np, fVector& k, const FP& timeStep, FP& time);
