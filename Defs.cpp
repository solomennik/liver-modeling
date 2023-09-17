#include "Defs.h"

int nDoses = 1;

FP addDoseMmh(const FP& time)	
{
	FP addDoseMmh = M_mh, timeNew = time / (4.48 * pow(10., -5.)) / 60.;
	//for (int i = 0; i < nDoses; ++i)
	//	if (timeNew >= (585. + i * 885.) && timeNew <= (1380. + i * 885.))	//time in minutes
	//		addDoseMmh = M_mh;

	return addDoseMmh;
}

FP addDoseMmr(const FP& time)	
{
	FP addDoseMmr = M_mr, timeNew = time / (4.48 * pow(10., -5.)) / 60.;
	//for (int i = 0; i < nDoses; ++i)
	//	if (timeNew >= (585. + i * 885.) && timeNew <= (1380. + i * 885.))	//time in minutes
	//		addDoseMmr = M_mr;

	return addDoseMmr;
}

FP addDoseF(const FP& time)	
{
	FP addDoseF = f_base, timeNew = time / (4.48 * pow(10., -5.)) / 60.;
	//for (int i = 0; i < nDoses; ++i)
	//	if (timeNew >= (585. + i * 885.) && timeNew <= (1380. + i * 885.))	//time in minutes
	//		addDoseF = f_max;

	return addDoseF;
}

FP addDoseWe(const FP& time)	
{
	FP addDoseWe = We_base, timeNew = time / (4.48 * pow(10., -5.)) / 60.;
	//for (int i = 0; i < nDoses; ++i)
	//	if (timeNew >= (6000. + i * 885.) && timeNew <= (6010. + i * 885.))	//time in minutes
	//		addDoseWe = We_max;

	return addDoseWe;
}

FP addDoseWa(const FP& time)	
{
	FP addDoseWa = Wa_base, timeNew = time / (4.48 * pow(10., -5.)) / 60.;
	//for (int i = 0; i < nDoses; ++i)
	//	if (timeNew >= (12000. + i * 885.) && timeNew <= (12010. + i * 885.))	//time in minutes
	//		addDoseWa = Wa_max;

	return addDoseWa;
}

FP addDoseWs(const FP& time)	
{
	FP addDoseWs = Ws_base, timeNew = time / (4.48 * pow(10., -5.)) / 60.;
	//for (int i = 0; i < nDoses; ++i)
	//	if (timeNew >= (12000. + i * 885.) && timeNew <= (12010. + i * 885.))	//time in minutes
	//		addDoseWs = Ws_max;

	return addDoseWs;
}

FP addDoseWn(const FP& time)	
{
	FP addDoseWn = Wn_base, timeNew = time / (4.48 * pow(10., -5.)) / 60.;
	//for (int i = 0; i < nDoses; ++i)
	//	if (timeNew >= (0. + i * 885.) && timeNew <= (10. + i * 885.))	//time in minutes
	//		addDoseWn = Wn_max;

	return addDoseWn;
}

void stepRungeKutt(FP& m_h_0, FP& m_h, FP& m_r_0, FP& m_r, FP& h_0, FP& h, FP& l_E_0, FP& l_E, FP& l_RB_0, FP& l_RB, FP& l_I_0, FP& l_I, FP& v_E_0, FP& v_E, FP& v_RB_0, FP& v_RB, FP& v_I_0, FP& v_I, FP& r_f_0, FP& r_f, FP& r_I_0, FP& r_I, FP& c_0, FP& c, FP& m_p_0, FP& m_p, FP& p_I_0, FP& p_I, FP& p_E_0, FP& p_E, FP& p_RB_0, FP& p_RB, FP& A_E_0, FP& A_E, FP& p_AB_0, FP& p_AB, FP& S_E_0, FP& S_E, FP& S_i_0, FP& S_i, FP& S_ih_0, FP& S_ih, FP& N_E_0, FP& N_E, FP& N_I_0, FP& N_I, FP& m_N_0, FP& m_N, FP& m_Np_0, FP& m_Np, fVector& k, const FP& timeStep, FP& time)
{
	k[0] = ((1. / J) * (addDoseMmh(time) / (1. + pow(k_mh * (1. + pow(c_0 / k_c, x_c)), x_h)) - del_mh * m_h_0)) * timeStep;
	k[4] = ((1. / J) * (addDoseMmr(time) / (1. + pow(k_mr * (1. + pow(c_0 / k_c, x_c)), x_r)) - del_mr * m_r_0)) * timeStep;
	k[8] = (M_h * m_h_0 - del_h * h_0 - e_S * S_i_0 * SE0_ * h_0 + e__S * S_ih_0 * SE0_ * (1. / S0_)) * timeStep;
	k[12] = ((1. / W) * (-alf_L * r_f_0 * l_E_0 + alf__L * l_RB_0 + W * h_v * ro_v * v_E_0)) * timeStep;
	k[16] = (alf_L * r_f_0 * l_E_0 - alf__L * l_RB_0 - b_L * l_RB_0) * timeStep;
	k[20] = (b_L * l_RB_0 - gam_L * l_I_0) * timeStep;
	k[24] = ((1. / W) * (-alf_v * r_f_0 * v_E_0 + alf__v * v_RB_0 - W * h_v * v_E_0 + W * w)) * timeStep;
	k[28] = (alf_v * r_f_0 * v_E_0 - alf__v * v_RB_0 - b_v * v_RB_0) * timeStep;
	k[32] = (b_v * v_RB_0 - gam_v * v_I_0) * timeStep;
	k[36] = (gam_r * r_I_0 + (m_l / v_l) * (-b_L * (r_f_0 * l_RB_0 / (1. - r_f_0)) - alf_L * l_E_0 * r_f_0 + alf__L * l_RB_0) - b_0 * pow(r_f_0, P) + (m_v / v_v) * (-b_v * (r_f_0 * v_RB_0 / (1. - r_f_0)) - alf_v * v_E_0 * r_f_0 + alf__v * v_RB_0) + (m_p_const / v_p) * (-b_p * (r_f_0 * p_RB_0 / (1. - r_f_0)) - alf_p * p_E_0 * r_f_0 + alf__p * p_RB_0)) * timeStep;
	k[40] = (M_r * m_r_0 - gam_r * r_I_0 + addDoseF(time) * b_0 * pow(r_f_0, P) + addDoseF(time) * (m_l / v_l) * b_L * l_RB_0 * (1. / (1. - r_f_0)) + addDoseF(time) * (m_v / v_v) * b_v * v_RB_0 * (1. / (1. - r_f_0))) * timeStep;
	k[44] = (R_chol_L * sig_l * gam_L * l_I_0 + R_chol_v * sig_v * gam_v * v_I_0 + M_c * h_0 - del_c * c_0) * timeStep;
	k[48] = ((1. / J) * (M_mp / (1. + pow(k_mp * (1. + pow(c_0 / k_c, x_c)), x_p)) - del_mp * m_p_0 - alf_N * m_N_0 * m_p_0 * S0_ + alf__N * m_Np_0)) * timeStep;
	k[52] = (M_p * m_p_0 - gam_p * p_I_0) * timeStep;
	k[56] = ((1. / W) * (-alf_p * r_f_0 * p_E_0 + alf__p * p_RB_0 + gam_p * p_I_0 + W * addDoseWe(time) - e_p * A_E_0 * p_E_0 + e__p * p_AB_0)) * timeStep;
	k[60] = (alf_p * r_f_0 * p_E_0 - alf__p * p_RB_0 - b_p * p_RB_0) * timeStep;
	k[64] = ((1. / W) * (-e_p * A_E_0 * p_E_0 + e__p * p_AB_0 + W * addDoseWa(time))) * timeStep;
	k[68] = (e_p * A_E_0 * p_E_0 - e__p * p_AB_0) * timeStep;
	k[72] = ((1. / W) * (-CLs * S_E_0 + addDoseWs(time))) * timeStep;
	k[76] = (CLs * S_E_0 - e_S * S_i_0 * S0_ * h_0 + e__S * S_ih_0) * timeStep;
	k[80] = (e_S * S_i_0 * S0_ * h_0 - e__S * S_ih_0) * timeStep;
	k[84] = ((1. / W) * (-CLn * N_E_0 + W * addDoseWn(time))) * timeStep;
	k[88] = (CLn * N_E_0 - (gam_N + del_N) * N_I_0) * timeStep;
	k[92] = (gam_N * N_I_0 * (N_E0_ / S0_) - alf_N * m_N_0 * m_p_0 * S0_ + alf__N * m_Np_0 - del_mN * m_N_0) * timeStep;
	k[96] = (alf_N * m_N_0 * m_p_0 * S0_ - (alf__N + del_Np) * m_Np_0) * timeStep;

	k[1] = ((1. / J) * (addDoseMmh(time + (timeStep / 2.)) / (1. + pow(k_mh * (1. + pow((c_0 + 0.5 * k[44]) / k_c, x_c)), x_h)) - del_mh * (m_h_0 + 0.5 * k[0]))) * timeStep;
	k[5] = ((1. / J) * (addDoseMmr(time + (timeStep / 2.)) / (1. + pow(k_mr * (1. + pow((c_0 + 0.5 * k[44]) / k_c, x_c)), x_r)) - del_mr * (m_r_0 + 0.5 * k[4]))) * timeStep;
	k[9] = (M_h * (m_h_0 + 0.5 * k[0]) - del_h * (h_0 + 0.5 * k[8]) - e_S * (S_i_0 + 0.5 * k[76]) * SE0_ * (h_0 + 0.5 * k[8]) + e__S * (S_ih_0 + 0.5 * k[80]) * SE0_ * (1. / S0_)) * timeStep;
	k[13] = ((1. / W) * (-alf_L * (r_f_0 + 0.5 * k[36]) * (l_E_0 + 0.5 * k[12]) + alf__L * (l_RB_0 + 0.5 * k[16]) + W * h_v * ro_v * (v_E_0 + 0.5 * k[24]))) * timeStep;
	k[17] = (alf_L * (r_f_0 + 0.5 * k[36]) * (l_E_0 + 0.5 * k[12]) - alf__L * (l_RB_0 + 0.5 * k[16]) - b_L * (l_RB_0 + 0.5 * k[16])) * timeStep;
	k[21] = (b_L * (l_RB_0 + 0.5 * k[16]) - gam_L * (l_I_0 + 0.5 * k[20])) * timeStep;
	k[25] = ((1. / W) * (-alf_v * (r_f_0 + 0.5 * k[36]) * (v_E_0 + 0.5 * k[24]) + alf__v * (v_RB_0 + 0.5 * k[28]) - W * h_v * (v_E_0 + 0.5 * k[24]) + W * w)) * timeStep;
	k[29] = (alf_v * (r_f_0 + 0.5 * k[36]) * (v_E_0 + 0.5 * k[24]) - alf__v * (v_RB_0 + 0.5 * k[28]) - b_v * (v_RB_0 + 0.5 * k[28])) * timeStep;
	k[33] = (b_v * (v_RB_0 + 0.5 * k[28]) - gam_v * (v_I_0 + 0.5 * k[32])) * timeStep;
	k[37] = (gam_r * (r_I_0 + 0.5 * k[40]) + (m_l / v_l) * (-b_L * ((r_f_0 + 0.5 * k[36]) * (l_RB_0 + 0.5 * k[16]) / (1. - (r_f_0 + 0.5 * k[36]))) - alf_L * (l_E_0 + 0.5 * k[12]) * (r_f_0 + 0.5 * k[36]) + alf__L * (l_RB_0 + 0.5 * k[16])) - b_0 * pow((r_f_0 + 0.5 * k[36]), P) + (m_v / v_v) * (-b_v * ((r_f_0 + 0.5 * k[36]) * (v_RB_0 + 0.5 * k[28]) / (1. - (r_f_0 + 0.5 * k[36]))) - alf_v * (v_E_0 + 0.5 * k[24]) * (r_f_0 + 0.5 * k[36]) + alf__v * (v_RB_0 + 0.5 * k[28])) + (m_p_const / v_p) * (-b_p * ((r_f_0 + 0.5 * k[36]) * (p_RB_0 + 0.5 * k[60]) / (1. - (r_f_0 + 0.5 * k[36]))) - alf_p * (p_E_0 + 0.5 * k[56]) * (r_f_0 + 0.5 * k[36]) + alf__p * (p_RB_0 + 0.5 * k[60]))) * timeStep;
	k[41] = (M_r * (m_r_0 + 0.5 * k[4]) - gam_r * (r_I_0 + 0.5 * k[40]) + addDoseF(time + (timeStep / 2.)) * b_0 * pow((r_f_0 + 0.5 * k[36]), P) + addDoseF(time + (timeStep / 2.)) * (m_l / v_l) * b_L * (l_RB_0 + 0.5 * k[16]) * (1. / (1. - (r_f_0 + 0.5 * k[36]))) + addDoseF(time + (timeStep / 2.)) * (m_v / v_v) * b_v * (v_RB_0 + 0.5 * k[28]) * (1. / (1. - (r_f_0 + 0.5 * k[36])))) * timeStep;
	k[45] = (R_chol_L * sig_l * gam_L * (l_I_0 + 0.5 * k[20]) + R_chol_v * sig_v * gam_v * (v_I_0 + 0.5 * k[32]) + M_c * (h_0 + 0.5 * k[8]) - del_c * (c_0 + 0.5 * k[44])) * timeStep;
	k[49] = ((1. / J) * (M_mp / (1. + pow(k_mp * (1. + pow((c_0 + 0.5 * k[44]) / k_c, x_c)), x_p)) - del_mp * (m_p_0 + 0.5 * k[48]) - alf_N * (m_N_0 + 0.5 * k[92]) * (m_p_0 + 0.5 * k[48]) * S0_ + alf__N * (m_Np_0 + 0.5 * k[96]))) * timeStep;
	k[53] = (M_p * (m_p_0 + 0.5 * k[48]) - gam_p * (p_I_0 + 0.5 * k[52])) * timeStep;
	k[57] = ((1. / W) * (-alf_p * (r_f_0 + 0.5 * k[36]) * (p_E_0 + 0.5 * k[56]) + alf__p * (p_RB_0 + 0.5 * k[60]) + gam_p * (p_I_0 + 0.5 * k[52]) + W * addDoseWe(time + (timeStep / 2.)) - e_p * (A_E_0 + 0.5 * k[64]) * (p_E_0 + 0.5 * k[56]) + e__p * (p_AB_0 + 0.5 * k[68]))) * timeStep;
	k[61] = (alf_p * (r_f_0 + 0.5 * k[36]) * (p_E_0 + 0.5 * k[56]) - alf__p * (p_RB_0 + 0.5 * k[60]) - b_p * (p_RB_0 + 0.5 * k[60])) * timeStep;
	k[65] = ((1. / W) * (-e_p * (A_E_0 + 0.5 * k[64]) * (p_E_0 + 0.5 * k[56]) + e__p * (p_AB_0 + 0.5 * k[68]) + W * addDoseWa(time + (timeStep / 2.)))) * timeStep;
	k[69] = (e_p * (A_E_0 + 0.5 * k[64]) * (p_E_0 + 0.5 * k[56]) - e__p * (p_AB_0 + 0.5 * k[68])) * timeStep;
	k[73] = ((1. / W) * (-CLs * (S_E_0 + 0.5 * k[72]) + addDoseWs(time + (timeStep / 2.)))) * timeStep;
	k[77] = (CLs * (S_E_0 + 0.5 * k[72]) - e_S * (S_i_0 + 0.5 * k[76]) * S0_ * (h_0 + 0.5 * k[8]) + e__S * (S_ih_0 + 0.5 * k[80])) * timeStep;
	k[81] = (e_S * (S_i_0 + 0.5 * k[76]) * S0_ * (h_0 + 0.5 * k[8]) - e__S * (S_ih_0 + 0.5 * k[80])) * timeStep;
	k[85] = ((1. / W) * (-CLn * (N_E_0 + 0.5 * k[84]) + W * addDoseWn(time + (timeStep / 2.)))) * timeStep;
	k[89] = (CLn * (N_E_0 + 0.5 * k[84]) - (gam_N + del_N) * (N_I_0 + 0.5 * k[88])) * timeStep;
	k[93] = (gam_N * (N_I_0 + 0.5 * k[88]) * (N_E0_ / S0_) - alf_N * (m_N_0 + 0.5 * k[92]) * (m_p_0 + 0.5 * k[48]) * S0_ + alf__N * (m_Np_0 + 0.5 * k[96]) - del_mN * (m_N_0 + 0.5 * k[92])) * timeStep;
	k[97] = (alf_N * (m_N_0 + 0.5 * k[92]) * (m_p_0 + 0.5 * k[48]) * S0_ - (alf__N + del_Np) * (m_Np_0 + 0.5 * k[96])) * timeStep;

	k[2] = ((1. / J) * (addDoseMmh(time + (timeStep / 2.)) / (1. + pow(k_mh * (1. + pow((c_0 + 0.5 * k[45]) / k_c, x_c)), x_h)) - del_mh * (m_h_0 + 0.5 * k[1]))) * timeStep;
	k[6] = ((1. / J) * (addDoseMmr(time + (timeStep / 2.)) / (1. + pow(k_mr * (1. + pow((c_0 + 0.5 * k[45]) / k_c, x_c)), x_r)) - del_mr * (m_r_0 + 0.5 * k[5]))) * timeStep;
	k[10] = (M_h * (m_h_0 + 0.5 * k[1]) - del_h * (h_0 + 0.5 * k[9]) - e_S * (S_i_0 + 0.5 * k[77]) * SE0_ * (h_0 + 0.5 * k[9]) + e__S * (S_ih_0 + 0.5 * k[81]) * SE0_ * (1. / S0_)) * timeStep;
	k[14] = ((1. / W) * (-alf_L * (r_f_0 + 0.5 * k[37]) * (l_E_0 + 0.5 * k[13]) + alf__L * (l_RB_0 + 0.5 * k[17]) + W * h_v * ro_v * (v_E_0 + 0.5 * k[25]))) * timeStep;
	k[18] = (alf_L * (r_f_0 + 0.5 * k[37]) * (l_E_0 + 0.5 * k[13]) - alf__L * (l_RB_0 + 0.5 * k[17]) - b_L * (l_RB_0 + 0.5 * k[17])) * timeStep;
	k[22] = (b_L * (l_RB_0 + 0.5 * k[17]) - gam_L * (l_I_0 + 0.5 * k[21])) * timeStep;
	k[26] = ((1. / W) * (-alf_v * (r_f_0 + 0.5 * k[37]) * (v_E_0 + 0.5 * k[25]) + alf__v * (v_RB_0 + 0.5 * k[29]) - W * h_v * (v_E_0 + 0.5 * k[25]) + W * w)) * timeStep;
	k[30] = (alf_v * (r_f_0 + 0.5 * k[37]) * (v_E_0 + 0.5 * k[25]) - alf__v * (v_RB_0 + 0.5 * k[29]) - b_v * (v_RB_0 + 0.5 * k[29])) * timeStep;
	k[34] = (b_v * (v_RB_0 + 0.5 * k[29]) - gam_v * (v_I_0 + 0.5 * k[33])) * timeStep;
	k[38] = (gam_r * (r_I_0 + 0.5 * k[41]) + (m_l / v_l) * (-b_L * ((r_f_0 + 0.5 * k[37]) * (l_RB_0 + 0.5 * k[17]) / (1. - (r_f_0 + 0.5 * k[37]))) - alf_L * (l_E_0 + 0.5 * k[13]) * (r_f_0 + 0.5 * k[37]) + alf__L * (l_RB_0 + 0.5 * k[17])) - b_0 * pow((r_f_0 + 0.5 * k[37]), P) + (m_v / v_v) * (-b_v * ((r_f_0 + 0.5 * k[37]) * (v_RB_0 + 0.5 * k[29]) / (1. - (r_f_0 + 0.5 * k[37]))) - alf_v * (v_E_0 + 0.5 * k[25]) * (r_f_0 + 0.5 * k[37]) + alf__v * (v_RB_0 + 0.5 * k[29])) + (m_p_const / v_p) * (-b_p * ((r_f_0 + 0.5 * k[37]) * (p_RB_0 + 0.5 * k[61]) / (1. - (r_f_0 + 0.5 * k[37]))) - alf_p * (p_E_0 + 0.5 * k[57]) * (r_f_0 + 0.5 * k[37]) + alf__p * (p_RB_0 + 0.5 * k[61]))) * timeStep;
	k[42] = (M_r * (m_r_0 + 0.5 * k[5]) - gam_r * (r_I_0 + 0.5 * k[41]) + addDoseF(time + (timeStep / 2.)) * b_0 * pow((r_f_0 + 0.5 * k[37]), P) + addDoseF(time + (timeStep / 2.)) * (m_l / v_l) * b_L * (l_RB_0 + 0.5 * k[17]) * (1. / (1. - (r_f_0 + 0.5 * k[37]))) + addDoseF(time + (timeStep / 2.)) * (m_v / v_v) * b_v * (v_RB_0 + 0.5 * k[29]) * (1. / (1. - (r_f_0 + 0.5 * k[37])))) * timeStep;
	k[46] = (R_chol_L * sig_l * gam_L * (l_I_0 + 0.5 * k[21]) + R_chol_v * sig_v * gam_v * (v_I_0 + 0.5 * k[33]) + M_c * (h_0 + 0.5 * k[9]) - del_c * (c_0 + 0.5 * k[45])) * timeStep;
	k[50] = ((1. / J) * (M_mp / (1. + pow(k_mp * (1. + pow((c_0 + 0.5 * k[45]) / k_c, x_c)), x_p)) - del_mp * (m_p_0 + 0.5 * k[49]) - alf_N * (m_N_0 + 0.5 * k[93]) * (m_p_0 + 0.5 * k[49]) * S0_ + alf__N * (m_Np_0 + 0.5 * k[97]))) * timeStep;
	k[54] = (M_p * (m_p_0 + 0.5 * k[49]) - gam_p * (p_I_0 + 0.5 * k[53])) * timeStep;
	k[58] = ((1. / W) * (-alf_p * (r_f_0 + 0.5 * k[37]) * (p_E_0 + 0.5 * k[57]) + alf__p * (p_RB_0 + 0.5 * k[61]) + gam_p * (p_I_0 + 0.5 * k[53]) + W * addDoseWe(time + (timeStep / 2.)) - e_p * (A_E_0 + 0.5 * k[65]) * (p_E_0 + 0.5 * k[57]) + e__p * (p_AB_0 + 0.5 * k[69]))) * timeStep;
	k[62] = (alf_p * (r_f_0 + 0.5 * k[37]) * (p_E_0 + 0.5 * k[57]) - alf__p * (p_RB_0 + 0.5 * k[61]) - b_p * (p_RB_0 + 0.5 * k[61])) * timeStep;
	k[66] = ((1. / W) * (-e_p * (A_E_0 + 0.5 * k[65]) * (p_E_0 + 0.5 * k[57]) + e__p * (p_AB_0 + 0.5 * k[69]) + W * addDoseWa(time + (timeStep / 2.)))) * timeStep;
	k[70] = (e_p * (A_E_0 + 0.5 * k[65]) * (p_E_0 + 0.5 * k[57]) - e__p * (p_AB_0 + 0.5 * k[69])) * timeStep;
	k[74] = ((1. / W) * (-CLs * (S_E_0 + 0.5 * k[73]) + addDoseWs(time + (timeStep / 2.)))) * timeStep;
	k[78] = (CLs * (S_E_0 + 0.5 * k[73]) - e_S * (S_i_0 + 0.5 * k[77]) * S0_ * (h_0 + 0.5 * k[9]) + e__S * (S_ih_0 + 0.5 * k[81])) * timeStep;
	k[82] = (e_S * (S_i_0 + 0.5 * k[77]) * S0_ * (h_0 + 0.5 * k[9]) - e__S * (S_ih_0 + 0.5 * k[81])) * timeStep;
	k[86] = ((1. / W) * (-CLn * (N_E_0 + 0.5 * k[85]) + W * addDoseWn(time + (timeStep / 2.)))) * timeStep;
	k[90] = (CLn * (N_E_0 + 0.5 * k[85]) - (gam_N + del_N) * (N_I_0 + 0.5 * k[89])) * timeStep;
	k[94] = (gam_N * (N_I_0 + 0.5 * k[89]) * (N_E0_ / S0_) - alf_N * (m_N_0 + 0.5 * k[93]) * (m_p_0 + 0.5 * k[49]) * S0_ + alf__N * (m_Np_0 + 0.5 * k[97]) - del_mN * (m_N_0 + 0.5 * k[93])) * timeStep;
	k[98] = (alf_N * (m_N_0 + 0.5 * k[93]) * (m_p_0 + 0.5 * k[49]) * S0_ - (alf__N + del_Np) * (m_Np_0 + 0.5 * k[97])) * timeStep;

	k[3] = ((1. / J) * (addDoseMmh(time + timeStep) / (1. + pow(k_mh * (1. + pow((c_0 + k[46]) / k_c, x_c)), x_h)) - del_mh * (m_h_0 + k[2]))) * timeStep;
	k[7] = ((1. / J) * (addDoseMmr(time + timeStep) / (1. + pow(k_mr * (1. + pow((c_0 + k[46]) / k_c, x_c)), x_r)) - del_mr * (m_r_0 + k[6]))) * timeStep;
	k[11] = (M_h * (m_h_0 + k[2]) - del_h * (h_0 + k[10]) - e_S * (S_i_0 + k[78]) * SE0_ * (h_0 + k[10]) + e__S * (S_ih_0 + k[82]) * SE0_ * (1. / S0_)) * timeStep;
	k[15] = ((1. / W) * (-alf_L * (r_f_0 + k[38]) * (l_E_0 + k[14]) + alf__L * (l_RB_0 + k[18]) + W * h_v * ro_v * (v_E_0 + k[26]))) * timeStep;
	k[19] = (alf_L * (r_f_0 + k[38]) * (l_E_0 + k[14]) - alf__L * (l_RB_0 + k[18]) - b_L * (l_RB_0 + k[18])) * timeStep;
	k[23] = (b_L * (l_RB_0 + k[18]) - gam_L * (l_I_0 + k[22])) * timeStep;
	k[27] = ((1. / W) * (-alf_v * (r_f_0 + k[38]) * (v_E_0 + k[26]) + alf__v * (v_RB_0 + k[30]) - W * h_v * (v_E_0 + k[26]) + W * w)) * timeStep;
	k[31] = (alf_v * (r_f_0 + k[38]) * (v_E_0 + k[26]) - alf__v * (v_RB_0 + k[30]) - b_v * (v_RB_0 + k[30])) * timeStep;
	k[35] = (b_v * (v_RB_0 + k[30]) - gam_v * (v_I_0 + k[34])) * timeStep;
	k[39] = (gam_r * (r_I_0 + k[42]) + (m_l / v_l) * (-b_L * ((r_f_0 + k[38]) * (l_RB_0 + k[18]) / (1. - (r_f_0 + k[38]))) - alf_L * (l_E_0 + k[14]) * (r_f_0 + k[38]) + alf__L * (l_RB_0 + k[18])) - b_0 * pow((r_f_0 + k[38]), P) + (m_v / v_v) * (-b_v * ((r_f_0 + k[38]) * (v_RB_0 + k[30]) / (1. - (r_f_0 + k[38]))) - alf_v * (v_E_0 + k[26]) * (r_f_0 + k[38]) + alf__v * (v_RB_0 + k[30])) + (m_p_const / v_p) * (-b_p * ((r_f_0 + k[38]) * (p_RB_0 + k[62]) / (1. - (r_f_0 + k[38]))) - alf_p * (p_E_0 + k[58]) * (r_f_0 + k[38]) + alf__p * (p_RB_0 + k[62]))) * timeStep;
	k[43] = (M_r * (m_r_0 + k[6]) - gam_r * (r_I_0 + k[42]) + addDoseF(time + timeStep) * b_0 * pow((r_f_0 + k[38]), P) + addDoseF(time + timeStep) * (m_l / v_l) * b_L * (l_RB_0 + k[18]) * (1. / (1. - (r_f_0 + k[38]))) + addDoseF(time + timeStep) * (m_v / v_v) * b_v * (v_RB_0 + k[30]) * (1. / (1. - (r_f_0 + k[38])))) * timeStep;
	k[47] = (R_chol_L * sig_l * gam_L * (l_I_0 + k[22]) + R_chol_v * sig_v * gam_v * (v_I_0 + k[34]) + M_c * (h_0 + k[10]) - del_c * (c_0 + k[46])) * timeStep;
	k[51] = ((1. / J) * (M_mp / (1. + pow(k_mp * (1. + pow((c_0 + k[46]) / k_c, x_c)), x_p)) - del_mp * (m_p_0 + k[50]) - alf_N * (m_N_0 + k[94]) * (m_p_0 + k[50]) * S0_ + alf__N * (m_Np_0 + k[98]))) * timeStep;
	k[55] = (M_p * (m_p_0 + k[50]) - gam_p * (p_I_0 + k[54])) * timeStep;
	k[59] = ((1. / W) * (-alf_p * (r_f_0 + k[38]) * (p_E_0 + k[58]) + alf__p * (p_RB_0 + k[62]) + gam_p * (p_I_0 + k[54]) + W * addDoseWe(time + timeStep) - e_p * (A_E_0 + k[66]) * (p_E_0 + k[58]) + e__p * (p_AB_0 + k[70]))) * timeStep;
	k[63] = (alf_p * (r_f_0 + k[38]) * (p_E_0 + k[58]) - alf__p * (p_RB_0 + k[62]) - b_p * (p_RB_0 + k[62])) * timeStep;
	k[67] = ((1. / W) * (-e_p * (A_E_0 + k[66]) * (p_E_0 + k[58]) + e__p * (p_AB_0 + k[70]) + W * addDoseWa(time + timeStep))) * timeStep;
	k[71] = (e_p * (A_E_0 + k[66]) * (p_E_0 + k[58]) - e__p * (p_AB_0 + k[70])) * timeStep;
	k[75] = ((1. / W) * (-CLs * (S_E_0 + k[74]) + addDoseWs(time + timeStep))) * timeStep;
	k[79] = (CLs * (S_E_0 + k[74]) - e_S * (S_i_0 + k[78]) * S0_ * (h_0 + k[10]) + e__S * (S_ih_0 + k[82])) * timeStep;
	k[83] = (e_S * (S_i_0 + k[78]) * S0_ * (h_0 + k[10]) - e__S * (S_ih_0 + k[82])) * timeStep;
	k[87] = ((1. / W) * (-CLn * (N_E_0 + k[86]) + W * addDoseWn(time + timeStep))) * timeStep;
	k[91] = (CLn * (N_E_0 + k[86]) - (gam_N + del_N) * (N_I_0 + k[90])) * timeStep;
	k[95] = (gam_N * (N_I_0 + k[90]) * (N_E0_ / S0_) - alf_N * (m_N_0 + k[94]) * (m_p_0 + k[50]) * S0_ + alf__N * (m_Np_0 + k[98]) - del_mN * (m_N_0 + k[94])) * timeStep;
	k[99] = (alf_N * (m_N_0 + k[94]) * (m_p_0 + k[50]) * S0_ - (alf__N + del_Np) * (m_Np_0 + k[98])) * timeStep;

	m_h = m_h_0 + (1. / 6.) * (k[0] + 2. * k[1] + 2. * k[2] + k[3]);
	m_r = m_r_0 + (1. / 6.) * (k[4] + 2. * k[5] + 2. * k[6] + k[7]);
	h = h_0 + (1. / 6.) * (k[8] + 2. * k[9] + 2. * k[10] + k[11]);
	l_E = l_E_0 + (1. / 6.) * (k[12] + 2. * k[13] + 2. * k[14] + k[15]);
	l_RB = l_RB_0 + (1. / 6.) * (k[16] + 2. * k[17] + 2. * k[18] + k[19]);
	l_I = l_I_0 + (1. / 6.) * (k[20] + 2. * k[21] + 2. * k[22] + k[23]);
	v_E = v_E_0 + (1. / 6.) * (k[24] + 2. * k[25] + 2. * k[26] + k[27]);
	v_RB = v_RB_0 + (1. / 6.) * (k[28] + 2. * k[29] + 2. * k[30] + k[31]);
	v_I = v_I_0 + (1. / 6.) * (k[32] + 2. * k[33] + 2. * k[34] + k[35]);
	r_f = r_f_0 + (1. / 6.) * (k[36] + 2. * k[37] + 2. * k[38] + k[39]);
	r_I = r_I_0 + (1. / 6.) * (k[40] + 2. * k[41] + 2. * k[42] + k[43]);
	c = c_0 + (1. / 6.) * (k[44] + 2. * k[45] + 2. * k[46] + k[47]);
	m_p = m_p_0 + (1. / 6.) * (k[48] + 2. * k[49] + 2. * k[50] + k[51]);
	p_I = p_I_0 + (1. / 6.) * (k[52] + 2. * k[53] + 2. * k[54] + k[55]);
	p_E = p_E_0 + (1. / 6.) * (k[56] + 2. * k[57] + 2. * k[58] + k[59]);
	p_RB = p_RB_0 + (1. / 6.) * (k[60] + 2. * k[61] + 2. * k[62] + k[63]);
	A_E = A_E_0 + (1. / 6.) * (k[64] + 2. * k[65] + 2. * k[66] + k[67]);
	p_AB = p_AB_0 + (1. / 6.) * (k[68] + 2. * k[69] + 2. * k[70] + k[71]);
	S_E = S_E_0 + (1. / 6.) * (k[72] + 2. * k[73] + 2. * k[74] + k[75]);
	S_i = S_i_0 + (1. / 6.) * (k[76] + 2. * k[77] + 2. * k[78] + k[79]);
	S_ih = S_ih_0 + (1. / 6.) * (k[80] + 2. * k[81] + 2. * k[82] + k[83]);
	N_E = N_E_0 + (1. / 6.) * (k[84] + 2. * k[85] + 2. * k[86] + k[87]);
	N_I = N_I_0 + (1. / 6.) * (k[88] + 2. * k[89] + 2. * k[90] + k[91]);
	m_N = m_N_0 + (1. / 6.) * (k[92] + 2. * k[93] + 2. * k[94] + k[95]);
	m_Np = m_Np_0 + (1. / 6.) * (k[96] + 2. * k[97] + 2. * k[98] + k[99]);

	m_h_0 = m_h;
	m_r_0 = m_r;
	h_0 = h;
	l_E_0 = l_E;
	l_RB_0 = l_RB;
	l_I_0 = l_I;
	v_E_0 = v_E;
	v_RB_0 = v_RB;
	v_I_0 = v_I;
	r_f_0 = r_f;
	r_I_0 = r_I;
	c_0 = c;
	m_p_0 = m_p;
	p_I_0 = p_I;
	p_E_0 = p_E;
	p_RB_0 = p_RB;
	A_E_0 = A_E;
	p_AB_0 = p_AB;
	S_E_0 = S_E;
	S_i_0 = S_i;
	S_ih_0 = S_ih;
	N_E_0 = N_E;
	N_I_0 = N_I;
	m_N_0 = m_N;
	m_Np_0 = m_Np;
}