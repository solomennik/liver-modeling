#include "Defs.h"


int main()
{
	fVector	k(100);

	FP m_h, m_r, h, l_E, l_RB, l_I, v_E, v_RB, v_I, r_f, r_I, c, m_p, p_I, p_E, p_RB, A_E, p_AB, S_E, S_i, S_ih, N_E, N_I, m_N, m_Np;

	// Initial values
	FP m_h_0 = 3.65 * pow(10., -8.), 
		m_r_0 = 3.862 * pow(10., -7.),
		h_0 = 0.615 * pow(10., -5.),
		l_E_0 = 0.05,
		l_RB_0 = 0.021, 
		l_I_0 = 0.017,
		v_E_0 = 0.1535, 
		v_RB_0 = 1.056, 
		v_I_0 = 0.855,
		r_f_0 = 0.678, 
		r_I_0 = 0.1659, 
		c_0 = 0.4945,
		m_p_0 = 3.65 * pow(10., -8.), 
		p_I_0 = 0., 
		p_E_0 = 1., 
		p_RB_0 = 0.,
		A_E_0 = 0.,
		p_AB_0 = 0.,
		S_E_0 = 0.,
		S_i_0 = 0.,
		S_ih_0 = 0.,
		N_E_0 = 0.,		
		N_I_0 = 0.,
		m_N_0 = 0.,
		m_Np_0 = 0.;

	FP simDuration = 61., timeStep = 0.00005, time;

	int N = int(ceil(simDuration / timeStep)), iStep;

	ofstream of_m_h("output/m_h.txt");
	ofstream of_m_r("output/m_r.txt");
	ofstream of_h("output/h.txt");
	ofstream of_l_E("output/l_E.txt");
	ofstream of_l_RB("output/l_RB.txt");
	ofstream of_l_I("output/l_I.txt");
	ofstream of_v_E("output/v_E.txt");
	ofstream of_v_RB("output/v_RB.txt");
	ofstream of_v_I("output/v_I.txt");
	ofstream of_r_f("output/r_f.txt");
	ofstream of_r_I("output/r_I.txt");
	ofstream of_c("output/c.txt");
	ofstream of_m_p("output/m_p.txt");
	ofstream of_p_I("output/p_I.txt");
	ofstream of_p_E("output/p_E.txt");
	ofstream of_p_RB("output/p_RB.txt");
	ofstream of_A_E("output/A_E.txt");
	ofstream of_p_AB("output/p_AB.txt");
	ofstream of_S_E("output/S_E.txt");
	ofstream of_S_i("output/S_i.txt");
	ofstream of_S_ih("output/S_ih.txt");
	ofstream of_N_E("output/N_E.txt");
	ofstream of_N_I("output/N_I.txt");
	ofstream of_m_N("output/m_N.txt");
	ofstream of_m_Np("output/m_Np.txt");

	for (iStep = 0; iStep < N; iStep += 1)
	{
		time = FP(iStep * timeStep);

		stepRungeKutt(m_h_0, m_h, m_r_0, m_r, h_0, h, l_E_0, l_E, l_RB_0, l_RB, l_I_0, l_I, v_E_0, v_E, v_RB_0, v_RB, v_I_0, v_I, r_f_0, r_f, r_I_0, r_I, c_0, c, m_p_0, m_p, p_I_0, p_I, p_E_0, p_E, p_RB_0, p_RB, A_E_0, A_E, p_AB_0, p_AB, S_E_0, S_E, S_i_0, S_i, S_ih_0, S_ih, N_E_0, N_E, N_I_0, N_I, m_N_0, m_N, m_Np_0, m_Np, k, timeStep, time);

		of_m_h << setprecision(8) << time << "\t" << m_h << endl;
		of_m_r << setprecision(8) << time << "\t" << m_r << endl;
		of_h << setprecision(8) << time << "\t" << h << endl;
		of_l_E << setprecision(8) << time << "\t" << l_E << endl;
		of_l_RB << setprecision(8) << time << "\t" << l_RB << endl;
		of_l_I << setprecision(8) << time << "\t" << l_I << endl;
		of_v_E << setprecision(8) << time << "\t" << v_E << endl;
		of_v_RB << setprecision(8) << time << "\t" << v_RB << endl;
		of_v_I << setprecision(8) << time << "\t" << v_I << endl;
		of_r_f << setprecision(8) << time << "\t" << r_f << endl;
		of_r_I << setprecision(8) << time << "\t" << r_I << endl;
		of_c << setprecision(8) << time << "\t" << c << endl;
		of_m_p << setprecision(8) << time << "\t" << m_p << endl;
		of_p_I << setprecision(8) << time << "\t" << p_I << endl;
		of_p_E << setprecision(8) << time << "\t" << p_E << endl;
		of_p_RB << setprecision(8) << time << "\t" << p_RB << endl;
		of_A_E << setprecision(8) << time << "\t" << A_E << endl;
		of_p_AB << setprecision(8) << time << "\t" << p_AB << endl;
		of_S_E << setprecision(8) << time << "\t" << S_E << endl;
		of_S_i << setprecision(8) << time << "\t" << S_i << endl;
		of_S_ih << setprecision(8) << time << "\t" << S_ih << endl;
		of_N_E << setprecision(8) << time << "\t" << N_E << endl;
		of_N_I << setprecision(8) << time << "\t" << N_I << endl;
		of_m_N << setprecision(8) << time << "\t" << m_N << endl;
		of_m_Np << setprecision(8) << time << "\t" << m_Np << endl;
	}

	of_m_h.close();
	of_m_r.close();
	of_h.close();
	of_l_E.close();
	of_l_RB.close();
	of_l_I.close();
	of_v_E.close();
	of_v_RB.close();
	of_v_I.close();
	of_r_f.close();
	of_r_I.close();
	of_c.close();
	of_m_p.close();
	of_p_I.close();
	of_p_E.close();
	of_p_RB.close();
	of_A_E.close();
	of_p_AB.close();
	of_S_E.close();
	of_S_i.close();
	of_S_ih.close();
	of_N_E.close();
	of_N_I.close();
	of_m_N.close();
	of_m_Np.close();

	system("pause");
	return 0;
}