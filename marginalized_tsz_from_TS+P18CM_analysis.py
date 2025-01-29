# This is the script for final figure with best-fit
# run with:
# $ python marginalized_tsz_from_TS+P18CM_analysis.py
import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pathlib
from scipy.interpolate import interp1d


TCMB = 2.726 #Kelvin
TCMB_uK = 2.726e6 #micro-Kelvin


hplanck=6.626068e-34 #MKS
kboltz=1.3806503e-23 #MKS
clight=299792458.0 #MKS
m_elec = 510.999 #keV

l_t21 = np.asarray([
68.5,
89.5,
117.0,
152.5,
198.0,
257.5,
335.5,
436.5,
567.5,
738.0,
959.5,
1247.5
])
dl_t21 = np.asarray([
0.048,
0.068, 0.101, 0.127, 0.161, 0.208, 0.261, 0.318, 0.380, 0.472, 0.550, 0.686
])
sigma_dl_t21 = np.asarray([
0.103,
0.055, 0.046, 0.038, 0.033, 0.037, 0.029, 0.039, 0.060, 0.097, 0.162, 0.285
])
# tsz function
def fsz(nu_ghz):
    nu = 1.e9*np.asarray(nu_ghz).astype(float)
    X = hplanck*nu/(kboltz*TCMB)
    resp = (X / np.tanh(X/2.0) - 4.0) * TCMB_uK #put explicitly into uK_CMB units, so that output ILC map is in Compton-y
#     resp[np.where(nu_ghz == None)] = 0. #this case is appropriate for HI or other maps that contain no CMB-relevant signals (and also no CIB); they're assumed to be denoted by None in nu_ghz
    return resp



path_to_files = str(pathlib.Path(__file__).parent.resolve())+'/'


def run(args):
    #



    if args.fg_from_P18CMB is not None:
        D = np.loadtxt(path_to_files+'szpowespectrum_measurement_urc_snr6_p18cmb_bf_fg_from_TSZ+P18_l_clyy_sigclyy_cib_ir_rs_cn.txt')
    else:
        D = np.loadtxt(path_to_files+'szpowespectrum_measurement_urc_snr6_p18cmb_bf_fg_from_TSZonly_l_clyy_sigclyy_cib_ir_rs_cn.txt')

    multipoles = D[:,0]
    marg_sz = D[:,1]
    marg_sz_yerr = D[:,2]
    Cl_cib =  D[:,3]
    Cl_ir = D[:,4]
    Cl_rs = D[:,5]
    Cl_cn =  D[:,6]



    R = np.loadtxt(path_to_files + 'szpowespectrum_measurement_urc_snr6_p18cmb_best_fit_curve_l_cl1h_cl2h_cl1h+2h.txt')
    multipoles_bf = R[:,0]
    cl_1h = R[:,1]
    cl_2h = R[:,2]

    label_size = 12
    title_size = 15
    legend_size = 13
    handle_length = 2
    if args.R21_full is not None:
        y_min = 5.e-4
        y_max = 4
    else:
        y_min = 1.e-2
        y_max = 2
    y_max = 10

    fig, ax1 = plt.subplots(1,1,figsize=(7,5))

    ax = ax1
    ax.tick_params(axis = 'x',which='both',length=5,direction='in', pad=10)
    ax.tick_params(axis = 'y',which='both',length=5,direction='in', pad=5)
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=label_size)
    plt.setp(ax.get_xticklabels(), fontsize=label_size)
    ax.grid( b=True, which="both", alpha=0.3, linestyle='--')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\ell$',size=title_size)

    ax.set_ylabel(r'$10^{12}\ell(\ell+1)\mathrm{C^{yy}_\ell/2\pi}$',size=title_size)


    if args.SZ_from_templates is None:
        PSZ = np.loadtxt(path_to_files+'tsz_143_eps0.50.dat')
        B12 = np.loadtxt(path_to_files+'cl_tsz_150_bat.dat')

        pow_planck = interp1d(PSZ[:,0],PSZ[:,1])
        tsz_planck = PSZ[:,1]
        tsz_planck = 1e12*tsz_planck/fsz(143)/fsz(143)/pow_planck(3000)
        A_sz_planck_mean = 5.
        A_sz_planck_sigma = 2.
        A_sz_planck_best_fit = 7.

        plt.plot(PSZ[:,0],A_sz_planck_best_fit*tsz_planck,label='P18 MFLIKE best-fit (EM12 template)',c='k',ls=':')
        plt.fill_between(PSZ[:,0],(A_sz_planck_mean-A_sz_planck_sigma)*tsz_planck,(A_sz_planck_mean+A_sz_planck_sigma)*tsz_planck,alpha=0.3,label = r'1-$\sigma$ P18')


        tsz_B12 = B12[:,1]
        pow_b12 = interp1d(B12[:,0],tsz_B12)
        tsz_B12 = 1e12*tsz_B12/pow_b12(3000)/fsz(148)/fsz(148)
        # tsz_B12 = 1e12*tsz_B12/pow_b12(3000)/(2.67e6)**2
        pow_tsz_b12 = interp1d(B12[:,0],tsz_B12)
        A_sz_act_mean = 5.29
        A_sz_act_sigma = 0.66
        A_sz_act_best_fit = 5.29
        # A_sz_act_best_fit = 1.

        plt.plot(B12[:,0][600:],A_sz_act_best_fit*tsz_B12[600:],label='ACT MFLIKE best-fit (B12 template)',c='r',ls='--')
        plt.fill_between(B12[:,0][600:],(A_sz_act_mean-A_sz_act_sigma)*tsz_B12[600:],(A_sz_act_mean+A_sz_act_sigma)*tsz_B12[600:],alpha=0.3,label = r'1-$\sigma$ C20')



        # path_to_chains = "/Users/boris/Work/CLASS-SZ/ATSZ_ACMB/analysis/atsz_analysis/"
        path_to_chains = "/Users/boris/Work/CLASS-SZ/ATSZ_ACMB/margestats_files/"
        # chain_root = 'act_extended_act_plus_planck_wbp_less353545_regcov1p13_mtsz_logspaced'
        # chain_root = "act_plus_planck_mcmc_d1p2_w1p4"
        # chain_root = "act_plus_planck_mcmc_d1p2_w1p4_18sept23-more_szbins_onlyACTcovfac"
        # chain_root = "act_plus_planck_mcmc_d1p2_w1p4_20sept23-not-more_szbins_APcovfac"
        # chain_root = "act_plus_planck_mcmc_d1p2_w1p4_20sept23-not-more_szbins_APcovfac.mcmc_2023-11-13_11.42.03"
        chain_root = "act_plus_planck_mcmc_d1p2_w1p4_Aonlycovfac_NOACT_cutlleq1500.mcmc_2023-11-15_09.25.00"
        # chain_root = "act_plus_planck_mcmc_d1p2_w1p4_Aonlycovfac.mcmc_2023-11-14_08.51.00"
        df_mtsz = open(path_to_chains + chain_root + '.margestats','r')
        lines = df_mtsz.readlines()
        for index, line in enumerate(lines[:10]):
              if line[0:9] == 'parameter':
                    idl = index
                    break
        print(idl)
        for index, line in enumerate(lines):
        #       print(lines)
        #       print(line)
        #       print(line.strip())
              lines[index] = line.strip()
              ls = lines[index].split(' ')
        #       print(ls)
              ls = [l for l in ls if l != '']

              lines[index] = ls[:3]
        lines = lines[idl:]


        all_rows = lines[1:]
        # print(lines)
        df_mtsz = pd.DataFrame(all_rows,columns=lines[0])

        id3000_ppp = np.where(df_mtsz['parameter']=='a_tsz_3000_ppp')[0]
        id3000_pp = np.where(df_mtsz['parameter']=='a_tsz_3000_pp')[0]
        id3000_p = np.where(df_mtsz['parameter']=='a_tsz_3000_p')[0]

        id3000 = np.where(df_mtsz['parameter']=='a_tsz_3000')[0]
        id2500 = np.where(df_mtsz['parameter']=='a_tsz_2500')[0]
        id2000 = np.where(df_mtsz['parameter']=='a_tsz_2000')[0]
        id1500 = np.where(df_mtsz['parameter']=='a_tsz_1500')[0]
        id1000 = np.where(df_mtsz['parameter']=='a_tsz_1000')[0]
        # df_mtsz['mean'][id3000]
        # df_mtsz['sddev'][id3000]

        atsz_mf = []
        atsz_mfe = []
        l_atsz_mf = []

        # atsz_mf.append(df_mtsz['mean'][id3000_ppp])
        # atsz_mfe.append(df_mtsz['sddev'][id3000_ppp])

        # atsz_mf.append(df_mtsz['mean'][id3000_pp])
        # atsz_mfe.append(df_mtsz['sddev'][id3000_pp])
        #
        # atsz_mf.append(df_mtsz['mean'][id3000_p])
        # atsz_mfe.append(df_mtsz['sddev'][id3000_p])

        atsz_mf.append(df_mtsz['mean'][id3000])
        atsz_mfe.append(df_mtsz['sddev'][id3000])
        l_atsz_mf.append(3050)
        atsz_mf.append(df_mtsz['mean'][id2500])
        atsz_mfe.append(df_mtsz['sddev'][id2500])
        l_atsz_mf.append(2500)
        atsz_mf.append(df_mtsz['mean'][id2000])
        atsz_mfe.append(df_mtsz['sddev'][id2000])
        l_atsz_mf.append(2000)
        atsz_mf.append(df_mtsz['mean'][id1500])
        atsz_mfe.append(df_mtsz['sddev'][id1500])
        l_atsz_mf.append(1500)
        atsz_mf.append(df_mtsz['mean'][id1000])
        atsz_mfe.append(df_mtsz['sddev'][id1000])
        l_atsz_mf.append(1000)

        atsz_mf = np.asarray(atsz_mf,dtype=np.float64)
        print("atsz_mf",atsz_mf)
        atsz_mfe = np.asarray(atsz_mfe,dtype=np.float64)
        l_atsz_mf = np.asarray(l_atsz_mf,dtype=np.float64)
        l_atsz_mf = [float(a) for a in l_atsz_mf]

        bin_mid = np.linspace(np.log(3e3),np.log(1e3),5)
        # bin_mid = np.linspace(np.log(1e3),np.log(1e4),8)[:-1]
        # bin_mid = bin_mid[::-1] ## more tsz bins ## 13 nov 23: shouldnt be inverted !!!,we fill the atsz starting from high-ell's !

        l_atsz_mf = np.exp(bin_mid)
        # l_atsz_mf = [3000,2500,2000,1500,1000]
        for a,l in zip(atsz_mf,l_atsz_mf):
            print(l,a)
        # atsz_mf = [5.29,5.29,5.29,5.29,5.29]
        atsz_mf = np.asarray([float(a)*pow_tsz_b12(l) for (a,l) in zip(atsz_mf,l_atsz_mf)])
        atsz_mfe = np.asarray([float(a)*pow_tsz_b12(l)for (a,l) in zip(atsz_mfe,l_atsz_mf)])

        print('tot SN:',np.sqrt(np.sum((atsz_mf/atsz_mfe)**2)))


        ax.errorbar(l_atsz_mf[3:],atsz_mf[3:],atsz_mfe[3:],color='k',markeredgecolor='r',
        alpha=1.,ls='None',marker='o',markersize=6,capsize=5,elinewidth=2,markeredgewidth=2,
        label = r'tSZ from multi-freq lkl')


    if args.R21_full is not None:

        ax.plot(multipoles_bf,cl_1h,color='k',ls='--',lw=2,alpha = 1.)
        ax.plot(multipoles_bf,cl_2h,color='k',ls='-.',lw=2,alpha = 1.)
        ax.plot(multipoles_bf,cl_1h+cl_2h,color='r',ls='-',lw=2,alpha = 1.)
        ax.plot(multipoles,Cl_cib,color='b',ls='--',alpha = .5)
        ir, = ax.plot(multipoles,Cl_ir,color='orange',ls='-.',alpha = .8)
        ir.set_dashes([4, 2, 1, 2, 1, 2, 1, 2])
        ax.plot(multipoles,Cl_rs,color='b',ls=':',alpha = .5)
        cn, = ax.plot(multipoles,Cl_cn,color='b',ls=':',alpha = .7)
        cn.set_dashes([4, 2, 1, 2, 1, 2])
        ax.errorbar(multipoles,marg_sz,yerr=[marg_sz_yerr,marg_sz_yerr],color='k',ls='None',alpha = 1.,label = r'R21 [uRC, $q_\mathrm{cut}=6$]',marker='o',markersize=3,capsize=5,elinewidth=2,markeredgewidth=2)

    else:
        # ax.plot(multipoles,marg_sz,color='k',ls='None',alpha = 1.,label = 'R21',marker='o',markersize=3,markeredgewidth=2)
        plt.fill_between(multipoles,marg_sz+marg_sz_yerr,marg_sz-marg_sz_yerr,color='k',alpha=0.2,label = r'1-$\sigma$ R21 [uRC, $q_\mathrm{cut}=6$]')

    ax.fill_between(l_t21,dl_t21-sigma_dl_t21,dl_t21+sigma_dl_t21,color='orange',alpha=0.2,label='Tani 21')
    ax.plot(l_t21,dl_t21,marker='o',c='k',markersize=2.,label='T21')

    ACTCellnew = 5.29/(2.67)**2
    ACTCellnewerr =	0.66/(2.67)**2

    # actpol + mbac
    ACTCellnew = 5.0/(2.67)**2
    ACTCellnewerr =	0.7/(2.67)**2

    SPTCellnew = 3.42/(2.84)**2
    SPTCellnewerr = 0.54/(2.84)**2

    ax.errorbar([3000],[SPTCellnew],
                yerr=[SPTCellnewerr],capthick=1,capsize=3,elinewidth=1,\
                 fmt='>',mfc='k',markersize=8,color = 'r')

    ax.errorbar([3000],[ACTCellnew],
                yerr=[ACTCellnewerr],capthick=1,capsize=3,elinewidth=2,\
                 fmt='<',mfc='k',markersize=8,color = 'r',zorder=100)
    ax.text(3.1e3, 3.e-1, 'SPT (2020)', fontsize=9,rotation = 0,
    verticalalignment='center',horizontalalignment='left')
    ax.text(3.1e3, 1.01e0, 'ACT (2020)', fontsize=9,rotation = 0,
    verticalalignment='center',horizontalalignment='left')



    P = np.loadtxt(path_to_files+'planck_sz_1712_00788v1.txt')
    # ax.errorbar(P[:,0][0:18],P[:,1][0:18],yerr = P[:,2][0:18],color='green',linestyle="None",alpha = .5,label = 'B18',marker='o',markersize=3,capsize=5,elinewidth=2,markeredgewidth=2)
    ax.plot(P[:,0][0:18],P[:,1][0:18],color='green',linestyle="None",alpha = .5,label = 'B18',marker='o',markersize=3,markeredgewidth=2)
    ax.fill_between(P[:,0][0:18],P[:,1][0:18]-P[:,2][0:18],P[:,1][0:18]+P[:,2][0:18],
    color='green',alpha=0.1,label=r'1-$\sigma$ B18')

    P = np.loadtxt(path_to_files+'Planck2015.txt')
    # ax.errorbar(P[:,0][0:18],P[:,1][0:18],yerr = P[:,4][0:18],color='purple',linestyle="None",alpha = .2,label = 'P15',marker='o',markersize=3,capsize=5,elinewidth=2,markeredgewidth=2)
    ax.plot(P[:,0][0:18],P[:,1][0:18],color='purple',linestyle="None",alpha = .2,label = 'P15',marker='o',markersize=3,markeredgewidth=2)
    ax.fill_between(P[:,0][0:18],P[:,1][0:18]-P[:,4][0:18],P[:,1][0:18]+P[:,4][0:18],
    color='purple',alpha=0.1,label=r'1-$\sigma$ P15')


    ax.set_ylim(y_min,y_max)
    ax.legend(loc=2,ncol=1,fontsize=10,frameon=False)

    plt.xlim(8,1.1e4)


    if args.R21_full is not None:
        ax.text(1.1e2, 2.4e-2, '1-halo', fontsize=11,rotation = 45,
        verticalalignment='bottom',horizontalalignment='left')
        ax.text(2e3, 6.5e-3, '2-halo', fontsize=11,rotation = -30,
        verticalalignment='bottom',horizontalalignment='left')
        ax.text(1.e3, 1.2, 'BEST-FIT TSZ', fontsize=9,rotation = 18,
        verticalalignment='bottom',horizontalalignment='left')

        ax.text(3.e2, 7e-2, 'IR Res.', fontsize=9,rotation = 58,
        verticalalignment='center',horizontalalignment='left')

        ax.text(4.7e2, 5e-3, 'Noise Res.', fontsize=9,rotation = 58,
        verticalalignment='center',horizontalalignment='left')

        ax.text(1e3, 4e-2, 'Radio Res.', fontsize=9,rotation = 0,
        verticalalignment='center',horizontalalignment='left')
        ax.text(1e3, 3.5e-3, 'CIB Res.', fontsize=9,rotation = 0,
        verticalalignment='center',horizontalalignment='left')
    fig.tight_layout()

    plt.savefig(path_to_files+'my_marginalized_tsz_fig_mftsz.pdf')
    plt.show(block=True)


def main():
    parser=argparse.ArgumentParser(description="plotting results")
    parser.add_argument("-R21_full",help="R21_full" ,dest="R21_full", type=str, required=False)
    parser.add_argument("-fg_from_P18CMB",help="fg_from_P18CMB" ,dest="fg_from_P18CMB", type=str, required=False)
    parser.add_argument("-SZ_from_templates",help="SZ_from_templates" ,dest="SZ_from_templates", type=str, required=False)
    parser.set_defaults(func=run)
    args=parser.parse_args()
    args.func(args)

if __name__=="__main__":
	main()
