# -*- coding: utf-8 -*-
"""
pw_eib_ring_beam_parallel.py

Solver for right-hand modes in electron-ion-beam (e-i-b) plasma

Yasuhito Narita and Uwe Motschmann
y.narita@tu-braunschweig.de

May 2026
License: MIT License

Usage:
    Refer to the README for detailed instructions.
    If you use this code, please cite it as follows:
    DOI: TBD.

    Narita, Y., and Motschmann, U.:
    Non-resonant beam instability driven by the ring-beam ion distribution,
    AIP Advances 16, 000000 (2026)
"""

# MIT License (Short form)
# Copyright (c) 2026 Y. Narita and U. Motschmann
#
# This software is released under the MIT License.
# http://opensource.org

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator


def disp_eib():

    # ================ #
    # input parameters #
    # ================ #

    # --------------------------------------------#
    # frequency parameters (in units of omega_ci) #
    # --------------------------------------------#
    omega_ce = -1836
    omega_ci = 1
    omega_cb = 1

    omega_pi = 1000
    omega_pb = 100

    # ---------------------------------------------#
    # omega_pe determined by the charge neutrality #
    # ---------------------------------------------#
    omega_pe = np.sqrt(
        (-omega_pi**2 / omega_ci - omega_pb**2 / omega_cb) * omega_ce
    )

    # --------------------------------------------------#
    # Alfven velocity Va (in units of speed of light c) #
    # --------------------------------------------------#
    v_a = 1 / np.sqrt(
        (omega_pe / omega_ce) ** 2
        + (omega_pi / omega_ci) ** 2
        + (omega_pb / omega_cb) ** 2
    )

    # ---------------------------------------------#
    # beam velocities (in units of speed of light) #
    # ---------------------------------------------#
    v_e = 0
    v_b = 20 * v_a

    #---------------------------------#
    # ring velocity (beam population) #
    #---------------------------------#
    v_perp = 10 * v_a
    vph = (v_perp**2) / 2

    # -----------------------------------------------------------------------#
    # wavenumber array (in units of ion inertial wavenumber as k*c/omega_pi) #
    # -----------------------------------------------------------------------#
    kcw01_min, kcw01_max, kcw01_tics = -0.2, 0.2, 4000
    karr = np.linspace(kcw01_min, kcw01_max, num=kcw01_tics, endpoint=False)


    # ============================= #
    # further variables, file names #
    # ============================= #

    #----------------------------------------------------------------------#
    # ion bulk velocity from current neutrality in electron rest frame v_e #
    #----------------------------------------------------------------------#
    v_i = -(omega_pb / omega_pi) ** 2 * (omega_ci / omega_cb) * v_b

    # ---------------------------#
    # output file specifications #
    # ---------------------------#
#    outfile = "fig_vpara_20_vperp_20"
    v_b_out = v_b/v_a
    v_perp_out = v_perp/v_a
    outfile = f"fig_vpara_{v_b_out:02.0f}_vperp_{v_perp_out:02.0f}"
    pdffile = outfile + ".pdf"

    # ---------------------------#
    # validation parameter codes #
    # ---------------------------#

    # ------------------------#
    # charge neutrality check #
    # ------------------------#
#    laneu = (
#        omega_pe**2 / omega_ce + omega_pi**2 / omega_ci + omega_pb**2 / omega_cb
#    )

    # ----------------------#
    # current density check #
    # ----------------------#
#    strom = (
#        (omega_pe**2 / omega_ce) * v_e
#        + (omega_pi**2 / omega_ci) * v_i
#        + (omega_pb**2 / omega_cb) * v_b
#    )

    # ---------------#
    # total momentum #
    # ---------------#
#    impuls = (
#        (omega_pe**2 / omega_ce**2) * v_e
#        + (omega_pi**2 / omega_ci**2) * v_i
#        + (omega_pb**2 / omega_cb**2) * v_b
#    )

    # -------------------#
    # firehose parameter #
    # -------------------#
#    fhp = (
#        (v_e * omega_pe / omega_ce) ** 2
#        + (v_i * omega_pi / omega_ci) ** 2
#        + (v_b * omega_pb / omega_cb) ** 2
#        - 0.5 * (v_perp * omega_pb / omega_cb) ** 2
#    )

    # -------------------#
    # console messages   #
    # -------------------#
    print(f'\n{" Frequency parameters (in units of omega_ci) ":=^40}')
    print(
        f"omega_pe = {omega_pe:.4f}\nomega_pi = {omega_pi}\nomega_pb = {omega_pb}"
    )
    print(f"omega_ce = {omega_ce}\nomega_ci = {omega_ci}\nomega_cb = {omega_cb}")
    print(f'\n{" Velocities (in units of speed of light c) ":=^40}')
    print(f"v_e/v_a  = {v_e}\nv_i      = {v_i:.4f}\nv_b      = {v_b:.4f}")
    print(f"v_perp   = {v_perp:.4f}\nAlfven = {v_a:.6f}\n")


    # ================ #
    # plot preparation #
    # ================ #

    # ---------------------------#
    # plot layout and parameters #
    # ---------------------------#
    cm = 1 / 2.54
    fig = plt.figure(figsize=(8.5 * cm, 6.0 * cm), facecolor="white")
    ax1 = fig.add_subplot(211, xlabel="")
    ax2 = fig.add_subplot(212, xlabel="")
    fig.subplots_adjust(bottom=0.2, left=0.2, top=0.88, right=0.9, hspace=0.1)

    # Storage arrays for unstable mode tracks
    k_keep, w_keep, g_keep = [], [], []
    all_roots_k, all_roots_real = [], []

    # ------------------------------------------------ #
    # loop over wavenumbers, root-finding and plotting #
    # ------------------------------------------------ #

    # --------------------#
    # root solving engine #
    # --------------------#
    for n, k_norm in enumerate(karr):

        #----------------------------------------------#
        # k is changed into kc by multiplying omega_pi #
        #----------------------------------------------#
        k = k_norm * omega_pi

        # ----------------------------------------- #
        # frequency factors for dispersion equation #
        # ----------------------------------------- #
        pk = np.poly1d([-k**2])

        pde = np.poly1d([1, -k * v_e])
        pdec = pde + omega_ce

        pdi = np.poly1d([1, -k * v_i])
        pdic = pdi + omega_ci

        pdb = np.poly1d([1, -k * v_b])
        pdbc = pdb + omega_cb

        poe = np.poly1d([omega_pe**2])
        poi = np.poly1d([omega_pi**2])
        pob = np.poly1d([omega_pb**2])

        # --------------------------------------------------- #
        # determinant of dispersion matrix in polynomial form #
        # --------------------------------------------------- #
        p = (
            pk * pdec * pdic * pdbc**2
            - poe * pde * pdic * pdbc**2
            - poi * pdi * pdec * pdbc**2
            - pob * (pdb * pdec * pdic * pdbc + vph * k**2 * pdec * pdic)
        )

        # ------------------ #
        # numpy root-finding #
        # ------------------ #
        nstellen = np.roots(p)

        # -------------------------#
        # map current step outputs #
        #--------------------------#
        for z in nstellen:
            all_roots_k.append(k / omega_pi)
            all_roots_real.append(z.real)
            if z.imag > 1.0e-4:
                k_keep.append(k / omega_pi)
                w_keep.append(z.real)
                g_keep.append(z.imag)

        # ------------- #
        # end_of_k_loop #
        # ------------- #

    # -------------------#
    # generate plot rows #
    # -------------------#
    ax1.plot(all_roots_k, all_roots_real, marker=".", color="0.6", ls="", ms=0.8)
    ax1.plot(k_keep, w_keep, marker=".", color="k", ls="", ms=1.2)
    ax2.plot(k_keep, g_keep, marker=".", color="k", ls="", ms=1.2)

    # --------------------#
    # axes bounds control #
    # --------------------#
    xmin, xmax = kcw01_min, kcw01_max
    ax1.set_xlim(xmin, xmax)
#    ax1.set_ylim(-1.5, 2.5)
    ax2.set_xlim(xmin, xmax)
#    ax2.set_ylim(0.0, 0.50)

#    xticks = [-0.2, 0, 0.2]
#    ax1.set_xticks(xticks)
#    ax1.set_yticks([-1, 0, 1, 2])
#    ax2.set_xticks(xticks)
#    ax2.set_yticks([0.0, 0.2, 0.4])
    ax1.set_xticklabels([])

    # Ticks density intervals
    ax1.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax2.xaxis.set_minor_locator(MultipleLocator(0.02))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax2.yaxis.set_minor_locator(MultipleLocator(0.05))

    for ax in (ax1, ax2):
        ax.tick_params(
            axis="x",
            which="both",
            length=4 if ax == ax2 else 2,
            labelbottom=ax == ax2,
            bottom=True,
            top=True,
            direction="in",
        )
        ax.tick_params(
            axis="y", which="both", left=True, right=True, direction="in"
        )

    # -----------------------#
    # labels and designations#
    # -----------------------#
#    title_text = (
#        r"$v_{\|\mathrm{b}}=20\,V_\mathrm{A}, v_{\perp\mathrm{b}}=20\,V_\mathrm{A}$"
#    )

    title_text = (
        rf"$v_{{\|\mathrm{{b}}}}={v_b_out:02.0f}\,V_\mathrm{{A}},  "
        rf"v_{{\perp\mathrm{{b}}}}={v_perp_out:02.0f}\,V_\mathrm{{A}}$"
    )

    ax1.set_title(title_text, fontsize=10.5)
    ax1.set_ylabel(r"$\omega/\Omega_\mathrm{i}$")
    ax2.set_xlabel(r"$kc/\omega_\mathrm{pi}$")
    ax2.set_ylabel(r"$\gamma/\Omega_\mathrm{i}$")

    # -----------------------------#
    # numpy save arrays in npz format #
    # -----------------------------#
    np.savez(
        outfile,
        wavenum=np.array(all_roots_k),
        freq_real=np.array(all_roots_real),
        unstable_k=np.array(k_keep),
        unstable_real=np.array(w_keep),
        unstable_imag=np.array(g_keep),
    )

    plt.tight_layout()
    plt.savefig(pdffile)
    print(f"Data saved to {outfile}.npz and plot generated as {pdffile}")


if __name__ == "__main__":
    disp_eib()

