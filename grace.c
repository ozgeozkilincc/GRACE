/*
 * grace.c
 * =======
 * Software: GRACE - Graetz Rapid Analysis for Channel Extraction
 * Version 1.0
 *
 * CITATION
 * --------
 * If you use this software for research purposes, please cite:
 *
 *   O. Ozkilinc, S. de-la-Huerta, V. Diez-Cabanes, A. Gutierrez-Vega,
 *   P. A. Marcos-Villa, A. Bol, M. Atilhan, S. Aparicio,
 *   Microfluidic Extraction of PFAS Using Natural Deep Eutectic Solvent
 *   Thin Films: Theory, Design, and Sustainability,
 *   Chemical Engineering Journal Advances (2026) 101461.
 *   https://doi.org/10.1016/j.ceja.2026.101461
 *
 * Contact: sapar@ubu.es
 *
 * LICENSE
 * -------
 * MIT License
 * Copyright (c) 2026 Ozge Ozkilinc (oox@ubu.es)
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
 * BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
 * ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
 * CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * DISCLAIMER
 * ----------
 * All outputs are idealised upper-bound estimates of single-pass
 * performance, computed under idealised assumptions (T = 25 C, dilute
 * loading, uniform film, fully developed laminar flow, single-component
 * PFAS). Validate experimentally before device fabrication.
 *
 * VERIFICATION (reference design, see paper)
 * ------------------------------------------
 * d = 50 um, L = 5 cm, Q = 0.07 uL/s, delta = 3 um, PFOS, Men:DecA
 *   Re     = 2.27
 *   dP     = 33.9 kPa
 *   Bi     = 2448
 *   Sh_eff = 3.654
 *   eta    = 96.5 %
 *   S      = 1298  (TBAC:DecA at same conditions: S = 224)
 *
 * COMPILE
 * -------
 *   cc grace.c -o grace -lm -O3
 *
 * USAGE
 * -----
 *   ./grace
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/* ── Physical constants ───────────────────────────────────────────────────── */
#define RHO_W   997.1        /* kg/m3  water density at 25 C                   */
#define MU_W    0.890e-3     /* Pa·s   water dynamic viscosity at 25 C         */
#define KM      1.86         /* Lévêque prefactor (Sh = KM*Gz^(1/3), Gz>10)    */
#define SH_INF  3.66         /* fully-developed Sherwood number asymptote      */
#define PI      3.14159265358979323846

/* ── Data structures ──────────────────────────────────────────────────────── */

typedef struct {
    char   name[64];
    double Dw;       /* aqueous diffusivity (m2/s)        */
    double K;        /* DES/water partition coefficient   */
} GRACECompound;

typedef struct {
    char   name[64];
    double R_D;      /* viscosity ratio mu_DES/mu_water   */
    double gamma;    /* interfacial tension (N/m)         */
} GRACESolvent;

typedef struct {
    /* ── inputs ── */
    double d;        /* channel diameter (m)              */
    double L;        /* channel length (m)                */
    double Q;        /* flow rate (m3/s)                  */
    double delta;    /* DES film thickness (m)            */
    int    N;        /* number of parallel channels       */
    double Dw;       /* aqueous diffusivity (m2/s)        */
    double K;        /* partition coefficient             */
    double R_D;      /* viscosity ratio                   */
    double gamma;    /* interfacial tension (N/m)         */
    /* ── outputs ── */
    double r_eff;    /* effective aqueous radius (m)      */
    double u;        /* mean velocity (m/s)               */
    double Re;       /* Reynolds number                   */
    double dP_kPa;   /* pressure drop (kPa)               */
    double tau_w;    /* wall shear stress (Pa)            */
    double tau_res;  /* residence time (s)                */
    double Gz;       /* Graetz number                     */
    double Sh;       /* mean Sherwood number              */
    double Bi;       /* Biot number                       */
    double Sh_eff;   /* effective Sherwood number         */
    double NTU;      /* number of transfer units          */
    double eta;      /* extraction efficiency (%)         */
    double S;        /* film safety factor                */
    double tau_crit; /* critical shear stress (Pa)        */
    double V_DES;    /* DES volume per channel (m3)       */
    double Q_tot;    /* total array throughput (L/h)      */
    double energy;   /* pumping energy (kWh/m3)           */
    int    stages_90;/* stages required for eta >= 90%   */
    int    stages_99;/* stages required for eta >= 99%   */
} GRACEResult;

/* ── Built-in compound library ────────────────────────────────────────────── */
/* Dw: Stokes-Einstein estimate at 25 C                                       */
/* K:  DES/water partition coefficient from Fan et al. (2024)                 */

static const GRACECompound GRACE_PFAS[] = {
    {"PFBS  (C4,  sulfonate)",   6.79e-10,  600.0},
    {"PFHxS (C6,  sulfonate)",   5.46e-10, 1800.0},
    {"PFOS  (C8,  sulfonate)",   5.27e-10, 2100.0},
    {"PFDS  (C10, sulfonate)",   4.73e-10,  200.0},
    {"PFDoA (C12, carboxylate)", 4.16e-10,   39.0},
};
#define GRACE_N_PFAS 5

/* ── Built-in solvent library ─────────────────────────────────────────────── */
/* R_D:   mu_DES/mu_water at 25 C measured by rheometry                       */
/* gamma: DES-water interfacial tension (N/m) by pendant-drop tensiometry     */

static const GRACESolvent GRACE_NADES[] = {
    {"Men:Thy   (menthol:thymol)",                39.6,  30e-3},
    {"Men:DecA  (menthol:decanoic acid)",          14.3,  29e-3},
    {"Thy:DecA  (thymol:decanoic acid)",           12.5,  25e-3},
    {"TBAC:DecA (tetrabutylammonium:decanoic acid)",298.0, 5e-3},
};
#define GRACE_N_NADES 4

/* ── Physics functions ────────────────────────────────────────────────────── */

/* Effective aqueous-core radius (m). Eq. 1. */
double grace_r_eff(double d, double delta)
{
    return d / 2.0 - delta;
}

/* Mean axial velocity for Poiseuille flow (m/s). Eq. 2. */
double grace_u_bar(double Q, double r_eff)
{
    return Q / (PI * r_eff * r_eff);
}

/* Reynolds number based on hydraulic diameter 2*r_eff. Eq. 3. */
double grace_reynolds(double u, double r_eff)
{
    return RHO_W * u * 2.0 * r_eff / MU_W;
}

/* Hagen-Poiseuille pressure drop (kPa). Eq. 4. */
double grace_dP_kPa(double Q, double L, double r_eff)
{
    return 8.0 * MU_W * Q * L / (PI * pow(r_eff, 4)) / 1e3;
}

/* Wall shear stress at the DES-water interface (Pa). Eq. 5. */
double grace_tau_wall(double Q, double r_eff)
{
    return 4.0 * MU_W * Q / (PI * pow(r_eff, 3));
}

/* Graetz number (dimensionless). Eq. 6. */
double grace_graetz(double r_eff, double u, double Dw, double L)
{
    return r_eff * r_eff * u / (Dw * L);
}

/*
 * Mean Sherwood number. Eq. 7.
 * Gz > 10: Lévêque entrance-region, Sh = 1.86 * Gz^(1/3)
 * Gz <= 10: fully developed asymptote, Sh = 3.66
 */
double grace_sherwood(double Gz)
{
    return (Gz > 10.0) ? KM * pow(Gz, 1.0/3.0) : SH_INF;
}

/*
 * Biot number (dimensionless). Eq. 10.
 * Bi = K * d / (R_D * delta)
 * Bi >> Sh_inf: aqueous-side limited (film resistance negligible)
 * Bi << Sh_inf: DES-film limited (film diffusion is the bottleneck)
 */
double grace_biot(double K, double d, double R_D, double delta)
{
    return K * d / (R_D * delta);
}

/*
 * Effective Sherwood number accounting for DES-film resistance. Eq. 11.
 * Sh_eff = Sh * Bi / (Sh + Bi)
 */
double grace_sh_eff(double Sh, double Bi)
{
    return Sh * Bi / (Sh + Bi);
}

/* Number of transfer units. Eq. 13. */
double grace_NTU(double Sh_eff, double Dw, double L, double d, double u)
{
    return 4.0 * Sh_eff * Dw * L / (d * d * u);
}

/* Single-pass extraction efficiency (%). Eq. 14. */
double grace_eta(double ntu)
{
    return (1.0 - exp(-ntu)) * 100.0;
}

/*
 * Cumulative cascade efficiency (%). Eq. 23.
 * eta_n = 1 - (1 - eta1/100)^n
 * Assumes fresh DES and identical geometry per stage.
 */
double grace_cascade(double eta1, int n)
{
    return (1.0 - pow(1.0 - eta1 / 100.0, n)) * 100.0;
}

/* Critical shear stress for DES film displacement (Pa). Eq. 17. */
double grace_tau_crit(double gamma_iw, double delta)
{
    return gamma_iw / delta;
}

/*
 * Film stability safety factor (dimensionless). Eq. 18.
 * S = tau_crit / tau_w
 * S >> 1: film stable; S < 1: film at risk of displacement.
 */
double grace_safety(double gamma_iw, double delta, double tau_w)
{
    return grace_tau_crit(gamma_iw, delta) / tau_w;
}

/* Mass-transfer regime classification from Bi. */
const char *grace_regime(double Bi)
{
    if (Bi > 100.0) return "Aqueous-side limited  (Bi >> Sh_inf)";
    if (Bi > 10.0)  return "Weakly Bi-sensitive   (Bi ~ 10-100) ";
    if (Bi > 1.0)   return "Bi-sensitive          (Bi ~ 1-10)   ";
    return                  "DES-film limited      (Bi << Sh_inf)";
}

/* ── Compute all outputs ──────────────────────────────────────────────────── */

int grace_compute(GRACEResult *r)
{
    r->r_eff = grace_r_eff(r->d, r->delta);
    if (r->r_eff <= 0.0) {
        fprintf(stderr,
            "\n  GRACE ERROR: delta=%.1f um >= d/2=%.1f um.\n"
            "  Aqueous core fully occluded. Choose smaller delta or larger d.\n",
            r->delta * 1e6, r->d * 1e6 / 2.0);
        return -1;
    }

    r->u        = grace_u_bar(r->Q, r->r_eff);
    r->Re       = grace_reynolds(r->u, r->r_eff);
    r->dP_kPa   = grace_dP_kPa(r->Q, r->L, r->r_eff);
    r->tau_w    = grace_tau_wall(r->Q, r->r_eff);
    r->tau_res  = r->L / r->u;
    r->Gz       = grace_graetz(r->r_eff, r->u, r->Dw, r->L);
    r->Sh       = grace_sherwood(r->Gz);
    r->Bi       = grace_biot(r->K, r->d, r->R_D, r->delta);
    r->Sh_eff   = grace_sh_eff(r->Sh, r->Bi);
    r->NTU      = grace_NTU(r->Sh_eff, r->Dw, r->L, r->d, r->u);
    r->eta      = grace_eta(r->NTU);
    r->tau_crit = grace_tau_crit(r->gamma, r->delta);
    r->S        = grace_safety(r->gamma, r->delta, r->tau_w);
    r->V_DES    = PI * r->d * r->delta * r->L;
    r->Q_tot    = r->N * r->Q * 3600.0 * 1000.0;
    r->energy   = r->dP_kPa / 3600.0;

    r->stages_90 = 1;
    r->stages_99 = 1;
    if (r->eta < 100.0) {
        for (int n = 1; n <= 50; n++) {
            double en = grace_cascade(r->eta, n);
            if (en >= 90.0 && r->stages_90 == 1) r->stages_90 = n;
            if (en >= 99.0) { r->stages_99 = n; break; }
        }
    }
    return 0;
}

/* ── Display helpers ──────────────────────────────────────────────────────── */

#define GRACE_LINE  "────────────────────────────────────────────────────────────\n"
#define GRACE_DLINE "════════════════════════════════════════════════════════════\n"

void grace_section(const char *title)
{
    printf("\n  ▸ %s\n", title);
    printf("  " GRACE_LINE);
}

void grace_row(const char *label, const char *value,
               const char *unit,  const char *note)
{
    printf("  %-32s %-14s %-12s %s\n", label, value, unit, note);
}

void grace_eta_bar(double eta_val)
{
    int filled = (int)(eta_val / 2.5);
    if (filled > 40) filled = 40;
    printf("  [");
    for (int i = 0; i < filled; i++) printf("█");
    for (int i = filled; i < 40; i++) printf("░");
    printf("]  %.2f%%\n", eta_val);
}

/* ── Input helpers ────────────────────────────────────────────────────────── */

double grace_get_double(const char *prompt, double def)
{
    char   buf[64];
    double val;
    printf("  → %s [%.4g]: ", prompt, def);
    fflush(stdout);
    if (fgets(buf, sizeof(buf), stdin) == NULL) return def;
    if (buf[0] == '\n') return def;
    if (sscanf(buf, "%lf", &val) == 1) return val;
    return def;
}

int grace_get_int(const char *prompt, int def)
{
    char buf[64];
    int  val;
    printf("  → %s [%d]: ", prompt, def);
    fflush(stdout);
    if (fgets(buf, sizeof(buf), stdin) == NULL) return def;
    if (buf[0] == '\n') return def;
    if (sscanf(buf, "%d", &val) == 1) return val;
    return def;
}

int grace_menu(int n_options)
{
    char buf[16];
    int  val;
    printf("  → Select [1-%d]: ", n_options);
    fflush(stdout);
    if (fgets(buf, sizeof(buf), stdin) == NULL) return 1;
    if (sscanf(buf, "%d", &val) == 1 &&
        val >= 1 && val <= n_options) return val;
    return 1;
}

/* ── Save results to file ─────────────────────────────────────────────────── */

void grace_save(const GRACEResult *r,
                const char *pfas_name,
                const char *nades_name)
{
    char fname[128];
    char pfas_short[16] = {0};
    sscanf(pfas_name, "%15s", pfas_short);
    snprintf(fname, sizeof(fname), "grace_result_%s.txt", pfas_short);

    FILE *f = fopen(fname, "w");
    if (!f) {
        printf("  GRACE ERROR: Cannot open %s for writing.\n", fname);
        return;
    }
    fprintf(f, "===============================================================\n");
    fprintf(f, "GRACE v1.0 -- Graetz Rapid Analysis for Channel Extraction\n");
    fprintf(f, "===============================================================\n");
    fprintf(f, "\n");
    fprintf(f, "CITATION\n");
    fprintf(f, "---------------------------------------------------------------\n");
    fprintf(f, "If you use GRACE for research purposes, please cite:\n");
    fprintf(f, "\n");
    fprintf(f, "  Ozkilinc, O.; de-la-Huerta, S.; Diez-Cabanes, V.;\n");
    fprintf(f, "  Gutierrez-Vega, A.; Marcos-Villa, P. A.; Bol, A.;\n");
    fprintf(f, "  Atilhan, M.; Aparicio, S. Microfluidic Extraction of PFAS\n");
    fprintf(f, "  Using Natural Deep Eutectic Solvent Thin Films: Theory,\n");
    fprintf(f, "  Design, and Sustainability. Chemical Engineering Journal\n");
    fprintf(f, "  Advances (2026) 101461.\n");
    fprintf(f, "  https://doi.org/10.1016/j.ceja.2026.101461\n");
    fprintf(f, "---------------------------------------------------------------\n");
    fprintf(f, "  Lead software author: Ozge Ozkilinc\n");
    fprintf(f, "  Corresponding author (paper): Santiago Aparicio, sapar@ubu.es\n");
    fprintf(f, "===============================================================\n");
    fprintf(f, "PFAS:    %s\n",   pfas_name);
    fprintf(f, "NADES:   %s\n",   nades_name);
    fprintf(f, "d        = %.1f um\n",    r->d      * 1e6);
    fprintf(f, "L        = %.1f cm\n",    r->L      * 1e2);
    fprintf(f, "Q        = %.4f uL/s\n",  r->Q      * 1e9);
    fprintf(f, "delta    = %.1f um\n",    r->delta  * 1e6);
    fprintf(f, "N        = %d channels\n\n", r->N);
    fprintf(f, "Dw       = %.3e m2/s\n",  r->Dw);
    fprintf(f, "K        = %.0f\n",        r->K);
    fprintf(f, "R_D      = %.1f\n",        r->R_D);
    fprintf(f, "r_eff    = %.2f um\n",     r->r_eff  * 1e6);
    fprintf(f, "u_bar    = %.4f cm/s\n",   r->u      * 1e2);
    fprintf(f, "Re       = %.2f\n",         r->Re);
    fprintf(f, "dP       = %.3f kPa\n",    r->dP_kPa);
    fprintf(f, "tau_w    = %.4f Pa\n",     r->tau_w);
    fprintf(f, "tau_res  = %.3f s\n",      r->tau_res);
    fprintf(f, "Gz       = %.4f\n",         r->Gz);
    fprintf(f, "Sh_bar   = %.4f\n",         r->Sh);
    fprintf(f, "Bi       = %.2f\n",         r->Bi);
    fprintf(f, "Sh_eff   = %.4f\n",         r->Sh_eff);
    fprintf(f, "NTU      = %.4f\n",         r->NTU);
    fprintf(f, "eta      = %.3f %%\n",     r->eta);
    fprintf(f, "S        = %.1f\n",         r->S);
    fprintf(f, "tau_crit = %.4f Pa\n",     r->tau_crit);
    fprintf(f, "Regime   = %s\n",           grace_regime(r->Bi));
    fprintf(f, "Stages for eta>=90%% = %d\n",  r->stages_90);
    fprintf(f, "Stages for eta>=99%% = %d\n",  r->stages_99);
    fprintf(f, "dP for eta>=99%%     = %.2f kPa\n",r->dP_kPa * r->stages_99);
    fprintf(f, "Q_tot    = %.4f L/h\n",    r->Q_tot);
    fprintf(f, "Energy   = %.4f x1e-3 kWh/m3\n", r->energy * 1e3);
    fprintf(f, "V_DES/ch = %.3f fL\n",     r->V_DES  * 1e18);
    fprintf(f, "\nNOTE: eta is an idealised upper bound on single-pass\n");
    fprintf(f, "performance; energy, cost, and waste are lower-bound burden\n");
    fprintf(f, "estimates. Assumptions: T=25C, dilute loading, uniform film,\n");
    fprintf(f, "fully developed laminar flow, single-component PFAS, no\n");
    fprintf(f, "prototype fabricated. Validate experimentally before device\n");
    fprintf(f, "fabrication.\n");

    fclose(f);
    printf("  Saved to %s\n", fname);
}

/* ── Main session ─────────────────────────────────────────────────────────── */

void grace_run(void)
{
    GRACEResult r;
    char pfas_name[64];
    char nades_name[64];

    printf("\n");
    printf("  " GRACE_DLINE);
    printf("  GRACE v1.0\n");
    printf("  Graetz Rapid Analysis for Channel Extraction\n");
    printf("  Biot-modified Graetz model\n");
    printf("  \n");
    printf("  Supporting Information for:\n");
    printf("  Ozkilinc, O.; de-la-Huerta, S.; Diez-Cabanes, V.;\n");
    printf("  Gutierrez-Vega, A.; Marcos-Villa, P. A.; Bol, A.; Atilhan, M.;\n");
    printf("  Aparicio, S. Microfluidic Extraction of PFAS Using Natural\n");
    printf("  Deep Eutectic Solvent Thin Films: Theory, Design, and\n");
    printf("  Sustainability. Chemical Engineering Journal Advances (2026) 101461.\n");
    printf("  https://doi.org/10.1016/j.ceja.2026.101461\n");
    printf("  \n");
    printf("  Lead software author: Ozge Ozkilinc  |  University of Burgos\n");
    printf("  Corresponding author (paper): Santiago Aparicio, sapar@ubu.es\n");
    printf("  " GRACE_DLINE);
    printf("\n");
    printf("  LICENCE:  MIT — free to use, modify, and redistribute.\n");
    printf("  CITATION: Ozkilinc et al., Chem. Eng. J. Adv. 2026, 101461.\n");
    printf("\n");
    printf("  MODEL ASSUMPTIONS:\n");
    printf("  - Dilute loading: C_DES << K*C0\n");
    printf("  - Uniform DES film of thickness delta\n");
    printf("  - Isothermal operation at T = 25 C\n");
    printf("  - Fully developed laminar flow (Re << 2100)\n");
    printf("  - Single-component PFAS (no competitive partitioning)\n");
    printf("  - eta is an idealised UPPER bound; energy/cost/waste are\n");
    printf("    LOWER-bound burden estimates. No prototype fabricated.\n");
    printf("\n");
    printf("  VERIFICATION (reference design, PFOS + Men:DecA):\n");
    printf("  d=50um  L=5cm  Q=0.07uL/s  delta=3um\n");
    printf("  Re=2.27 | dP=33.9kPa | Bi=2448 | eta=96.5%% | S=1298\n");
    printf("  (TBAC:DecA worst-case film stability at same conditions: S=224)\n");
    printf("\n");
    printf("  Press Enter at any prompt to accept the default value.\n");
    printf("\n");

    /* PFAS selection */
    printf("  PFAS Compound:\n");
    for (int i = 0; i < GRACE_N_PFAS; i++)
        printf("  [%d]  %s\n", i+1, GRACE_PFAS[i].name);
    printf("  [%d]  Custom compound\n", GRACE_N_PFAS+1);
    int pfas_idx = grace_menu(GRACE_N_PFAS + 1);

    if (pfas_idx <= GRACE_N_PFAS) {
        r.Dw = GRACE_PFAS[pfas_idx-1].Dw;
        r.K  = GRACE_PFAS[pfas_idx-1].K;
        strncpy(pfas_name, GRACE_PFAS[pfas_idx-1].name,
                sizeof(pfas_name)-1);
        pfas_name[sizeof(pfas_name)-1] = '\0';
    } else {
        printf("\n  Enter custom PFAS properties:\n");
        printf("  (Dw from Stokes-Einstein; K from batch DES/water extraction)\n");
        double Dw_in = grace_get_double(
            "Aqueous diffusivity Dw (x1e-10 m2/s)", 5.27);
        r.Dw = Dw_in * 1e-10;
        r.K  = grace_get_double("Partition coefficient K", 2100.0);
        strncpy(pfas_name, "Custom", sizeof(pfas_name)-1);
    }

    /* NADES selection */
    printf("\n  NADES System:\n");
    for (int i = 0; i < GRACE_N_NADES; i++)
        printf("  [%d]  %s  (R_D=%.1f)\n",
               i+1, GRACE_NADES[i].name, GRACE_NADES[i].R_D);
    printf("  [%d]  Custom NADES\n", GRACE_N_NADES+1);
    int nades_idx = grace_menu(GRACE_N_NADES + 1);

    if (nades_idx <= GRACE_N_NADES) {
        r.R_D   = GRACE_NADES[nades_idx-1].R_D;
        r.gamma = GRACE_NADES[nades_idx-1].gamma;
        strncpy(nades_name, GRACE_NADES[nades_idx-1].name,
                sizeof(nades_name)-1);
        nades_name[sizeof(nades_name)-1] = '\0';
    } else {
        printf("\n  Enter custom NADES properties:\n");
        printf("  (R_D = mu_DES/mu_water measured by rheometry at 25 C)\n");
        printf("  (gamma measured by pendant-drop tensiometry)\n");
        r.R_D   = grace_get_double("Viscosity ratio R_D", 14.3);
        double g = grace_get_double(
            "Interfacial tension gamma_iw (mN/m)", 29.0);
        r.gamma = g * 1e-3;
        strncpy(nades_name, "Custom", sizeof(nades_name)-1);
    }

    /* Geometry and flow */
    printf("\n  " GRACE_DLINE);
    printf("  Channel Geometry and Operating Conditions\n");
    printf("  Practical ranges: d=15-200um | L=0.3-10cm |\n");
    printf("  Q=0.01-0.5uL/s | delta=0.5-10um | dP < 50kPa\n");
    printf("  " GRACE_DLINE);
    double d_um     = grace_get_double("Channel diameter d (um)",      50.0);
    double L_cm     = grace_get_double("Channel length L (cm)",         5.0);
    double Q_uls    = grace_get_double("Flow rate Q (uL/s)",            0.07);
    double delta_um = grace_get_double("DES film thickness delta (um)", 3.0);
    r.N             = grace_get_int(   "Number of parallel channels N",    1);

    r.d     = d_um     * 1e-6;
    r.L     = L_cm     * 1e-2;
    r.Q     = Q_uls    * 1e-9;
    r.delta = delta_um * 1e-6;

    if (grace_compute(&r) != 0) return;

    if (r.delta >= 0.35 * (r.d / 2.0)) {
        printf("\n  GRACE WARNING: delta/r = %.2f > 0.35.\n",
               r.delta / (r.d / 2.0));
        printf("  Film occupies >35%% of channel radius.\n");
        printf("  Model validity may be reduced at this film fraction.\n\n");
    }

    /* Results */
    printf("\n");
    printf("  " GRACE_DLINE);
    printf("  GRACE v1.0 — Results\n");
    printf("  " GRACE_DLINE);

    char buf[64];

    grace_section("Selected System");
    grace_row("PFAS compound",     pfas_name,  "", "");
    grace_row("NADES system",      nades_name, "", "");
    snprintf(buf, sizeof(buf), "%.2f", r.Dw * 1e10);
    grace_row("Dw", buf, "x1e-10 m2/s", "aqueous diffusivity at 25 C");
    snprintf(buf, sizeof(buf), "%.0f", r.K);
    grace_row("K (DES/water)", buf, "--", "partition coefficient");
    snprintf(buf, sizeof(buf), "%.1f", r.R_D);
    grace_row("R_D (mu_DES/mu_w)", buf, "--", "viscosity ratio");
    snprintf(buf, sizeof(buf), "%.1f", r.gamma * 1e3);
    grace_row("gamma_iw", buf, "mN/m", "interfacial tension");

    grace_section("Geometry and Hydrodynamics  [Eqs. 1-5]");
    snprintf(buf, sizeof(buf), "%.2f", r.r_eff * 1e6);
    grace_row("R_eff", buf, "um", "= d/2 - delta");
    snprintf(buf, sizeof(buf), "%.4f", r.u * 1e2);
    grace_row("Mean velocity u_bar", buf, "cm/s", "");
    snprintf(buf, sizeof(buf), "%.2f", r.Re);
    grace_row("Reynolds number Re", buf, "--",
              r.Re < 2100 ? "laminar OK" : "TURBULENT -- model invalid!");
    snprintf(buf, sizeof(buf), "%.3f", r.dP_kPa);
    grace_row("Pressure drop dP", buf, "kPa",
              r.dP_kPa < 50.0 ? "within 50 kPa limit" :
                                 "EXCEEDS 50 kPa ceiling");
    snprintf(buf, sizeof(buf), "%.4f", r.tau_w);
    grace_row("Wall shear stress tau_w", buf, "Pa", "");
    snprintf(buf, sizeof(buf), "%.3f", r.tau_res);
    grace_row("Residence time tau", buf, "s", "= L/u_bar");

    grace_section("Aqueous-Side Mass Transfer  [Eqs. 6-8]");
    snprintf(buf, sizeof(buf), "%.4f", r.Gz);
    grace_row("Graetz number Gz", buf, "--",
              r.Gz <= 10.0 ? "Gz<=10  fully developed  Sh=3.66" :
                             "Gz>10   entrance region  Sh=1.86*Gz^(1/3)");
    snprintf(buf, sizeof(buf), "%.4f", r.Sh);
    grace_row("Mean Sherwood Sh_bar", buf, "--", "");

    grace_section("DES-Film Mass Transfer  [Eqs. 9-12]");
    snprintf(buf, sizeof(buf), "%.2f", r.Bi);
    grace_row("Biot number Bi", buf, "--", grace_regime(r.Bi));
    snprintf(buf, sizeof(buf), "%.4f", r.Sh_eff);
    grace_row("Sh_eff", buf, "--", "= Sh_bar*Bi/(Sh_bar+Bi)");

    grace_section("Overall Extraction Performance  [Eqs. 13-15]");
    snprintf(buf, sizeof(buf), "%.4f", r.NTU);
    grace_row("NTU", buf, "--", "= 4*Sh_eff*Dw*L/(d^2*u)");
    printf("\n  Extraction efficiency eta:\n  ");
    grace_eta_bar(r.eta);
    printf("\n");
    snprintf(buf, sizeof(buf), "eta = %.3f%%", r.eta);
    grace_row("", buf, "",
              r.eta >= 90.0 ? "Above 90%% target" :
              r.eta >= 80.0 ? "Near target -- consider longer L or smaller d" :
                              "Below target -- increase L, reduce d or Q");

    grace_section("DES Film Stability  [Eqs. 17-18]");
    snprintf(buf, sizeof(buf), "%.4f", r.tau_crit);
    grace_row("tau_crit = gamma_iw/delta", buf, "Pa",
              "threshold for film displacement");
    snprintf(buf, sizeof(buf), "%.1f", r.S);
    grace_row("Safety factor S", buf, "--",
              r.S >= 10.0 ? "Film stable" :
              r.S >= 1.0  ? "Marginal -- verify experimentally" :
                            "Film at risk of displacement!");

    grace_section("Array Performance");
    snprintf(buf, sizeof(buf), "%.4f", r.Q_tot);
    grace_row("Total throughput", buf, "L/h", "");
    snprintf(buf, sizeof(buf), "%.4f", r.energy * 1e3);
    grace_row("Energy consumption", buf, "x1e-3 kWh/m3",
              "cf. GAC=0.5, NF/RO=3.5, AOPs=15 kWh/m3");
    snprintf(buf, sizeof(buf), "%.3f", r.V_DES * 1e18);
    grace_row("DES volume/channel", buf, "fL", "= pi*d*delta*L");
    snprintf(buf, sizeof(buf), "%.4f", r.V_DES * r.N * 1e15);
    grace_row("DES volume total", buf, "pL", "");

    grace_section("Multi-Stage Cascade Guidance  [Eq. 23]");
    printf("  eta_n = 1 - (1-eta1)^n  |  fresh DES assumed per stage\n");
    printf("  dP_tot = n*dP_stage  |  repressurisation if dP_tot > 50 kPa\n\n");
    snprintf(buf, sizeof(buf), "%d", r.stages_90);
    grace_row("Stages for eta >= 90%%", buf, "stage(s)", "");
    snprintf(buf, sizeof(buf), "%d", r.stages_99);
    grace_row("Stages for eta >= 99%%", buf, "stage(s)", "");
    snprintf(buf, sizeof(buf), "%.2f", r.dP_kPa);
    grace_row("dP per stage", buf, "kPa", "");
    snprintf(buf, sizeof(buf), "%.2f", r.dP_kPa * r.stages_99);
    grace_row("dP for eta >= 99%%", buf, "kPa",
              r.dP_kPa * r.stages_99 < 50.0 ?
              "single pump sufficient" :
              "repressurisation required between stages");

    /* Summary */
    printf("\n  " GRACE_DLINE);
    printf("  GRACE SUMMARY\n");
    printf("  %s | %s\n", pfas_name, nades_name);
    printf("  d=%.0fum  L=%.1fcm  Q=%.3fuL/s  delta=%.1fum  N=%d\n",
           d_um, L_cm, Q_uls, delta_um, r.N);
    printf("  Re=%.2f  Gz=%.4f  Bi=%.1f  eta=%.2f%%  dP=%.2fkPa  S=%.0f\n",
           r.Re, r.Gz, r.Bi, r.eta, r.dP_kPa, r.S);
    printf("  eta is an idealised upper bound; energy/cost/waste are lower-\n");
    printf("  bound estimates (T=25C, dilute loading, uniform film).\n");
    printf("  No prototype fabricated -- validate experimentally.\n");
    printf("  " GRACE_DLINE);

    /* Save */
    printf("\n  → Save results to file? [y/N]: ");
    fflush(stdout);
    char choice[8];
    if (fgets(choice, sizeof(choice), stdin) &&
        (choice[0] == 'y' || choice[0] == 'Y'))
        grace_save(&r, pfas_name, nades_name);
}

/* ── Entry point ──────────────────────────────────────────────────────────── */

int main(void)
{
    char again[8];
    do {
        grace_run();
        printf("\n  → Run another calculation? [y/N]: ");
        fflush(stdout);
        if (!fgets(again, sizeof(again), stdin)) break;
    } while (again[0] == 'y' || again[0] == 'Y');

    printf("\n  GRACE session ended.\n");
    printf("  Please cite: Ozkilinc et al., Chem. Eng. J. Adv. 2026, 101461.\n");
    printf("  DOI: 10.1016/j.ceja.2026.101461 | Contact: sapar@ubu.es\n\n");
    return 0;
}
