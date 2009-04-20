/* 2004-2005 Futoshi Shimizu */
#define CUI_GLOBAL
#include "cui.h"
#ifdef USE_P3D
#   include "p3dp.h"
#endif
#include "A.h"

#include <sys/fcntl.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/stat.h>
#include <sys/un.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <readline/readline.h>
#include <readline/history.h>

static char buf[CUI_LINEMAX], cmdline[CUI_LINEMAX];
static char fname[FILENAME_MAX] = "";
static int quit = 0;
static int frontend = -1;
#define nonvacant(iw) (AX_cid[iw]!=0)

#ifdef ATOMEYE_LIB
static void (*atomeyelib_on_click_atom)(int atom);
static void (*atomeyelib_on_close)();
static void (*atomeyelib_on_advance)(char *instr);
#endif

double cui_wtime(void)
{
    struct timeval now;
    gettimeofday(&now, NULL);
    return now.tv_sec + now.tv_usec * 1.e-6;
}


static char *cui_stripspace(char *line)
{
    char *s, *t;

    s = line;
    while (*s && isspace(*s)) s++;
    if (!*s) return s;

    t = s + strlen(s) - 1;
    while (t > s && isspace(*t)) t--;
    *++t = 0;

    return s;
}

static char cui_show_syntax[] = "syntax: ";

static void cui_send(char *line, FILE *fp)
{
    char *s, sendbuf[CUI_LINEMAX];
    if (fp) {
        if (fileno(fp) == STDIN_FILENO) fp = stdout;
        strncpy(sendbuf, line, sizeof(sendbuf)-1);
        s = cui_stripspace(sendbuf);
        fputs(strcat(s, "\n"), fp);
        fflush(fp);
    }
}


static int cui_send_result(char *line, FILE *fp)
{
    int fd = fileno(fp);
    if (!line || !*line) {
        if (fd != frontend)
            cui_send(CUI_PROTOCOL_OK, fp);
        return 0;
    }
    else {
        int len, len0;
        char *s, sendbuf[CUI_LINEMAX];
        strncpy(sendbuf, line, sizeof(sendbuf)-1);
        s = cui_stripspace(sendbuf);
        len = strlen(s);
        len0 = len - 1;
        if (len0 > strlen(CUI_PROTOCOL_OK)&&
                strcmp(s+len0-strlen(CUI_PROTOCOL_OK),"\n"CUI_PROTOCOL_OK)==0) {
            if (fd == frontend)
                *(s+len-strlen(CUI_PROTOCOL_OK)) = 0;
            cui_send(s, fp);
            return 0;
        }
        if (len0 > strlen(CUI_PROTOCOL_NG)&&
                strcmp(s+len0-strlen(CUI_PROTOCOL_NG),"\n"CUI_PROTOCOL_NG)==0) {
            if (fd == frontend)
                *(s+len-strlen(CUI_PROTOCOL_NG)) = 0;
            cui_send(s, fp);
            return 1;
        }
        cui_send(s, fp);
        if (fd != frontend)
            cui_send(CUI_PROTOCOL_NG, fp);
        return 1;
    }
}


static int cui_recv_print(FILE *fp, FILE *out)
{
    char recvbuf[CUI_LINEMAX];
    if (fp) {
        while (fgets(recvbuf, sizeof(recvbuf), fp)) {
            if (!*recvbuf || strcmp(recvbuf, CUI_PROTOCOL_QUIT"\n") == 0)
                return -1;
            else if (strcmp(recvbuf, CUI_PROTOCOL_OK"\n") == 0)
                return 0;
            else if (strcmp(recvbuf, CUI_PROTOCOL_NG"\n") == 0)
                return 1;
            if (out) {
                fputs(recvbuf, out);
                fflush(out);
            }
        }
    }
    return -1;
}


static void stdin_gets(char *line)
{
    if (line) {
        strncpy(cmdline, line, sizeof(cmdline));
        free(line);
    }
}

static char *gui_readline_gets(int iw, char *prompt, char *danswer)
{
    static char answer[CUI_LINEMAX] = "";
    if (AX_display[iw]) {
        if (IS_MANAGER) {
            xterm_get_focus(iw);
            if (frontend >= 0)
                rl_callback_handler_remove();
            clear_stdin_buffer();
            strncpy(answer, readline_gets(prompt, danswer), sizeof(answer));
            if (frontend >= 0)
                rl_callback_handler_install(CUI_PROMPT, stdin_gets);
        }
#ifdef USE_P3D
    }
    if (p3dp_enabled) {
        int len = strlen(answer);
        p3d_bcast(p3dp_cell, &len, 1, MPI_INT, 0);
        if (len > 0)
            p3d_bcast(p3dp_cell, answer, len+1, MPI_CHAR, 0);
    }
    if (AX_display[iw]) {
#endif
        if (IS_MANAGER)
            xterm_release_focus(iw);
    }
    return answer;
}

/*******************************************************************/

struct aec {
    char *name;
    bool (*proc)(int iw, char *instr, char **outstr);
    char *instr;
    char *syntax;
};

struct gkt {
    KeySym keysym;
    struct aec *normal, *shift, *ctrl, *ctrl_shift, *meta, *meta_shift;
    char *n, *s, *c, *cs, *m, *ms;
};
#define RPT2(name) name, name
#define RPT3(name) name, name, name
#define RPT4(name) name, name, name, name
#define RPT6(name) name, name, name, name, name, name
#define GKT(keysym) keysym, RPT6(NULL)
static struct gkt gui_key_table[] = {

/* TTY Functions */
    {GKT(XK_BackSpace),},
    {GKT(XK_Tab),       RPT6("toggle_parallel_projection")},
    {GKT(XK_Linefeed),},
    {GKT(XK_Clear),},
    {GKT(XK_Return),    RPT6("finish_rcut_patch")},
    {GKT(XK_Pause),},
    {GKT(XK_Scroll_Lock),},
    {GKT(XK_Sys_Req),},
    {GKT(XK_Escape),    RPT6("isoatomic_reference_imprint")},
    {GKT(XK_Delete),    "load_config_forward", NULL, "load_config_last"},

/* Cursor control & motion */
    {GKT(XK_Home),      "change_bond_radius_inc",
                        "change_view_angle_amplification_inc",
                        "rcut_patch_inc"},
    {GKT(XK_End),       "change_bond_radius_dec",
                        "change_view_angle_amplification_dec",
                        "rcut_patch_dec"},
    {GKT(XK_Right),     "rotate_1_dec",             "shift_cutting_plane_inc",
                RPT2(   "translate_0_inc"), RPT2(   "shift_xtal_0_dec"      )},
    {GKT(XK_Left),      "rotate_1_inc",             "shift_cutting_plane_dec",
                RPT2(   "translate_0_dec"), RPT2(   "shift_xtal_0_inc"      )},
    {GKT(XK_Up),        "rotate_0_dec",             "rotate_2_inc",
                        "translate_1_dec",          "translate_2_inc",
                        "shift_xtal_1_inc",         "shift_xtal_2_dec"      },
    {GKT(XK_Down),      "rotate_0_inc",             "rotate_2_dec",
                        "translate_1_inc",          "translate_2_dec",
                        "shift_xtal_1_dec",         "shift_xtal_2_inc"      },
    {GKT(XK_Page_Up),   "change_atom_r_ratio_inc",
                        "change_aux_property_threshold_upper_inc",
                        "change_aux_property_threshold_lower_inc",
                        NULL, "advance_dec"},
    {GKT(XK_Page_Down), "change_atom_r_ratio_dec",
                        "change_aux_property_threshold_upper_dec",
                        "change_aux_property_threshold_lower_dec",
                        NULL, "advance_inc"},

/* Misc Functions */
    {GKT(XK_Select),},
    {GKT(XK_Print),},
    {GKT(XK_Execute),},
    {GKT(XK_Insert),    "load_config_backward", NULL, "load_config_first"},
    {GKT(XK_Undo),},
    {GKT(XK_Redo),},
    {GKT(XK_Menu),},
    {GKT(XK_Find),},
    {GKT(XK_Cancel),},
    {GKT(XK_Help),},
    {GKT(XK_Break),},
    {GKT(XK_Mode_switch),},

/* Keypad Functions, keypad numbers */
    {GKT(XK_KP_Home),   "change_bond_radius_inc"},
    {GKT(XK_KP_End),    "change_bond_radius_dec"},
    {GKT(XK_KP_Right),  "rotate_1_dec",             "shift_cutting_plane_inc",
                RPT2(   "translate_0_inc"), RPT2(   "shift_xtal_0_dec"      )},
    {GKT(XK_KP_Left),   "rotate_1_inc",             "shift_cutting_plane_dec",
                RPT2(   "translate_0_dec"), RPT2(   "shift_xtal_0_inc"      )},
    {GKT(XK_KP_Up),     "rotate_0_dec",             "rotate_2_inc",
                        "translate_1_dec",          "translate_2_inc",
                        "shift_xtal_1_inc",         "shift_xtal_2_dec"      },
    {GKT(XK_KP_Down),   "rotate_0_inc",             "rotate_2_dec",
                        "translate_1_inc",          "translate_2_dec",
                        "shift_xtal_1_dec",         "shift_xtal_2_inc"      },
    {GKT(XK_KP_Page_Up), "change_atom_r_ratio_inc",
                        "change_aux_property_threshold_upper_inc",
                        "change_aux_property_threshold_lower_inc"},
    {GKT(XK_KP_Page_Down), "change_atom_r_ratio_dec",
                        "change_aux_property_threshold_upper_dec",
                        "change_aux_property_threshold_lower_dec"},
    {GKT(XK_KP_Insert),},
    {GKT(XK_KP_Delete),},
    {GKT(XK_KP_Enter),},
    {GKT(XK_KP_Multiply),},
    {GKT(XK_KP_Add),},
    {GKT(XK_KP_Subtract),},
    {GKT(XK_KP_Decimal),},
    {GKT(XK_KP_Divide),},
    {GKT(XK_KP_0),},
    {GKT(XK_KP_1),},
    {GKT(XK_KP_2),},
    {GKT(XK_KP_3),},
    {GKT(XK_KP_4),},
    {GKT(XK_KP_5),},
    {GKT(XK_KP_6),},
    {GKT(XK_KP_7),},
    {GKT(XK_KP_8),},
    {GKT(XK_KP_9),},

/* Auxilliary Functions */
    {GKT(XK_F1),},
    {GKT(XK_F2),    RPT2("reset_scratch"), RPT4("free_scratch")},
    {GKT(XK_F3),    RPT6("scratch_coloring")},
    {GKT(XK_F4),},
    {GKT(XK_F5),},
    {GKT(XK_F6),},
    {GKT(XK_F7),    RPT6("capture_eps")},
    {GKT(XK_F8),    RPT6("find_atom")},
    {GKT(XK_F9),    RPT6("load_config")},
    {GKT(XK_F10),   RPT6("reload_config")},
    {GKT(XK_F11),   RPT6("load_aux")},
    {GKT(XK_F12),   RPT6("load_atom_color")},

/* Latin 1 */
    {GKT(XK_space),},
    {GKT(XK_exclam),},
    {GKT(XK_numbersign),},
    {GKT(XK_dollar),},
    {GKT(XK_percent),},
    {GKT(XK_parenright),},
    {GKT(XK_asterisk),},
    {GKT(XK_plus),},
    {GKT(XK_comma), "print_atom_info_pair"},
    {GKT(XK_period),"print_atom_info_triplet"},
    {GKT(XK_slash), "print_atom_info_quartet"},
    {GKT(XK_0), "select_gear_0",            "cutting_plane_0",
                "aux_property_coloring_20", "shift_cutting_plane_to_anchor_0",
                "aux_property_coloring_0",  "delete_cutting_plane_0"},
    {GKT(XK_1), "select_gear_1",            "cutting_plane_1",
                "aux_property_coloring_21", "shift_cutting_plane_to_anchor_1",
                "aux_property_coloring_1",  "delete_cutting_plane_1"},
    {GKT(XK_2), "select_gear_2",            "cutting_plane_2",
                "aux_property_coloring_22", "shift_cutting_plane_to_anchor_2",
                "aux_property_coloring_2",  "delete_cutting_plane_2"},
    {GKT(XK_3), "select_gear_3",            "cutting_plane_3",
                "aux_property_coloring_23", "shift_cutting_plane_to_anchor_3",
                "aux_property_coloring_3",  "delete_cutting_plane_3"},
    {GKT(XK_4), "select_gear_4",            "cutting_plane_4",
                "aux_property_coloring_24", "shift_cutting_plane_to_anchor_4",
                "aux_property_coloring_4",  "delete_cutting_plane_4"},
    {GKT(XK_5), "select_gear_5",            "cutting_plane_5",
                "aux_property_coloring_25", "shift_cutting_plane_to_anchor_5",
                "aux_property_coloring_5",  "delete_cutting_plane_5"},
    {GKT(XK_6), "select_gear_6",            "cutting_plane_6",
                "aux_property_coloring_26", "shift_cutting_plane_to_anchor_6",
                "aux_property_coloring_6",  "delete_cutting_plane_6"},
    {GKT(XK_7), "select_gear_7",            "cutting_plane_7",
                "aux_property_coloring_27", "shift_cutting_plane_to_anchor_7",
                "aux_property_coloring_7",  "delete_cutting_plane_7"},
    {GKT(XK_8), "select_gear_8",            "cutting_plane_8",
                "aux_property_coloring_28", "shift_cutting_plane_to_anchor_8",
                "aux_property_coloring_8",  "delete_cutting_plane_8"},
    {GKT(XK_9), "select_gear_9",            "cutting_plane_9",
                "aux_property_coloring_29", "shift_cutting_plane_to_anchor_9",
                "aux_property_coloring_9",  "delete_cutting_plane_9"},
    {GKT(XK_agrave),    "select_gear_0",    "cutting_plane_0",
                NULL,                       "shift_cutting_plane_to_anchor_0",
                "aux_property_coloring_0",  "delete_cutting_plane_0"},
    {GKT(XK_ampersand), "select_gear_1",    "cutting_plane_1",
                NULL,                       "shift_cutting_plane_to_anchor_1",
                "aux_property_coloring_1",  "delete_cutting_plane_1"},
    {GKT(XK_eacute),    "select_gear_2",    "cutting_plane_2",
                NULL,                       "shift_cutting_plane_to_anchor_2",
                "aux_property_coloring_2",  "delete_cutting_plane_2"},
    {GKT(XK_quotedbl),  "select_gear_3",    "cutting_plane_3",
                NULL,                       "shift_cutting_plane_to_anchor_3",
                "aux_property_coloring_3",  "delete_cutting_plane_3"},
    {GKT(XK_apostrophe),    "select_gear_4","cutting_plane_4",
                NULL,                       "shift_cutting_plane_to_anchor_4",
                "aux_property_coloring_4",  "delete_cutting_plane_4"},
    {GKT(XK_parenleft), "select_gear_5",    "cutting_plane_5",
                NULL,                       "shift_cutting_plane_to_anchor_5",
                "aux_property_coloring_5",  "delete_cutting_plane_5"},
    {GKT(XK_minus), "select_gear_6",        "cutting_plane_6",
                NULL,                       "shift_cutting_plane_to_anchor_6",
                RPT2("change_aux_colormap_prev")},
    {GKT(XK_egrave),    "select_gear_7",    "cutting_plane_7",
                NULL,                       "shift_cutting_plane_to_anchor_7",
                "aux_property_coloring_7",  "delete_cutting_plane_7"},
    {GKT(XK_underscore),    "select_gear_8","cutting_plane_8",
                NULL,                       "shift_cutting_plane_to_anchor_8",
                "aux_property_coloring_8",  "delete_cutting_plane_8"},
    {GKT(XK_ccedilla),  "select_gear_9",    "cutting_plane_9",
                NULL,                       "shift_cutting_plane_to_anchor_9",
                "aux_property_coloring_9",  "delete_cutting_plane_9"},

    {GKT(XK_colon),     "save_atom_indices"},
    {GKT(XK_semicolon), "save_atom_indices"},
    {GKT(XK_less),},
    {GKT(XK_equal), RPT4(NULL), RPT2("change_aux_colormap_next")},
    {GKT(XK_greater),},
    {GKT(XK_question),},
    {GKT(XK_at),},
    {GKT(XK_bracketleft),},
    {GKT(XK_backslash),},
    {GKT(XK_bracketright),},
    {GKT(XK_asciicircum),},
    {GKT(XK_grave),},
    {GKT(XK_a), "look_at_the_anchor",       "cutting_plane_a",
                "toggle_aux_property_thresholds_saturation",
                "shift_cutting_plane_to_anchor_a",
                "aux_property_coloring_10",  "delete_cutting_plane_a"},
    {GKT(XK_b), "toggle_bond_mode",         "cutting_plane_b",
                "toggle_bond_mode",         "shift_cutting_plane_to_anchor_b",
                "aux_property_coloring_11", "delete_cutting_plane_b"},
    {GKT(XK_c), "clone",                    "cutting_plane_c",
                "close",                    "shift_cutting_plane_to_anchor_c",
                "aux_property_coloring_12", "delete_cutting_plane_c"},
    {GKT(XK_d), "change_bgcolor",           "cutting_plane_d",
                "change_bgcolor",           "shift_cutting_plane_to_anchor_d",
                "aux_property_coloring_13", "delete_cutting_plane_d"},
    {GKT(XK_e), "capture_eps",              "cutting_plane_e",
                NULL,                       "shift_cutting_plane_to_anchor_e",
                "aux_property_coloring_14", "delete_cutting_plane_e"},
    {GKT(XK_f), "find_atom",                "cutting_plane_f",
                NULL,                       "shift_cutting_plane_to_anchor_f",
                "aux_property_coloring_15", "delete_cutting_plane_f"},
    {GKT(XK_g), "observer_goto",            "change_shear_strain_subtract_mean",
                "aux_property_coloring_16", NULL,
                "shear_strain_coloring",  NULL},
    {GKT(XK_h), NULL,                       "change_central_symm_neighbormax",
                "aux_property_coloring_17", NULL,
                "central_symmetry_coloring",  NULL},
    {GKT(XK_i), "change_wireframe_mode", "change_cutting_plane_wireframe_mode",
                "aux_property_coloring_18", NULL},
    {GKT(XK_j), "capture_jpg",              NULL,
                "aux_property_coloring_19", NULL},
    {GKT(XK_k), "toggle_coordination_coloring", NULL,
                "aux_property_coloring_30"},
    {GKT(XK_l), "normal_coloring",          NULL,
                "aux_property_coloring_31", NULL, "load_script"},
    {GKT(XK_m), RPT4(NULL), "change_aux_colormap", NULL},
    {GKT(XK_n), "new"},
    {GKT(XK_o), "original_normal_coloring"},
    {GKT(XK_p), "capture_png",              "flip_cutting_plane"},
    {GKT(XK_q), "close"},
    {GKT(XK_r), "toggle_rcut_patch_mode", NULL,"reset_aux_property_thresholds"},
    {GKT(XK_s), NULL, NULL, "resize", "resize", "save"},
    {GKT(XK_t), "xtal_origin_goto", NULL,
                "toggle_aux_property_thresholds_rigid"},
    {GKT(XK_u), "viewframe_upright"},
    {GKT(XK_v), "toggle_shell_viewer_mode"},
    {GKT(XK_w), "reset_anchor"},
    {GKT(XK_x), "toggle_xtal_mode"},
    {GKT(XK_y), "script_animate"},
    {GKT(XK_z), "xtal_origin_zero", "xtal_origin_half", },
    {GKT(XK_braceleft),},
    {GKT(XK_bar),},
    {GKT(XK_braceright),},
    {GKT(XK_asciitilde),},

    {GKT(XK_nobreakspace),},
    {GKT(XK_exclamdown),},
    {GKT(XK_cent),},
    {GKT(XK_sterling),},
    {GKT(XK_currency),},
    {GKT(XK_yen),},
    {GKT(XK_brokenbar),},
    {GKT(XK_section),},
    {GKT(XK_diaeresis),},
    {GKT(XK_copyright),},
    {GKT(XK_ordfeminine),},
    {GKT(XK_guillemotleft),},
    {GKT(XK_notsign),},
    {GKT(XK_hyphen),},
    {GKT(XK_registered),},
    {GKT(XK_macron),},
    {GKT(XK_degree),},
    {GKT(XK_plusminus),},
    {GKT(XK_twosuperior),},
    {GKT(XK_threesuperior),},
    {GKT(XK_acute),},
    {GKT(XK_mu),},
    {GKT(XK_paragraph),},
    {GKT(XK_periodcentered),},
    {GKT(XK_cedilla),},
    {GKT(XK_onesuperior),},
    {GKT(XK_masculine),},
    {GKT(XK_guillemotright),},
    {GKT(XK_onequarter),},
    {GKT(XK_onehalf),},
    {GKT(XK_threequarters),},
    {GKT(XK_questiondown),},
    {GKT(XK_Agrave),},
    {GKT(XK_Aacute),},
    {GKT(XK_Acircumflex),},
    {GKT(XK_Atilde),},
    {GKT(XK_Adiaeresis),},
    {GKT(XK_Aring),},
    {GKT(XK_AE),},
    {GKT(XK_Ccedilla),},
    {GKT(XK_Egrave),},
    {GKT(XK_Eacute),},
    {GKT(XK_Ecircumflex),},
    {GKT(XK_Ediaeresis),},
    {GKT(XK_Igrave),},
    {GKT(XK_Iacute),},
    {GKT(XK_Icircumflex),},
    {GKT(XK_Idiaeresis),},
    {GKT(XK_ETH),},
    {GKT(XK_Ntilde),},
    {GKT(XK_Ograve),},
    {GKT(XK_Oacute),},
    {GKT(XK_Ocircumflex),},
    {GKT(XK_Otilde),},
    {GKT(XK_Odiaeresis),},
    {GKT(XK_multiply),},
    {GKT(XK_Ooblique),},
    {GKT(XK_Ugrave),},
    {GKT(XK_Uacute),},
    {GKT(XK_Ucircumflex),},
    {GKT(XK_Udiaeresis),},
    {GKT(XK_Yacute),},
    {GKT(XK_THORN),},
    {GKT(XK_ssharp),},
    {GKT(XK_aacute),},
    {GKT(XK_acircumflex),},
    {GKT(XK_atilde),},
    {GKT(XK_adiaeresis),},
    {GKT(XK_aring),},
    {GKT(XK_ae),},
    {GKT(XK_ecircumflex),},
    {GKT(XK_ediaeresis),},
    {GKT(XK_igrave),},
    {GKT(XK_iacute),},
    {GKT(XK_icircumflex),},
    {GKT(XK_idiaeresis),},
    {GKT(XK_eth),},
    {GKT(XK_ntilde),},
    {GKT(XK_ograve),},
    {GKT(XK_oacute),},
    {GKT(XK_ocircumflex),},
    {GKT(XK_otilde),},
    {GKT(XK_odiaeresis),},
    {GKT(XK_division),},
    {GKT(XK_oslash),},
    {GKT(XK_ugrave),},
    {GKT(XK_uacute),},
    {GKT(XK_ucircumflex),},
    {GKT(XK_udiaeresis),},
    {GKT(XK_yacute),},
    {GKT(XK_thorn),},
    {GKT(XK_ydiaeresis),},

/* Japanese keyboard */
    {GKT(XK_Muhenkan),},
    {GKT(XK_Henkan),},
    {GKT(XK_Hiragana_Katakana),},
    {GKT(XK_Zenkaku_Hankaku),},

    {GKT(NoSymbol)}
};

static struct gkt *gktp_bykeysym(KeySym keysym)
{
    struct gkt *gktp;
    for (gktp = gui_key_table; gktp->keysym != NoSymbol; gktp++)
        if (gktp->keysym == keysym)
            return gktp;
    return NULL;
}

static struct aec **aecpp_bykeyname(char *keyname)
{
    char *s, buf2[CUI_LINEMAX];
    KeySym keysym;
    struct gkt *gktp;
    int shift = 0, ctrl = 0, meta = 0;

    if (!keyname || !*keyname) return NULL;
    strncpy(buf2, keyname, sizeof(buf2));

    if ((s = strchr(buf2, '+'))) {
        *s++ = 0;
        if (*s) {
            switch (toupper(*s++)) {
            case 'S': shift = 1;    break;
            case 'C': ctrl  = 1;    break;
            case 'M': meta  = 1;    break;
            }
        }
        if (*s) {
            switch (toupper(*s++)) {
            case 'S': shift = 1;    break;
            case 'C': ctrl  = 1;    break;
            case 'M': meta  = 1;    break;
            }
        }
    }

    if ((keysym = XStringToKeysym(buf2)) != NoSymbol &&
                                (gktp = gktp_bykeysym(keysym))) {
        if (ctrl)
            if (shift)  return &gktp->ctrl_shift;
            else        return &gktp->ctrl;
        else if (meta)
            if (shift)  return &gktp->meta_shift;
            else        return &gktp->meta;
        else
            if (shift)  return &gktp->shift;
            else        return &gktp->normal;
    }

    return NULL;
}

static int rgb_url_printed = 0;
static int maxfd_plus_1 = 0, listen_fd = -1;
static fd_set allset;

static struct {
    FILE *fp;
    int fd;
    double pause;
    int opener;
    int writable;
} cmdf[CUI_N_FILE] = {{NULL, -1, 0.0, 0, 0}};

static int cmdfdopen(int fd, char *mode, int opener)
{
    int i;

    for (i = 0; i < CUI_N_FILE; i++)
        if (!cmdf[i].fp) break;
    if (i == CUI_N_FILE)
        return -1;

    if (fd == -1 || !(cmdf[i].fp = fdopen(fd, mode)))
        return -1;

    cmdf[i].fd = fd;
    cmdf[i].pause = 0.0;
    cmdf[i].opener = opener;
    cmdf[i].writable = (strcmp(mode, "r+") == 0 || fd == STDIN_FILENO) ? 1 : 0;
    setbuf(cmdf[i].fp, NULL);
    
    FD_SET(cmdf[i].fd, &allset);
    if (cmdf[i].fd >= maxfd_plus_1)
        maxfd_plus_1 = cmdf[i].fd + 1;
    if (opener >= 0)
        cmdf[opener].pause = -1.0;

    return i;
}


static void cmdfclose(int k)
{
    int i, k0 = k, k1 = k + 1;

    if (k < 0) {
        k0 = 0;
        k1 = CUI_N_FILE;
    }

    for (i = k0; i < k1; i++) {
        if (!cmdf[i].fp) continue;

        if (cmdf[i].fd == STDIN_FILENO)
            quit = 1;

        if (cmdf[i].fd == frontend)
            rl_callback_handler_remove();
        else if (cmdf[i].writable) {
            cui_send(CUI_PROTOCOL_QUIT, cmdf[i].fp);
        }

        FD_CLR(cmdf[i].fd, &allset);
        if (maxfd_plus_1 == cmdf[i].fd + 1) {
            int j;
            maxfd_plus_1 = 0;
            if (listen_fd >= maxfd_plus_1) maxfd_plus_1 = listen_fd + 1;
            for (j = 0; j < CUI_N_FILE; j++) {
                if (cmdf[j].fp && cmdf[j].fd >= maxfd_plus_1)
                    maxfd_plus_1 = cmdf[j].fd + 1;
            }
        }
        cmdf[i].fd = -1;
        if (cmdf[i].opener >= 0)
            cmdf[cmdf[i].opener].pause = 0.0;

        fclose(cmdf[i].fp);
        cmdf[i].fp = NULL;
        cmdf[i].writable = 0;
    }
}

enum isv_type { ISV_int, ISV_double, ISV_v2, ISV_v3, ISV_v4};

#define NAVI(type,mem)\
            {ISV_##type, "n->"#mem, offsetof(Navigator,mem)}
            /*{ISV_##type, "n->"#mem, (size_t) &((Navigator *)0)->mem}*/
#define NAVI_POINTER(iw,isvtp)  (((char*)&n[(iw)])+(isvtp)->offset)

#define AX3D(type,mem)\
            {ISV_##type, "AX_3D->"#mem, offsetof(AX_3D_Viewport,mem)}
            /*{ISV_##type, "AX_3D->"#mem, (size_t)&((AX_3D_Viewport *)0)->mem}*/
#define AX3D_POINTER(iw,isvtp)  (((char*)&AX_3D[(iw)])+(isvtp)->offset)

#define ACC(num) \
            {ISV_v3,"CC["#num"]", (size_t) (\
                                num*sizeof(Atom_coordination_color)+\
                                offsetof(Atom_coordination_color,r))}
                                /*(size_t)&((Atom_coordination_color *)0)->r)}*/
#define ACC_POINTER(isvtp)      (((char*)CoordColor)+(isvtp)->offset)

static struct isvt {
    enum isv_type type;
    char *name;
    size_t offset;
} internal_state_variable_table[] = {

    NAVI(v3, bgcolor), /* background color r,g,b */
  /*NAVI(int, lx),*//* last pointer position */
  /*NAVI(int, ly),*/
  /*NAVI(long, last_button_press_time), for capturing double click */
    NAVI(int, anchor), /* -1: hook[] as anchor, >=0: certain atom */
/*
#if (ATOM_STACK_SIZE == 4)
    NAVI(int, atom_stack[0]),*//* keep track of selected atoms *//*
    NAVI(int, atom_stack[1]), NAVI(int, atom_stack[2]),NAVI(int, atom_stack[3]),
#else
#   error
#endif*/
    NAVI(v3, hook), /* real space coordinates */
    NAVI(int, suppress_printout), /* term printout */
    NAVI(double, delta), /* rate of change (gear-box) */
    NAVI(double, mgs_radius), /* magic sphere for pointer rotation */
    NAVI(double, atom_r_ratio), /* compared to ATOM_RADIUS */
    NAVI(int, color_mode), /* NORMAL, COORD, AUXILIARY or SCRATCH */
    NAVI(int, last_color_mode),
    NAVI(int, wireframe_mode), /* draw H-box */
    NAVI(int, bond_mode), /* draw bonds */
    NAVI(double, bond_radius), /* in Angstrom */
  /*NAVI(int, bond_atom_color_need_update), delayed update due to atom color */
    NAVI(int, shell_viewer_mode), /* xv and ghostview */
    NAVI(int, xtal_mode),
        /* to allow for crystal translation and shear-encoding */
    NAVI(v3, xtal_origin), /* shifted xtal origin */
  /*NAVI(int, bond_xtal_origin_need_update), delayed bond update due to shift */
    NAVI(int, last_surface_id), /* for if the pointer is not pointing at H box*/
  /*NAVI(int, mitosis),*//* for spawning children thread */
    NAVI(int, auxiliary_idx),
        /* which auxiliary property should we render */
    NAVI(int, auxiliary_cmap), /* which colormap should we use */
#if (CONFIG_MAX_AUXILIARY+MAX_GEO_MEASURES == 50) /* 18 */
#   define NAVI_AUX_THR(i)  NAVI(v2, auxiliary_threshold[i])
    NAVI_AUX_THR(0), NAVI_AUX_THR(1), NAVI_AUX_THR(2), NAVI_AUX_THR(3),
    NAVI_AUX_THR(4), NAVI_AUX_THR(5), NAVI_AUX_THR(6), NAVI_AUX_THR(7),
    NAVI_AUX_THR(8), NAVI_AUX_THR(9), NAVI_AUX_THR(10),NAVI_AUX_THR(11),
    NAVI_AUX_THR(12),NAVI_AUX_THR(13),NAVI_AUX_THR(14),NAVI_AUX_THR(15),
    NAVI_AUX_THR(16),NAVI_AUX_THR(17),NAVI_AUX_THR(18),NAVI_AUX_THR(19),
    NAVI_AUX_THR(20),NAVI_AUX_THR(21),NAVI_AUX_THR(22),NAVI_AUX_THR(23),
    NAVI_AUX_THR(24),NAVI_AUX_THR(25),NAVI_AUX_THR(26),NAVI_AUX_THR(27),
    NAVI_AUX_THR(28),NAVI_AUX_THR(29),NAVI_AUX_THR(30),NAVI_AUX_THR(31),
    NAVI_AUX_THR(32),NAVI_AUX_THR(33),NAVI_AUX_THR(34),NAVI_AUX_THR(35),
    NAVI_AUX_THR(36),NAVI_AUX_THR(37),NAVI_AUX_THR(38),NAVI_AUX_THR(39),
    NAVI_AUX_THR(40),NAVI_AUX_THR(41),NAVI_AUX_THR(42),NAVI_AUX_THR(43),
    NAVI_AUX_THR(44),NAVI_AUX_THR(45),NAVI_AUX_THR(46),NAVI_AUX_THR(47),
#else
#   error
#endif
    NAVI(int, auxiliary_thresholds_saturation),
        /* 0:invisible if outside thresholds */
    NAVI(int, auxiliary_thresholds_rigid), /* 0: floating thresholds */
    NAVI(int, parallel_projection),
        /* parallel / perspective projections */
    NAVI(int, glob_advance), /* how fast the file list advances */
#if (AX_3D_MAX_FILTER_PLANE == 16)
#   define NAVI_FP(i) \
            NAVI(int, fp[i].wireframe_mode), NAVI(v3, fp[i].s0),\
            NAVI(v4, fp[i].dx_input), NAVI(v4, fp[i].dx_cache)
    NAVI_FP(0), NAVI_FP(1), NAVI_FP(2), NAVI_FP(3), NAVI_FP(4), NAVI_FP(5),
    NAVI_FP(6), NAVI_FP(7), NAVI_FP(8), NAVI_FP(9), NAVI_FP(10),NAVI_FP(11),
    NAVI_FP(12),NAVI_FP(13),NAVI_FP(14),NAVI_FP(15),
#else
#   error
#endif
    NAVI(int, just_activated_fp),  /* index of the fp under focus */

    AX3D(v3, x), /* coordinates of the viewpoint */
    AX3D(double, k), /* conversion factor from radian to window pixels*/
    AX3D(v3, V[0]), /* V[i][0-2] is the normalized ith axis of viewport*/
    AX3D(v3, V[1]),
    AX3D(v3, V[2]),
    AX3D(double, zcut), /* (0,zcut] of the view frustum */
    AX3D(double, wx),/* viewpoint coordinates in window frame (pixels)*/
    AX3D(double, wy),
#if (AX_3D_MAX_FILTER_PLANE == 16)
#   define AX3D_FP(i)   AX3D(double, fp[i].d0),     AX3D(v3, fp[i].dx)
    AX3D_FP(0), AX3D_FP(1), AX3D_FP(2), AX3D_FP(3), AX3D_FP(4), AX3D_FP(5),
    AX3D_FP(6), AX3D_FP(7), AX3D_FP(8), AX3D_FP(9), AX3D_FP(10),AX3D_FP(11),
    AX3D_FP(12),AX3D_FP(13),AX3D_FP(14),AX3D_FP(15), /* filter planes */
#else
#   error
#endif
#if (ATOM_COORDINATION_MAX+1 == 25)
    ACC(0),  ACC(1),  ACC(2),  ACC(3),  ACC(4),
    ACC(5),  ACC(6),  ACC(7),  ACC(8),  ACC(9),
    ACC(10), ACC(11), ACC(12), ACC(13), ACC(14),
    ACC(15), ACC(16), ACC(17), ACC(18), ACC(19),
    ACC(20), ACC(21), ACC(22), ACC(23), ACC(24),
#else
#   error
#endif

    {0, NULL, 0}
};

static struct aec *aecp_byname(char *name);
static bool proc_set(int iw, char *instr, char **outstr)
{
    char name[CUI_LINEMAX] = "";

    *outstr = NULL;
    buf[0] = 0;
    if (sscanf(instr, " %[^ ] %[^\n]", name, buf) < 2)
        goto error;

    if (strncmp(name, "key->", strlen("key->")) == 0) {
        struct aec **aecpp, *aecp;
        if (!(aecpp = aecpp_bykeyname(name + strlen("key->")))) {
            goto error;
        }
        else if (!(aecp = aecp_byname(buf))) {
            goto error;
        }
        *aecpp = aecp;
        return FALSE;
    }
    else {
        struct isvt *isvtp;
        for (isvtp = internal_state_variable_table; isvtp->name; isvtp++) {
            if (strcmp(isvtp->name, name) == 0) {
                void *pointer;
                switch (*isvtp->name) {
                case 'n': pointer = NAVI_POINTER(iw, isvtp); break;
                case 'A': pointer = AX3D_POINTER(iw, isvtp); break;
                case 'C': pointer =  ACC_POINTER(    isvtp); break;
                default: goto error;
                }
                switch (isvtp->type) {
                case ISV_int:
                    if (sscanf(buf, "%d",  (int*)   pointer) == 1)
                        return FALSE;
                    else
                        goto error;
                case ISV_double:
                    if (sscanf(buf, "%lf", (double*)pointer) == 1)
                        return FALSE;
                    else
                        goto error;
                case ISV_v2:
                    if (sscanf(buf, "%lf %lf",
                                (double*)pointer, (double*)pointer+1) > 0)
                        return FALSE;
                    else
                        goto error;
                case ISV_v3:
                    if (sscanf(buf, "%lf %lf %lf", V3e((double*)pointer)) > 0)
                        return FALSE;
                    else
                        goto error;
                case ISV_v4:
                    if (sscanf(buf, "%lf %lf %lf %lf",
                                (double*)pointer,   (double*)pointer+1,
                                (double*)pointer+2, (double*)pointer+3) > 0)
                        return FALSE;
                    else
                        goto error;
                }
            }
        }
    }
error:
    *outstr = cui_show_syntax;
    return FALSE;
}


#define SETKEYPRINT(fp, keystr, modstr, entry)\
        if (gktp->entry && gktp->entry->name && *gktp->entry->name)\
            fprintf(fp, "set key->%s%s %s\n", keystr, modstr, gktp->entry->name)

static bool proc_save(int iw, char *instr, char **outstr)
{
    char fname[FILENAME_MAX] = CUI_SCRIPT_DEFAULT;
    int isv = 1, key = 0;

    if (!instr)
        goto error;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(
            gui_readline_gets(iw, "\nSave script_file", CUI_SCRIPT_DEFAULT));

    *outstr = NULL;
    buf[0] = 0;

    if (!IS_MANAGER) return FALSE;

    if (sscanf(instr, " %[^ ] %[^ ]", fname, buf) == 2) {
        if (strcmp(buf, "both") == 0)
            key = 1;
        else if (strcmp(buf, "key") == 0) {
            key = 1;
            isv = 0;
        }
    }

    if (fname[0]) {
        FILE *fp = fopen(fname, "w");
        if (!fp)
            goto error;
        if (isv) {
            struct isvt *isvtp;
            void *pointer;

            for (isvtp = internal_state_variable_table; isvtp->name; isvtp++) {
                switch (*isvtp->name) {
                case 'n': pointer = NAVI_POINTER(iw, isvtp); break;
                case 'A': pointer = AX3D_POINTER(iw, isvtp); break;
                case 'C': pointer =  ACC_POINTER(    isvtp); break;
                default: goto error;
                }
                switch (isvtp->type) {
                case ISV_int:
                    fprintf(fp, "set %s %d\n",
                        isvtp->name, *(int*)   pointer);
                    break;
                case ISV_double:
                    fprintf(fp, "set %s %.15g\n",
                        isvtp->name, *(double*)pointer);
                    break;
                case ISV_v2:
                    fprintf(fp, "set %s %.15g %.15g\n",
                        isvtp->name, *(double*)pointer, *((double*)pointer+1));
                    break;
                case ISV_v3:
                    fprintf(fp, "set %s %.15g %.15g %.15g\n",
                        isvtp->name, V3E((double*)pointer));
                    break;
                case ISV_v4:
                    fprintf(fp, "set %s %.15g %.15g %.15g %.15g\n",
                        isvtp->name,*(double*)pointer  , *((double*)pointer+1),
                                    *(double*)pointer+2, *((double*)pointer+3));
                    break;
                }
            }
            fprintf(fp, "resize %d %d\n", AX_size[iw].width,AX_size[iw].height);
            fprintf(fp, "redraw\n");
        }

        if (key) {
            struct gkt *gktp;
            char *keystr;
            for (gktp = gui_key_table; gktp->keysym != NoSymbol; gktp++) {
                keystr = XKeysymToString(gktp->keysym);
                if (!keystr) continue;
                SETKEYPRINT(fp, keystr, ""   , normal    );
                SETKEYPRINT(fp, keystr, "+S" , shift     );
                SETKEYPRINT(fp, keystr, "+C" , ctrl      );
                SETKEYPRINT(fp, keystr, "+CS", ctrl_shift);
                SETKEYPRINT(fp, keystr, "+M" , meta      );
                SETKEYPRINT(fp, keystr, "+MS", meta_shift);
            }
        }

        fclose(fp);
        return FALSE;
    }
error:
    *outstr = cui_show_syntax;
    return FALSE;
}


static bool proc_nop(int iw, char *instr, char **outstr)
{
    *outstr = NULL;
    return FALSE;
}

static bool proc_internal_error(int iw, char *instr, char **outstr)
{
    *outstr = "internal error";
    return FALSE;
}

static bool proc_load_script(int iw, char *instr, char **outstr)
{
    int i = -1;
    *outstr = cui_show_syntax;
    if (!instr)
        return FALSE;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(
            gui_readline_gets(iw, "\nLoad script_file", CUI_SCRIPT_DEFAULT));

    *outstr = NULL;

    if (!IS_MANAGER) return FALSE;
    strncpy(fname, CUI_SCRIPT_DEFAULT, sizeof(fname));
    sscanf(instr, " %[^ ]%d", fname, &i);

    if (cmdfdopen(open(fname, O_RDONLY), "r", i) == -1)
        *outstr = "open failed";

    return FALSE;
}

static bool proc_next(int iw, char *instr, char **outstr)
{
    int jw;
#ifdef USE_P3D
    if (p3dp_enabled)
        return FALSE;
#endif
    *outstr = NULL;
    for (jw = iw + 1; jw < AX_MAXWIN; jw++) {
        if (nonvacant(jw)) {
            cui_iw = jw;
            AXSetName(iw);
            AXSetName(jw);
            return FALSE;
        }
    }
    for (jw = 0; jw < iw; jw++) {
        if (nonvacant(jw)) {
            cui_iw = jw;
            AXSetName(iw);
            AXSetName(jw);
            return FALSE;
        }
    }
    cui_iw = -1;
    return FALSE;
}

static bool proc_close(int iw, char *instr, char **outstr)
{
#ifdef ATOMEYE_LIB
  (*atomeyelib_on_close)();
#endif
    if (iw == cui_iw)
        proc_next(iw, NULL, outstr);
    if (cui_iw < 0)
        cmdfclose(-1);
    AX_closewindow(iw);
#ifdef USE_P3D
    if (p3dp_enabled) longjmp(quit_env, 1);
    else
#endif
    pthread_exit ((void *)NULL);
    exit(-1);
}

static bool proc_quit(int iw, char *instr, char **outstr)
{
    quit = 1;
    proc_close(iw, NULL, outstr);
    exit (-1);
}


static bool proc_new(int iw, char *instr, char **outstr)
{
    pthread_t tid;

#ifdef ATOMEYE_LIB
    *outstr = "new window creation disabled";
    return FALSE;
#endif
#ifdef USE_P3D
    if (p3dp_enabled) {
        *outstr = "not supported";
        return FALSE;
    }
#endif
    if (instr && strcmp(instr, "clone") == 0)
        n[iw].mitosis = iw; /* clone */
    else
        n[iw].mitosis = -1;
    
    *outstr = NULL;
    pthread_create(&tid, NULL,  (void *(*)(void *)) thread_start,
                                (void *)(&(n[iw].mitosis)));
    return FALSE;
}


static bool proc_key(int iw, char *instr, char **outstr)
{
    struct aec **aecpp;

    *outstr = NULL;
    if (instr && (aecpp = aecpp_bykeyname(instr)) && *aecpp) {
        if (strcmp((*aecpp)->name, "load_script") == 0)
            goto error;
        else if ((*aecpp)->instr && strcmp((*aecpp)->instr,CUI_ARG_REQUIRED)!=0)
            return (*(*aecpp)->proc)(iw, (*aecpp)->instr, outstr);
        else
            return (*(*aecpp)->proc)(iw, NULL, outstr);
    }
error:
    *outstr = cui_show_syntax;
    return FALSE;
}


static bool proc_free_scratch(int iw, char *instr, char **outstr)
{
    extern AX_Float *scratch_r, *scratch_g, *scratch_b;
    extern int scratch_np;
    int i = ((scratch_r!=NULL)+(scratch_g!=NULL)+(scratch_b!=NULL));

    *outstr = NULL;
    Free(scratch_r);
    Free(scratch_g);
    Free(scratch_b);
    sprintf(buf, "Previous scratch colors freed (%s bytes).",
                    strmem((long)i*scratch_np*sizeof(AX_Float)));
    *outstr = buf;
    scratch_np = 0;
    return FALSE;
}


static bool proc_scratch_coloring(int iw, char *instr, char **outstr)
{
    extern AX_Float *scratch_r, *scratch_g, *scratch_b;
    extern int scratch_np;
    int i;

#ifdef USE_P3D
    if (p3dp_enabled) {
        *outstr = "not supported";
        return FALSE;
    }
#endif
    *outstr = NULL;
    if (np != scratch_np) {
        sprintf(buf,"scratch_color: current np=%d not equals to that of\n"
                    "old scratch=%d, Please rescratch.", np, scratch_np);
        *outstr = buf;
        return FALSE;
    }
    strcpy(AX_title[iw], str4(fbasename, " (", "scratch", ")"));
    AXSetName(iw);
    XStoreName(AX_display[iw],xterm_win,AX_title[iw]);
    XSetIconName(AX_display[iw],xterm_win,AX_title[iw]);
    for (i=0; i<np; i++)
        AX_3D_AssignRGB(B->BALL[i],scratch_r[i],scratch_g[i],scratch_b[i]);
    n[iw].color_mode = COLOR_MODE_SCRATCH;
    if (n[iw].bond_mode) {
        bond_xtal_origin_update(iw);
        bond_atom_color_update(iw);
    }
    else {
        n[iw].bond_xtal_origin_need_update = TRUE;
        n[iw].bond_atom_color_need_update = TRUE;
    }
    return TRUE;
}


static bool proc_reset_scratch(int iw, char *instr, char **outstr)
{
    extern AX_Float *scratch_r, *scratch_g, *scratch_b;
    extern int scratch_np, scratch_n[3];
    char *c;
    int i, j;
    double ds[3], new_s[3];

#ifdef USE_P3D
    if (p3dp_enabled) {
        *outstr = "not supported";
        return FALSE;
    }
#endif
    *outstr = NULL;
    if (instr) {
        c = cui_stripspace(strncpy(buf, instr,sizeof(buf)));
        if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
            instr = gui_readline_gets(iw,
                            "\nScratch into n1xn2xn3 color blocks", "2 2 2");
        sscanf(instr, " %d %d %d\n", scratch_n, scratch_n+1, scratch_n+2);
    }
    else
        scratch_n[0] = scratch_n[1] = scratch_n[2] = 2;

    REALLOC(rescratch, scratch_r, np, AX_Float);
    REALLOC(rescratch, scratch_g, np, AX_Float);
    REALLOC(rescratch, scratch_b, np, AX_Float);
    scratch_np = np;
    if (scratch_n[0] < 1) scratch_n[0] = 2;
    if (scratch_n[1] < 1) scratch_n[1] = 2;
    if (scratch_n[2] < 1) scratch_n[2] = 2;
    V3mM3 (n[iw].xtal_origin, HI, ds);
    V3TRIM (ds, ds);
    for (i=0; i<np; i++) {
        V3SUB ( &(s[DIMENSION*i]), ds, new_s );
        V3TriM ( new_s );
        j = INT(new_s[0] * scratch_n[0]) +
            INT(new_s[1] * scratch_n[1]) +
            INT(new_s[2] * scratch_n[2]);
        if (ISEVEN(j)) { /* firebrick4 */
            scratch_r[i] = 159 / 255.;
            scratch_g[i] =  26 / 255.;
            scratch_b[i] =  26 / 255.;
        }
        else { /* CadetBlue3 */
            scratch_r[i] = 122 / 255.;
            scratch_g[i] = 197 / 255.;
            scratch_b[i] = 205 / 255.;
        }
    }
    return proc_scratch_coloring(iw, instr, outstr);
}


static bool proc_toggle_coordination_coloring(int iw,char*instr,char**outstr)
{
    if (n[iw].color_mode != COLOR_MODE_COORD)
    {
        if (n[iw].xtal_mode)
        {
            n[iw].last_color_mode = n[iw].color_mode;
            return(assign_coordination_color(iw));
        }
        else return (FALSE);
    }
    else if (n[iw].last_color_mode == COLOR_MODE_NORMAL)
        return(assign_normal_color(iw));
    else if (n[iw].last_color_mode == COLOR_MODE_AUXILIARY)
        return(color_encode_auxiliary(iw));
    else if (n[iw].color_mode == COLOR_MODE_SCRATCH)
        /*return(scratch_color(iw));*/
        return(proc_scratch_coloring(iw, instr, outstr));
    else return (FALSE);
}


static bool proc_redraw(int iw, char *instr, char **outstr)
{
    switch (n[iw].color_mode) {
    case COLOR_MODE_COORD:
        assign_coordination_color(iw);
        break;
    case COLOR_MODE_NORMAL:
        assign_normal_color(iw);
        break;
    case COLOR_MODE_AUXILIARY:
        color_encode_auxiliary(iw);
        break;
    case COLOR_MODE_SCRATCH:
        /*scratch_color(iw);*/
        proc_scratch_coloring(iw, instr, outstr);
        break;
    }
    bond_xtal_origin_update(iw);
    bond_atom_color_update(iw);
    cui_diligence = TRUE;
    return TRUE;
}


static bool proc_translate(int iw, char *instr, char **outstr)
{
    int axis = 0;
    double d = 0.0;

    if (!instr)
        goto error;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(gui_readline_gets(iw, "\nTranslate", "0 0.0"));

    *outstr = NULL;
    buf[0] = 0;
    switch (sscanf(instr, " %d%lf %[^ ]", &axis, &d, buf))
    {
    case 3:
        if (strcmp(buf, "delta") == 0)
            d *= n[iw].delta;
        else
            goto error;
    case 2:
        if (axis < 0 || axis > 2)
            goto error;
        else
            return translate(iw, axis, d);
    default:
        /* error */
        break;
    }
error:
    *outstr = cui_show_syntax;
    return FALSE;
}


static bool proc_shift_xtal(int iw, char *instr, char **outstr)
{
    int axis = 0;
    double d = 0.0;

    if (!n[iw].xtal_mode)
        return FALSE;

    if (!instr)
        goto error;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(gui_readline_gets(iw, "\nXtal shift", "0 0.0"));

    *outstr = NULL;
    buf[0] = 0;
    switch (sscanf(instr, " %d%lf %[^ ]", &axis, &d, buf))
    {
    case 3:
        if (strcmp(buf, "delta") == 0)
            d *= n[iw].delta;
        else
            goto error;
    case 2:
        if (axis < 0 || axis > 2)
            goto error;
        else
            return xtal_shift(iw, axis, d);
    default:
        /* error */
        break;
    }
error:
    *outstr = cui_show_syntax;
    return FALSE;
}


static bool proc_rotate(int iw, char *instr, char **outstr)
{
    int axis = 0;
    double theta = 0.0;

    if (!instr)
        goto error;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(gui_readline_gets(iw, "\nRotate axis", "0 0.0"));

    *outstr = NULL;
    buf[0] = 0;
    switch (sscanf(instr, " %d%lf %[^ ]", &axis, &theta, buf))
    {
    case 3:
        if (strcmp(buf, "delta") == 0)
            theta *= n[iw].delta;
        else
            goto error;
    case 2:
        if (axis < 0 || axis > 2)
            goto error;
        else
            return axis_rotate(iw, axis, theta*PI);
    default:
        /* error */
        break;
    }
error:
    *outstr = cui_show_syntax;
    return FALSE;
}


static bool proc_advance(int iw, char *instr, char **outstr)
{
    double delta = 0.0;
    double tmp[3];

    *outstr = cui_show_syntax;
    if (!instr)
        return FALSE;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        if (AX_display[iw])
            printf("\nChange the distance between viewpoint and anchor:\n"
                    "if delta >=0, d -> d/(1+delta), else d -> d * (1-delta).");
        instr = cui_stripspace(gui_readline_gets(iw, "\ndelta", "0.0"));
    }

    buf[0] = 0;
    if (sscanf(instr, " %lf %s", &delta, buf) < 1)
        return FALSE;

    *outstr = NULL;
    if (strcmp(buf, "delta") == 0)
        delta *= n[iw].delta;

    if (delta >= 0) delta = 1/(1+delta);
    else delta = (1-delta);
    if (n[iw].anchor >= 0) {
#ifdef USE_P3D
        if (p3dp_enabled) return FALSE;
#endif
        V3SUB (AX_3D[iw].x, B->BALL[n[iw].anchor].x, tmp);
        V3ADDMUL (B->BALL[n[iw].anchor].x, delta, tmp, AX_3D[iw].x);
    }
    else {
        V3SUB (AX_3D[iw].x, n[iw].hook, tmp);
        V3ADDMUL (n[iw].hook, delta, tmp, AX_3D[iw].x);
    }
    return TRUE;
}


static bool proc_shift_cutting_plane(int iw, char *instr, char **outstr)
{
    double d;

    if (!instr)
        goto error;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(
                        gui_readline_gets(iw, "\nShift filter plane", "0.0"));

    *outstr = NULL;
    buf[0] = 0;
    if (sscanf(instr, " %lf %[^ ]", &d, buf) == 2)
    {
        if (strcmp(buf, "delta") == 0)
            d *= n[iw].delta;
        else
            goto error;
    }
    return shift_filter_plane(iw, d);
error:
    *outstr = "parameter error";
    return FALSE;
}


static bool proc_change_bgcolor(int iw, char *instr, char **outstr)
{
    double c[3];

    *outstr = cui_show_syntax;
    if (!instr)
        return FALSE;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        if (AX_display[iw] && !rgb_url_printed) {
            printf("\nColor choices at:\n%s\n", RGB_URL);
            rgb_url_printed = 1;
        }
        sprintf(buf, "%.3f %.3f %.3f",
                n[iw].bgcolor[0], n[iw].bgcolor[1], n[iw].bgcolor[2]);
        instr = cui_stripspace(
                    gui_readline_gets(iw,"\nChange background color",buf));
    }

    V3EQV(n[iw].bgcolor, c);
    if (sscanf(instr, "%lf %lf %lf", c, c+1, c+2) < 1)
        return FALSE;
    *outstr = NULL;
    if (c[0]>1) c[0]/=255;
    if (c[1]>1) c[1]/=255;
    if (c[2]>1) c[2]/=255;
    if ( V3NE(c, n[iw].bgcolor) ) {
        V3EQV(c, n[iw].bgcolor);
        return TRUE;
    }
    else return FALSE;
}

static bool proc_change_atom_r_ratio(int iw, char *instr, char **outstr)
{
    double c = 0.0;
    if (!instr)
        goto error;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(
                    gui_readline_gets(iw, "\nChange atom r ratio", "0.0"));

    *outstr = NULL;
    buf[0] = 0;
    if (sscanf(instr, " %lf %[^ ]", &c, buf) == 2)
    {
        if (strcmp(buf, "delta") == 0)
            c *= n[iw].delta;
        else
            goto error;
    }

    if (c >= 0)
        change_atom_r_ratio(n[iw].atom_r_ratio*=(1+c));
    else
        change_atom_r_ratio(n[iw].atom_r_ratio/=(1-c));
    return TRUE;

error:
    *outstr = "parameter error";
    return FALSE;
}


static bool proc_change_bond_radius(int iw, char *instr, char **outstr)
{
    double c = 0.0;
    if (!n[iw].bond_mode) {
        *outstr = NULL;
        return FALSE;
    }

    if (!instr)
        goto error;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(
                    gui_readline_gets(iw, "\nChange bond radius", "0.0"));

    *outstr = NULL;
    buf[0] = 0;
    if (sscanf(instr, " %lf %[^ ]", &c, buf) == 2)
    {
        if (strcmp(buf, "delta") == 0)
            c *= n[iw].delta;
        else
            goto error;
    }

    if (c >= 0)
        change_bond_radius(n[iw].bond_radius*=(1+c));
    else
        change_bond_radius(n[iw].bond_radius/=(1-c));
    return TRUE;

error:
    *outstr = "parameter error";
    return FALSE;
}


static bool proc_change_view_angle_amplification(int iw,char*instr,char**outstr)
{
    double c = 0.0;
    if (!instr || !*instr)
        goto error;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(gui_readline_gets(iw,
                            "\nChange view angle amplification", "0.0"));

    *outstr = NULL;
    buf[0] = 0;
    if (sscanf(instr, " %lf %[^ ]", &c, buf) == 2)
    {
        if (strcmp(buf, "delta") == 0)
            c *= n[iw].delta;
        else
            goto error;
    }

    if (c >= 0)
        AX_3D[iw].k *= (1+c);
    else
        AX_3D[iw].k /= (1-c);
    return TRUE;

error:
    *outstr = "parameter error";
    return FALSE;
}


static bool proc_toggle_parallel_projection(int iw, char *instr, char **outstr)
{
    if (n[iw].parallel_projection)
        return(parallel_to_perspective(iw));
    else
        return(perspective_to_parallel(iw));
}


static bool proc_toggle_bond_mode(int iw, char *instr, char **outstr)
{
    n[iw].bond_mode = !n[iw].bond_mode;
    if (n[iw].bond_mode)
    {
        if (n[iw].bond_xtal_origin_need_update)
            bond_xtal_origin_update (iw);
        if (n[iw].bond_atom_color_need_update)
            bond_atom_color_update (iw);
    }
    return (TRUE);
}



static bool proc_normal_coloring(int iw, char *instr, char **outstr)
{
    *outstr = NULL;
    if (instr && strcmp(instr, "original") == 0)
        return assign_original_normal_color(iw);
    else
        return assign_normal_color(iw);
}


bool proc_change_atom_color(int iw, char *instr, char **outstr)
{
    int items, i;
    double c[4];
    char question[MAX_FILENAME_SIZE];

    if (!instr || sscanf(instr, " %d %[^\n]", &i, buf) != 2) {
        *outstr = cui_show_syntax;
        return FALSE;
    }
    instr = buf;

    *outstr = NULL;
#ifdef USE_P3D
    if (!p3dp_enabled || p3dp_rank_grab == p3d_rank(p3dp_cell)) {
#endif
    c[0] = B->BALL[i].r;
    c[1] = B->BALL[i].g;
    c[2] = B->BALL[i].b;
    c[3] = B->BALL[i].radius / n[iw].atom_r_ratio;
#ifdef USE_P3D
    }
    if (p3dp_enabled)
        p3d_bcast(p3dp_cell, c, 4, MPI_DOUBLE, p3dp_rank_grab);
#endif
    if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        if (!rgb_url_printed) {
            if (AX_display[iw])
                printf("\nColor choices at:\n%s\n", RGB_URL);
            rgb_url_printed = 1;
        }
#ifdef USE_P3D
        if (p3dp_enabled) {
            if (p3dp_rank_grab == p3d_rank(p3dp_cell))
                sprintf(question, "\nChange color [radius] of atom (\"%s\")",
                    COMPACT_SYMBOL(Symbol(symbol,i)));
            p3d_bcast(p3dp_cell, question,
                          sizeof(question), MPI_CHAR, p3dp_rank_grab);
        }
        else
#endif
        sprintf(question, "\nChange color [radius] of atom-%d (\"%s\")",
                i, COMPACT_SYMBOL(Symbol(symbol,i)));
        sprintf(buf, "%.3f %.3f %.3f [%.3f]", c[0], c[1], c[2], c[3]);
        instr = cui_stripspace(gui_readline_gets(iw, question, buf));
    }
#ifdef USE_P3D
    if (!p3dp_enabled || p3dp_rank_grab == p3d_rank(p3dp_cell)) {
#endif
    items = sscanf(instr, "%lf %lf %lf %lf", c, c+1, c+2, c+3);
    if (c[0] < 0 || c[1] < 0 || c[2] < 0) {
        c[0] = c[1] = c[2] = -1;
        c[3] = B->BALL[i].radius / n[iw].atom_r_ratio;
    }
    else if (items == 1) {
        c[3] = c[0];
        c[0] = B->BALL[i].r;
    }
    else {
        if (c[0]>1) c[0]/=255;
        if (c[1]>1) c[1]/=255;
        if (c[2]>1) c[2]/=255;
    }
    AX_3D_AssignRGB(B->BALL[i], c[0],c[1],c[2]);
    B->BALL[i].radius = c[3] * n[iw].atom_r_ratio;
#ifdef USE_P3D
    }
#endif
    bond_atom_color_update(iw);
    return TRUE;
}


bool proc_change_bond_color(int iw, char *instr, char **outstr)
{
    int items, i;
    double c[4];
    char question[MAX_FILENAME_SIZE];

    if (!instr || sscanf(instr, " %d %[^\n]", &i, buf) != 2) {
        *outstr = cui_show_syntax;
        return FALSE;
    }
    instr = buf;

    *outstr = NULL;
#ifdef USE_P3D
    if (!p3dp_enabled || p3dp_rank_grab == p3d_rank(p3dp_cell)) {
#endif
    c[0] = C->CYLINDER[i].r;
    c[1] = C->CYLINDER[i].g;
    c[2] = C->CYLINDER[i].b;
    c[3] = C->CYLINDER[i].radius;
#ifdef USE_P3D
    }
    if (p3dp_enabled)
        p3d_bcast(p3dp_cell, c, 4, MPI_DOUBLE, p3dp_rank_grab);
#endif
    if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        if (!rgb_url_printed) {
            if (AX_display[iw])
                printf("\nColor choices at:\n%s\n", RGB_URL);
            rgb_url_printed = 1;
        }
#ifdef USE_P3D
        if (p3dp_enabled)
            strncpy(question,
                        "\nChange color [radius] of bond", sizeof(question));
        else
#endif
        sprintf(question, "\nChange color [radius] of bond-%d", i);
        sprintf(buf, "%.3f %.3f %.3f [%.3f]", c[0], c[1], c[2], c[3]);
        instr = cui_stripspace(gui_readline_gets(iw, question, buf));
    }
#ifdef USE_P3D
    if (!p3dp_enabled || p3dp_rank_grab == p3d_rank(p3dp_cell)) {
#endif
    items = sscanf(instr, "%lf %lf %lf %lf", c, c+1, c+2, c+3);
    if (c[0] < 0 || c[1] < 0 || c[2] < 0) {
        c[0] = c[3] = -1;
    }
    else if (items == 1) {
        c[3] = c[0];
        c[0] = C->CYLINDER[i].r;
    }
    else {
        if (c[0]>1) c[0]/=255;
        if (c[1]>1) c[1]/=255;
        if (c[2]>1) c[2]/=255;
    }
    C->CYLINDER[i].r = c[0];
    C->CYLINDER[i].g = c[1];
    C->CYLINDER[i].b = c[2];
    C->CYLINDER[i].radius = c[3];
#ifdef USE_P3D
    }
#endif
    return TRUE;
}


bool proc_change_normal_color(int iw, char *instr, char **outstr)
{
    int t, items, z;
    double c[4];
    char question[MAX_FILENAME_SIZE] = "";
    /* selection means there must be something there already */
    if (n[iw].color_mode != COLOR_MODE_NORMAL) {
        *outstr = "You need to be in NORMAL color mode first.";
        return FALSE;
    }
    else if (!instr || (sscanf(instr, " %d %[^\n]", &z, buf) != 2 &&
                        sscanf(instr, " %s %[^\n]", question, buf) != 2)) {
        *outstr = cui_show_syntax;
        return FALSE;
    }
    if (question[0]) z = search_atom_by_symbol(question);
    instr = buf;

    *outstr = NULL;
    for (t = 0; t < ct->t; t++)
        if (ct->Z[t] == z) break;
    if (t == ct->t) {
        t = -1;
        c[0] = c[1] = c[2] = c[3] = 0;
    }
    else {
        c[0] = ATOM_Color_R(ct->Z[t]);
        c[1] = ATOM_Color_G(ct->Z[t]);
        c[2] = ATOM_Color_B(ct->Z[t]);
        c[3] = ATOM_Radius(ct->Z[t]);
    }
    if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        if (!rgb_url_printed) {
            if (AX_display[iw])
                printf("\nColor choices at:\n%s\n", RGB_URL);
            rgb_url_printed = 1;
        }
        sprintf(question, "\nChange color [radius] of type-%d (\"%s\") atoms",
                t, atom_symbol(z));
                /*t, COMPACT_SYMBOL(Symbol(symbol,ct->first[t])));*/
        sprintf(buf, "%.3f %.3f %.3f [%.3f]", c[0], c[1], c[2], c[3]);
        instr = cui_stripspace(gui_readline_gets(iw, question, buf));
    }
    items = sscanf(instr, "%lf %lf %lf %lf", c, c+1, c+2, c+3);

    if (t >= 0) {
        if (c[0] < 0 || c[1] < 0 || c[2] < 0) {
            c[0] = c[1] = c[2] = -1;
            c[3] = ATOM_Radius(ct->Z[t]);
        }
        else if (items == 1) {
            c[3] = c[0];
            c[0] = ATOM_Color_R(ct->Z[t]);
        }
        else {
            if (c[0]>1) c[0]/=255;
            if (c[1]>1) c[1]/=255;
            if (c[2]>1) c[2]/=255;
        }
        ATOM_Color_R(ct->Z[t]) = c[0];
        ATOM_Color_G(ct->Z[t]) = c[1];
        ATOM_Color_B(ct->Z[t]) = c[2];
        ATOM_Radius(ct->Z[t])  = c[3];
    }
    return assign_normal_color(iw);
}


bool proc_change_coordination_color(int iw, char *instr, char **outstr)
{
    int coord;
    char question[MAX_FILENAME_SIZE];
    double c[3];

    if (!n[iw].xtal_mode) {
        *outstr = "Coordination color change works only in Xtal mode.";
        goto error;
    }
    else if (n[iw].color_mode != COLOR_MODE_COORD) {
        *outstr = "You need to be in COORDINATION color mode first.";
        goto error;
    }
    else if (!instr || sscanf(instr, " %d %[^\n]", &coord, buf) != 2) {
        *outstr = cui_show_syntax;
        return FALSE;
    }
    instr = buf;

    c[0] = CoordColor[coord].r;
    c[1] = CoordColor[coord].g;
    c[2] = CoordColor[coord].b;

    if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        if (!rgb_url_printed) {
            if (AX_display[iw])
                printf("\nColor choices at:\n%s\n", RGB_URL);
            rgb_url_printed = 1;
        }
        sprintf(question, "\nChange color of %d-coordinated atoms", coord);
        sprintf(buf, "%.3f %.3f %.3f", c[0], c[1], c[2]);
        instr = cui_stripspace(gui_readline_gets(iw, question, buf));
    }
    sscanf(instr, " %lf %lf %lf", c, c+1, c+2);

    if (c[0] < 0 || c[1] < 0 || c[2] < 0)
        c[0] = c[1] = c[2] = -1;
    else {
        if (c[0] > 1) c[0]/=255;
        if (c[1] > 1) c[1]/=255;
        if (c[2] > 1) c[2]/=255;
    }

    CoordColor[coord].r = c[0];
    CoordColor[coord].g = c[1];
    CoordColor[coord].b = c[2];
    return assign_coordination_color(iw);

error:
    if (AX_display[iw] && !cui_xterm_win)
        printf("%s\n", *outstr);
    return FALSE;
}

static int get_number(int iw, char *instr)
{
    int number = -1;
    char c;

    if (!instr || !*instr)
        return -1;

    if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(gui_readline_gets(iw, "\nNumber", "0"));

    if (sscanf(instr, "%d", &number) != 1 &&
        sscanf(instr, "%c", &c) == 1 && isalpha(c))
        number = c - 'a' + 10;

    return number;
}

static bool proc_aux_property_coloring(int iw, char *instr, char **outstr)
{
    int i, number;

    // Allow aux props to be looked up by name as well as by number
    if (sscanf(instr, "%d", &number) != 1) {
      number = -1;
      for (i=0; i < CONFIG_num_auxiliary; i++) {
	if (strcmp(CONFIG_auxiliary_name[i], instr) == 0) {
	  number = i;
	  break;
	}
      }
    } 

    if (XIN(number, 0, CONFIG_MAX_AUXILIARY+MAX_GEO_MEASURES)) {
        n[iw].auxiliary_idx = number;
        if ( OUW(n[iw].auxiliary_idx,CONFIG_num_auxiliary) &&
             OUW(n[iw].auxiliary_idx-CONFIG_MAX_AUXILIARY,MAX_GEO_MEASURES) ) {
            sprintf(buf,"Auxiliary %d is not available.\n",n[iw].auxiliary_idx);
            *outstr = buf;
            return (FALSE);
        }
        *outstr = NULL;
        n[iw].color_mode = COLOR_MODE_AUXILIARY;
        return color_encode_auxiliary(iw);
    }

    *outstr = cui_show_syntax;
    return FALSE;
}

static bool proc_central_symmetry_coloring(int iw, char *instr, char **outstr)
{
    *outstr = NULL;
    n[iw].auxiliary_idx = 49;
    n[iw].color_mode = COLOR_MODE_AUXILIARY;
    n[iw].auxiliary_threshold[n[iw].auxiliary_idx][0] = 0.00376;
    n[iw].auxiliary_thresholds_rigid = 1;
    color_encode_auxiliary(iw);
    return TRUE ;
}

static bool proc_change_aux_property_threshold(int iw,char *instr,char **outstr)
{
    int idx = n[iw].auxiliary_idx, ul = 0;
    double value;

    *outstr = NULL;
    if (n[iw].color_mode != COLOR_MODE_AUXILIARY) 
        return FALSE;

    if (!instr)
        goto error;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        value = n[iw].auxiliary_threshold[idx][ul];
        sprintf(buf, "%d %.3f", ul, value);
        instr = cui_stripspace(gui_readline_gets(iw,
                        "\nChange uxiliary property threshold", buf));
    }

    buf[0] = 0;
    if (sscanf(instr, " %d%lf", &ul, &value) == 2 ||
        sscanf(instr, " %s%lf", buf, &value) == 2) {

        if (buf[0]) {
            if (strcmp(buf, "lower") == 0)
                ul = 0;
            else if (strcmp(buf, "upper") == 0)
                ul = 1;
            else
                goto error;
        }

        if (strstr(instr, "delta"))
            n[iw].auxiliary_threshold[idx][ul]
                    += ( n[iw].auxiliary_threshold[idx][1] -
                         n[iw].auxiliary_threshold[idx][0] )
                        * value * n[iw].delta;
        else if ((ul == 0 && value <= n[iw].auxiliary_threshold[idx][1]) ||
                 (ul == 1 && value >= n[iw].auxiliary_threshold[idx][0]))
            n[iw].auxiliary_threshold[idx][ul] = value;
        else
            goto error;

        if (AX_display[iw]) {
            if (idx < CONFIG_num_auxiliary)
                printf("Thresholds of %s = [%g, %g]%s%s.\n",
                       CONFIG_auxiliary_name[idx],
                       n[iw].auxiliary_threshold[idx][0],
                       n[iw].auxiliary_threshold[idx][1],
                       (*CONFIG_auxiliary_unit[idx]==EOS)?"":" ",
                         CONFIG_auxiliary_unit[idx]);
            else
                printf("Thresholds of %s = [%g, %g].\n",
                       geolist[idx-CONFIG_MAX_AUXILIARY].token, 
                       n[iw].auxiliary_threshold[idx][0],
                       n[iw].auxiliary_threshold[idx][1]);
        }

        return color_encode_auxiliary(iw);
    }

error:
    *outstr = "parameter error";
    return FALSE;
}


static bool proc_reset_aux_property_thresholds(int iw,char *instr,char **outstr)
{
    *outstr = NULL;
    if (n[iw].color_mode == COLOR_MODE_AUXILIARY) {
        int idx = n[iw].auxiliary_idx;
        reset_auxiliary_threshold(iw, idx);
        if (AX_display[iw])
            printf("\nAuxiliary[%d] = %s [%s]'s thresholds have been reset:\n"
                "Thresholds now = [%g, %g].\n\n",
                idx, CONFIG_auxiliary_name[idx], CONFIG_auxiliary_unit[idx],
                n[iw].auxiliary_threshold[idx][0],
                n[iw].auxiliary_threshold[idx][1]);
        return color_encode_auxiliary(iw);
    }
    return FALSE;
}

static bool proc_toggle_aux_property_thresholds_saturation(
                    int iw, char *instr, char **outstr)
{
    *outstr = NULL;
    if (n[iw].color_mode == COLOR_MODE_AUXILIARY) {
        n[iw].auxiliary_thresholds_saturation =
                        !n[iw].auxiliary_thresholds_saturation;
        return color_encode_auxiliary(iw);
    }
    else
        return FALSE;
}

static bool proc_toggle_aux_property_thresholds_rigid(
                    int iw,char *instr,char **outstr)
{
    *outstr = NULL;
    if (n[iw].color_mode == COLOR_MODE_AUXILIARY) {
        if (AX_display[iw])
            printf("floating auxiliary thresholds flag = %s.\n",
                        (n[iw].auxiliary_thresholds_rigid =
                        !n[iw].auxiliary_thresholds_rigid) ? "OFF" : "ON");
    }
    return FALSE;
}


static bool proc_rcut_patch(int iw, char *instr, char **outstr)
{
    if (!instr || !*instr) {
        *outstr = cui_show_syntax;
        return FALSE;
    }

    *outstr = NULL;
    buf[0] = 0;
    sscanf(instr, " %s", buf);
    if (strcmp(buf, "start") == 0 || strcmp(buf, "toggle") == 0) {
        char ini[3], inj[3], *q;
        int i, j, patch_Zi, patch_Zj, k, *count;

        if (rcut_patching) {
            if (strcmp(buf, "toggle") == 0)
                goto finish;
            else {
                if (AX_display[iw])
                    printf(
                        "You haven't yet finished patching the last rcut.\n"
                        "Press R to submit the last item.\n");
                return FALSE;
            }
        }
        if (rcut_patch_pairname[0] == EOS) { /* popularity contest */
            CALLOC (start_rcut_patch, count, SQUARE(ct->t), int);
            for (i=0; i<np; i++)
                for (j=N->idx[i]; j<N->idx[i+1]; j++)
                    count[ct->t*((int)tp[i])+((int)tp[N->list[j]])]++;
#ifdef USE_P3D
            if (p3dp_enabled) {
                int *c_l;
                CALLOC(start_rcut_patch, c_l, SQUARE(ct->t), int);
                for (i = 0; i < SQUARE(ct->t); i++)
                    c_l[i] = count[i];
                p3d_reduce(p3dp_cell, c_l, count,SQUARE(ct->t),MPI_INT,MPI_SUM);
                free(c_l);
            }
#endif
            for (i=0; i<ct->t; i++)
                for (j=i+1; j<ct->t; j++)
                    count[ct->t*i+j] = count[ct->t*j+i] =
                        count[ct->t*i+j] + count[ct->t*j+i];
            for (i=j=k=0; i<SQUARE(ct->t); i++)
                if (count[i] > j)
                {
                    j = count[i];
                    k = i;
                }
            free (count);
            i = k / ct->t;
            j = k % ct->t;
            sprintf (rcut_patch_pairname, "%s %s",
                     COMPACT_SYMBOL(ATOM_SYMBOL(ct->Z[i])),
                     COMPACT_SYMBOL(ATOM_SYMBOL(ct->Z[j])));
        }
        if (ct->t > 1) {
            if (strcmp(buf, "start") == 0)
                strncpy(buf,strstr(instr,"start") +strlen("start"),sizeof(buf));
            else
                strncpy(buf,strstr(instr,"toggle")+strlen("toggle"),
                                                                   sizeof(buf));
            instr = cui_stripspace(buf);
            if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
                strncpy(buf, rcut_patch_pairname, sizeof(buf));
                instr = cui_stripspace(gui_readline_gets(iw,
                    "\nStart patching neighbor distance cutoff between", buf));
            }
            q = (strlen(instr)) ? instr : rcut_patch_pairname;
            sscanf (q, "%2s %2s", ini, inj);
            SAFE_SYMBOL(ini);
            SAFE_SYMBOL(inj);
            patch_Zi = search_atom_by_symbol(ini);
            patch_Zj = search_atom_by_symbol(inj);
        }
        else {
            patch_Zi = patch_Zj = ct->Z[0];
            if (AX_display[iw])
                printf("\n"
                    "Start patching neighbor distance cutoff between %s %s:\n",
                    COMPACT_SYMBOL(ATOM_SYMBOL(patch_Zi)),
                    COMPACT_SYMBOL(ATOM_SYMBOL(patch_Zj)));
        }
        if (AX_display[iw])
            sprintf(rcut_patch_pairname, "%s %s",
                COMPACT_SYMBOL(ATOM_SYMBOL(patch_Zi)),
                COMPACT_SYMBOL(ATOM_SYMBOL(patch_Zj)));
        for (k=0; k<rcut_patch_top; k++)
            if ( ( (rcut_patch[k].Zi == patch_Zi) &&
                   (rcut_patch[k].Zj == patch_Zj) ) ||
                 ( (rcut_patch[k].Zi == patch_Zj) &&
                   (rcut_patch[k].Zj == patch_Zi) ) ) break;
        if (k == rcut_patch_top) { /* new one: */
            rcut_patch_item = rcut_patch_top++;
            rcut_patch[rcut_patch_item].Zi = MIN(patch_Zi,patch_Zj);
            rcut_patch[rcut_patch_item].Zj = MAX(patch_Zi,patch_Zj);
            rcut_patch[rcut_patch_item].rcut = NEIGHBORLIST_RCUT_RATIO *
                ( ATOM_RADIUS_IN_A(patch_Zi) + ATOM_RADIUS_IN_A(patch_Zj) );
        }
        else { /* the pair is already in patch list */
            rcut_patch_item = k;
        }
        if (rcut_patch_top >= RCUT_PATCH_MAX) {
            if (AX_display[iw]) {
                printf ("RCUT_PATCH_MAX = %d reached.\n", RCUT_PATCH_MAX);
                printf ("Cannot add rcut patch no more.\n\n");
            }
            rcut_patch_top--;
            return FALSE;
        }
        else rcut_patching = 1;
        if (AX_display[iw])
            printf ("RCUT(%s) = %g.\n", rcut_patch_pairname,
                rcut_patch[rcut_patch_item].rcut);
        return FALSE;
    }
    else if (strcmp(buf, "finish") == 0) {
finish:
        if (rcut_patching) {
            if (AX_display[iw]) {
                printf ("RCUT(%s) = %g.\n\n",
                    rcut_patch_pairname, rcut_patch[rcut_patch_item].rcut);
                Sfpr(stdout, "rlist = %M (reduced)\n ", N->rcut, ct->t,ct->t);
            print_coordination_histogram(); cr();
            }
            rcut_patching = 0;
        }
        return FALSE;
    }
    else if (rcut_patching) {
        double rcut, c = 0.0;
        if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
            instr = cui_stripspace(
                            gui_readline_gets(iw, "\nChange cutoff", "0.0"));
        buf[0] = 0;
        if (sscanf(instr, " %lf %[^ ]", &c, buf) == 2) {
            if (strcmp(buf, "delta") == 0)
                c *= n[iw].delta;
            else {
                *outstr = cui_show_syntax;
                return FALSE;
            }
        }
        rcut = rcut_patch[rcut_patch_item].rcut + c;
        if (rcut < 0)
            rcut = 0;
        else if (c > 0) {
            int i;
            double thickness[DIMENSION];
            M3rowthicknesses(H, thickness);
            for (i = 0; i < DIMENSION; i++) {
                if (thickness[i] < 2 * rcut) {
                    sprintf(buf, "losing r<rcut images may happen "
                            "in direction %d, request ignored", i);
                    *outstr = buf;
                    return FALSE;
                }
            }
        }
        rcut_patch[rcut_patch_item].rcut = rcut;
        if (AX_display[iw])
            printf("rcut(%s) = %g.\n", rcut_patch_pairname,
                                rcut_patch[rcut_patch_item].rcut);
        return apply_rcut_patch(iw);
    }
    return FALSE;
}

static bool proc_start_rcut_patch(int iw, char *instr, char **outstr)
{
    char instr2[CUI_LINEMAX] = "start ";
    if (instr)
        strncat(instr2, instr, sizeof(instr2));
    return proc_rcut_patch(iw, instr2, outstr);
}

static bool proc_finish_rcut_patch(int iw, char *instr, char **outstr)
{
    return proc_rcut_patch(iw, "finish", outstr);
}

static bool proc_toggle_rcut_patch_mode(int iw, char *instr, char **outstr)
{
    char instr2[CUI_LINEMAX] = "toggle ";
    if (instr)
        strncat(instr2, instr, sizeof(instr2));
    return proc_rcut_patch(iw, instr2, outstr);
}


static bool proc_select_gear(int iw, char *instr, char **outstr)
{
    int number = get_number(iw, instr);
    if (XIN(number, 0, 9)) {
        n[iw].delta = gearbox[number];
        *outstr = NULL;
        return FALSE;
    }
    *outstr = cui_show_syntax;
    return FALSE;
}


static bool new_plane(int iw, char *instr, char **outstr, int number)
{
    double tmp[3], dxkj[4], dxij[4];
    char danswer[CUI_LINEMAX];

    if ( V3EQZERO(n[iw].fp[number].dx_input) )
        V3ASSIGN(1,1,1,n[iw].fp[number].dx_input);
    V3mM3 (n[iw].hook, HI, n[iw].fp[number].s0);

    if (n[iw].atom_stack[0]>=0&&n[iw].atom_stack[1]>=0&&n[iw].atom_stack[2]>=0){
#if USE_P3D
      if (p3dp_enabled) {
        double sk[3], sj[3], si[3];
        if (            n[iw].atom_stack[2] !=      n[iw].atom_stack[0] ||
                   p3dp_n[iw].rank_stack[2] != p3dp_n[iw].rank_stack[0]) {
            if (        n[iw].atom_stack[0] ==      n[iw].atom_stack[1] &&
                   p3dp_n[iw].rank_stack[0] == p3dp_n[iw].rank_stack[1]) {
                p3dp_s(sk, n[iw].atom_stack[0], p3dp_n[iw].rank_stack[0]);
                p3dp_s(si, n[iw].atom_stack[2], p3dp_n[iw].rank_stack[2]);
                p3dp_atom_pair_s(sk, si, n[iw].fp[number].dx_input);
                V3ADDMULMUL(0.5, sk, 0.5, si, n[iw].fp[number].s0);
            }
            else if (   n[iw].atom_stack[1] !=      n[iw].atom_stack[2] ||
                   p3dp_n[iw].rank_stack[1] != p3dp_n[iw].rank_stack[2]) {
                p3dp_s(sk, n[iw].atom_stack[0], p3dp_n[iw].rank_stack[0]);
                p3dp_s(sj, n[iw].atom_stack[1], p3dp_n[iw].rank_stack[1]);
                p3dp_s(si, n[iw].atom_stack[2], p3dp_n[iw].rank_stack[2]);
                p3dp_atom_triplet_s(sk, sj, si, dxkj, dxij);
                V3CROSS(dxkj, dxij, n[iw].fp[number].dx_input);
                V3EQV(sk, n[iw].fp[number].s0);
            }
        }
      }
      else {
#endif
        if (    n[iw].atom_stack[2] != n[iw].atom_stack[0]) {
            if (n[iw].atom_stack[0] == n[iw].atom_stack[1]) {
                atom_pair (n[iw].atom_stack[0], n[iw].atom_stack[2],
                           n[iw].fp[number].dx_input);
                V3ADDMULMUL (0.5, &s[DIMENSION*n[iw].atom_stack[0]],
                             0.5, &s[DIMENSION*n[iw].atom_stack[2]],
                             n[iw].fp[number].s0);
            }
            else if (n[iw].atom_stack[1] != n[iw].atom_stack[2]){
                atom_triplet (n[iw].atom_stack[0], n[iw].atom_stack[1],
                              n[iw].atom_stack[2], dxkj, dxij);
                V3CROSS (dxkj, dxij, n[iw].fp[number].dx_input);
                V3EQV (&s[DIMENSION*n[iw].atom_stack[0]], n[iw].fp[number].s0);
            }
        }
#if USE_P3D
      }
#endif
    }

    if (!instr || !*instr || strstr(instr, CUI_ARG_REQUIRED)) {
        sprintf (danswer, "%g %g %g %g %g %g",
                 V3E(n[iw].fp[number].dx_input), V3E(n[iw].fp[number].s0));
        sprintf(buf, "\nCutting plane %d's dx dy dz s0 s1 s2", number);
        instr = gui_readline_gets(iw, buf, danswer);
    }

    if (sscanf(instr, "%lf %lf %lf %lf %lf %lf",
               V3e(AX_3D[iw].fp[number].dx), V3e(n[iw].fp[number].s0)) < 3) {
        *outstr = "parameter error";
        return FALSE;
    }
    else if (V3NEZERO(AX_3D[iw].fp[number].dx)) {
        V3EQV( AX_3D[iw].fp[number].dx, n[iw].fp[number].dx_input );
        V3normalize (AX_3D[iw].fp[number].dx);
        AX_3D[iw].fp[number].d0 =
            V3ADOT(n[iw].fp[number].s0, H, AX_3D[iw].fp[number].dx, tmp);
        *outstr = NULL;
        return TRUE;

    }
    else {
        if (AX_display[iw])
            printf("(%g %g %g) is unacceptable as a plane normal!\n",
                   V3E(AX_3D[iw].fp[number].dx));
        *outstr = danswer;
        return FALSE;
    }
}


static bool proc_cutting_plane(int iw, char *instr, char **outstr)
{
    int number = get_number(iw, instr), j;

    if (!XIN(number, 0, 15)) {
        *outstr = "parameter error";
        return FALSE;
    }
    else
        *outstr = NULL;

    if (n[iw].anchor >= 0)
        V3EQV(B->BALL[n[iw].anchor].x, n[iw].hook);

    if (V3EQZERO(AX_3D[iw].fp[number].dx)) {
    /* not active but will be activated and focused */
        n[iw].just_activated_fp = number;
        if (V3EQZERO(n[iw].fp[number].dx_cache)) { /* no cached plane normal */
            char *s = NULL;
            if (instr && strcmp(instr, CUI_ARG_REQUIRED) != 0) {
                s = instr;
                while (*s && isspace(*s)) s++;
                if (*s) {
                    s++;
                    while (*s && !isspace(*s)) s++;
                }
            }
            return new_plane(iw, s, outstr, number);
        }
        else { /* there is cached plane normal */
            double tmp[3];
            V3EQV( n[iw].fp[number].dx_cache, AX_3D[iw].fp[number].dx );
            AX_3D[iw].fp[number].d0 =
                V3ADOT(n[iw].fp[number].s0, H, AX_3D[iw].fp[number].dx,
                       tmp);
            return TRUE;
        }
    }
    else if (n[iw].just_activated_fp != number) {
    /* activated but not the focus */
        n[iw].just_activated_fp = number;
        /* just gain focus */
        V3EQV( AX_3D[iw].fp[number].dx, n[iw].fp[number].dx_cache );
        V3ZERO( AX_3D[iw].fp[number].dx );
        for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
            if ( V3NEZERO(AX_3D[iw].fp[j].dx) )
                n[iw].just_activated_fp = j;
        return TRUE;
    }
    else {
    /* activated and is the focus */
        V3EQV( AX_3D[iw].fp[number].dx, n[iw].fp[number].dx_cache );
        V3ZERO( AX_3D[iw].fp[number].dx );
        for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
            if ( V3NEZERO(AX_3D[iw].fp[j].dx) )
                n[iw].just_activated_fp = j;
        return TRUE;
    }
    return FALSE;
}


static bool proc_shift_cutting_plane_to_anchor(int iw,char *instr,char **outstr)
{
    int number = get_number(iw, instr);
    double tmp[3];

    if (!XIN(number, 0, 15)) {
        *outstr = "parameter error";
        return FALSE;
    }
    else
        *outstr = NULL;

    if (n[iw].anchor >= 0)
        V3EQV(B->BALL[n[iw].anchor].x, n[iw].hook);

    if (V3EQZERO(AX_3D[iw].fp[number].dx)) {
    /* not active but will be activated and focused */
        n[iw].just_activated_fp = number;
        if (V3EQZERO(n[iw].fp[number].dx_cache)) { /* no cached plane normal */
            return FALSE;
        }
        else { /* there is cached plane normal */
            V3mM3 (n[iw].hook, HI, n[iw].fp[number].s0);
            V3EQV( n[iw].fp[number].dx_cache, AX_3D[iw].fp[number].dx );
            AX_3D[iw].fp[number].d0 =
                V3ADOT(n[iw].fp[number].s0, H, AX_3D[iw].fp[number].dx, tmp);
            return TRUE;
        }
    }
    else if (n[iw].just_activated_fp != number) {
    /* activated but not the focus */
        n[iw].just_activated_fp = number;
        V3mM3 (n[iw].hook, HI, n[iw].fp[number].s0);
        AX_3D[iw].fp[number].d0 =
            V3ADOT(n[iw].fp[number].s0, H, AX_3D[iw].fp[number].dx, tmp);
        return TRUE;
    }
    else {
    /* activated and is the focus */
        V3mM3 (n[iw].hook, HI, n[iw].fp[number].s0);
        AX_3D[iw].fp[number].d0 =
            V3ADOT(n[iw].fp[number].s0, H, AX_3D[iw].fp[number].dx, tmp);
        return TRUE;
    }
    return FALSE;
}


static bool proc_delete_cutting_plane(int iw, char *instr, char **outstr)
{
    int number = get_number(iw, instr), j;

    if (!XIN(number, 0, 15)) {
        *outstr = "parameter error";
        return FALSE;
    }
    else
        *outstr = NULL;

    if (n[iw].anchor >= 0)
        V3EQV(B->BALL[n[iw].anchor].x, n[iw].hook);

    if (V3EQZERO(AX_3D[iw].fp[number].dx)) {
    /* not active but will be activated and focused */
        n[iw].just_activated_fp = number;
        if (V3EQZERO(n[iw].fp[number].dx_cache)) { /* no cached plane normal */
            return FALSE;
        }
        else { /* there is cached plane normal */
            V3ZERO( n[iw].fp[number].dx_cache );
            for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
                if ( V3NEZERO(AX_3D[iw].fp[j].dx) ) n[iw].just_activated_fp = j;
            return FALSE;
        }
    }
    else if (n[iw].just_activated_fp != number) {
    /* activated but not the focus */
        n[iw].just_activated_fp = number;
        V3ZERO( n[iw].fp[number].dx_cache );
        V3ZERO( AX_3D[iw].fp[number].dx );
        for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
            if ( V3NEZERO(AX_3D[iw].fp[j].dx) )
                n[iw].just_activated_fp = j;
        return TRUE;
    }
    else {
    /* activated and is the focus */
        V3ZERO( n[iw].fp[number].dx_cache );
        V3ZERO( AX_3D[iw].fp[number].dx );
        for (j=0; j<AX_3D_MAX_FILTER_PLANE; j++)
            if ( V3NEZERO(AX_3D[iw].fp[j].dx) )
                n[iw].just_activated_fp = j;
        return TRUE;
    }
    return FALSE;
}


static bool proc_flip_cutting_plane(int iw, char *instr, char **outstr)
{
    int i = n[iw].just_activated_fp;
    *outstr = NULL;
    if (V3NEZERO(AX_3D[iw].fp[i].dx))
    {
        AX_3D[iw].fp[i].d0 = -AX_3D[iw].fp[i].d0;
        V3NeG ( AX_3D[iw].fp[i].dx );
        return (TRUE);
    }
    return (FALSE);
}


static bool proc_capture(int iw, char *instr, char **outstr)
{
    int new_resolution, writable;
    int new_res;
    char danswer[MAX_FILENAME_SIZE] = "", *answer,
        buffer[MAX_FILENAME_SIZE]={' ', ' ', ' ', ' '}, *fname;
    struct CS {
        Pixmap new_pixmap;
        Drawable old_drawable;
        GC old_gc;
        AXSize old_size;
        double scale;
    } cs;
    extern void (*cui_CaptureResize)(int, int, struct CS *);
    extern void (*cui_CaptureRecover)(int, struct CS *);

    *outstr = cui_show_syntax;
    buf[0] = '.';
    if (!instr || sscanf(instr, " %[^ ] %[^\n]", buf+1, danswer) <= 0 ||
                strlen(buf) != 4)
        return FALSE;
    buf[1] = tolower(buf[1]);
    buf[2] = tolower(buf[2]);
    buf[3] = tolower(buf[3]);
    if (!(strcmp(buf,".png")==0||strcmp(buf,".jpg")==0||strcmp(buf,".eps") ==0))
        return FALSE;

    if (buf[1] == 'e') 
        new_resolution = (CUI_EPS_RESOLUTION_DEFAULT > AX_MAXWIDTH) ?
                            AX_MAXWIDTH : CUI_EPS_RESOLUTION_DEFAULT;
    else
        new_resolution = MAX(AX_size[iw].width, AX_size[iw].height);
    new_res = new_resolution;
    fname = buffer+4;
    if (strcmp(danswer, CUI_ARG_REQUIRED) == 0) {
        sprintf (danswer, "%s%s %d", fbasename, buf, new_res);
        answer = gui_readline_gets(iw, "\nSave screen on", danswer);
        sscanf(answer, "%s %d", fname, &new_res);
    }
    else
        sscanf(danswer, "%s %d", fname, &new_res);
    if (new_res > 1) 
        new_resolution = new_res;
    else if (new_res < -1) 
        new_resolution = -new_res;
    if (new_resolution <= 1 || !*fname)
        return FALSE;

    *outstr = NULL;
    if (new_res > 0) {
        if (strcasecmp(eos(fname)-4, buf)) strcat(fname, buf);
        if (n[iw].color_mode == COLOR_MODE_AUXILIARY)
            save_auxiliary_colormap(iw,str2(fname, ".cmap.eps"));
    }
    if (IS_MANAGER)
        writable = (tested_to_be_writable(fname));
#ifdef USE_P3D
    if (p3dp_enabled)
        p3d_bcast(p3dp_cell, &writable, 1, MPI_INT, 0);
#endif
    if (writable)
    {
        cui_CaptureResize(iw, new_resolution, &cs);
#ifdef USE_P3D
        if (p3dp_enabled)
                MPI_Barrier(p3d_comm(p3dp_cell));
#endif
        if (IS_MANAGER) {
            switch (buf[1]) {
            case 'p':
                sprintf(danswer, "image saved on \"%s\" (%d bytes)\n", fname,
                                            AX_save_pixmap_as_PNG (iw,fname));
                break;
            case 'j':
                sprintf(danswer, "image saved on \"%s\" (%d bytes)\n", fname,
                                            AX_save_pixmap_as_JPG (iw,fname));
                break;
            case 'e':
                sprintf(danswer, "image saved on \"%s\" (%d bytes)\n", fname,
                                            AX_save_pixmap_as_EPS (iw,fname));
                sprintf (fname, "-geometry %dx%d+%d+%d %s",
                         2*cs.old_size.width, 2*cs.old_size.height,
                         ZERO_RAMP(AX_root_x[iw]-cs.old_size.width),
                         ZERO_RAMP(AX_root_y[iw]-cs.old_size.height/2),
                         absolute_pathname(fname));
                break;
            }
            if (AX_display[iw] && !cui_xterm_win)
                printf(danswer);
            strncat(danswer, CUI_PROTOCOL_OK, sizeof(danswer));
            *outstr = danswer;
        }
        AXQueryPointer (iw);
        cui_CaptureRecover(iw, &cs);
        if (n[iw].shell_viewer_mode && IS_MANAGER) {
            switch (buf[1]) {
            case 'p':
            case 'j':
                sprintf (fname, "-geometry +%d+%d %s",
                     ZERO_RAMP(AX_root_x[iw]-AX_size[iw].width/2),
                     ZERO_RAMP(AX_root_y[iw]-AX_size[iw].height/4),
                     absolute_pathname(fname));
                try_to_runbg (NUMBER_RASTER_VIEWERS, raster_viewers, fname);
                break;
            case 'e':
                try_to_runbg(NUMBER_POSTSCRIPT_VIEWERS,
                        postscript_viewers, fname);
                break;
            }
        }
    }
    else
    {
        if (AX_display[iw])
        printf ("\n** %s: **\n""** This file is unwritable! **\n", fname);
        *outstr = "file is unwritable";
    }

    cui_diligence = TRUE;
    return FALSE;
}

static bool proc_change_wireframe_mode(int iw, char *instr, char **outstr)
{
    int number = get_number(iw, instr);

    if (XIN(number, 0, 4)) {
        n[iw].wireframe_mode = number;
    }
    else
        n[iw].wireframe_mode = (n[iw].wireframe_mode+1) % 5;

    *outstr = NULL;
    return TRUE;
}


static bool proc_change_cutting_plane_wireframe_mode
(int iw, char *instr, char **outstr)
{
    int i = n[iw].just_activated_fp;
    int number = get_number(iw, instr);

    *outstr = NULL;
    if (!V3NEZERO(AX_3D[iw].fp[i].dx))
        return FALSE;

    
    if (XIN(number, 0, FILTER_PLANE_WIREFRAME_GRADES))
        n[iw].fp[i].wireframe_mode = number;
    else
        n[iw].fp[i].wireframe_mode = ( n[iw].fp[i].wireframe_mode + 1 ) %
                                   (FILTER_PLANE_WIREFRAME_GRADES + 1 );

    return TRUE;
}


static bool proc_load_config(int iw, char *instr, char **outstr)
{
    int i, j, k, old_np;
    char fname[MAX_FILENAME_SIZE], oldfname[MAX_FILENAME_SIZE]; 
    V3 hook_s, tmp, dx;
    char *old_symbol=NULL;
    bool incompatible_config;

    *outstr = NULL;
    if (!instr)
        goto error;
    else if (!*instr)
        instr = config_fname;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
        instr = cui_stripspace(
            gui_readline_gets(iw, "\nLoad configuration", config_fname));

    strcpy(oldfname, config_fname);
    if (n[iw].anchor >= 0) {
        /* the new configuration may not even have the atom */
        V3EQV (B->BALL[n[iw].anchor].x, n[iw].hook);
        n[iw].anchor = -1;
    }
    /* hook_s[] is what is kept invariant */
    V3mM3 (n[iw].hook, HI, hook_s);

    strncpy(fname, instr, sizeof(fname));
    strcpy(config_fname,fname);

#ifdef USE_P3D
    if (!p3dp_enabled) {
#endif
      if (!strstr(config_fname, ".nc") && !strstr(config_fname, ".xyz")) {
    if (!Fexists(config_fname)) {
        if (AX_display[iw])
            printf("\n** %s: **\n** There is no such file! **\n", config_fname);
        strcpy(config_fname, oldfname);
        return FALSE;
    }
    if (!Freadable(config_fname)) {
        if (AX_display[iw])
            printf("\n** %s: **\n** This file is unreadable! **\n",
                                                                config_fname);
        strcpy(config_fname, oldfname);
        return FALSE;
    }
      }
#ifdef USE_P3D
    }
#endif
#ifndef ATOMEYE_LIB
    if (AX_display[iw]) printf("\n");
#endif

    old_np = np;
    CLONE(symbol, SYMBOL_SIZE*np, char, old_symbol);
    i = CONFIG_LOAD(config_fname, Config_Aapp_to_Alib);

    for (k=0; k<CONFIG_num_auxiliary; k++)
        if (*blank_advance(CONFIG_auxiliary_name[k])==EOS)
            sprintf(CONFIG_auxiliary_name[k], "auxiliary%d", k);
#ifndef ATOMEYE_LIB
    rebind_CT (Config_Aapp_to_Alib, "", ct, &tp); cr();
#else
    rebind_ct (Config_Aapp_to_Alib, "", ct, &tp, NULL); cr();
#endif
    Neighborlist_Recreate_Form (Config_Aapp_to_Alib, ct, N);
    if (i == CONFIG_CFG_LOADED)
        N->s_overflow_err_handler =
            NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_FOLD_INTO_PBC;
    else
        N->s_overflow_err_handler =
            NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_BOUNDING_BOX;
    N->small_cell_err_handler = NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_MULTIPLY;
    for (i=0; i<ct->t; i++)
        for (j=i; j<ct->t; j++)
            for (k=0; k<rcut_patch_top; k++)
                if ( ( ( (rcut_patch[k].Zi == ct->Z[i]) &&
                         (rcut_patch[k].Zj == ct->Z[j]) ) ||
                       ( (rcut_patch[k].Zi == ct->Z[j]) &&
                         (rcut_patch[k].Zj == ct->Z[i]) ) ) )
                    NEIGHBOR_TABLE(N->rcut,ct,i,j) =
                        NEIGHBOR_TABLE(N->rcut,ct,j,i) = rcut_patch[k].rcut;
    Neighborlist_Recreate (Config_Aapp_to_Alib, stdout, ct, &tp, N);
    V3mM3 (hook_s, H, tmp);
    V3SUB (tmp, n[iw].hook, dx);
    V3EQV (tmp, n[iw].hook);
    V3AdD (dx, AX_3D[iw].x);
    M3InV (H, HI, volume);
    lengthscale = cbrt(volume);
    V3ASSIGN (0.5,0.5,0.5,tmp);
    V3mM3 (tmp, H, cm);
    geo_clear_has_evaluated_flags();
    evaluate_geo_measures(); Free(s1); Free(mass);

    if  ( ComputeLeastSquareStrain )
    {
        if (ConfigChecksum(Config_Aapp_to_Alib) != ref->checksum) {
            static char buf2[CUI_LINEMAX];
            snprintf(buf2, sizeof(buf2),
                    "This configuration is not isoatomic with the imprinted "
                    "reference\n%s. Least-square strain NOT calculated.\n",
                    ref_fbasename);
            *outstr = buf2;
            fprintf(stderr, "%s", buf2);
        }
        else LeastSquareStrain_Append();
    }

    incompatible_config = (np != old_np) ||
        memcmp(symbol, old_symbol, SYMBOL_SIZE*MIN(np,old_np));
    Free(old_symbol);

    if (incompatible_config)
        Config_to_3D_Balls (n[iw].atom_r_ratio);
    else for (i=0; i<np; i++) V3mM3 ( &(s[DIMENSION*i]), H, B->BALL[i].x );

    atom_xtal_origin (n[iw].xtal_origin);
    if (!n[iw].auxiliary_thresholds_rigid) {
        for (i=0; i<CONFIG_num_auxiliary; i++)
            reset_auxiliary_threshold(iw,i);
        for (i=0; i<MAX_GEO_MEASURES; i++)
            if (geolist[i].has_evaluated)
                reset_auxiliary_threshold(iw,CONFIG_MAX_AUXILIARY+i);
    }
    if (!temporary_disable_bond) Config_to_3D_Bonds (n[iw].bond_radius);
    select_fbasename (config_fname);
    if ((n[iw].xtal_mode) && (n[iw].color_mode == COLOR_MODE_COORD))
        assign_coordination_color(iw);
    else if (n[iw].color_mode == COLOR_MODE_AUXILIARY)
        color_encode_auxiliary(iw);
    else if (n[iw].color_mode == COLOR_MODE_SCRATCH)
        /*scratch_color (iw);*/
        proc_scratch_coloring(iw, instr, outstr);
    else {
        strcpy (AX_title[iw],fbasename);
        AXSetName (iw);
        XStoreName(AX_display[iw],xterm_win,AX_title[iw]);
        XSetIconName(AX_display[iw],xterm_win,AX_title[iw]);
        if (!temporary_disable_bond) {
            bond_xtal_origin_update (iw);
            bond_atom_color_update (iw);
        }
    }
    return TRUE;

error:
    *outstr = "parameter error";
    return FALSE;
}


static bool proc_load_config_advance(int iw, char *instr, char **outstr)
{
  char *p;

#ifdef ATOMEYE_LIB
  (*atomeyelib_on_advance)(instr);
  return FALSE;
#endif
    *outstr = "parameter error";
    if (strstr(config_fname, ".nc") || strstr(config_fname, ".xyz")) {
      // Remove trailing frame number from current filename 
      // and append instruction "forward", "backward", "first" or "last"
      strcpy (buf, config_fname);
      p = buf;
      strsep(&p, ":");
      strcat(buf, ":");
      strcat(buf, instr);
      if(strcmp(instr,"forward") == 0 || strcmp(instr,"backward") == 0)
	sprintf(buf, "%s:%d", buf, n[iw].glob_advance);
      return proc_load_config(iw, buf, outstr);
    } else {
      if (instr) {
        int how_much;
        strncpy (buf, config_fname, sizeof(buf));
        if (sscanf(instr, " %d", &how_much) == 1)
            Numerically_sorted_glob_advance(buf, how_much);
        else if (strstr(instr, "forward")) 
            Numerically_sorted_glob_advance(buf,  n[iw].glob_advance);
        else if (strstr(instr, "backward")) 
            Numerically_sorted_glob_advance(buf, -n[iw].glob_advance);
        else if (strstr(instr, "first")) 
            Numerically_sorted_glob_first(buf);
        else if (strstr(instr, "last")) 
            Numerically_sorted_glob_last(buf);
	else if (strstr(instr, "reload"))
	  ; // do nothing - reload same file
        else
            return FALSE;

        *outstr = NULL;
        return proc_load_config(iw, buf, outstr);
      }
    }
    return FALSE;
}

#define SCRIPT_LINE_SIZE  512
#define SCRIPT_LINE_CHAR (SCRIPT_LINE_SIZE-1)
/* Use script to produce jpeg frames */
static bool proc_script_animate(int iw, char *instr, char **outstr)
{
    char output[SCRIPT_LINE_SIZE],*ptr = NULL;
    int i,hascontent,quality,frames;
    FILE *fp = NULL;
    glob_t globbuf;
    
    *outstr = "parameter error";
    if (!instr)
        return FALSE;
    else if (!*instr || strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        strncpy(buf, "scr_anim", sizeof(buf));
        i = 1;
    }
    else {
        strncpy(buf, instr, sizeof(buf));
        i = 0;
    }

    if (IS_MANAGER) {
        if (!(fp=ropen(buf)) && i) {
            if (AX_display[iw]) {
                printf("\nAnimation script \"%s\" does not exist,\n", buf);
                xterm_get_focus(iw); clear_stdin_buffer();
                if (!strcasecmp("y", readline_gets
                            ("Do you want a default one created (y/n)?","y"))) {
                    numerically_sorted_glob (config_fname, &globbuf);
                    fp=wopen(buf);
                    fprintf (fp, "%d\n", AX_JPG_DEF_QUALITY);
                    for (i=0; i<globbuf.gl_pathc; i++)
                        fprintf (fp,"%s Jpg/%05d.jpg\n",globbuf.gl_pathv[i],i);
                    globfree (&globbuf);
                    fclose(fp);
                    fp = ropen(buf);
                }
                else {
                    xterm_release_focus(iw);
                    *outstr = "Animation script does not exist.";
                }
            }
        }
        if (fp && !(ptr=fgets(buf,SCRIPT_LINE_SIZE,fp))) {
            *outstr = "There is nothing in animation script.";
            if (AX_display[iw]) printf("\n%s\n", *outstr);
            fclose(fp);
            fp = NULL;
        }
        if (fp) {
            for (hascontent=i=0; (buf[i]!=EOS) &&
                    (ISDIGIT(buf[i]) || ISBLANK(buf[i]) || (buf[i]=='\n')); i++)
                if (ISALNUM(buf[i])) hascontent=1;
            if (!hascontent) {
                *outstr = "There is no content in animation script.";
                fclose(fp);
                fp = NULL;
            }
        }
        if (fp) {
            if (buf[i] == EOS) {
                sscanf (buf, "%d", &quality);
                if ((quality<0) || (quality>100))
                {
                    *outstr = "quality is out of valid range ([0,100]).";
                    fclose(fp);
                    fp = NULL;
                }
                else if (AX_display[iw]) printf ("\nquality = %d\n", quality);
                ptr = fgets(buf,SCRIPT_LINE_SIZE,fp);
            }
            else
                quality = AX_JPG_DEF_QUALITY;
        }
        if (!fp)
            quality = -1;
    }
#ifdef USE_P3D
    if (p3dp_enabled) p3d_bcast(p3dp_cell, &quality, 1, MPI_INT, 0);
#endif
    if (quality < 0)
        return FALSE;

    *outstr = NULL;
    frames = 0;
    /* If bonds are not on now, there is no need to refresh */
    temporary_disable_bond = !n[iw].bond_mode;
    /* cylinder data structure during the rendering.        */
    while (1) {
        int count;
        if (IS_MANAGER)
            count = (ptr) ? strlen(buf) + 1 : -1;
#ifdef USE_P3D
        if (p3dp_enabled) p3d_bcast(p3dp_cell, &count, 1, MPI_INT, 0);
#endif
        if (count < 0) break;
#ifdef USE_P3D
        if (p3dp_enabled) p3d_bcast(p3dp_cell, buf, count, MPI_CHAR, 0);
#endif
        buf[SCRIPT_LINE_CHAR] = EOS;
        sscanf (buf, "%s %s", config_fname, output);

        proc_load_config(iw, config_fname, outstr);

        paint_scene(iw);
        AX_dump(iw); AX_show(iw);
        if (str_caseend_with(output,".png"))
            AX_save_pixmap_as_png
                (iw,AX_JPG_QUALITY_TO_PNG_LEVEL(quality),output);
        else if (str_caseend_with(output,".eps"))
            AX_save_pixmap_as_eps(iw,quality,output);
        else AX_save_pixmap_as_jpg(iw,quality,output);
        frames++;
        if (IS_MANAGER)
            ptr = fgets(buf,SCRIPT_LINE_SIZE,fp);
    }
    if (IS_MANAGER) {
        fclose(fp);
        if (AX_display[i])
            printf("%d frames saved.\n\n", frames);
    }
    if (temporary_disable_bond) {
        Config_to_3D_Bonds (n[iw].bond_radius);
        if (n[iw].bond_mode) {
            bond_xtal_origin_update (iw);
            bond_atom_color_update(iw);
        }
        else {
            n[iw].bond_xtal_origin_need_update = TRUE;
            n[iw].bond_atom_color_need_update = TRUE;
        }
        temporary_disable_bond = 0;
    }
    return FALSE;
}
#undef SCRIPT_LINE_CHAR
#undef SCRIPT_LINE_SIZE


/* Load color/radii file for atoms */
static bool proc_load_atom_color(int iw, char *instr, char **outstr)
{
    char fname[MAX_FILENAME_SIZE];
    FILE *fp;
    int m, j, k, items;
    double c[4];

#ifdef USE_P3D
    if (p3dp_enabled) {
        *outstr = "not supported";
        return FALSE;
    }
#endif
    *outstr = NULL;
    if (!instr) {
        *outstr = "parameter error";
        return FALSE;
    }
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        sprintf (buf, "%s.clr", fbasename);
        instr = gui_readline_gets(iw, "\nLoad color properties from", buf);
    }
    strncpy(fname, instr, sizeof(fname));
    if (Freadable(fname)) {
        fp = ROpen(fname);
        for (m=0; ; m++)
            if (fgets(buf,sizeof(buf),fp)) {
                if (m >= np) {
                    sprintf(buf, "** %s has more rows than atoms **", fname);
                    *outstr = buf;
                    return FALSE;
                }
                items = sscanf (buf, "%lf %lf %lf %lf", c, c+1, c+2, c+3);
                if (items == 1) {
                    c[3] = c[0];
                    c[0] = B->BALL[m].r;
                    c[1] = B->BALL[m].g;
                    c[2] = B->BALL[m].b;
                }
                else if (items == 3)
                    c[3] = B->BALL[m].radius / n[iw].atom_r_ratio;
                if (c[0]>1) c[0]/=255;
                if (c[1]>1) c[1]/=255;
                if (c[2]>1) c[2]/=255;
                AX_3D_AssignRGB(B->BALL[m], c[0], c[1], c[2]);
                B->BALL[m].radius = c[3] * n[iw].atom_r_ratio;
            }
            else if (m < np) {
                if (m == 0) {
                    sprintf(buf, "** %s has no data **", fname);
                    *outstr = buf;
                    return FALSE;
                }
                else if ( ISFACTOR (m, np) ) { /* make a bold guess */
                    for (j=m; j<np; j++) {
                        k = j % m;
                        AX_3D_AssignRGB (B->BALL[j], B->BALL[k].r,
                                         B->BALL[k].g, B->BALL[k].b);
                        B->BALL[j].radius = B->BALL[k].radius;
                    }
                }
                else {
                    sprintf(buf, "** premature ending of %s **", fname);
                    *outstr = buf;
                }
                bond_atom_color_update (iw);
                return TRUE;
            }
            else break;
    }
    else
    {
        sprintf(buf, "** %s: **\n** This file is unreadable! **", fname);
        *outstr = buf;
        return FALSE;
    }
    bond_atom_color_update(iw);
    return TRUE;
}


static bool proc_load_aux(int iw, char *instr, char **outstr)
{
    char fname[MAX_FILENAME_SIZE];
    FILE *fp;
    int i,j,k,m;
    SimpleStatistics ss;

#ifdef USE_P3D
    if (p3dp_enabled) {
        *outstr = "not supported";
        return FALSE;
    }
#endif
    *outstr = NULL;
    if (!instr) {
        *outstr = "parameter error";
        return FALSE;
    }
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        sprintf (buf, "%s.aux", fbasename);
        instr = gui_readline_gets(iw, "\nLoad auxiliary properties from", buf);
    }
    strncpy(fname, instr, sizeof(fname));

    if (Freadable(fname)) {
        fp = ropen(fname);
        CONFIG_num_auxiliary = CONFIG_MAX_AUXILIARY;
        for (k=0; k<CONFIG_num_auxiliary; k++)
            REALLOC (load_auxiliary_from_file,
                     CONFIG_auxiliary[k], np, double);
        for (m=0; ; m++)
            if (fgets(buf,sizeof(buf),fp)) {
                if (m >= np) {
                    sprintf(buf, "** %s has more rows than atoms **", fname);
                    *outstr = buf;
                    goto safe_exit;
                }
                for (i=j=0; ;) {
                    while ( ISNOTDIGIT(buf[i]) && (buf[i]!='.') &&
                            (buf[i]!=EOS) ) i++;
                    if (buf[i] == EOS) {
                        if (m==0) {
                            CONFIG_num_auxiliary = j;
                            if (AX_display[iw])
                                printf("number of auxiliaries found = %d.\n",j);
                            if (j==0) goto safe_exit;
                            for (k=j; k<CONFIG_MAX_AUXILIARY; k++)
                                Free(CONFIG_auxiliary[k]);
                            for (k=0; k<CONFIG_num_auxiliary; k++) {
                                sprintf(CONFIG_auxiliary_name[k],
                                        "auxiliary%d", k);
                                CONFIG_auxiliary_unit[k][0] = EOS;
                            }
                        }
                        else if (j != CONFIG_num_auxiliary) {
                            sprintf(buf, "** %s corrupted **", fname);
                            *outstr = buf;
                            goto safe_exit;
                        }
                        break;
                    }
                    for (k=i; (buf[k]!='\n') & (buf[k]!=EOS) & (buf[k]!=' ');
                         k++);
                    if (k >= sizeof(buf)-1) {
                        sprintf(buf, "load_auxiliary_from_file: %s"
                                     " line too long", fname);
                        *outstr = buf;
                        goto safe_exit;
                    }
                    if (j >= CONFIG_MAX_AUXILIARY) {
                        sprintf(buf, "load_auxiliary_from_file: number of "
                                     "entries > CONFIG_MAX_AUXILIARY=%d", 
                                     CONFIG_MAX_AUXILIARY);
                        *outstr = buf;
                        goto safe_exit;
                    }
                    if (buf[k]==EOS) {
                        CONFIG_auxiliary[j++][m] = atof(&buf[i]);
                        i = k;
                    }
                    else {
                        buf[k] = EOS;
                        CONFIG_auxiliary[j++][m] = atof(&buf[i]);
                        i = k+1;
                    }
                }
            }
            else if (m < np) {
                sprintf(buf, "** premature ending of %s **", fname);
                *outstr = buf;
                goto safe_exit;
            }
            else break;
    }
    else {
        sprintf(buf, "** %s: **\n** This file is unreadable! **", fname);
        *outstr = buf;
        return FALSE;
    }
    for (i=0; i<CONFIG_num_auxiliary; i++) {
        CalculateSimpleStatistics
            (np, CHARP(CONFIG_auxiliary[i]), sizeof(double),
             IOVAL_DOUBLE, &ss);
        n[iw].auxiliary_threshold[i][0] = ss.min;
        n[iw].auxiliary_threshold[i][1] = ss.max;
    }
    n[iw].color_mode = COLOR_MODE_AUXILIARY;
    if (OUW(n[iw].auxiliary_idx,CONFIG_num_auxiliary))
        n[iw].auxiliary_idx=0;
    return (color_encode_auxiliary(iw));
  safe_exit:
    if (AX_display[iw])
        printf ("** all auxiliary properties freed **\n");
    Config_free_auxiliary();
    return FALSE;
}


static bool proc_look_at_the_anchor(int iw, char *instr, char **outstr)
{
    *outstr = NULL;
    if (instr && sscanf(instr, " %s", buf) == 1) {
        if (strcmp(buf, "upright") == 0) {
            M3IDENTITY(AX_3D[iw].V);
        }
        else if (strcmp(buf, "reset") == 0) {
            n[iw].anchor = -1;
            V3EQV(cm, n[iw].hook);
        }
    }
    return look_at_the_anchor(iw);
}


static bool proc_observer_goto(int iw, char *instr, char **outstr)
{
    double s[3];

    V3mM3 (AX_3D[iw].x, HI, s);
    *outstr = cui_show_syntax;
    if (!instr)
        return FALSE;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        sprintf(buf, "%g %g %g", s[0],s[1],s[2]);
        instr = cui_stripspace(
            gui_readline_gets(iw, "\nObserver goto", buf));
    }

    sscanf(instr, "%lf %lf %lf", s, s+1, s+2);
    V3mM3 (s, H, AX_3D[iw].x);

    *outstr = NULL;
    return TRUE;
}

static bool proc_xtal_origin_goto(int iw, char *instr, char **outstr)
{
    double old_s[3], s[3];

    if (!n[iw].xtal_mode) {
        *outstr = "Crystal translation is only available under Xtal mode";
        return FALSE;
    }

    V3mM3(AX_3D[iw].x, HI, old_s);
    *outstr = cui_show_syntax;
    if (!instr || !*instr)
        return FALSE;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        sprintf(buf, "%g %g %g", old_s[0],old_s[1],old_s[2]);
        instr = cui_stripspace(
            gui_readline_gets(iw, "\nCrystal origin s0,s1,s2", buf));
    }

    *outstr = NULL;
    sscanf(instr, "%lf %lf %lf", s, s+1, s+2);
    V3TRIM (s,s);
    if (V3EQ(old_s,s)) return FALSE;
    V3mM3 (s, H, n[iw].xtal_origin);
    atom_xtal_origin (n[iw].xtal_origin);
    if (n[iw].bond_mode) bond_xtal_origin_update(iw);
    else n[iw].bond_xtal_origin_need_update = TRUE;
    return TRUE;
}

static bool proc_find_atom(int iw, char *instr, char **outstr)
{
    static int last_atom = 0;
    int i;
    char question[MAX_FILENAME_SIZE];

#if USE_P3D
    if (p3dp_enabled) {
        *outstr = "not supported in parallel mode";
        return FALSE;
    }
#endif
    *outstr = cui_show_syntax;
    if (!instr || !*instr)
        return FALSE;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        sprintf(question, "\nFind atom [0-%d]", np-1);
        sprintf(buf, "%d", last_atom);
        instr = cui_stripspace(gui_readline_gets(iw, question, buf));
    }

    sscanf (instr, "%d", &i);
    if ((i < 0) || (i >= np))
        *outstr = "find_atom: illegal index";
    else {
        n[iw].anchor = i;
        print_atom(iw,i);
        last_atom = i;
        *outstr = NULL;
    }
    return FALSE;
}


static bool proc_resize(int iw, char *instr, char **outstr)
{
    int new_width = AX_size[iw].width, new_height = AX_size[iw].height;

    *outstr = cui_show_syntax;
    if (!instr)
        return FALSE;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        sprintf(buf, "%d %d", new_width, new_height);
        instr = cui_stripspace(
            gui_readline_gets(iw, "\nInput new window width height", buf));
    }

    if (sscanf(instr, " %d %d", &new_width, &new_height) < 1)
        return FALSE;

    *outstr = NULL;
    if ( (new_width == 0) || (new_width > AX_MAXWIDTH) ) {
        if (AX_display[iw])
            printf("width = %d, set to AX_MAXWIDTH = %d\n",
                    new_width, AX_MAXWIDTH);
        new_width = AX_MAXWIDTH;
    }
    if ( (new_height==0) || (new_height > AX_MAXHEIGHT) ) {
        if (AX_display[iw])
            printf("height = %d, set to AX_MAXHEIGHT = %d\n",
                    new_height, AX_MAXHEIGHT);
        new_height = AX_MAXHEIGHT;
    }
    if ( (new_width < 0) || (new_height < 0) ||
         ( (new_width ==  AX_size[iw].width) &&
           (new_height == AX_size[iw].height) ) ) return FALSE;
    AX_size[iw].width = new_width;
    AX_size[iw].height = new_height;
    AX_resizewindow(iw, TRUE);
    cui_diligence = TRUE;
    return TRUE;
}


static bool proc_change_aux_colormap(int iw, char *instr, char **outstr)
{
    int i;

    *outstr = NULL;
    if (!instr || !*instr)
        goto error;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        if (AX_display[iw]) {
        printf("\nPlease choose a colormap index from:\n");
            for (i=0; i<AX_MAX_CMAP; i++)
                printf ("%2d: %s: %s\n",
                        i, AX_cmap_funs[i].name, AX_cmap_funs[i].description);
        }
        sprintf (buf, "%d", n[iw].auxiliary_cmap);
        instr = cui_stripspace(gui_readline_gets(iw, "Colormap index", buf));
    }

    if (strcmp(instr, "next") == 0)
        i = (n[iw].auxiliary_cmap + 1) % AX_MAX_CMAP;
    else if (strcmp(instr, "prev") == 0)
        i = (n[iw].auxiliary_cmap + AX_MAX_CMAP - 1) % AX_MAX_CMAP;
    else {
        i = -1;
        i = atoi(instr);
        if (OUW(i,AX_MAX_CMAP))
            goto error;
    }

    if ((i != n[iw].auxiliary_cmap) ||
            (n[iw].color_mode != COLOR_MODE_AUXILIARY)) {
        n[iw].auxiliary_cmap = i;
        n[iw].color_mode = COLOR_MODE_AUXILIARY;
        return color_encode_auxiliary(iw);
    }
    else
        return FALSE;

error:
    *outstr = cui_show_syntax;
    return FALSE;
}


static bool proc_print_atom_info(int iw, char *instr, char **outstr)
{
    if (instr) {
        *outstr = NULL;

        if (n[iw].atom_stack[0] < 0 || n[iw].atom_stack[1] < 0)
            return FALSE;
        else if (strstr(instr, "pair")) {
            if (AX_display[iw]) {
                printf ("\n");
#ifdef USE_P3D
                if (p3dp_enabled)
                    p3dp_print_atom_pair_info(iw,   n[iw].atom_stack[0],
                                                    n[iw].atom_stack[1],
                                                    p3dp_n[iw].rank_stack[0],
                                                    p3dp_n[iw].rank_stack[1]);
                else
#endif
                print_atom_pair_info(iw,            n[iw].atom_stack[0],
                                                    n[iw].atom_stack[1]);
            }
            return FALSE;
        }

        if (n[iw].atom_stack[2] < 0)
            return FALSE;
        else if (strstr(instr, "triplet")) {
            if (AX_display[iw]) {
                printf ("\n");
#ifdef USE_P3D
                if (p3dp_enabled)
                p3dp_print_atom_triplet_info(iw,    n[iw].atom_stack[0],
                                                    n[iw].atom_stack[1],
                                                    n[iw].atom_stack[2],
                                                    p3dp_n[iw].rank_stack[0],
                                                    p3dp_n[iw].rank_stack[1],
                                                    p3dp_n[iw].rank_stack[2]);
                else
#endif
                print_atom_triplet_info(iw,         n[iw].atom_stack[0],
                                                    n[iw].atom_stack[1],
                                                    n[iw].atom_stack[2]);
            }
            return FALSE;
        }

        if (n[iw].atom_stack[3] < 0)
            return FALSE;
        else if (strstr(instr, "quartet")) {
            if (AX_display[iw]) {
                printf ("\n");
#ifdef USE_P3D
                if (p3dp_enabled)
                p3dp_print_atom_quartet_info(iw,    n[iw].atom_stack[0],
                                                    n[iw].atom_stack[1],
                                                    n[iw].atom_stack[2],
                                                    n[iw].atom_stack[3],
                                                    p3dp_n[iw].rank_stack[0],
                                                    p3dp_n[iw].rank_stack[1],
                                                    p3dp_n[iw].rank_stack[2],
                                                    p3dp_n[iw].rank_stack[3]);
                else
#endif
                print_atom_quartet_info(iw,         n[iw].atom_stack[0],
                                                    n[iw].atom_stack[1],
                                                    n[iw].atom_stack[2],
                                                    n[iw].atom_stack[3]);
            }
            return FALSE;
        }
    }
    *outstr = cui_show_syntax;
    return FALSE;
}

static bool proc_save_atom_indices
(int iw, char *instr, char **outstr)
{
#if USE_P3D
    if (p3dp_enabled) {
        *outstr = "not supported in parallel mode";
        return FALSE;
    }
#endif
    if ( (n[iw].atom_stack[0] >= 0) &&
         (n[iw].atom_stack[1] >= 0) &&
         (n[iw].atom_stack[2] >= 0) ) {
        *outstr = NULL;
        printf ("\n");
        save_atoms_in_monoclinic_filter (iw);
    }
    else
        *outstr = cui_show_syntax;
    return FALSE;
}

static bool proc_change_central_symm_neighbormax
(int iw, char *instr, char **outstr)
{
    int new_central_symm_neighbormax;
    extern int central_symm_neighbormax;

    *outstr = cui_show_syntax;
    if (!instr)
        return FALSE;
    else if (strcmp(instr, CUI_ARG_REQUIRED) == 0) {
        sprintf(buf, "%d", central_symm_neighbormax);
        instr = cui_stripspace(
            gui_readline_gets(iw, "\ncentral_symm_neighbormax", buf));
    }

    sscanf(instr, "%d", &new_central_symm_neighbormax);
    if (new_central_symm_neighbormax <= 0)
        new_central_symm_neighbormax = coordination_crystal;
    /* it must be an even number */
    new_central_symm_neighbormax = new_central_symm_neighbormax / 2 * 2;
    if ( new_central_symm_neighbormax == central_symm_neighbormax )
        return (FALSE);
    central_symm_neighbormax = new_central_symm_neighbormax;
    evaluate_central_symm (geo+GEO_CENTRAL_SYMM);
    n[iw].color_mode = COLOR_MODE_AUXILIARY;
    n[iw].auxiliary_idx = CONFIG_MAX_AUXILIARY + GEO_CENTRAL_SYMM;
    return color_encode_auxiliary(iw);
}


static bool proc_timer(int iw, char *instr, char **outstr)
{
    double time;
    char title[CUI_LINEMAX] = "", *t = title;
    int reset = 0;

    buf[0] = 0;
    if (instr) {
        strncpy(title, instr, sizeof title);
        if (title[0] && sscanf(title, "%s", buf) && strcmp(buf, "reset") == 0) {
            reset = 1;
            t = cui_stripspace(strstr(title, "reset")+strlen("reset"));
        }
    }

    time = cui_wtime();
    if (*t)
        sprintf(buf, "%s: %.3f\n"CUI_PROTOCOL_OK, t, time - cui_time);
    else
        sprintf(buf, "%.3f\n"CUI_PROTOCOL_OK, time - cui_time);
    *outstr = buf;
    if (reset)
        cui_time = time;

    return FALSE;
}


static bool proc_isoatomic_reference_imprint(int iw, char *instr, char **outstr)
{
#ifdef USE_P3D
    if (p3dp_enabled) {
        *outstr = "not supported";
        return FALSE;
    }
#endif
    if  ( (!ComputeLeastSquareStrain) ||
          (strcmp(ref_fbasename, fbasename)!=0) ) {
        ComputeLeastSquareStrain = 1;
        strcpy (ref_fbasename, fbasename);
        IsoAtomicReferenceReImprint (Config_Aapp_to_Alib, ref);
        snprintf(buf, sizeof(buf), "\"%s\" is now reference for least-square "
               "strain calculation.\n", ref_fbasename);
    }
    else {
        ComputeLeastSquareStrain = 0;
        IsoAtomicReferenceFree (ref);
        LeastSquareStrain_Free();
        strncpy(buf, "Least-square strain calculation "
               "turned OFF.\n", sizeof(buf));
    }
    *outstr = strncat(buf, CUI_PROTOCOL_OK, sizeof(buf));
    return FALSE;
}


static bool proc_toggle_shell_viewer_mode(int iw, char *instr, char **outstr)
{
    n[iw].shell_viewer_mode  = !n[iw].shell_viewer_mode;
    *outstr = (n[iw].shell_viewer_mode) ?
            "Shell viewer auto-invoke is turned ON.\n"  CUI_PROTOCOL_OK :
            "Shell viewer auto-invoke is turned OFF.\n" CUI_PROTOCOL_OK ;
    return FALSE;
}


static bool proc_toggle_xtal_mode(int iw, char *instr, char **outstr)
{
    *outstr = NULL;
    n[iw].xtal_mode  = !n[iw].xtal_mode;
    return FALSE;
}


static bool proc_lazydraw(int iw, char *instr, char **outstr)
{
    *outstr = NULL;
    cui_diligence = FALSE;
    return FALSE;
}


static bool proc_change_shear_strain_subtract_mean
(int iw, char *instr, char **outstr)
{
    *outstr = NULL;
    return change_shear_strain_subtract_mean(iw);
}


static bool proc_open_port(int iw, char *instr, char **outstr)
{
    *outstr = NULL;
    if (listen_fd < 0 && (listen_fd=socket(PF_INET, SOCK_STREAM, 0)) != -1) {
        int i, port_from = CUI_PORT_FROM, port_to = CUI_PORT_TO;
        struct sockaddr_in inet_address;
        if (*instr == ':')
            instr++;
        if (sscanf(instr, " %d-%d", &port_from, &port_to) == 1)
            port_to = port_from;
        
        memset(&inet_address, 0, sizeof(inet_address));
        inet_address.sin_family = PF_INET;
        inet_address.sin_addr.s_addr = htonl(INADDR_ANY);
        
        for (i = port_from; i <= port_to; i++) {
            inet_address.sin_port = htons(i);
            if (bind(listen_fd, (struct sockaddr*)&inet_address,
                                            sizeof(inet_address)) >= 0) {
                if (listen(listen_fd, 5) < 0) {
                    close(listen_fd);
                    listen_fd = -1;
                }
                else {
                    FD_SET(listen_fd, &allset);
                    if (listen_fd >= maxfd_plus_1)
                        maxfd_plus_1 = listen_fd + 1;
                }
                sprintf(buf, "port: %d\n"CUI_PROTOCOL_OK, i);
                *outstr = buf;
                break;
            }
        }
    }
    return FALSE;
}


static bool proc_close_port(int iw, char *instr, char **outstr)
{
    if (listen_fd >= 0) {
        int j;
        close(listen_fd);
        listen_fd = -1;
        maxfd_plus_1 = 0;
        for (j = 0; j < CUI_N_FILE; j++) {
            if (cmdf[j].fp && cmdf[j].fd >= maxfd_plus_1)
                maxfd_plus_1 = cmdf[j].fd + 1;
        }
    }
    return FALSE;
}


static bool proc_syntax(int iw, char *instr, char **outstr)
{
    if (!instr || !*instr)
        *outstr = cui_show_syntax;
    else {
        struct aec *aecp;
        char name[CUI_LINEMAX] = "", *s;
        strncpy(buf, instr, sizeof(buf));
        s = cui_stripspace(buf);
        sscanf(s, " %[^ ]", name);
        if ((aecp = aecp_byname(name))) {
            if (aecp->syntax)
                sprintf(buf, "%s%s %s", cui_show_syntax, name, aecp->syntax);
            else
                sprintf(buf, "%s%s",    cui_show_syntax, name);
        }
        else {
            sprintf(buf, "%sunkown command \"%s\"", cui_show_syntax, name);
        }
        *outstr = buf;
    }
    return FALSE;
}

/*
static bool proc_subcommand(int iw, char *instr, char **outstr);
*/

#define AEC_NP(name) #name, proc_##name
#define AEC_EP(name, ext) #name"_"#ext, proc_##name
static struct aec atomeye_commands[] = {
    {AEC_NP(nop), NULL, NULL},
    {AEC_NP(next), NULL, NULL},
    {AEC_NP(close), NULL, NULL}, {"close_window", proc_close, NULL, NULL},
    {CUI_PROTOCOL_QUIT, proc_quit, NULL, NULL},
    {AEC_NP(new), NULL, "[clone]"},
    {"clone", proc_new, "clone", NULL},
    {AEC_NP(resize), CUI_ARG_REQUIRED, "width [height]"},
    {AEC_NP(set), CUI_ARG_REQUIRED, "variable_name value"},
    {AEC_NP(save), CUI_ARG_REQUIRED, "file_name [key|both]"},
    {AEC_NP(redraw), NULL, NULL},
    {AEC_NP(lazydraw), NULL, NULL},
    {AEC_NP(key), CUI_ARG_REQUIRED, "key_string"},
    {AEC_NP(change_bgcolor), CUI_ARG_REQUIRED, "R G B"},
    {AEC_NP(toggle_parallel_projection), NULL, NULL},
    {AEC_NP(toggle_bond_mode), NULL, NULL},
    {AEC_NP(toggle_coordination_coloring), NULL, NULL},
    {AEC_NP(normal_coloring), CUI_ARG_REQUIRED, "[original]"},
    {"original_normal_coloring", proc_normal_coloring, "original", NULL},
    {AEC_NP(change_atom_color),     CUI_ARG_REQUIRED,   "i [r g b] [R]"},
    {AEC_NP(change_normal_color),   CUI_ARG_REQUIRED,   "Z|Symbol [r g b] [R]"},
    {AEC_NP(change_coordination_color), CUI_ARG_REQUIRED,"c r g b"},
    {AEC_NP(change_bond_color),     CUI_ARG_REQUIRED,   "i [r g b] [R]"},
    {AEC_NP(translate), CUI_ARG_REQUIRED, "axis delta"},
    {AEC_EP(translate, 0_inc), "0  1 delta", NULL},
    {AEC_EP(translate, 0_dec), "0 -1 delta", NULL},
    {AEC_EP(translate, 1_inc), "1  1 delta", NULL},
    {AEC_EP(translate, 1_dec), "1 -1 delta", NULL},
    {AEC_EP(translate, 2_inc), "2  1 delta", NULL},
    {AEC_EP(translate, 2_dec), "2 -1 delta", NULL},
    {AEC_NP(shift_xtal), CUI_ARG_REQUIRED, "axis delta"},
    {AEC_EP(shift_xtal, 0_inc), "0  1 delta", NULL},
    {AEC_EP(shift_xtal, 0_dec), "0 -1 delta", NULL},
    {AEC_EP(shift_xtal, 1_inc), "1  1 delta", NULL},
    {AEC_EP(shift_xtal, 1_dec), "1 -1 delta", NULL},
    {AEC_EP(shift_xtal, 2_inc), "2  1 delta", NULL},
    {AEC_EP(shift_xtal, 2_dec), "2 -1 delta", NULL},
    {AEC_NP(rotate), CUI_ARG_REQUIRED, "axis theta"},
    {AEC_EP(rotate, 0_inc), "0  1 delta", NULL},
    {AEC_EP(rotate, 0_dec), "0 -1 delta", NULL},
    {AEC_EP(rotate, 1_inc), "1  1 delta", NULL},
    {AEC_EP(rotate, 1_dec), "1 -1 delta", NULL},
    {AEC_EP(rotate, 2_inc), "2  1 delta", NULL},
    {AEC_EP(rotate, 2_dec), "2 -1 delta", NULL},
    {AEC_NP(advance), CUI_ARG_REQUIRED, "delta"},
    {AEC_EP(advance, inc), " 1 delta", NULL},
    {AEC_EP(advance, dec), "-1 delta", NULL},
    {AEC_NP(shift_cutting_plane), CUI_ARG_REQUIRED, NULL},
    {AEC_EP(shift_cutting_plane, inc), " .33333333 delta", NULL},
    {AEC_EP(shift_cutting_plane, dec), "-.33333333 delta", NULL},
    {AEC_NP(change_view_angle_amplification), CUI_ARG_REQUIRED, NULL},
    {AEC_EP(change_view_angle_amplification, inc), " 1 delta", NULL},
    {AEC_EP(change_view_angle_amplification, dec), "-1 delta", NULL},
    {AEC_NP(change_atom_r_ratio), CUI_ARG_REQUIRED, NULL},
    {AEC_EP(change_atom_r_ratio, inc), " 1 delta", NULL},
    {AEC_EP(change_atom_r_ratio, dec), "-1 delta", NULL},
    {AEC_NP(change_bond_radius), CUI_ARG_REQUIRED, NULL},
    {AEC_EP(change_bond_radius, inc), " 1 delta", NULL},
    {AEC_EP(change_bond_radius, dec), "-1 delta", NULL},
    {AEC_NP(change_aux_colormap), CUI_ARG_REQUIRED, "number|prev|next"},
    {AEC_EP(change_aux_colormap, prev), "prev", NULL},
    {AEC_EP(change_aux_colormap, next), "next", NULL},
    {AEC_NP(aux_property_coloring), CUI_ARG_REQUIRED, "number"},
    {AEC_EP(aux_property_coloring,  0), " 0", NULL},
    {AEC_EP(aux_property_coloring,  1), " 1", NULL},
    {AEC_EP(aux_property_coloring,  2), " 2", NULL},
    {AEC_EP(aux_property_coloring,  3), " 3", NULL},
    {AEC_EP(aux_property_coloring,  4), " 4", NULL},
    {AEC_EP(aux_property_coloring,  5), " 5", NULL},
    {AEC_EP(aux_property_coloring,  6), " 6", NULL},
    {AEC_EP(aux_property_coloring,  7), " 7", NULL},
    {AEC_EP(aux_property_coloring,  8), " 8", NULL},
    {AEC_EP(aux_property_coloring,  9), " 9", NULL},
    {AEC_EP(aux_property_coloring, 10), "10", NULL},
    {AEC_EP(aux_property_coloring, 11), "11", NULL},
    {AEC_EP(aux_property_coloring, 12), "12", NULL},
    {AEC_EP(aux_property_coloring, 13), "13", NULL},
    {AEC_EP(aux_property_coloring, 14), "14", NULL},
    {AEC_EP(aux_property_coloring, 15), "15", NULL},
    {AEC_EP(aux_property_coloring, 16), "16", NULL},
    {AEC_EP(aux_property_coloring, 17), "17", NULL},
    {AEC_EP(aux_property_coloring, 18), "18", NULL},
    {AEC_EP(aux_property_coloring, 19), "19", NULL},
    {AEC_EP(aux_property_coloring, 20), "20", NULL},
    {AEC_EP(aux_property_coloring, 21), "21", NULL},
    {AEC_EP(aux_property_coloring, 22), "22", NULL},
    {AEC_EP(aux_property_coloring, 23), "23", NULL},
    {AEC_EP(aux_property_coloring, 24), "24", NULL},
    {AEC_EP(aux_property_coloring, 25), "25", NULL},
    {AEC_EP(aux_property_coloring, 26), "26", NULL},
    {AEC_EP(aux_property_coloring, 27), "27", NULL},
    {AEC_EP(aux_property_coloring, 28), "28", NULL},
    {AEC_EP(aux_property_coloring, 29), "29", NULL},
    {AEC_EP(aux_property_coloring, 30), "30", NULL},
    {AEC_EP(aux_property_coloring, 31), "31", NULL},
    {"shear_strain_coloring",     proc_aux_property_coloring, "48", NULL},
    {AEC_NP(central_symmetry_coloring), NULL, NULL},
    {AEC_EP(aux_property_coloring,  a), "10", NULL}, /* obsolete */
    {AEC_EP(aux_property_coloring,  b), "11", NULL}, /* obsolete */
    {AEC_EP(aux_property_coloring,  c), "12", NULL}, /* obsolete */
    {AEC_EP(aux_property_coloring,  d), "13", NULL}, /* obsolete */
    {AEC_EP(aux_property_coloring,  e), "14", NULL}, /* obsolete */
    {AEC_EP(aux_property_coloring,  f), "15", NULL}, /* obsolete */
    {AEC_EP(aux_property_coloring,  g), "32", NULL}, /* obsolete */
    {AEC_EP(aux_property_coloring,  h), "33", NULL}, /* obsolete */
    {AEC_NP(change_aux_property_threshold), CUI_ARG_REQUIRED, NULL},
    {AEC_EP(change_aux_property_threshold, lower_inc), "0  1 delta", NULL},
    {AEC_EP(change_aux_property_threshold, lower_dec), "0 -1 delta", NULL},
    {AEC_EP(change_aux_property_threshold, upper_inc), "1  1 delta", NULL},
    {AEC_EP(change_aux_property_threshold, upper_dec), "1 -1 delta", NULL},
    {AEC_NP(reset_aux_property_thresholds), NULL, NULL},
    {AEC_NP(toggle_aux_property_thresholds_saturation), NULL, NULL},
    {AEC_NP(toggle_aux_property_thresholds_rigid), NULL, NULL},
  /*{AEC_NP(select_gear), CUI_ARG_REQUIRED, "number"},*/
    {AEC_EP(select_gear, 0), "0", NULL},
    {AEC_EP(select_gear, 1), "1", NULL},
    {AEC_EP(select_gear, 2), "2", NULL},
    {AEC_EP(select_gear, 3), "3", NULL},
    {AEC_EP(select_gear, 4), "4", NULL},
    {AEC_EP(select_gear, 5), "5", NULL},
    {AEC_EP(select_gear, 6), "6", NULL},
    {AEC_EP(select_gear, 7), "7", NULL},
    {AEC_EP(select_gear, 8), "8", NULL},
    {AEC_EP(select_gear, 9), "9", NULL},
    {AEC_NP(cutting_plane), CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, 0), " 0 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, 1), " 1 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, 2), " 2 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, 3), " 3 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, 4), " 4 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, 5), " 5 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, 6), " 6 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, 7), " 7 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, 8), " 8 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, 9), " 9 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, a), "10 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, b), "11 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, c), "12 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, d), "13 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, e), "14 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(cutting_plane, f), "15 " CUI_ARG_REQUIRED, NULL},
  /*{AEC_NP(shift_cutting_plane_to_anchor), CUI_ARG_REQUIRED},*/
    {AEC_EP(shift_cutting_plane_to_anchor, 0), " 0", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, 1), " 1", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, 2), " 2", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, 3), " 3", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, 4), " 4", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, 5), " 5", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, 6), " 6", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, 7), " 7", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, 8), " 8", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, 9), " 9", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, a), "10", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, b), "11", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, c), "12", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, d), "13", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, e), "14", NULL},
    {AEC_EP(shift_cutting_plane_to_anchor, f), "15", NULL},
  /*{AEC_NP(delete_cutting_plane), CUI_ARG_REQUIRED, NULL},*/
    {AEC_EP(delete_cutting_plane, 0), " 0 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, 1), " 1 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, 2), " 2 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, 3), " 3 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, 4), " 4 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, 5), " 5 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, 6), " 6 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, 7), " 7 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, 8), " 8 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, 9), " 9 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, a), "10 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, b), "11 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, c), "12 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, d), "13 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, e), "14 " CUI_ARG_REQUIRED, NULL},
    {AEC_EP(delete_cutting_plane, f), "15 " CUI_ARG_REQUIRED, NULL},
    {AEC_NP(flip_cutting_plane), NULL, NULL},
    {AEC_NP(change_cutting_plane_wireframe_mode), NULL, NULL},
    {AEC_NP(change_wireframe_mode), NULL, NULL},
    {AEC_NP(capture), CUI_ARG_REQUIRED, "image_format [file_name]"},
    {"capture_png", proc_capture, "png " CUI_ARG_REQUIRED, "[file_name]"},
    {"capture_jpg", proc_capture, "jpg " CUI_ARG_REQUIRED, "[file_name]"},
    {"capture_eps", proc_capture, "eps " CUI_ARG_REQUIRED, "[file_name]"},

    {AEC_NP(load_config),                CUI_ARG_REQUIRED, NULL},
    {"reload_config",       proc_load_config_advance,    "reload", NULL},
  /*{AEC_NP(load_config_advance),        CUI_ARG_REQUIRED, NULL},*/
    {"load_config_forward", proc_load_config_advance,   "forward", NULL},
    {"load_config_backward",proc_load_config_advance,   "backward", NULL},
    {"load_config_first",   proc_load_config_advance,   "first", NULL},
    {"load_config_last",    proc_load_config_advance,   "last", NULL},
    {AEC_NP(script_animate), CUI_ARG_REQUIRED, NULL},
    {AEC_NP(load_atom_color),CUI_ARG_REQUIRED, NULL},
    {AEC_NP(load_aux), CUI_ARG_REQUIRED, NULL},

    {AEC_NP(look_at_the_anchor), NULL, NULL},
    {"reset_anchor",        proc_look_at_the_anchor,    "reset", NULL},
    {"viewframe_upright",   proc_look_at_the_anchor,    "upright", NULL},
    {AEC_NP(observer_goto),     CUI_ARG_REQUIRED, "s0 s1 s2"},
    {AEC_NP(xtal_origin_goto),  CUI_ARG_REQUIRED, "s0 s1 s2"},
    {"xtal_origin_zero",    proc_xtal_origin_goto, "0 0 0", NULL},
    {"xtal_origin_half",    proc_xtal_origin_goto, "0.5 0.5 0.5", NULL},
    {AEC_NP(find_atom),     CUI_ARG_REQUIRED, "i"},

    {AEC_NP(rcut_patch), CUI_ARG_REQUIRED,  "start symbol symbol"
                                            "|finish|value [delta]"},
    {AEC_EP(rcut_patch, inc),       " 1 delta", NULL},
    {AEC_EP(rcut_patch, dec),       "-1 delta", NULL},
    {AEC_NP(start_rcut_patch),      CUI_ARG_REQUIRED, "symbol symbol"},
    {AEC_NP(toggle_rcut_patch_mode),CUI_ARG_REQUIRED, NULL},
    {AEC_NP(finish_rcut_patch),     NULL, NULL},

  /*{AEC_NP(print_atom_info),   CUI_ARG_REQUIRED,   "pair|triplet|quartet"},*/
    {AEC_EP(print_atom_info, pair),     "pair",     NULL},
    {AEC_EP(print_atom_info, triplet),  "triplet",  NULL},
    {AEC_EP(print_atom_info, quartet),  "quartet",  NULL},

    {AEC_NP(save_atom_indices), NULL, NULL},

    {AEC_NP(timer), CUI_ARG_REQUIRED, "[reset] title"},
    {"reset_timer", proc_timer, "reset " CUI_ARG_REQUIRED, "title"},

    {AEC_NP(toggle_shell_viewer_mode), NULL, NULL},
    {AEC_NP(toggle_xtal_mode), NULL, NULL},

    {AEC_NP(syntax), CUI_ARG_REQUIRED, "commandname"},

    {AEC_NP(change_shear_strain_subtract_mean), NULL, NULL},
    {AEC_NP(change_central_symm_neighbormax), CUI_ARG_REQUIRED, "num"},

    {AEC_NP(isoatomic_reference_imprint), NULL, NULL},

    {AEC_NP(scratch_coloring), NULL, NULL},
    {AEC_NP(reset_scratch), CUI_ARG_REQUIRED, "n1 n2 n3"},
    {AEC_NP(free_scratch), NULL, NULL},

    {AEC_NP(load_script),  CUI_ARG_REQUIRED, "scriptfile"},
    {"disconnect",  proc_internal_error, NULL, NULL},
    {"end",         proc_internal_error, NULL, NULL},
    {"pause",       proc_internal_error, NULL, NULL},

    {AEC_NP(open_port),  CUI_ARG_REQUIRED, "[port]"},
    {AEC_NP(close_port), NULL, NULL},

    {NULL, NULL, NULL, NULL},
};



static struct aec *aecp_byname(char *name)
{
    struct aec *aecp;
    if (name && *name)
        for (aecp = atomeye_commands; aecp->name; aecp++)
            if (strcmp(name, aecp->name) == 0)
                return aecp;
    return NULL;
}

#define DEFAULTKEY(gktp, str, entry) \
    do { struct aec *aecp;\
        if (gktp->str && *gktp->str) \
            if ((aecp = aecp_byname(gktp->str)))\
                gktp->entry = aecp;\
            else {\
                fprintf(stderr, "unknown name \"%s\"\n", gktp->str);\
                exit(1);\
            }\
        else\
            if ((aecp = aecp_byname("nop")))\
                gktp->entry = aecp;\
            else {\
                fprintf(stderr, "unknown name \"nop\"\n");\
                exit(1);\
            }\
    } while (0)

static void gui_defaultkey(void)
{
    struct gkt *gktp;
    for (gktp = gui_key_table; gktp->keysym != NoSymbol; gktp++) {
        DEFAULTKEY(gktp, n,  normal);
        DEFAULTKEY(gktp, s,  shift);
        DEFAULTKEY(gktp, c,  ctrl);
        DEFAULTKEY(gktp, cs, ctrl_shift);
        DEFAULTKEY(gktp, m,  meta);
        DEFAULTKEY(gktp, ms, meta_shift);
    }
}



bool cui_treatevent(int iw)
{
    bool result = FALSE;
    static int firsttime = 1;

    if (firsttime && IS_MANAGER) {
        firsttime = 0;
        if (frontend >= 0) {
            rl_initialize();
            rl_redisplay();
        }
        else if (cui_stdin)
            cui_send(CUI_PROTOCOL_OK, stdout);
    }

    if (quit) {
        char *outstr;
        proc_close(iw, NULL, &outstr);
    }
    else if (iw == cui_iw) {
        char *line = cmdline;
        int len, i = 0;
        struct timeval timeout = CUI_TIMEOUT;
        fd_set rfds = allset;
        if (IS_MANAGER) {
            for (i = 0; i < CUI_N_FILE; i++) {
                if (cmdf[i].fp && cmdf[i].pause) {
                    if (cmdf[i].pause > 0.0) {
                        if (cmdf[i].pause < cui_wtime()) {
                            cmdf[i].pause = 0.0;
                            continue;
                        }
                    }
                    FD_CLR(cmdf[i].fd, &rfds);
                }
            }

            select(maxfd_plus_1, &rfds, NULL, NULL, &timeout);

            len = 0;
            for (i = 0; i < CUI_N_FILE; i++) {
                if (cmdf[i].fp && FD_ISSET(cmdf[i].fd, &rfds)) {
                    cmdline[0] = 0;
                    if (cmdf[i].fd == frontend)
                        rl_callback_read_char();
                    else if (!fgets(cmdline, sizeof(cmdline), cmdf[i].fp)) {
                        cmdfclose(i);
                        cmdline[0] = 0;
                    }
                    if ((len = strlen(line = cui_stripspace(cmdline))))
                        break;
                }
            }
            if (len > 0 && cmdf[i].fd == frontend)
                fputc('\r', stdout);
        }
#ifdef USE_P3D
        if (p3dp_enabled) {
            if (quit)
                len = -1;
            p3d_bcast(p3dp_cell, &len, 1, MPI_INT, 0);
            if (len > 0)
                p3d_bcast(p3dp_cell, line, len+1, MPI_CHAR, 0);
            else if (len < 0)
                quit = 1;
        }
#endif
        if (len > 0) {
            struct aec *aecp;
            char *outstr = "parameter error";

            if (!line || !*line) return FALSE;

            outstr = NULL;
            buf[0] = 0;
            sscanf(line, "%s", buf);

            if (buf[0] == '#') {
            }
            else if (strcmp(buf,"end") == 0 || strcmp(buf,"disconnect") == 0) {
                if (IS_MANAGER)
                    cmdfclose(i);
#ifdef USE_P3D
                if (p3dp_enabled)
                    p3d_bcast(p3dp_cell, &quit, 1, MPI_INT, 0);
#endif
            }
            else if (strncmp(buf, "pause", strlen("pause")) == 0) {
                if (IS_MANAGER) {
                    double seconds = 0;
                    if (sscanf(line,"pause %lf",&seconds) == 1 && seconds >= 0)
                        cmdf[i].pause = cui_wtime() + seconds;
                    else
                        outstr = "parameter error";
                }
            }
            else if ((aecp = aecp_byname(buf))) {
                char args[CUI_LINEMAX] = "", *instr = aecp->instr;
                if (strcmp(aecp->name, "load_script") == 0) {
                    strncpy(fname, CUI_SCRIPT_DEFAULT, sizeof(fname));
                    sscanf(line + strlen(aecp->name), " %s", fname);
                    sprintf(args, "%s %d", fname, i);
                    instr = args;
                }
                else if (instr) {
                    if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
                        instr = cui_stripspace(line + strlen(aecp->name));
                    else if (strstr(instr, CUI_ARG_REQUIRED)) {
                        strncpy(args, instr, sizeof(args));
                        *strstr(args, CUI_ARG_REQUIRED) = 0;
                        instr = strncat(args,
                                cui_stripspace(line + strlen(aecp->name)),
                                sizeof(args));
                    }
                }
                result = (*aecp->proc)(cui_iw, instr, &outstr);
                if (outstr) {
                    if (!*outstr)
                        outstr = NULL;
                    else if (outstr == cui_show_syntax && aecp->syntax) {
                        sprintf(buf, "%s%s %s",
                                cui_show_syntax, aecp->name, aecp->syntax);
                        outstr = buf;
                    }
                }
            }
            else
                outstr = "unknown command";

            if (IS_MANAGER && cmdf[i].writable) {
                int j = cui_send_result(outstr, cmdf[i].fp);
                char *s;
                if (cmdf[i].fd == frontend) {
                    if (j == 0 && *(s = cui_stripspace(cmdline)))
                        add_history(cmdline);
                    rl_initialize();
                    rl_redisplay();
                }
            }
        }
        if (IS_MANAGER && listen_fd >= 0 && FD_ISSET(listen_fd, &rfds)) {
            struct sockaddr_in inet_address;
            int j, fd;
            socklen_t l;
            if ((fd=accept(listen_fd,(struct sockaddr*)&inet_address,&l)) == -1)
                fprintf(stderr, "inet: accept\n");
            else if ((j = cmdfdopen(fd, "r+", -1)) == -1) {
                strncpy(buf, "open failed\n", sizeof(buf));
                write(fd, buf, strlen(buf));
                close(fd);
            }
            else {
                cui_send(CUI_GREETINGS,     cmdf[j].fp);
                cui_send(CUI_PROTOCOL_OK,   cmdf[j].fp);
            }
        }
    }
    else {
        struct timeval timeout = CUI_TIMEOUT;
        select(0, NULL, NULL, NULL, &timeout);
    }
    return result;
}


bool gui_treatevent(int iw)
{
    int i=0;
    bool pointer_in_window;
    AXSize newsize;
    static int WindowUnMapped[AX_MAXWIN] = {0};
    pointer_in_window = AXQueryPointer(iw);
    switch (AX_event[iw].type)
    {
        case UnmapNotify:
            for (i=0; i<AX_MAXWIN; i++)
                if ((AX_cid[i]!=0) && (i!=iw) && (!WindowUnMapped[i])) break;
            if (i == AX_MAXWIN)
                XIconifyWindow(AX_display[iw], xterm_win, AX_screen[iw]);
            WindowUnMapped[iw] = 1;
            return (FALSE);
        case MapNotify:
            XMapWindow(AX_display[iw], xterm_win);
            AXRaiseWindow(iw);
            WindowUnMapped[iw] = 0;
            return (FALSE);
        case KeyPress:
	    if (!pointer_in_window)
                return (FALSE);
            else {
                struct gkt *gktp;
                if ((gktp = gktp_bykeysym(AXKeysym0(iw)))) {
		  
                    char *outstr = NULL;
                    struct aec *aecp;
                    if (AXCTRL(iw))
                        if (AXShft(iw)) aecp = gktp->ctrl_shift;
                        else            aecp = gktp->ctrl;
                    else if (AXLOCK(iw) || AXMETA(iw))
                        if (AXShft(iw)) aecp = gktp->meta_shift;
                        else            aecp = gktp->meta;
                    else
                        if (AXShft(iw)) aecp = gktp->shift;
                        else            aecp = gktp->normal;
                    if (aecp->proc && *aecp->proc) {
                        bool result = (*aecp->proc)(iw, aecp->instr, &outstr);
                        if (IS_MANAGER && outstr && *outstr && frontend >= 0)
                            cui_send_result(outstr, stdin);
                        return result;
                    }
                }
                return FALSE;
            }
        case ButtonPress:
#ifdef USE_P3D
            p3dp_bcast_int2(&AX_win_x[iw], &AX_win_y[iw]);
#endif
            if (AXPRESSBTN(iw) == 4) {
                if (AXSHFT(iw) && n[iw].xtal_mode)
                    if (AXCTRL(iw)) return(xtal_shift(iw, 2, -3*n[iw].delta));
                    else return(xtal_shift(iw, 2, -n[iw].delta));
                else
                    if (AXCTRL(iw)) return(foo_advance(iw, 5*n[iw].delta));
                    else return(foo_advance(iw, n[iw].delta));
            }
            else if (AXPRESSBTN(iw) == 5) {
                if (AXSHFT(iw) && n[iw].xtal_mode)
                    if (AXCTRL(iw)) return(xtal_shift(iw, 2, 3*n[iw].delta));
                    else return(xtal_shift(iw, 2, n[iw].delta));
                else
                    if (AXCTRL(iw)) return(foo_advance(iw, -5*n[iw].delta));
                    else return(foo_advance(iw, -n[iw].delta));
            }
            else if ( (AXPRESSBTN(iw) == 2) || (AXPRESSBTN(iw) == 3) ||
                      (AXBTNTIME(iw) - n[iw].last_button_press_time <
                       DEFAULT_DBLECLICK_IN_MS) || AXCTRL(iw) ||
                      AXLOCK(iw) || AXMETA(iw) ) {
                char instr[CUI_LINEMAX], *outstr;
#ifdef USE_P3D
                if (p3dp_enabled)
                    i = p3dp_AX_3D_Balls_Zgrab(
                            iw, B, AX_win_x[iw], AX_win_y[iw], &p3dp_rank_grab);
                else
#endif
                i = AX_3D_Balls_Zgrab(iw, B, AX_win_x[iw], AX_win_y[iw]);
                if (i >= 0) {
#ifdef ATOMEYE_LIB
		  (*atomeyelib_on_click_atom)(i);
#endif
#ifdef USE_P3D
                    if (p3dp_enabled)
                        p3dp_atom_stack_insert(iw, i, p3dp_rank_grab);
                    else
#endif
                    atom_stack_insert (iw, i);
                    if (AXCTRL(iw) && AXSHFT(iw)) {
                        if (n[iw].color_mode == COLOR_MODE_NORMAL) {
                            int z;
#ifdef USE_P3D
                            if (!p3dp_enabled ||
                                    p3dp_rank_grab == p3d_rank(p3dp_cell))
#endif
                            z = ct->Z[(int)tp[i]];
#ifdef USE_P3D
                            if (p3dp_enabled)
                                p3d_bcast(p3dp_cell, &z, 1, MPI_INT,
                                                                p3dp_rank_grab);
#endif

                            sprintf(instr, "%d %s", z,
                                    (AXPRESSBTN(iw)!=1)?"-1":CUI_ARG_REQUIRED);
                            return proc_change_normal_color(iw, instr, &outstr);
                        }
                        else if (n[iw].color_mode == COLOR_MODE_COORD) {
                            int coord;
#ifdef USE_P3D
                            if (!p3dp_enabled ||
                                    p3dp_rank_grab == p3d_rank(p3dp_cell))
#endif
                            coord = coordination[i];
#ifdef USE_P3D
                            if (p3dp_enabled)
                                p3d_bcast(p3dp_cell, &coord, 1, MPI_INT,
                                                                p3dp_rank_grab);
#endif
                            sprintf(instr, "%d %s", coord,
                                    (AXPRESSBTN(iw)!=1)?"-1":CUI_ARG_REQUIRED);
                            return proc_change_coordination_color(iw, instr,
                                                    &outstr);
                        }
                    }
                    else if ( AXLOCK(iw) || AXMETA(iw) ) {
                        sprintf(instr, "%d %s", i,
                                    (AXPRESSBTN(iw)!=1)?"-1":CUI_ARG_REQUIRED);
                        return proc_change_atom_color(iw, instr, &outstr);
                    }
                    else if ( (AXPRESSBTN(iw) != 1) || AXCTRL(iw) ) {
#ifdef USE_P3D
                      if (p3dp_enabled) {
                        if (p3d_rank(p3dp_cell) == p3dp_rank_grab)
                            V3EQV(B->BALL[i].x, n[iw].hook);
                        p3d_bcast(p3dp_cell,
                                    n[iw].hook, 3, MPI_DOUBLE, p3dp_rank_grab);
                      }
                      else
#endif
                        n[iw].anchor = i;
/*
#ifdef USE_P3D
                        if (p3dp_enabled)
                            p3dp_atom_stack_insert(iw, i, p3dp_rank_grab);
                        else
#endif
                        atom_stack_insert (iw, i);
*/
                    }
#ifdef USE_P3D
                    if (!p3dp_enabled || p3d_rank(p3dp_cell) == p3dp_rank_grab)
#endif
		      print_atom (iw,i);
                }
                else if (n[iw].bond_mode)
                {
#ifdef USE_P3D
                    if (p3dp_enabled)
                        i = p3dp_AX_3D_Cylinders_Zgrab(
                            iw, C, AX_win_x[iw],AX_win_y[iw], &p3dp_rank_grab);
                    else
#endif
                    i = AX_3D_Cylinders_Zgrab(iw, C, AX_win_x[iw],AX_win_y[iw]);
                    if (i >= 0)
                    {
                        if ( AXLOCK(iw) || AXMETA(iw) ) {
                            sprintf(instr, "%d %s", i,
                                    (AXPRESSBTN(iw)!=1)?"-1":CUI_ARG_REQUIRED);
                            return proc_change_bond_color(iw, instr, &outstr);
                        }
#ifdef USE_P3D
                        if (!p3dp_enabled||p3d_rank(p3dp_cell)==p3dp_rank_grab)
#endif
                        print_bond(iw,i);
                        if ( (AXPRESSBTN(iw) != 1) || AXCTRL(iw) )
                        {
                            n[iw].anchor = -1;
#ifdef USE_P3D
                            if (!p3dp_enabled ||
                                    p3d_rank(p3dp_cell) == p3dp_rank_grab)
#endif
                            hook_to_cylinder(i, n[iw].hook);
#ifdef USE_P3D
                            if (p3dp_enabled)
                                p3d_bcast(p3dp_cell,
                                    n[iw].hook, 3, MPI_DOUBLE, p3dp_rank_grab);
#endif
                        }
                    }
                }
                n[iw].lx = AX_win_x[iw];
                n[iw].ly = AX_win_y[iw];
                return (FALSE);
            }
            else
            {
                n[iw].lx = AX_win_x[iw];
                n[iw].ly = AX_win_y[iw];
                n[iw].last_button_press_time = AXBTNTIME(iw);
                return (FALSE);
            }
        case MotionNotify:
#ifdef USE_P3D
            p3dp_bcast_int2(&AX_win_x[iw], &AX_win_y[iw]);
#endif
            if (AXMOTIONBTN2(iw) || AXMOTIONBTN3(iw))
                return( pointer_advance(iw, AX_win_x[iw], AX_win_y[iw]) );
            else if (AXMOTIONBTN1(iw))
            {
                if (AXSHFT(iw) && n[iw].xtal_mode)
                    return( pointer_xtal_shift
                            (iw, AX_win_x[iw], AX_win_y[iw]) );
                else if (AXCTRL(iw) || AXSHFT(iw))
                    return( pointer_translate
                            (iw, AX_win_x[iw], AX_win_y[iw]) );
                else return( pointer_rotate(iw, AX_win_x[iw], AX_win_y[iw]) );
            }
        case Expose:
            AXGetGeometry (iw, newsize);
            if ( (newsize.width  != AX_size[iw].width) ||
                 (newsize.height != AX_size[iw].height) )
            {
                /* there is a 4-byte alignment requirement for shm pixmap */
                /* in 16-bit graphics, therefore width must be even. */
                AX_size[iw] = newsize;
                if (AX_size[iw].width & 1) AX_size[iw].width++;
                AX_resizewindow(iw,False);
                n[iw].mgs_radius = DEFAULT_MGS_RATIO *
                    sqrt((double)AX_size[iw].width*AX_size[iw].height/PI);
                return (TRUE);
            }
            return (TRUE);
        default: return (FALSE);
    }
} /* end treatevent() */


int cui_config_load_A3(Alib_Declare_Config)
{
    int i;
    double mul = 1.15;
    M3 h = {{28.83722875,0.,0.}, {0.,23.069783,0.}, {0.,0.,9.0782}};
    struct cfg {
        char sym[3];
        V3 s;
    } sample[] = {
        {"Au",{0.1000,0.3250,0.3631}}, {"Au",{0.1000,0.4500,0.3631}},
        {"Au",{0.1000,0.5750,0.3631}}, {"Au",{0.1000,0.7000,0.3631}},
        {"Au",{0.1000,0.8250,0.3631}}, {"Au",{0.2000,0.2000,0.3631}},
        {"Au",{0.2000,0.5750,0.3631}}, {"Au",{0.3000,0.2000,0.3631}},
        {"Au",{0.3000,0.5750,0.3631}}, {"Au",{0.4000,0.3250,0.3631}},
        {"Au",{0.4000,0.4500,0.3631}}, {"Au",{0.4000,0.5750,0.3631}},
        {"Au",{0.4000,0.7000,0.3631}}, {"Au",{0.4000,0.8250,0.3631}},
        {"Au",{0.6000,0.2000,0.3631}}, {"Au",{0.6000,0.7000,0.3631}},
        {"Au",{0.7000,0.2000,0.3631}}, {"Au",{0.7000,0.4500,0.3631}},
        {"Au",{0.7000,0.8250,0.3631}}, {"Au",{0.8000,0.2000,0.3631}},
        {"Au",{0.8000,0.4500,0.3631}}, {"Au",{0.8000,0.8250,0.3631}},
        {"Au",{0.9000,0.2000,0.3631}}, {"Au",{0.9000,0.3250,0.3631}},
        {"Au",{0.9000,0.5750,0.3631}}, {"Au",{0.9000,0.7000,0.3631}},
                                       
        {"Ag",{0.0500,0.1375,0.5877}}, {"Ag",{0.0500,0.2625,0.5877}},
        {"Ag",{0.0500,0.3875,0.5877}}, {"Ag",{0.0500,0.5125,0.5877}},
        {"Ag",{0.0500,0.6375,0.5877}}, {"Ag",{0.0500,0.7625,0.5877}},
        {"Ag",{0.0500,0.8875,0.5877}}, {"Ag",{0.1500,0.1375,0.5877}},
        {"Ag",{0.1500,0.2625,0.5877}}, {"Ag",{0.1500,0.3875,0.5877}},
        {"Ag",{0.1500,0.5125,0.5877}}, {"Ag",{0.1500,0.6375,0.5877}},
        {"Ag",{0.1500,0.7625,0.5877}}, {"Ag",{0.1500,0.8875,0.5877}},
        {"Ag",{0.2500,0.1375,0.5877}}, {"Ag",{0.2500,0.2625,0.5877}},
        {"Ag",{0.2500,0.3875,0.5877}}, {"Ag",{0.2500,0.5125,0.5877}},
        {"Ag",{0.2500,0.6375,0.5877}}, {"Ag",{0.2500,0.7625,0.5877}},
        {"Ag",{0.2500,0.8875,0.5877}}, {"Ag",{0.3500,0.1375,0.5877}},
        {"Ag",{0.3500,0.2625,0.5877}}, {"Ag",{0.3500,0.3875,0.5877}},
        {"Ag",{0.3500,0.5125,0.5877}}, {"Ag",{0.3500,0.6375,0.5877}},
        {"Ag",{0.3500,0.7625,0.5877}}, {"Ag",{0.3500,0.8875,0.5877}},
        {"Ag",{0.4500,0.1375,0.5877}}, {"Ag",{0.4500,0.2625,0.5877}},
        {"Ag",{0.4500,0.3875,0.5877}}, {"Ag",{0.4500,0.5125,0.5877}},
        {"Ag",{0.4500,0.6375,0.5877}}, {"Ag",{0.4500,0.7625,0.5877}},
        {"Ag",{0.4500,0.8875,0.5877}}, {"Ag",{0.5500,0.1375,0.5877}},
        {"Ag",{0.5500,0.2625,0.5877}}, {"Ag",{0.5500,0.3875,0.5877}},
        {"Ag",{0.5500,0.5125,0.5877}}, {"Ag",{0.5500,0.6375,0.5877}},
        {"Ag",{0.5500,0.7625,0.5877}}, {"Ag",{0.5500,0.8875,0.5877}},
        {"Ag",{0.6500,0.1375,0.5877}}, {"Ag",{0.6500,0.2625,0.5877}},
        {"Ag",{0.6500,0.3875,0.5877}}, {"Ag",{0.6500,0.5125,0.5877}},
        {"Ag",{0.6500,0.6375,0.5877}}, {"Ag",{0.6500,0.7625,0.5877}},
        {"Ag",{0.6500,0.8875,0.5877}}, {"Ag",{0.7500,0.1375,0.5877}},
        {"Ag",{0.7500,0.2625,0.5877}}, {"Ag",{0.7500,0.3875,0.5877}},
        {"Ag",{0.7500,0.5125,0.5877}}, {"Ag",{0.7500,0.6375,0.5877}},
        {"Ag",{0.7500,0.7625,0.5877}}, {"Ag",{0.7500,0.8875,0.5877}},
        {"Ag",{0.8500,0.1375,0.5877}}, {"Ag",{0.8500,0.2625,0.5877}},
        {"Ag",{0.8500,0.3875,0.5877}}, {"Ag",{0.8500,0.5125,0.5877}},
        {"Ag",{0.8500,0.6375,0.5877}}, {"Ag",{0.8500,0.7625,0.5877}},
        {"Ag",{0.8500,0.8875,0.5877}}, {"Ag",{0.9500,0.1375,0.5877}},
        {"Ag",{0.9500,0.2625,0.5877}}, {"Ag",{0.9500,0.3875,0.5877}},
        {"Ag",{0.9500,0.5125,0.5877}}, {"Ag",{0.9500,0.6375,0.5877}},
        {"Ag",{0.9500,0.7625,0.5877}}, {"Ag",{0.9500,0.8875,0.5877}},

        {"",{.0,.0,.0}},
    };

    M3MULTIPLY(mul, h, H);

    i = 0;
    while (*sample[i].sym)
        i++;
    
    *np = i;
    Config_realloc (Config_Alib_to_Alib);

    if (CONFIG_num_auxiliary) {
        for (i = 1; i < CONFIG_num_auxiliary; i++) {
            Free (CONFIG_auxiliary[i]);
            CONFIG_auxiliary_name[i][0] = EOS;
            CONFIG_auxiliary_unit[i][0] = EOS;
        }
    }
    CONFIG_num_auxiliary = 1;
    REALLOC (Config_load, CONFIG_auxiliary[0], *np, double);
    strncpy(CONFIG_auxiliary_name[0],"linear",sizeof(CONFIG_auxiliary_name[0]));
    CONFIG_auxiliary_unit[0][0] = EOS;

    i = 0;
    while (*sample[i].sym) {
        (*mass)[i] = 1.0;
        strncpy(SYMBOL(i), sample[i].sym, 3);
        (*s )[i * DIMENSION    ] = sample[i].s[0];
        (*s )[i * DIMENSION + 1] = sample[i].s[1];
        (*s )[i * DIMENSION + 2] = sample[i].s[2];
        CONFIG_auxiliary[0][i]   = sample[i].s[0]-sample[i].s[1];
        i++;
    }
    return CONFIG_CFG_LOADED;
}



static char *command_generator(const char *text, int state)
{
    static int list_index, len;
    char *name;

    if (!state) {
        list_index = 0;
        len = strlen(text);
    }

    while ((name = atomeye_commands[list_index].name)) {
        list_index++;
        if (strncmp(name, text, len) == 0) {
            char *r = malloc(strlen(name)+1);
            if (r)
                strcpy(r, name);
            return (r);
        }
    }

    return ((char *)NULL);
}

static char **atomeye_completion(const char *text, int start, int end)
{
    char **matches;
    extern char **rl_completion_matches();

    matches = (char **)NULL;

    if (start == 0)
        matches = rl_completion_matches (text, command_generator);

    return (matches);
}


static void cui_initialize_readline(void)
{
    rl_readline_name = "AtomEye";
    rl_attempted_completion_function = atomeye_completion;
}


static FILE *terminal_fp, *terminal_out = NULL;
static char *terminal_prompt = "";

static void terminal_exit(int i)
{
    if (terminal_out)
        fprintf(terminal_out, "\n");
    if (*terminal_prompt)
        rl_callback_handler_remove();
    exit(i);
}


static void terminal_treatline(char *line)
{
    if (line) {
        if (*line) {
            char *s = cui_stripspace(line);
            int status = -1;

            if (*s) {
                cui_send(s, terminal_fp);
                status = cui_recv_print(terminal_fp, terminal_out);
            }

            if (*terminal_prompt) {
                if (status < 0)
                    terminal_exit(0);
                else if (status == 0 && *line)
                    add_history(line);
                free(line);
            }
        }
    }
    else {
        cui_send("disconnect", terminal_fp);
        terminal_exit(0);
    }
}


static int cui_terminal(void)
{
    int fd, port = CUI_PORT_FROM;
    char host[CUI_LINEMAX] = "localhost";
    char buf[CUI_LINEMAX];
    struct sockaddr_in addr;
    struct hostent *hp;

    if (sscanf(cui_hostport, " %[^:]:%d", host, &port) == 1) {
        int tmp = 0;
        if (sscanf(host, " %d", &tmp) == 1 && strchr(host, '.') == NULL) {
            strncpy(host, "localhost", sizeof(host));
            port = tmp;
        }
    }
    if (port < CUI_PORT_FROM || port > CUI_PORT_TO) {
        fprintf(stderr, "error: (port < %d || prot > %d)\n",
                                CUI_PORT_FROM, CUI_PORT_TO);
        return 1;
    }

    if ((fd = socket(PF_INET, SOCK_STREAM, 0)) < 0)
        return 1;

    memset(&addr, 0, sizeof(addr));
    addr.sin_family = PF_INET;
    if ((hp = gethostbyname(host)))
        memcpy(&addr.sin_addr, hp->h_addr, hp->h_length);
    else if ((addr.sin_addr.s_addr = inet_addr(host)) == -1)
        return 1;
    addr.sin_port = htons(port);
    if (connect(fd, (struct sockaddr*)&addr, sizeof(addr)) != 0)
        return 1;

    if ((terminal_fp = fdopen(fd, "r+"))) {
        int maxfd_plus_1;
        fd_set allset;

        if (isatty(STDIN_FILENO)) {
            terminal_prompt = CUI_PROMPT;
            terminal_out = stderr;
            setbuf(terminal_fp, NULL);
        }

        cui_recv_print(terminal_fp, terminal_out);

        FD_ZERO(&allset);
        FD_SET(STDIN_FILENO, &allset);
        FD_SET(fd, &allset);
        maxfd_plus_1 = (fd > STDIN_FILENO) ? fd + 1 : STDIN_FILENO + 1;

        if (*terminal_prompt) {
            cui_initialize_readline();
            rl_callback_handler_install(terminal_prompt, terminal_treatline);
        }

        while (1) {
            fd_set rfds = allset;
            if (select(maxfd_plus_1, &rfds, NULL, NULL, NULL) == -1)
                terminal_exit(0);
            if (FD_ISSET(fd, &rfds) && fgets(buf, sizeof(buf), terminal_fp)) {
                char *s = cui_stripspace(buf);
                if (strcmp(s, CUI_PROTOCOL_QUIT) == 0) {
                    terminal_exit(0);
                }
                else {
                    cui_send(CUI_PROTOCOL_OK, terminal_fp); /* NG? */
                }
            }
            if (FD_ISSET(STDIN_FILENO, &rfds)) {
                if (*terminal_prompt)
                    rl_callback_read_char();
                else {
                    if (fgets(buf, sizeof(buf), stdin))
                        terminal_treatline(buf);
                    else
                        terminal_exit(0);
                }
            }
        }
    }

    return 0;
}

struct co {
    char *name;
    int *pointer, value;
    char **strp;
};

static char *cui_xterm;

static struct co cui_options[] = {
    {"v2",              &cui_enabled,   0,  NULL},
    {"nowindow",        &cui_enabled,  -1,  NULL},
    {"nostdin",         &cui_stdin,     0,  NULL},
    {"listen",          &cui_listen,    1,  NULL},
    {"listen=",         &cui_listen,    1,  &cui_hostport},
    {"c",               &cui_connect,   1,  NULL},
    {"connect",         &cui_connect,   1,  NULL},
    {"connect=",        &cui_connect,   1,  &cui_hostport},
    {"f=",              NULL,           0,  &cui_scrfname},
    {NULL, NULL, 0}
};

#ifdef ATOMEYE_LIB
int cui_init(int *argc, char ***argv,  void (*on_click)(int atom), void (*on_close)(), void (*on_advance)(char *instr))
#else
int cui_init(int *argc, char ***argv)
#endif
{
    int i, retval = 0;
    char *str, *display_str = NULL, *TTYname = "/dev/null";
    extern char config_fname[MAX_FILENAME_SIZE];
    extern char xterm_identifier[XTERM_IDENTIFIER_SIZE];

#ifdef ATOMEYE_LIB
    atomeyelib_on_click_atom = on_click;
    atomeyelib_on_close = on_close;
    atomeyelib_on_advance = on_advance;
#endif

    cui_time = cui_wtime();
    cui_xterm = CUI_XTERM_DEFAULT;
    cui_xterm_win = ((str = getenv("WINDOWID"))) ? atol(str) : 0;
    cui_enabled = 1;
    cui_listen = 0;
    cui_connect = 0;
    cui_stdin = 1;
    cui_hostport = "";
    cui_scrfname = "";
    cui_geometry = "";
    config_fname[0] = 0;
    xterm_identifier[0] = 0;
    cui_diligence = TRUE;

#ifdef USE_P3D
    if ((i = p3dp_init(argc, argv)) < 0)
        return i;
#endif

    if (argv && *argv) {
	for (i = 1; i < *argc; i++) {
            if (*(str = (*argv)[i]) != '-') {
                if (!config_fname[0])
                    strcpy(config_fname, str);
                else {
                    retval = 1;
                    if (!xterm_identifier[0])
                        strcpy(xterm_identifier, str);
                    else
                        TTYname = str;
                }
            }
            else if (strcmp(str, "-nofep") == 0 && cui_stdin) {
                cui_stdin = -1;
            }
            else if (strcmp(str, "-display") == 0 && ++i < *argc) {
                display_str = (*argv)[i];
            }
            else if (strcmp(str, "-geometry") == 0 && ++i < *argc) {
                    cui_geometry = (*argv)[i];
            }
            else if (strcmp(str, "-xterm") == 0) {
                    cui_xterm_win = 0;
            }
            else if (strncmp(str, "-xterm=", strlen("-xterm=")) == 0) {
                    cui_xterm_win = 0;
                    cui_xterm = str + strlen("-xterm=");
            }
            else if (strncmp(str, "-xterm_win=", strlen("-xterm_win=")) == 0) {
                    cui_xterm_win = atol(str + strlen("-xterm_win="));
            }
            else {
                struct co *cop;
                for (cop = cui_options; cop->name; cop++) {
                    if (strcmp(str+1, cop->name) == 0) {
                        if (cop->pointer)
                            *cop->pointer = cop->value;
                        break;
                    }
                    else if (cop->strp) {
                        int len = strlen(cop->name);
                        if (strncmp(str+1, cop->name, len) == 0) {
                            if (cop->pointer)
                                *cop->pointer = cop->value;
                            *cop->strp = str + len + 1;
                            break;
                        }
                    }
                }
                if (!cop->name) {
                    fprintf(stderr, "unknown option \"%s\"\n", str);
                    return -1;
                }
            }
        }
    }

    if (cui_connect) {
        if (cui_listen)
            return -1; /* error */
#ifdef USE_P3D
        if (p3dp_enabled)
            return -2; /* error */
#endif
        cui_enabled = 1;
        if (cui_terminal()) {
            fprintf(stderr, "%s: failed to connect.\n", (*argv)[0]);
            exit(1);
        } 
        exit(0);
    }

    if ((cui_listen || *cui_scrfname) && !cui_enabled)
        cui_enabled = 1;

    if (!config_fname[0]) {
        if (cui_enabled)
            strcpy(config_fname, "/dev/null");
        else
            return -1; /* error */
    }

    if (!cui_enabled || display_str)
        cui_xterm_win = 0;

    if (display_str) {
        char *display_env =
                    (char*)malloc(strlen("DISPLAY=")+strlen(display_str)+1);
        if (display_env) {
            strcpy(display_env, "DISPLAY=");
            strcat(display_env, display_str);
            putenv(display_env);
        }
    }

#ifdef ATOMEYE_LIB
    /* Pretend we've got an xterm even if we haven't */
    cui_xterm_win = -1;
#endif
    if (!retval) {
        char **p;

        cui_argc = 0;

#ifdef USE_P3D
        if (p3dp_enabled) {
            if (!p3dp_mpiexec_argc) goto renderer;
        }
        else
#endif
	  if (cui_enabled < 0 || (cui_enabled && cui_xterm_win))
            goto renderer;

#ifdef USE_P3D
        if (cui_enabled >= 0)
#endif
        if (!cui_xterm_win) {
            int len;
            strncpy(buf, cui_xterm, sizeof(buf));
            len = strlen(buf);
            for (i = 0; i < len; i++)
                if (isspace(buf[i]))
                    buf[i] = 0;
            i = 0;
            while (i < len) {
                while ((!buf[i] || isspace(buf[i])) && i < len)
                    i++;
                if (i < len) {
                    cui_argv[cui_argc++] = buf + i;
                    while (buf[i] && !isspace(buf[i]) && i < len)
                        i++;
                }
            }
            cui_argv[cui_argc++] = "-xrm";
            cui_argv[cui_argc++] = xterm_identifier;
            cui_argv[cui_argc++] = "-e";
        }

#ifdef USE_P3D
        if (p3dp_enabled) {
            for (p = p3dp_mpiexec_argv; *p; p++)
                cui_argv[cui_argc++] = *p;
        }
#endif
        TimeRandomize();
        RandomBase64String(XTERM_IDENTIFIER_SIZE, xterm_identifier);
        cui_argv[cui_argc++] = (*argv)[0];
        cui_argv[cui_argc++] = config_fname;
        cui_argv[cui_argc++] = xterm_identifier;
        cui_argv[cui_argc++] = ttyname(STDERR_FILENO);

        if (!cui_enabled)
            cui_argv[cui_argc++] = "-v2";
        else {
            if (*cui_geometry) {
                cui_argv[cui_argc++] = "-geometry";
                cui_argv[cui_argc++] = cui_geometry;
            }

            if (cui_enabled < 0)
                cui_argv[cui_argc++] = "-nowindow";

            if (*cui_scrfname) {
                strncpy(fname, "-f=", sizeof(fname));
                strncat(fname, cui_scrfname, sizeof(fname));
                cui_argv[cui_argc++] = fname;
            }

            if (*cui_hostport) {
                char *chp =
                    (char*)malloc(strlen("-listen=")+strlen(cui_hostport)+1);
                if (chp) {
                    strcpy(chp, "-listen=");
                    strcat(chp, cui_hostport);
                    cui_argv[cui_argc++] = chp;
                }
            }
            else if (cui_listen)
                cui_argv[cui_argc++] = "-listen";
        }

#ifdef USE_P3D
        if (p3dp_enabled) {
            cui_argv[cui_argc++] = p3dp_decomp;
            if (!cui_xterm_win)
                cui_argv[cui_argc++] = "-xterm";
        }
#endif

        cui_argv[cui_argc] = NULL;
        execvp(cui_argv[0], cui_argv);
        exit(-2);
    }

renderer:
    if (!cui_enabled) {
        redirect_stderr_to(wOpen(TTYname));
    } else
    {
        gui_defaultkey();
        if (*cui_geometry) {
            if (sscanf(cui_geometry, "%dx%d", &cui_startup_width,
                                            &cui_startup_height) < 1)
                sscanf(cui_geometry, "x%d", &cui_startup_height);
        }
    
        if (IS_MANAGER) {
            char *home;
            int i;
    
            for (i = 0; i < CUI_N_FILE; i++) cmdf[i].fp = NULL;
            FD_ZERO(&allset);
    
            /* startup file */
            fname[0] = 0;
            if (*cui_scrfname)
                strncpy(fname, cui_scrfname, sizeof(fname));
            else if ((home = getenv("HOME"))) {
                strncpy(fname, home, sizeof(fname));
                strncat(fname, "/", sizeof(fname));
                strncat(fname, CUI_SCRIPT_STARTUP, sizeof(fname));
            }
            if (fname[0])
                cmdfdopen(open(fname, O_RDONLY), "r", -1);
    
            /* greetings */
            cui_send(CUI_GREETINGS, stderr);

#ifdef USE_P3D
            if (p3dp_enabled) {
                sprintf(buf, "[parallel: %d]", p3d_nprocs(p3dp_cell));
                cui_send(buf, stderr);
            }
#endif
            /* inet */
            if (cui_listen) {
                char *outstr;
                proc_open_port(0, cui_hostport, &outstr);
                cui_send(outstr, stderr);
            }

            /* frontend */
            if (cui_stdin) {
                if ((i = cmdfdopen(STDIN_FILENO,
                            isatty(STDIN_FILENO) ? "r+" : "r", -1)) >= 0) {
                    if (i != 0) {
                        cmdf[i].pause = -1.0;
                        cmdf[0].opener = i;
                    }
                    if (cui_stdin > 0 && isatty(STDIN_FILENO)) {
                        frontend = cmdf[i].fd;
                        cui_initialize_readline();
                        rl_callback_handler_install(CUI_PROMPT, stdin_gets);
                    }
                }
            }
        }
    
        cui_iw = 0;
        cui_title[0] = '*';
    }

    return 0;
}

/* API functions to be called from python */
#ifdef ATOMEYE_LIB

#include <xyz_netcdf.h>
#include <atomeyelib.h>

int atomeyelib_redraw(int iw) {
  /* send an expose event to the window so that it gets redrawn */

  XExposeEvent expose;
  expose.type = Expose;
  expose.display = AX_display[iw];
  expose.window = AX_win[iw];
  expose.x = 0;
  expose.y = 0;
  expose.width = AX_size[iw].width;
  expose.height = AX_size[iw].height;
  expose.count = 0;

  proc_redraw(iw, NULL, NULL);
  return XSendEvent(AX_display[iw], AX_win[iw], FALSE, 0, (XEvent *) &expose);
}

int atomeyelib_close(int iw) {

  XKeyEvent event;
  event.type = KeyPress;
  event.display = AX_display[iw];
  event.window = AX_win[iw];
  event.time = CurrentTime;
  event.x = 1;
  event.y = 1;
  event.same_screen = TRUE;

  event.type = KeyPress;
  event.keycode = XKeysymToKeycode(AX_display[iw], XK_q);
  event.state = 0;

  return XSendEvent(AX_display[iw], AX_win[iw], FALSE, 0, (XEvent *) &event);
}

int atomeyelib_run_command(int iw, char *line, char **outstr) {
  int result, i;
  struct aec *aecp;

  *outstr = NULL;
  buf[0] = 0;
  sscanf(line, "%s", buf);

  if ((aecp = aecp_byname(buf))) {
    char args[CUI_LINEMAX] = "", *instr = aecp->instr;
    if (strcmp(aecp->name, "load_script") == 0) {
      strncpy(fname, CUI_SCRIPT_DEFAULT, sizeof(fname));
      sscanf(line + strlen(aecp->name), " %s", fname);
      sprintf(args, "%s %d", fname, i);
      instr = args;
    }
    else if (instr) {
      if (strcmp(instr, CUI_ARG_REQUIRED) == 0)
	instr = cui_stripspace(line + strlen(aecp->name));
      else if (strstr(instr, CUI_ARG_REQUIRED)) {
	strncpy(args, instr, sizeof(args));
	*strstr(args, CUI_ARG_REQUIRED) = 0;
	instr = strncat(args,
			cui_stripspace(line + strlen(aecp->name)),
			sizeof(args));
      }
    }
    result = (*aecp->proc)(iw, instr, outstr);
    if (*outstr) {
      if (*outstr == cui_show_syntax && aecp->syntax) {
	sprintf(buf, "%s%s %s",
		cui_show_syntax, aecp->name, aecp->syntax);
	*outstr = buf;
      }
    }
  }
  else
    *outstr = "unknown command";

  if (result) atomeyelib_redraw(iw);
  return result;
}

/* Copy data from Atoms C structure in memory */
int atomeyelib_load_libatoms(int iw, Atoms *atoms, char *title, char **outstr) 
{
    int i, j, k, old_np;
    V3 hook_s, tmp, dx;
    char *old_symbol=NULL;
    bool incompatible_config;
/*     static int firsttime = 1; */
    
    *outstr = NULL;

    n[iw].suppress_printout = 1;

/*     if (!firsttime) { */
/*       firsttime = 0; */
/*       if (n[iw].anchor >= 0) { */
/*         /\* the new configuration may not even have the atom *\/ */
/*         V3EQV (B->BALL[n[iw].anchor].x, n[iw].hook); */
/*         n[iw].anchor = -1; */
/*       } */
/*       /\* hook_s[] is what is kept invariant *\/ */
/*       V3mM3 (n[iw].hook, HI, hook_s); */
/*     } */

    old_np = np;
    CLONE(symbol, SYMBOL_SIZE*np, char, old_symbol);
    Config_load_libatoms(atoms, NULL, Config_Aapp_to_Alib);

    for (k=0; k<CONFIG_num_auxiliary; k++)
        if (*blank_advance(CONFIG_auxiliary_name[k])==EOS)
            sprintf(CONFIG_auxiliary_name[k], "auxiliary%d", k);
#ifndef ATOMEYE_LIB
    rebind_CT (Config_Aapp_to_Alib, "", ct, &tp); cr();
#else
    rebind_ct (Config_Aapp_to_Alib, "", ct, &tp, NULL); 
#endif
    Neighborlist_Recreate_Form (Config_Aapp_to_Alib, ct, N);
    N->s_overflow_err_handler =
      NEIGHBORLIST_S_OVERFLOW_ERR_HANDLER_FOLD_INTO_PBC;
    N->small_cell_err_handler = NEIGHBORLIST_SMALL_CELL_ERR_HANDLER_MULTIPLY;
    for (i=0; i<ct->t; i++)
        for (j=i; j<ct->t; j++)
            for (k=0; k<rcut_patch_top; k++)
                if ( ( ( (rcut_patch[k].Zi == ct->Z[i]) &&
                         (rcut_patch[k].Zj == ct->Z[j]) ) ||
                       ( (rcut_patch[k].Zi == ct->Z[j]) &&
                         (rcut_patch[k].Zj == ct->Z[i]) ) ) )
                    NEIGHBOR_TABLE(N->rcut,ct,i,j) =
                        NEIGHBOR_TABLE(N->rcut,ct,j,i) = rcut_patch[k].rcut;
    Neighborlist_Recreate (Config_Aapp_to_Alib, NULL, ct, &tp, N);
    V3mM3 (hook_s, H, tmp);
    V3SUB (tmp, n[iw].hook, dx);
    V3EQV (tmp, n[iw].hook);
    V3AdD (dx, AX_3D[iw].x);
    M3InV (H, HI, volume);
    lengthscale = cbrt(volume);
    V3ASSIGN (0.5,0.5,0.5,tmp);
    V3mM3 (tmp, H, cm);
    geo_clear_has_evaluated_flags();
    evaluate_geo_measures(); Free(s1); Free(mass);

    if  ( ComputeLeastSquareStrain )
    {
        if (ConfigChecksum(Config_Aapp_to_Alib) != ref->checksum) {
            static char buf2[CUI_LINEMAX];
            snprintf(buf2, sizeof(buf2),
                    "This configuration is not isoatomic with the imprinted "
                    "reference\n%s. Least-square strain NOT calculated.\n",
                    ref_fbasename);
            *outstr = buf2;
            fprintf(stderr, "%s", buf2);
        }
        else LeastSquareStrain_Append();
    }

    incompatible_config = (np != old_np) ||
        memcmp(symbol, old_symbol, SYMBOL_SIZE*MIN(np,old_np));
    Free(old_symbol);

    if (incompatible_config)
        Config_to_3D_Balls (n[iw].atom_r_ratio);
    else for (i=0; i<np; i++) V3mM3 ( &(s[DIMENSION*i]), H, B->BALL[i].x );

    atom_xtal_origin (n[iw].xtal_origin);
    if (!n[iw].auxiliary_thresholds_rigid) {
        for (i=0; i<CONFIG_num_auxiliary; i++)
            reset_auxiliary_threshold(iw,i);
        for (i=0; i<MAX_GEO_MEASURES; i++)
            if (geolist[i].has_evaluated)
                reset_auxiliary_threshold(iw,CONFIG_MAX_AUXILIARY+i);
    }
    if (!temporary_disable_bond) Config_to_3D_Bonds (n[iw].bond_radius);

    strcpy(fbasename,title);

    if ((n[iw].xtal_mode) && (n[iw].color_mode == COLOR_MODE_COORD))
        assign_coordination_color(iw);
    else if (n[iw].color_mode == COLOR_MODE_AUXILIARY)
        color_encode_auxiliary(iw);
    else if (n[iw].color_mode == COLOR_MODE_SCRATCH)
        /*scratch_color (iw);*/
        proc_scratch_coloring(iw, "", outstr);
    else {
      strcpy (AX_title[iw],fbasename);
      AXSetName (iw);
      if (!temporary_disable_bond) {
	bond_xtal_origin_update (iw);
	bond_atom_color_update (iw);
      }
    }
    return TRUE;

error:
    *outstr = "parameter error";
    return FALSE;  
}


bool atomeyelib_help(int iw, char *instr, char **outstr)
{
    if (!instr || !*instr)
        *outstr = cui_show_syntax;
    else {
        struct aec *aecp;
        char name[CUI_LINEMAX] = "", *s;
        strncpy(buf, instr, sizeof(buf));
        s = cui_stripspace(buf);
        sscanf(s, " %[^ ]", name);
        if ((aecp = aecp_byname(name))) {
            if (aecp->syntax)
                sprintf(buf, "%s%s %s", cui_show_syntax, name, aecp->syntax);
            else
                sprintf(buf, "%s%s",    cui_show_syntax, name);
        }
        else {
            sprintf(buf, "%sunknown command \"%s\"", cui_show_syntax, name);
        }
        *outstr = buf;
    }
    return FALSE;
}


/* int atomeyelib_set_output(int on_off) */
/* { */
/*   static int stdout_fd = -1; */
/*   static fpos_t stdout_pos;  */

/*   if (on_off) { */
/*     if (stdout_fd != -1) { */
/*       fflush(stdout); */
/*       dup2(stdout_fd, fileno(stdout)); */
/*       close(stdout_fd); */
/*       clearerr(stdout); */
/*       fsetpos(stdout, &stdout_pos);  */
/*     } */
/*   } else { */
/*     fflush(stdout); */
/*     fgetpos(stdout, &stdout_pos); */
/*     stdout_fd = dup(fileno(stdout)); */
/*     freopen("/dev/null", "w", stdout); */
/*   } */
/*   return TRUE; */
/* } */

#endif
