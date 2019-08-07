#include "x2b_A1B2Z2_C1D4_deg3_exp0_v1x.h" 
 

////////////////////////////////////////////////////////////////////////////////

namespace {

//----------------------------------------------------------------------------//

struct variable {
    double v_exp0(const double& r0, const double& k,
                 const double * p1, const double * p2 );
                 
    double v_exp(const double& k,
                 const double * p1, const double * p2 );

    double v_coul0(const double& r0, const double& k,
                  const double * p1, const double * p2 );
                  
    double v_coul(const double& k,
                  const double * p1, const double * p2 );
                  
    void grads(const double& gg, double * grd1, double * grd2,
               const double * p1, const double * p2);

    double g[3]; // diff(value, p1 - p2)
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_exp0(const double& r0, const double& k,
                       const double * p1, const double * p2)
{
    g[0] = p1[0] - p2[0];
    g[1] = p1[1] - p2[1];
    g[2] = p1[2] - p2[2];

    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);

    const double exp1 = std::exp(k*(r0 - r));
    const double gg = - k*exp1/r;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return exp1;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_exp(const double& k,
                       const double * p1, const double * p2)
{
    g[0] = p1[0] - p2[0];
    g[1] = p1[1] - p2[1];
    g[2] = p1[2] - p2[2];

    const double r = std::sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);

    const double exp1 = std::exp(k*(- r));
    const double gg = - k*exp1/r;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return exp1;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //


double variable::v_coul(const double& k,
                        const double * p1, const double * p2)
{
    g[0] = p1[0] - p2[0];
    g[1] = p1[1] - p2[1];
    g[2] = p1[2] - p2[2];

    const double rsq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
    const double r = std::sqrt(rsq);

    const double exp1 = std::exp(k*(-r));
    const double rinv = 1.0/r;
    const double val = exp1*rinv;

    const double gg = - (k + rinv)*val*rinv;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return val;
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

double variable::v_coul0(const double& r0, const double& k,
                        const double * p1, const double * p2)
{
    g[0] = p1[0] - p2[0];
    g[1] = p1[1] - p2[1];
    g[2] = p1[2] - p2[2];

    const double rsq = g[0]*g[0] + g[1]*g[1] + g[2]*g[2];
    const double r = std::sqrt(rsq);

    const double exp1 = std::exp(k*(r0 - r));
    const double rinv = 1.0/r;
    const double val = exp1*rinv;

    const double gg = - (k + rinv)*val*rinv;

    g[0] *= gg;
    g[1] *= gg;
    g[2] *= gg;

    return val;
}

//----------------------------------------------------------------------------//

void variable::grads(const double& gg, double * grd1, double * grd2, 
                     const double * p1, const double * p2) {
    for (size_t i = 0; i < 3 ; i++) {
        double d = gg*g[i];
        grd1[i] += d;
        grd2[i] -= d;
    }
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

struct monomer {
    double oh1[3];
    double oh2[3];

    void setup(const double* ohh,
               const double& in_plane_g, const double& out_of_plane_g,
               double x1[3], double x2[3]);

    void grads(const double* g1, const double* g2,
               const double& in_plane_g, const double& out_of_plane_g,
               double* grd) const;
};

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void monomer::setup(const double* ohh,
                    const double& in_plane_g, const double& out_of_plane_g,
                    double* x1, double* x2)
{
    for (int i = 0; i < 3; ++i) {
        oh1[i] = ohh[i + 3] - ohh[i];
        oh2[i] = ohh[i + 6] - ohh[i];
    }

    const double v[3] = {
        oh1[1]*oh2[2] - oh1[2]*oh2[1],
        oh1[2]*oh2[0] - oh1[0]*oh2[2],
        oh1[0]*oh2[1] - oh1[1]*oh2[0]
    };

    for (int i = 0; i < 3; ++i) {
        const double in_plane = ohh[i] + 0.5*in_plane_g*(oh1[i] + oh2[i]);
        const double out_of_plane = out_of_plane_g*v[i];

        x1[i] = in_plane + out_of_plane;
        x2[i] = in_plane - out_of_plane;
    }
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void monomer::grads(const double* g1, const double* g2,
                    const double& in_plane_g, const double& out_of_plane_g,
                    double* grd) const
{
    const double gm[3] = {
        g1[0] - g2[0],
        g1[1] - g2[1],
        g1[2] - g2[2]
    };

    const double t1[3] = {
        oh2[1]*gm[2] - oh2[2]*gm[1],
        oh2[2]*gm[0] - oh2[0]*gm[2],
        oh2[0]*gm[1] - oh2[1]*gm[0]
    };

    const double t2[3] = {
        oh1[1]*gm[2] - oh1[2]*gm[1],
        oh1[2]*gm[0] - oh1[0]*gm[2],
        oh1[0]*gm[1] - oh1[1]*gm[0]
    };

    for (int i = 0; i < 3; ++i) {
        const double gsum = g1[i] + g2[i];
        const double in_plane = 0.5*in_plane_g*gsum;

        const double gh1 = in_plane + out_of_plane_g*t1[i];
        const double gh2 = in_plane - out_of_plane_g*t2[i];

        grd[i + 0] += gsum - (gh1 + gh2); // O
        grd[i + 3] += gh1; // H1
        grd[i + 6] += gh2; // H2
    }
}

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
} // namespace

////////////////////////////////////////////////////////////////////////////////

namespace x2b_A1B2Z2_C1D4_deg3_exp0 {

//----------------------------------------------------------------------------//

x2b_A1B2Z2_C1D4_v1x::x2b_A1B2Z2_C1D4_v1x(std::string mon1, std::string mon2) {

    // =====>> SECTION CONSTRUCTOR <<=====
    // =>> PASTE RIGHT BELOW THIS LINE <==
    if (mon1 == "h2o" && mon2  == "ch4") {
        coefficients = std::vector<double> {
 1.291778713538232e+01, // 0
 5.638813587262306e+00, // 1
 1.442044516821717e+01, // 2
-3.279979118833674e+00, // 3
 2.219620240170074e+01, // 4
 6.206879225499147e+01, // 5
-3.591293720043276e+00, // 6
 5.168778008335407e+00, // 7
 1.394942565561287e+01, // 8
 3.149119643990327e+01, // 9
 1.930614473711115e+00, // 10
-3.308932365870695e+01, // 11
 2.340868906806053e+01, // 12
-6.763082357113222e+00, // 13
 2.492078467807481e+01, // 14
 5.493500393059817e-01, // 15
 1.198252977377918e+00, // 16
-2.823285357000332e+00, // 17
-5.633212084401434e-01, // 18
 1.080218138907881e-03, // 19
 1.870150247883740e+01, // 20
-1.050885624613773e+00, // 21
-8.649502272156409e+00, // 22
 5.316773441989630e+01, // 23
 2.511297152422869e-01, // 24
-1.001361865533858e-01, // 25
-1.846799942408567e+00, // 26
 1.479923093221109e+01, // 27
 1.252578349323374e+00, // 28
-2.271997050492584e-02, // 29
 1.030042286432924e+00, // 30
 1.099805632702580e+00, // 31
-7.341159250021915e-03, // 32
 1.040472942005260e+01, // 33
-1.307031950240526e-02, // 34
-4.921826506263285e-01, // 35
-2.223853883662865e+01, // 36
-1.530117182646826e-01, // 37
-3.533914109453805e+01, // 38
-1.704726388100217e+00, // 39
 9.517889162415166e-02, // 40
-2.450679964098376e-04, // 41
 5.433336393434935e-02, // 42
 1.843482356991758e-04, // 43
 9.797751570177250e-02, // 44
-1.546239657722971e-01, // 45
 7.481695254030254e+00, // 46
 4.962360523652224e-02, // 47
 6.136526210780423e-02, // 48
-3.578253005147145e+01, // 49
-6.723388851822471e+00, // 50
-3.543659183091473e+01, // 51
 1.064914085134314e-02, // 52
-5.194978242724287e-01, // 53
 5.827691342237258e-01, // 54
 1.882514127978259e+00, // 55
 2.063270277669162e+01, // 56
-6.540486377222112e+00, // 57
 3.557217915083224e-01, // 58
-4.640661800153709e-02, // 59
 4.822553579430539e+01, // 60
 7.021961639195605e-02, // 61
-5.954823425138384e-01, // 62
-4.866233077516360e+00, // 63
 1.990944471095904e+01, // 64
 5.443638001964342e-01, // 65
-8.893953853485369e+00, // 66
 4.230105416619289e-02, // 67
-6.752751344836645e+00, // 68
-6.556934805686359e+01, // 69
 9.296681645149385e-01, // 70
-1.965135176501090e+01, // 71
 6.062244806477967e-01, // 72
-2.234469955600580e-02, // 73
 1.559486360842538e-04, // 74
-4.691008741219837e+01, // 75
-5.180554494897914e-04, // 76
 7.564325233450248e-04, // 77
-3.798345352555770e-01, // 78
-1.205483499750835e-07, // 79
-3.684254979261697e-02, // 80
 5.462771013168652e+00, // 81
 6.226084395356756e-01, // 82
-3.265526401592174e-02, // 83
-4.808124817647002e-01, // 84
 7.040008701739083e-03, // 85
-4.148751940678502e-01, // 86
-6.007497948306673e-01, // 87
 2.363606844868535e+00, // 88
 2.409487636881835e+00, // 89
 2.913524356019750e-03, // 90
-3.633571497281363e+00, // 91
-9.870244704809719e+00, // 92
 1.871295503802080e-01, // 93
 7.126451376161566e-05, // 94
 1.553908282480471e-02, // 95
 8.718039238523223e-02, // 96
-5.058881456134242e-04, // 97
-4.658172631851469e-03, // 98
-3.012872596271748e+01, // 99
 4.885356951104531e-01, // 100
 7.759194971161875e-05, // 101
-1.807662615964473e-05, // 102
 5.796024312257224e+01, // 103
 7.850981672069157e-05, // 104
-9.865099008381700e-03, // 105
-3.799682443323078e+00, // 106
-1.709931124176518e+00, // 107
 2.415023403652249e-05, // 108
-1.086014439277443e+00, // 109
-7.471183430642950e-01, // 110
-2.042478157372797e-05, // 111
 5.269322799086015e+00, // 112
-4.767157034447156e-03, // 113
-5.049226574532259e+00, // 114
 9.201563030024811e+00, // 115
 1.220010479031038e-01, // 116
-6.280037874938206e-01, // 117
 2.182797154012973e+00, // 118
-2.223466950604664e+00, // 119
-5.986593393883997e-01, // 120
 7.441993420380273e-03, // 121
 1.361586573471242e-01, // 122
 1.846984833453719e-03, // 123
 5.089726485771395e-03, // 124
-4.182446557247895e-02, // 125
-8.086938950149123e-06, // 126
-4.294773245534452e-04, // 127
 2.530731351351520e+00, // 128
 1.744105463363929e+00, // 129
 1.076005249785383e+02, // 130
 1.153747333346181e+00, // 131
 3.576166902449290e-03, // 132
 3.010751107832063e-03, // 133
 8.510783598880891e-06, // 134
-2.555637535831878e-05, // 135
-9.847440107684545e-06, // 136
-7.463705963727423e-02, // 137
-1.319971791092909e-02, // 138
-2.695456607292037e-03, // 139
-5.040895706626520e-01, // 140
 2.414780034912544e+00, // 141
 1.626987586013982e-01, // 142
-2.055743828718932e-01, // 143
 6.861637939779593e+00, // 144
 3.288384493965721e-06, // 145
-4.657552227383827e-01, // 146
-4.539259393733526e+00, // 147
-1.822483892911952e-02, // 148
 1.243475553646200e-02, // 149
 5.273634801336976e-03, // 150
 3.623503465438488e-06, // 151
 1.326359517269674e-02, // 152
 1.251983660306919e+00, // 153
 5.290578339660936e+00, // 154
 1.865865023279976e-01, // 155
-5.002045918166546e-02, // 156
-2.927090526363945e+01, // 157
 2.047692055936500e+01, // 158
 1.381153082533654e+00, // 159
 1.661123127894352e-01, // 160
 5.124058946436833e-02, // 161
 3.086502243985460e-02, // 162
 5.203664789838554e-03, // 163
 9.884650955481916e-03, // 164
-3.291368098459466e-01, // 165
 2.386810619293034e+01, // 166
 3.529448413704904e-03, // 167
-4.690260566020915e-05, // 168
 2.159807367618208e-02, // 169
 8.122358436900308e-01, // 170
 6.724286839367601e+00, // 171
 1.284560100139823e+00, // 172
-8.602787470651876e-01, // 173
 2.722003673107825e+00, // 174
 1.444884654336454e-01, // 175
-2.406408255547162e+01, // 176
 1.214484884624669e+00, // 177
-2.644477643352262e+00, // 178
-1.145727969854399e+01, // 179
-1.190044337744878e+00, // 180
-1.130257071907318e+00, // 181
-3.491103447574893e-03, // 182
-1.891676754639449e-06, // 183
 2.334690655180634e-02, // 184
 4.691096374960424e+00, // 185
 4.151956349366356e-03, // 186
-2.197615263229161e-01, // 187
-3.972394815080920e-05, // 188
-8.165374166529675e-04, // 189
 1.369504916505922e+00, // 190
-3.832810952556958e-01, // 191
 8.112628530389443e+00, // 192
-6.286825912092102e-02, // 193
-6.491772233115948e-04, // 194
 2.898555941913045e-01, // 195
-9.165108895935835e-01, // 196
 1.183577557387791e+01, // 197
-2.820746919129897e-03, // 198
 7.926437253468760e-02, // 199
 5.342396185655850e+00, // 200
-1.226777941534845e+00, // 201
 4.028813802028553e-01, // 202
 1.097396604723944e-01, // 203
 1.464182197848229e-03, // 204
-1.674766949616971e+01, // 205
 9.999577377230777e-03, // 206
 8.488854876693573e-07, // 207
 2.465279252083182e-02, // 208
 5.801469233050071e-05, // 209
-2.214789285060068e+00, // 210
 1.141628932904846e+00, // 211
 2.487511832548881e+01, // 212
-9.739532133490058e-05, // 213
 4.356370597433608e-02, // 214
-6.964575662281742e+00, // 215
 4.795570159076989e-02, // 216
 9.765347501700621e-02, // 217
-6.009715768741929e+00, // 218
 2.796092623932332e-03, // 219
 3.463377850936466e+01, // 220
 1.474276036000630e-03, // 221
 8.412218162061317e-04, // 222
 1.235259978999048e-01, // 223
 1.512046985535377e-01, // 224
 3.862092124064131e-02, // 225
 1.234959953077248e-05, // 226
 4.198574375580526e-02, // 227
-7.914213921212405e+00, // 228
-1.083745781169305e-02, // 229
-2.391981083193258e-05, // 230
-4.549600659453264e+00, // 231
 2.303095233054306e-01, // 232
 1.044359931056302e-03, // 233
 2.264754104020949e+00, // 234
 6.149189955440972e+00, // 235
 8.789523963319409e-04, // 236
 1.520153912471034e+00, // 237
 9.088999103028770e-04, // 238
 5.543898316649059e+00, // 239
-3.748517850743265e+00, // 240
 8.568030032699944e-01, // 241
 4.011122732739716e-02, // 242
 8.153803085134638e-01, // 243
 1.445445767972670e-05, // 244
-1.104046541823230e-03, // 245
-7.109633499725333e-04, // 246
 5.436774334378062e-04, // 247
 1.125721910135185e+00, // 248
-9.764505697285774e-01, // 249
 2.898410074783787e-02, // 250
 2.056254083921268e-03, // 251
 1.713968726047619e-02, // 252
 3.711573792550134e-02, // 253
-1.691778021356527e-01, // 254
 4.486665734298521e-01, // 255
-1.948716246092321e-02, // 256
 7.046123796204029e-07, // 257
 7.416854104847133e+00, // 258
-6.077663164028767e-02, // 259
-6.738370125293429e+01, // 260
 9.260530666047347e-03, // 261
 1.238709735456836e-01, // 262
-1.100073411512487e+00, // 263
 3.035299329853389e-04, // 264
 3.893344708803968e+00, // 265
 2.165070048747651e-05, // 266
-4.796270281789978e+00, // 267
 1.913685707344445e-01, // 268
 7.253038635229029e-04, // 269
-5.718033398855965e-01, // 270
-1.072889792017336e-03, // 271
-9.752626817025231e-01, // 272
 2.464515909472969e-01, // 273
-6.639023511023822e-01, // 274
 3.198908636254036e-02, // 275
 9.164855066006602e-06, // 276
 1.171191899907343e+01, // 277
-1.496479972944073e-03, // 278
 6.713858170424129e-03, // 279
 5.007774232193078e-01, // 280
-3.922171292728332e-03, // 281
-5.123680705059623e-05, // 282
-1.336199203791317e+00, // 283
-8.306979712196263e-03, // 284
 1.189843047045328e+00, // 285
 1.175588099791827e+00, // 286
-8.726307113926701e-04, // 287
 1.496514755599319e-01, // 288
 2.981511582397324e-06, // 289
-7.034387601002652e-02, // 290
-3.193949081782105e+00, // 291
 7.688929541389313e-03, // 292
 2.247983865037268e+01, // 293
-8.606773323450055e-03, // 294
 6.941561884134629e-01, // 295
-8.766035673280856e-01, // 296
 5.150310999104287e+00, // 297
 1.596435705166030e-02, // 298
-1.153353653023311e-03, // 299
 1.438426968327244e-02, // 300
-1.842021940460480e+01, // 301
-1.387911597983954e+00, // 302
-1.428952335112277e+01, // 303
-7.367425349463584e+00, // 304
 1.833582106448321e-01, // 305
 2.664274925831816e-04, // 306
 5.561315586921332e-02, // 307
-3.813783216369284e+00, // 308
 2.079026782033144e-01, // 309
-4.561341837964809e-01, // 310
 9.940267623941462e-03, // 311
-7.502690315246306e-06, // 312
-9.159513756796597e-03, // 313
 1.340200339731768e-02, // 314
-5.564307713399175e-02, // 315
-1.463410675279343e-02, // 316
-3.490845900559923e-04, // 317
-6.992994349581140e+00, // 318
-1.567263428943910e-01, // 319
 7.847028578975192e-02, // 320
 1.121495629984286e-01, // 321
 1.007378743191846e+00, // 322
 1.822835045953743e-01, // 323
-1.763213791239304e-02, // 324
 2.056865186576417e+02, // 325
-1.066400504075611e-03, // 326
 3.689459841792663e+00, // 327
-3.060169612853612e+00, // 328
-2.296980832624350e-02, // 329
 5.014904806605320e-03, // 330
 2.045964154411254e+00, // 331
-1.837011844182590e-03, // 332
-1.551391503880563e-05, // 333
-5.163957885359048e-03, // 334
 7.497254464073717e-02, // 335
-7.207322389328786e-01, // 336
 6.776114707645124e+00, // 337
 7.525704392760977e-01, // 338
 3.572212936731471e-02, // 339
-2.136041745413528e-02, // 340
 5.770674497340013e-03, // 341
 1.431025839664138e-03, // 342
-1.923226610131200e+00, // 343
 1.912886639673983e-05, // 344
 1.597142993166820e-04, // 345
-3.307680730911380e-03, // 346
-4.922540241434566e+00, // 347
 2.867601534922616e-03, // 348
-3.066087905321299e+00, // 349
-4.553664107460084e-02, // 350
 1.322138329633372e-01, // 351
 4.544369996152396e-01, // 352
-1.138155914195157e+00, // 353
-2.378170894399298e-02, // 354
 2.475719250669047e-05, // 355
-8.431507971337724e-07, // 356
-3.868140749368091e-03, // 357
 6.952410523048593e-03, // 358
 4.842105858598649e-03, // 359
-2.241035107018743e+01, // 360
-1.141533163467902e-03, // 361
-3.593958403712560e-07, // 362
-1.374049443628304e+00, // 363
 1.172674587725694e+01, // 364
-8.884406925585567e-01, // 365
 1.838407601066553e+00, // 366
 1.664518227667212e-01, // 367
-2.541498762465235e+01, // 368
 3.240086213597889e+00, // 369
-7.044640290700685e-02, // 370
-9.165728264023809e+01, // 371
 4.580315137054395e+00, // 372
-8.243876937387001e-05, // 373
-4.472563228714126e+00, // 374
-1.987414020213625e-01, // 375
-3.160514211233190e-02, // 376
-3.720748171169798e-02, // 377
-4.005809108256976e-01, // 378
-5.576436541157079e-01, // 379
-8.765807923468997e-04, // 380
 2.479244442985117e-02, // 381
-3.217546939478989e-01, // 382
 4.514172532468397e-01, // 383
 2.132579213249199e-01, // 384
 1.939312837453414e-03, // 385
-5.463453210355973e+00, // 386
-2.634624008523841e-03, // 387
-6.567750034220347e-01, // 388
 5.317765562767528e+01, // 389
 6.027016034849169e-01, // 390
 3.235619363869706e+01, // 391
-1.050874785876603e-02, // 392
-1.406901008052624e-03, // 393
-3.103046552343679e-03, // 394
 1.745786330352498e+01, // 395
-1.785998950055137e-02, // 396
 1.502223398926982e-01, // 397
-1.064871274875248e-02, // 398
-4.790562167412788e-03, // 399
 7.374816389204007e-03, // 400
-3.334577351470413e-04, // 401
-1.439319966305114e+00, // 402
 3.343781455112465e-03, // 403
-7.788010411715439e-02, // 404
 2.323108965286399e+00, // 405
 2.324978532853036e+00, // 406
-1.867143708396256e-03, // 407
-1.040769578752323e-03, // 408
-1.966392402085471e-04, // 409
 9.306144098400076e-06, // 410
-1.381675673807580e-03, // 411
-1.881411481872492e-01, // 412
-1.281369023464886e-04, // 413
-9.139300349432373e-04, // 414
 2.730698806261698e+00, // 415
 2.485016281024538e+01, // 416
 1.772258756554937e-03, // 417
-5.348379824855107e-02, // 418
-7.202328023909980e+00, // 419
-1.200100098458796e-04, // 420
-7.006352298818381e-01, // 421
-2.412569067892058e+00, // 422
-1.578625902526947e-03, // 423
 1.760951060497483e+00, // 424
 1.251853722141828e-03, // 425
 6.034176682070334e-04, // 426
 1.975662094186727e+00, // 427
-1.899231105928469e-01, // 428
 2.010391049064938e-02, // 429
 1.119775264967099e-01, // 430
 4.101414273976683e+01, // 431
 1.081869010720396e+01, // 432
 1.562332490627762e+00, // 433
 1.037584855255481e+00, // 434
 6.857548721808206e+01, // 435
-6.038099672566764e-06, // 436
 1.848541676300455e-03, // 437
-7.461843014065644e+00, // 438
 5.911239020163262e-03, // 439
 5.817765838483787e-05, // 440
-1.257356066634111e+00, // 441
 1.234288766276253e+00, // 442
 2.648096134285459e-03, // 443
-3.088485531771722e-02, // 444
-3.093293255537717e+01, // 445
-1.625080043276579e-02, // 446
-2.074668455913734e-02, // 447
 1.727720656554861e-04, // 448
-5.388885282926590e+00, // 449
-5.422785025351184e+00, // 450
 3.165198171634380e+01, // 451
 4.901266705778912e+01, // 452
-7.041641867334233e+01, // 453
-4.339782081773161e-04, // 454
 2.976363556926400e+00, // 455
-1.646538676958855e-04, // 456
-1.621159209349989e-01, // 457
-1.753199348276763e+00, // 458
 1.167922035183819e-04, // 459
 6.559834337367958e+00, // 460
-6.979408716406912e-03, // 461
 3.040961445244802e+00, // 462
-1.730390434062906e+01, // 463
 3.283683804285535e+00, // 464
 8.895417282254409e-02, // 465
-2.046619036235156e-02, // 466
 1.493709068731651e-03, // 467
 2.715187597621080e-03, // 468
-3.776108295094393e+00, // 469
-1.019586998141061e-02, // 470
 6.772336286711326e-05, // 471
-4.185923972890027e+00, // 472
-7.726406143939756e-03, // 473
-3.885359519926254e-04, // 474
 4.273996244781337e-01, // 475
-3.255981940284439e-02, // 476
 1.971168864214568e-01, // 477
-4.029782330014042e-03, // 478
-9.888003204109391e-04, // 479
-2.444059153727687e-02, // 480
-1.267440936003245e-02, // 481
 4.018709393125027e-01, // 482
 5.122366160930230e-01, // 483
 2.218807739442711e+00, // 484
-1.317120576584223e-01, // 485
-6.218558758822166e-03, // 486
 1.299540538730405e+00, // 487
-6.320307016091915e-02, // 488
-2.575348996611741e-06, // 489
-3.646043025706732e-03, // 490
-5.432249587926850e-01, // 491
 1.239846393057434e-06, // 492
-1.887236038997717e+01, // 493
 5.386340139203918e-02, // 494
-4.406833344923537e-03, // 495
-2.249774413725945e-02, // 496
-5.924626897193252e-05, // 497
 1.246898295878332e-01, // 498
-3.539501219785788e+00, // 499
 2.103497453044589e+01, // 500
-1.893443055672670e-02, // 501
-6.765487305740141e-02, // 502
 5.544831044632087e-05, // 503
 4.769255866722236e-01, // 504
-2.791286864147979e+00, // 505
-4.945027602291759e-02, // 506
 9.525325525526419e-01, // 507
 1.115128379502024e+01, // 508
 8.064169353362323e-03, // 509
-4.755015993052600e-01, // 510
-1.014447655422842e-03, // 511
 3.505607250438209e-08, // 512
-3.777042407093373e-01, // 513
-5.999929900271412e+00, // 514
 6.680946409696037e-02, // 515
-6.815689985910561e-04, // 516
 1.046789332176415e-07, // 517
 1.531206918440141e-05, // 518
-3.729845597740258e-01, // 519
-2.255548442090649e-02, // 520
-6.620843436537120e-01, // 521
 1.266668931920032e+02, // 522
-6.835540063213635e-03, // 523
-6.113986247930577e+02, // 524
-1.991876987617095e+01, // 525
-6.650588700499896e+00, // 526
-1.042775000216342e-07, // 527
 7.253322283801049e-02, // 528
 2.481983359416030e-03, // 529
-6.680683226774506e-02, // 530
-2.013213358643868e-02, // 531
 1.591816890364112e+01, // 532
 1.980383629117110e-03, // 533
-1.099527972774104e-11, // 534
 1.588389052011755e-05, // 535
 1.380395672878348e-03, // 536
-1.033892541237770e+00, // 537
 1.947973461641924e+01, // 538
 5.219866210230528e+01, // 539
-8.681027787568946e-04, // 540
 2.275804004317869e-05, // 541
-2.502845035203813e-03, // 542
-2.986516334024229e-02, // 543
 2.367887608848886e+00, // 544
-2.943722269959519e-02, // 545
-1.128681647859539e-06, // 546
-9.201276340960095e-03, // 547
 5.606940694566580e-02, // 548
-5.683633174933000e+01, // 549
 7.713508931733262e-01, // 550
 3.366817968373659e+00, // 551
-1.867520508285203e-01, // 552
-1.352984792928633e-02, // 553
-1.972043245101288e-02, // 554
 1.865797447586403e-03, // 555
 1.343640840956069e+01, // 556
-7.966090393641488e-03, // 557
-5.230381698102403e+00, // 558
-1.763769243114671e+00, // 559
-6.782476671830215e-06, // 560
-1.368769358215754e-08, // 561
 6.114825964247039e+00, // 562
-2.443668165152092e-05, // 563
 1.441091708080751e-06, // 564
-1.140358109396201e-02, // 565
-3.682579627579224e+00, // 566
-4.047617701328023e-02, // 567
 1.863857520919913e-02, // 568
-2.899940165812726e-03, // 569
-1.521723773106403e-05, // 570
-1.086743910286433e-06, // 571
 7.882328762649680e+00, // 572
 8.892024683962831e+00, // 573
 1.669453770716013e-01, // 574
 1.549998881784788e-06, // 575
-5.354190408152913e-01, // 576
 2.793356812339502e-02, // 577
 1.315645408234747e-01, // 578
 2.618245924649467e-01, // 579
-1.425567847673159e-01, // 580
-2.810390243685471e-05, // 581
 6.319140197671757e-04, // 582
-3.660757906563632e+01, // 583
-1.330030068937110e+00, // 584
-2.454265622994447e-02, // 585
-2.358013878054862e+01, // 586
-1.217673211232488e-03, // 587
-1.529960893841101e-02, // 588
-1.460639867614487e-02, // 589
-4.565692042275246e-01, // 590
-2.708708117710081e-02, // 591
 2.204392662958351e+01, // 592
-1.960938159061999e-04, // 593
 2.954850125528728e-01, // 594
 7.016316849858209e-04, // 595
-2.307425652319090e-02, // 596
-7.906065433844760e-04, // 597
 2.362753465306763e-04, // 598
-4.573466998083948e+01, // 599
 2.167156606533230e+00, // 600
 7.464735827593693e-04, // 601
 4.472242023673328e-03}; // 602

      m_k_intra_AB =  1.687639706175954e+00; // A^(-1))
      m_d_intra_AB =  4.767741146570285e+00; // A^(-1))
      m_k_intra_BB =  6.452329685729143e-01; // A^(-1))
      m_d_intra_BB =  2.779024709627167e+00; // A^(-1))
      m_k_intra_CD =  1.969180198641554e+00; // A^(-1))
      m_d_intra_CD =  1.529140402135649e+00; // A^(-1))
      m_k_intra_DD =  7.617646808139036e-01; // A^(-1))
      m_d_intra_DD =  5.418827218087009e+00; // A^(-1))
      m_k_AC =  1.421128400517191e+00; // A^(-1))
      m_d_AC =  3.059161965065690e+00; // A^(-1))
      m_k_AD =  2.278635677372154e+00; // A^(-1))
      m_d_AD =  4.416012049142289e+00; // A^(-1))
      m_k_BC =  1.403203969175851e+00; // A^(-1))
      m_d_BC =  3.607335264427631e+00; // A^(-1))
      m_k_BD =  8.069925525486774e-01; // A^(-1))
      m_d_BD =  1.809061271452679e+00; // A^(-1))
      m_k_CZ =  8.994857338690051e-01; // A^(-1))
      m_d_CZ =  6.999903347087452e+00; // A^(-1))
      m_k_DZ =  3.201841430649240e-01; // A^(-1))
      m_d_DZ =  4.745950304692959e+00; // A^(-1))
      m_r2i =  8.000000000000000e+00; // A
      m_r2f =  9.000000000000000e+00; // A
   }
}



double x2b_A1B2Z2_C1D4_v1x::f_switch(const double& r, double& g) const
{
    if (r > m_r2f) {
        g = 0.0;
        return 0.0;
    } else if (r > m_r2i) {
        const double t1 = M_PI/(m_r2f - m_r2i);
        const double x = (r - m_r2i)*t1;
        g = - std::sin(x)*t1/2.0;
        return (1.0 + std::cos(x))/2.0;
    } else {
        g = 0.0;
        return 1.0;
    }
}

//----------------------------------------------------------------------------//

double x2b_A1B2Z2_C1D4_v1x::eval(const double* xyz1, const double* xyz2, const size_t ndim) const
{

    std::vector<double> energies(ndim,0.0);

    for (size_t j = 0; j < ndim; j++) {
        double mon1[9];
        double mon2[15];

        std::copy(xyz1 + j * 9, xyz1 + (j+1) * 9, mon1);
        std::copy(xyz2 + j * 15, xyz2 + (j+1) * 15, mon2);


        const double d12[3] = {mon1[0] -  mon2[0],
                               mon1[1] -  mon2[1],
                               mon1[2] -  mon2[2]};
    
        const double r12sq = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
        const double r12 = std::sqrt(r12sq);
    
        if (r12 > m_r2f)
            return 0.0;
    
        double xcrd[24]; // coordinates of real sites ONLY
    
        std::copy(mon1, mon1 + 9, xcrd);
        std::copy(mon2, mon2 + 15, xcrd + 9);
        
        double v[38];
        
        double sw = 0.0;
        double gsw = 0.0;
        
        const double* A_1_a= xcrd + 0;
        const double* B_1_a= xcrd + 3;
        const double* B_2_a= xcrd + 6;
    
        const double* C_1_b= xcrd + 9;
        const double* D_1_b= xcrd + 12;
        const double* D_2_b= xcrd + 15;
        const double* D_3_b= xcrd + 18;
        const double* D_4_b= xcrd + 21;
    
        double Z_1_a[3];
        double Z_2_a[3];
    
    
    //    vsites virt;
        double w12 =     -9.721486914088159e-02;  //from MBpol
        double w13 =     -9.721486914088159e-02;
        double wcross =   9.859272078406150e-02;
    
        monomer m;
        
        m.setup(A_1_a, w12, wcross,
                 Z_1_a, Z_2_a);
                        
        variable vr[38];
        
        v[0]  = vr[0].v_exp0(m_d_intra_AB, m_k_intra_AB, A_1_a, B_1_a);
        v[1]  = vr[1].v_exp0(m_d_intra_AB, m_k_intra_AB, A_1_a, B_2_a);
        v[2]  = vr[2].v_exp0(m_d_intra_BB, m_k_intra_BB, B_1_a, B_2_a);
    
        v[3]  = vr[3].v_exp0(m_d_intra_CD, m_k_intra_CD, C_1_b, D_1_b);
        v[4]  = vr[4].v_exp0(m_d_intra_CD, m_k_intra_CD, C_1_b, D_2_b);
        v[5]  = vr[5].v_exp0(m_d_intra_CD, m_k_intra_CD, C_1_b, D_3_b);
        v[6]  = vr[6].v_exp0(m_d_intra_CD, m_k_intra_CD, C_1_b, D_4_b);
        v[7]  = vr[7].v_exp0(m_d_intra_DD, m_k_intra_DD, D_1_b, D_2_b);
        v[8]  = vr[8].v_exp0(m_d_intra_DD, m_k_intra_DD, D_1_b, D_3_b);
        v[9]  = vr[9].v_exp0(m_d_intra_DD, m_k_intra_DD, D_1_b, D_4_b);
        v[10]  = vr[10].v_exp0(m_d_intra_DD, m_k_intra_DD, D_2_b, D_3_b);
        v[11]  = vr[11].v_exp0(m_d_intra_DD, m_k_intra_DD, D_2_b, D_4_b);
        v[12]  = vr[12].v_exp0(m_d_intra_DD, m_k_intra_DD, D_3_b, D_4_b);
    
        v[13]  = vr[13].v_exp0(m_d_AC, m_k_AC, A_1_a, C_1_b);
        v[14]  = vr[14].v_exp0(m_d_AD, m_k_AD, A_1_a, D_1_b);
        v[15]  = vr[15].v_exp0(m_d_AD, m_k_AD, A_1_a, D_2_b);
        v[16]  = vr[16].v_exp0(m_d_AD, m_k_AD, A_1_a, D_3_b);
        v[17]  = vr[17].v_exp0(m_d_AD, m_k_AD, A_1_a, D_4_b);
    
        v[18]  = vr[18].v_exp0(m_d_BC, m_k_BC, B_1_a, C_1_b);
        v[19]  = vr[19].v_exp0(m_d_BD, m_k_BD, B_1_a, D_1_b);
        v[20]  = vr[20].v_exp0(m_d_BD, m_k_BD, B_1_a, D_2_b);
        v[21]  = vr[21].v_exp0(m_d_BD, m_k_BD, B_1_a, D_3_b);
        v[22]  = vr[22].v_exp0(m_d_BD, m_k_BD, B_1_a, D_4_b);
    
        v[23]  = vr[23].v_exp0(m_d_BC, m_k_BC, B_2_a, C_1_b);
        v[24]  = vr[24].v_exp0(m_d_BD, m_k_BD, B_2_a, D_1_b);
        v[25]  = vr[25].v_exp0(m_d_BD, m_k_BD, B_2_a, D_2_b);
        v[26]  = vr[26].v_exp0(m_d_BD, m_k_BD, B_2_a, D_3_b);
        v[27]  = vr[27].v_exp0(m_d_BD, m_k_BD, B_2_a, D_4_b);
    
        v[28]  = vr[28].v_coul0(m_d_CZ, m_k_CZ, Z_1_a, C_1_b);
        v[29]  = vr[29].v_coul0(m_d_DZ, m_k_DZ, Z_1_a, D_1_b);
        v[30]  = vr[30].v_coul0(m_d_DZ, m_k_DZ, Z_1_a, D_2_b);
        v[31]  = vr[31].v_coul0(m_d_DZ, m_k_DZ, Z_1_a, D_3_b);
        v[32]  = vr[32].v_coul0(m_d_DZ, m_k_DZ, Z_1_a, D_4_b);
    
        v[33]  = vr[33].v_coul0(m_d_CZ, m_k_CZ, Z_2_a, C_1_b);
        v[34]  = vr[34].v_coul0(m_d_DZ, m_k_DZ, Z_2_a, D_1_b);
        v[35]  = vr[35].v_coul0(m_d_DZ, m_k_DZ, Z_2_a, D_2_b);
        v[36]  = vr[36].v_coul0(m_d_DZ, m_k_DZ, Z_2_a, D_3_b);
        v[37]  = vr[37].v_coul0(m_d_DZ, m_k_DZ, Z_2_a, D_4_b);
    
         
        
        sw = f_switch(r12, gsw);
        
        energies[j] = sw*polynomial::eval(coefficients.data(), v);

    }

    double energy = 0.0;
    for (size_t i = 0; i < ndim; i++) {
      energy += energies[i];
    }

    return energy;
}


double x2b_A1B2Z2_C1D4_v1x::eval(const double* xyz1, const double* xyz2,
                double * grad1, double * grad2, const size_t ndim) const
{
    std::vector<double> energies(ndim,0.0);

    for (size_t j = 0; j < ndim; j++) {
        double mon1[9];
        double mon2[15];

        std::copy(xyz1 + j * 9, xyz1 + (j+1) * 9, mon1);
        std::copy(xyz2 + j * 15, xyz2 + (j+1) * 15, mon2);


        const double d12[3] = {mon1[0] -  mon2[0],
                               mon1[1] -  mon2[1],
                               mon1[2] -  mon2[2]};
    
        const double r12sq = d12[0]*d12[0] + d12[1]*d12[1] + d12[2]*d12[2];
        const double r12 = std::sqrt(r12sq);
    
        if (r12 > m_r2f)
            return 0.0;
    
        double xcrd[24]; // coordinates of real sites ONLY
    
        std::copy(mon1, mon1 + 9, xcrd);
        std::copy(mon2, mon2 + 15, xcrd + 9);
        
        double v[38];
        
        double sw = 0.0;
        double gsw = 0.0;
        
        const double* A_1_a= xcrd + 0;
        const double* B_1_a= xcrd + 3;
        const double* B_2_a= xcrd + 6;
    
        const double* C_1_b= xcrd + 9;
        const double* D_1_b= xcrd + 12;
        const double* D_2_b= xcrd + 15;
        const double* D_3_b= xcrd + 18;
        const double* D_4_b= xcrd + 21;
    
        double Z_1_a[3];
        double Z_2_a[3];
    
    
        //vsites virt;
        double w12 =     -9.721486914088159e-02;  //from MBpol
        double w13 =     -9.721486914088159e-02;
        double wcross =   9.859272078406150e-02;
    
        monomer m;
        
        m.setup(A_1_a, w12, wcross, 
                 Z_1_a, Z_2_a);
                        
        variable vr[38];
        
        v[0]  = vr[0].v_exp0(m_d_intra_AB, m_k_intra_AB, A_1_a, B_1_a);
        v[1]  = vr[1].v_exp0(m_d_intra_AB, m_k_intra_AB, A_1_a, B_2_a);
        v[2]  = vr[2].v_exp0(m_d_intra_BB, m_k_intra_BB, B_1_a, B_2_a);
    
        v[3]  = vr[3].v_exp0(m_d_intra_CD, m_k_intra_CD, C_1_b, D_1_b);
        v[4]  = vr[4].v_exp0(m_d_intra_CD, m_k_intra_CD, C_1_b, D_2_b);
        v[5]  = vr[5].v_exp0(m_d_intra_CD, m_k_intra_CD, C_1_b, D_3_b);
        v[6]  = vr[6].v_exp0(m_d_intra_CD, m_k_intra_CD, C_1_b, D_4_b);
        v[7]  = vr[7].v_exp0(m_d_intra_DD, m_k_intra_DD, D_1_b, D_2_b);
        v[8]  = vr[8].v_exp0(m_d_intra_DD, m_k_intra_DD, D_1_b, D_3_b);
        v[9]  = vr[9].v_exp0(m_d_intra_DD, m_k_intra_DD, D_1_b, D_4_b);
        v[10]  = vr[10].v_exp0(m_d_intra_DD, m_k_intra_DD, D_2_b, D_3_b);
        v[11]  = vr[11].v_exp0(m_d_intra_DD, m_k_intra_DD, D_2_b, D_4_b);
        v[12]  = vr[12].v_exp0(m_d_intra_DD, m_k_intra_DD, D_3_b, D_4_b);
    
        v[13]  = vr[13].v_exp0(m_d_AC, m_k_AC, A_1_a, C_1_b);
        v[14]  = vr[14].v_exp0(m_d_AD, m_k_AD, A_1_a, D_1_b);
        v[15]  = vr[15].v_exp0(m_d_AD, m_k_AD, A_1_a, D_2_b);
        v[16]  = vr[16].v_exp0(m_d_AD, m_k_AD, A_1_a, D_3_b);
        v[17]  = vr[17].v_exp0(m_d_AD, m_k_AD, A_1_a, D_4_b);
    
        v[18]  = vr[18].v_exp0(m_d_BC, m_k_BC, B_1_a, C_1_b);
        v[19]  = vr[19].v_exp0(m_d_BD, m_k_BD, B_1_a, D_1_b);
        v[20]  = vr[20].v_exp0(m_d_BD, m_k_BD, B_1_a, D_2_b);
        v[21]  = vr[21].v_exp0(m_d_BD, m_k_BD, B_1_a, D_3_b);
        v[22]  = vr[22].v_exp0(m_d_BD, m_k_BD, B_1_a, D_4_b);
    
        v[23]  = vr[23].v_exp0(m_d_BC, m_k_BC, B_2_a, C_1_b);
        v[24]  = vr[24].v_exp0(m_d_BD, m_k_BD, B_2_a, D_1_b);
        v[25]  = vr[25].v_exp0(m_d_BD, m_k_BD, B_2_a, D_2_b);
        v[26]  = vr[26].v_exp0(m_d_BD, m_k_BD, B_2_a, D_3_b);
        v[27]  = vr[27].v_exp0(m_d_BD, m_k_BD, B_2_a, D_4_b);
    
        v[28]  = vr[28].v_coul0(m_d_CZ, m_k_CZ, Z_1_a, C_1_b);
        v[29]  = vr[29].v_coul0(m_d_DZ, m_k_DZ, Z_1_a, D_1_b);
        v[30]  = vr[30].v_coul0(m_d_DZ, m_k_DZ, Z_1_a, D_2_b);
        v[31]  = vr[31].v_coul0(m_d_DZ, m_k_DZ, Z_1_a, D_3_b);
        v[32]  = vr[32].v_coul0(m_d_DZ, m_k_DZ, Z_1_a, D_4_b);
    
        v[33]  = vr[33].v_coul0(m_d_CZ, m_k_CZ, Z_2_a, C_1_b);
        v[34]  = vr[34].v_coul0(m_d_DZ, m_k_DZ, Z_2_a, D_1_b);
        v[35]  = vr[35].v_coul0(m_d_DZ, m_k_DZ, Z_2_a, D_2_b);
        v[36]  = vr[36].v_coul0(m_d_DZ, m_k_DZ, Z_2_a, D_3_b);
        v[37]  = vr[37].v_coul0(m_d_DZ, m_k_DZ, Z_2_a, D_4_b);
    
         
        
        double g[38];
        
        sw = f_switch(r12, gsw);

        energies[j] = polynomial::eval(coefficients.data(), v, g);
        
        double xgrd[30];
        std::fill(xgrd, xgrd + 30, 0.0);
    
        double* A_1_a_g = xgrd + 0;
        double* B_1_a_g = xgrd + 3;
        double* B_2_a_g = xgrd + 6;
    
        double* C_1_b_g = xgrd + 9;
        double* D_1_b_g = xgrd + 12;
        double* D_2_b_g = xgrd + 15;
        double* D_3_b_g = xgrd + 18;
        double* D_4_b_g = xgrd + 21;
    
        double* Z_1_a_g = xgrd + 24;
        double* Z_2_a_g = xgrd + 27;
    
        vr[0].grads(g[0], A_1_a_g, B_1_a_g, A_1_a, B_1_a);
        vr[1].grads(g[1], A_1_a_g, B_2_a_g, A_1_a, B_2_a);
        vr[2].grads(g[2], B_1_a_g, B_2_a_g, B_1_a, B_2_a);
    
        vr[3].grads(g[3], C_1_b_g, D_1_b_g, C_1_b, D_1_b);
        vr[4].grads(g[4], C_1_b_g, D_2_b_g, C_1_b, D_2_b);
        vr[5].grads(g[5], C_1_b_g, D_3_b_g, C_1_b, D_3_b);
        vr[6].grads(g[6], C_1_b_g, D_4_b_g, C_1_b, D_4_b);
        vr[7].grads(g[7], D_1_b_g, D_2_b_g, D_1_b, D_2_b);
        vr[8].grads(g[8], D_1_b_g, D_3_b_g, D_1_b, D_3_b);
        vr[9].grads(g[9], D_1_b_g, D_4_b_g, D_1_b, D_4_b);
        vr[10].grads(g[10], D_2_b_g, D_3_b_g, D_2_b, D_3_b);
        vr[11].grads(g[11], D_2_b_g, D_4_b_g, D_2_b, D_4_b);
        vr[12].grads(g[12], D_3_b_g, D_4_b_g, D_3_b, D_4_b);
    
        vr[13].grads(g[13], A_1_a_g, C_1_b_g, A_1_a, C_1_b);
        vr[14].grads(g[14], A_1_a_g, D_1_b_g, A_1_a, D_1_b);
        vr[15].grads(g[15], A_1_a_g, D_2_b_g, A_1_a, D_2_b);
        vr[16].grads(g[16], A_1_a_g, D_3_b_g, A_1_a, D_3_b);
        vr[17].grads(g[17], A_1_a_g, D_4_b_g, A_1_a, D_4_b);
    
        vr[18].grads(g[18], B_1_a_g, C_1_b_g, B_1_a, C_1_b);
        vr[19].grads(g[19], B_1_a_g, D_1_b_g, B_1_a, D_1_b);
        vr[20].grads(g[20], B_1_a_g, D_2_b_g, B_1_a, D_2_b);
        vr[21].grads(g[21], B_1_a_g, D_3_b_g, B_1_a, D_3_b);
        vr[22].grads(g[22], B_1_a_g, D_4_b_g, B_1_a, D_4_b);
    
        vr[23].grads(g[23], B_2_a_g, C_1_b_g, B_2_a, C_1_b);
        vr[24].grads(g[24], B_2_a_g, D_1_b_g, B_2_a, D_1_b);
        vr[25].grads(g[25], B_2_a_g, D_2_b_g, B_2_a, D_2_b);
        vr[26].grads(g[26], B_2_a_g, D_3_b_g, B_2_a, D_3_b);
        vr[27].grads(g[27], B_2_a_g, D_4_b_g, B_2_a, D_4_b);
    
        vr[28].grads(g[28], Z_1_a_g, C_1_b_g, Z_1_a, C_1_b);
        vr[29].grads(g[29], Z_1_a_g, D_1_b_g, Z_1_a, D_1_b);
        vr[30].grads(g[30], Z_1_a_g, D_2_b_g, Z_1_a, D_2_b);
        vr[31].grads(g[31], Z_1_a_g, D_3_b_g, Z_1_a, D_3_b);
        vr[32].grads(g[32], Z_1_a_g, D_4_b_g, Z_1_a, D_4_b);
    
        vr[33].grads(g[33], Z_2_a_g, C_1_b_g, Z_2_a, C_1_b);
        vr[34].grads(g[34], Z_2_a_g, D_1_b_g, Z_2_a, D_1_b);
        vr[35].grads(g[35], Z_2_a_g, D_2_b_g, Z_2_a, D_2_b);
        vr[36].grads(g[36], Z_2_a_g, D_3_b_g, Z_2_a, D_3_b);
        vr[37].grads(g[37], Z_2_a_g, D_4_b_g, Z_2_a, D_4_b);
    
    
    
        // ##DEFINE HERE## the redistribution of the gradients
        
    
        m.grads(Z_1_a_g, Z_2_a_g, 
                 w12, wcross, A_1_a_g);
        
        for (int i = 0; i < 9; ++i) {
            grad1[i+j*9] += sw*xgrd[i];
        }
    
        for (int i = 0; i < 15; ++i) {
            grad2[i+j*15] += sw*xgrd[i + 9];
        }
    
        // gradient of the switch
    
        gsw *= energies[j]/r12;
        for (int i = 0; i < 3; ++i) {
            const double d = gsw*d12[i];
            grad1[i+j*9] += d;
            grad2[i+j*15] -= d;
        }

    }

    double energy = 0.0;
    for (size_t i = 0; i < ndim; i++) {
      energy += energies[i];
    }

    return energy;

}

} // namespace x2b_A1B2Z2_C1D4

////////////////////////////////////////////////////////////////////////////////