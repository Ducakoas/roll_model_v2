function [dxdt,y,vx,r_t,Y,phi_t,phi_s] = vehicle_18dof_full_roll_road(x,u,mu_road)
%#codegen
% =========================================================================
% 18-STATE ARTICULATED TRACTOR-SEMITRAILER PLANT
% DYNAMIC ROLL OF BOTH TRACTOR AND SEMITRAILER
%
% This file extends the user's current 16-state plant by adding:
%
%   x(17) = phi_s       semitrailer roll angle [rad]
%   x(18) = dphi_s      semitrailer roll rate  [rad/s]
%
% The simultaneous algebraic unknown vector is now
%
%   q = [dvx dvy dr_t dr_s Hx Hy ddphi_t ddphi_s]^T
%
% and an 8x8 system is solved at every plant evaluation.
%
% State vector:
%   1  vx
%   2  vy
%   3  r_t
%   4  r_s
%   5  gamma
%   6  psi
%   7  Y
%   8  X
%   9  omega_fL
%  10  omega_fR
%  11  omega_rL
%  12  omega_rR
%  13  omega_sL
%  14  omega_sR
%  15  phi_t
%  16  dphi_t
%  17  phi_s
%  18  dphi_s
%
% Input vector remains 9 physical channels:
%   1  delta_f
%   2  T_drive (tractor drive axle D1; original symbol T_f_drive)
%   3  M_act
%   4  Tb_fL
%   5  Tb_fR
%   6  Tb_rL
%   7  Tb_rR
%   8  Tb_sL
%   9  Tb_sR
%
% IMPORTANT MODELING NOTE
% -----------------------
% The original 16-state model already contains tractor roll directly in the
% longitudinal/lateral/yaw equations.  This extension applies the analogous
% semitrailer sprung-mass and roll-yaw terms and adds an independent
% semitrailer roll equation.
%
% Tractor:
%
% mt*(dvx-r_t*vy) - mst*ht*dphi_t*r_t = Fx_t + Hx
%
% mt*(dvy+r_t*vx) + mst*ht*ddphi_t = Fy_t + Hy
%
% Izt*dr_t - Ixzt*ddphi_t = Mz_t - Lwt*Hy
%
% Semitrailer roll inertial vector is written first in the semitrailer frame:
%
%   c_s = [ -mss*hs*dphi_s*r_s ;
%            mss*hs*ddphi_s ]
%
% and rotated into the tractor frame before entering the translational
% equations.  The trailer yaw equation also contains -Ixzs*ddphi_s.
%
% Roll equations use an energy-consistent relative-roll coupling:
%
%   M_w = K_phi_w*(phi_t-phi_s) + C_phi_w*(dphi_t-dphi_s)
%
% Tractor receives -M_w and semitrailer receives +M_w.
%
% PARAMETER WARNING
% -----------------
% m_st and I_xzt were already scaled from Karimyan et al. (2023) ratios in
% the 16-state model.  The NEW trailer-roll parameters m_ss, I_xxs, I_xzs,
% K_phi_s, C_phi_s, K_phi_w and C_phi_w are also ratio-scaled estimates.
% Replace them with TruckSim-identified values for final validation.
%
% This is therefore a verification extension of the current model, not a
% claim that the current TruckSim semitrailer roll parameters have already
% been identified.
% =========================================================================

%% 1. MAIN VEHICLE PARAMETERS - from 14DOF_roll.par
% Lead unit sprung 4455 kg, front/rear unsprung axle groups 570/735 kg.
% Trailer sprung 6000 kg, trailer unsprung axle group 640 kg.
% NO EXTRA PAYLOAD is declared in the expanded loaded-combination dataset.
% Total CG yaw inertias below use the parallel-axis theorem with unsprung
% axle-group masses concentrated at axle centers; rotary axle inertias are
% not fully represented. These are approximations, not direct TruckSim Izz.
m_t   = 5760.000000;
m_s   = 6640.000000;
m_tot = m_t + m_s;
g     = 9.81;
I_zt  = 39504.373902;
I_zs  = 186213.725301;
% Sprung masses and their inertia products as given in TruckSim.
m_st  = 4455.0;
m_ss  = 6000.0;
I_xzt = 1626.0;   % sign convention still must be verified
I_xzs = 0.0;
% Estimated parallel-axis inertias about roll axes in current 8x8 solve.
% Tractor roll center interpolated between steer and drive axle;
% trailer roll center approximated from its axle (NOT identified 3D axis).
I_xxt = 3988.825584;
I_xxs = 10413.450000;
% Wheel spin inertias retained pending wheel/hub/tire combined identification.
I_wf = 40.0;
I_wr = 40.0;
I_ws = 50.0;

%% 2. GEOMETRY - TruckSim positions with TOTAL CG convention
% Tractor axle x: front 0.000 m, rear 3.500 m, hitch 2.750 m.
% Sprung CG x = 1.110 m; full CG includes both unsprung axle groups.
% Trailer x from hitch: sprung CG 5.220 m, axle 8.500 m.
L_ft = 1.305130208;  % full tractor CG -> front axle
L_rt = 2.194869792;  % full tractor CG -> rear axle
L_wt = 1.444869792;  % full tractor CG -> fifth wheel
L_t  = L_ft + L_rt;
L_rs = 8.500000000;  % fifth wheel -> semitrailer axle
L_cs = 5.536144578;  % fifth wheel -> FULL semitrailer CG
L_sa = L_rs - L_cs;
T_f = 2.030;         % steer axle left/right tire center track [m]
T_r = 1.863;         % drive axle, paired wheel center track [m]
T_s = 1.863;         % trailer axle, paired wheel center track [m]
R_f = 0.510;
R_r = 0.510;
R_s = 0.510;

%% 3. TIRE PARAMETERS - TruckSim FY_TIRE_CARPET / FX_TIRE_CARPET
% The six wheel states represent 2 steer tires and 4+4 dual tires.
% Per-contact stiffness is interpolated from TruckSim at local Fz.
% Small-angle Fy slope: 1 deg row; small-slip Fx secant: kappa=0.025 row.
% Actual TruckSim uses a full nonlinear 2D force carpet, tire relaxation,
% and compliance; this remains a LOAD-DEPENDENT DUGOFF approximation.
mu = min(max(mu_road,0.10),1.20); % road-aware reduced Dugoff model

%% 4. AERODYNAMIC / ROLLING RESISTANCE
C_D   = 0.56;  % FX_AERO_SHAPING_TABLE, 0deg
A_a   = 10.0;  % AREA_AERO
rho_a = 1.206;  % D_AIR
c_rr_f = 0.0041;  % RR_C, RR_V handled in function
c_rr_r = 0.0041;  % RR_C, RR_V handled in function
c_rr_s = 0.0041;  % RR_C, RR_V handled in function

%% 5. ROLL / LOAD-TRANSFER PARAMETERS - TruckSim reduced equivalent
% Sprung mass CG heights: tractor 1.175m, trailer 1.000m.
% Roll center heights are approximated from axle wheel-center heights and
% 53mm below (steer) / 195mm above (rear/trailer) axle roll centers.
% Trailer roll center along hitch axle is NOT fully determined by par file.
h_t = 0.618064315;  % tractor CG to modeled roll axis, estimated
h_s = 0.275000000;  % trailer CG to axle roll center, estimated
% Two springs/dampers per axle; lever arms: half left/right separation.
% Auxiliary roll moment is supplied in Nm/deg -> converted to Nm/rad.
K_phi_f = 586026.625861;
K_phi_r = 435943.669270;
K_phi_s = 521887.338539;
C_phi_f = 9075.000000;
C_phi_r = 15000.000000;
C_phi_s = 15000.000000;
K_phi = K_phi_f + K_phi_r;
C_phi = C_phi_f + C_phi_r;
eta_phi_f = K_phi_f / max(K_phi,1.0);
eta_phi_r = 1.0 - eta_phi_f;
% Hitch: OPT_ROLL_HITCH=0 and MX_HITCH_TABLE are provided by TruckSim.
% Table's exact role under this option needs solver/Echo confirmation.
% Set use_hitch_Mx_table=0 to isolate the hitch model experimentally.
use_hitch_Mx_table = 1.0;

%% 6. READ STATES
vx = x(1);
vy = x(2);
r_t = x(3);
r_s = x(4);
gamma = x(5);
psi   = x(6);
Y     = x(7);
X     = x(8);
omega_fL = x(9);
omega_fR = x(10);
omega_rL = x(11);
omega_rR = x(12);
omega_sL = x(13);
omega_sR = x(14);
phi_t  = x(15);
dphi_t = x(16);
% NEW independent semitrailer roll states
phi_s  = x(17);
dphi_s = x(18);

%% 7. READ INPUTS
delta_f   = u(1);
T_f_drive = u(2);
M_act     = u(3);
Tb_fL = max(u(4),0.0);
Tb_fR = max(u(5),0.0);
Tb_rL = max(u(6),0.0);
Tb_rR = max(u(7),0.0);
Tb_sL = max(u(8),0.0);
Tb_sR = max(u(9),0.0);

%% 8. STATIC FIFTH-WHEEL / AXLE LOADS
% Static equilibrium of the semitrailer:
F_h0   = m_s*g*L_sa/L_rs;
F_zs_0 = m_s*g*L_cs/L_rs;
% The tractor supports its own weight plus the downward fifth-wheel load.
F_zt_sup = m_t*g + F_h0;
F_zf_0 = (L_rt*F_zt_sup - L_wt*F_h0)/L_t;
F_zr_0 = F_zt_sup - F_zf_0;

%% 9. BASIC KINEMATICS
eps_v = 1.0e-4;
% Tractor front contact-point velocities in tractor frame
Vx_fL_t = vx - 0.5*T_f*r_t;
Vx_fR_t = vx + 0.5*T_f*r_t;
Vy_f_t  = vy + L_ft*r_t;
% Rotate front contact velocities into steered wheel frames
cDel = cos(delta_f);
sDel = sin(delta_f);
Vx_fL_w =  cDel*Vx_fL_t + sDel*Vy_f_t;
Vy_fL_w = -sDel*Vx_fL_t + cDel*Vy_f_t;
Vx_fR_w =  cDel*Vx_fR_t + sDel*Vy_f_t;
Vy_fR_w = -sDel*Vx_fR_t + cDel*Vy_f_t;
% Tractor rear wheel velocities
Vx_rL_w = vx - 0.5*T_r*r_t;
Vx_rR_w = vx + 0.5*T_r*r_t;
Vy_r_w  = vy - L_rt*r_t;
% Fifth-wheel velocity in tractor frame
Vx_h_t = vx;
Vy_h_t = vy - L_wt*r_t;
% Fifth-wheel velocity expressed in semitrailer frame
cg = cos(gamma);
sg = sin(gamma);
Vx_h_s =  cg*Vx_h_t - sg*Vy_h_t;
Vy_h_s =  sg*Vx_h_t + cg*Vy_h_t;
% Semitrailer axle wheel velocities
Vx_sL_w = Vx_h_s - 0.5*T_s*r_s;
Vx_sR_w = Vx_h_s + 0.5*T_s*r_s;
Vy_s_w  = Vy_h_s - L_rs*r_s;

%% 10. SLIP ANGLES
alpha_fL = -atan2(Vy_fL_w,abs(Vx_fL_w)+eps_v);
alpha_fR = -atan2(Vy_fR_w,abs(Vx_fR_w)+eps_v);
alpha_rL = -atan2(Vy_r_w, abs(Vx_rL_w)+eps_v);
alpha_rR = -atan2(Vy_r_w, abs(Vx_rR_w)+eps_v);
alpha_sL = -atan2(Vy_s_w, abs(Vx_sL_w)+eps_v);
alpha_sR = -atan2(Vy_s_w, abs(Vx_sR_w)+eps_v);

%% 11. LONGITUDINAL SLIP RATIOS
sigma_fL = clipval((R_f*omega_fL - Vx_fL_w)/(abs(Vx_fL_w)+eps_v),-0.95,5.0);
sigma_fR = clipval((R_f*omega_fR - Vx_fR_w)/(abs(Vx_fR_w)+eps_v),-0.95,5.0);
sigma_rL = clipval((R_r*omega_rL - Vx_rL_w)/(abs(Vx_rL_w)+eps_v),-0.95,5.0);
sigma_rR = clipval((R_r*omega_rR - Vx_rR_w)/(abs(Vx_rR_w)+eps_v),-0.95,5.0);
sigma_sL = clipval((R_s*omega_sL - Vx_sL_w)/(abs(Vx_sL_w)+eps_v),-0.95,5.0);
sigma_sR = clipval((R_s*omega_sR - Vx_sR_w)/(abs(Vx_sR_w)+eps_v),-0.95,5.0);

%% 12. AERODYNAMIC DRAG
c_vloss = 30.0;  % N/(m/s), identified lumped drivetrain/parasitic loss
F_aero = 0.5*rho_a*C_D*A_a*vx*abs(vx);
F_loss = c_vloss*vx;
F_a = F_aero + F_loss;

%% 13. INITIAL ACCELERATION ESTIMATES FOR LOAD ITERATION
F_rr_0 = c_rr_f*F_zf_0 + c_rr_r*F_zr_0 + c_rr_s*F_zs_0;
a_x_hat = ( ...
      T_f_drive/R_f ...
    - (Tb_fL+Tb_fR+Tb_rL+Tb_rR)/R_f ...
    - (Tb_sL+Tb_sR)/R_s ...
    - F_a ...
    - F_rr_0 )/m_tot;
a_y_t_hat = vx*r_t;
a_y_s_hat = Vx_h_s*r_s;

%% 14. INITIALIZE FINAL-SOLVE VARIABLES
dv_x = 0.0;
dv_y = 0.0;
dr_t = 0.0;
dr_s = 0.0;
H_x  = 0.0;
H_y  = 0.0;
ddphi_t = 0.0;
ddphi_s = 0.0;
% Force variables initialized for code generation
F_zf_L = 0.5*F_zf_0;
F_zf_R = 0.5*F_zf_0;
F_zr_L = 0.5*F_zr_0;
F_zr_R = 0.5*F_zr_0;
F_zs_L = 0.5*F_zs_0;
F_zs_R = 0.5*F_zs_0;
Fx_fL_tire = 0.0; Fy_fL_tire = 0.0;
Fx_fR_tire = 0.0; Fy_fR_tire = 0.0;
Fx_rL_tire = 0.0; Fy_rL_tire = 0.0;
Fx_rR_tire = 0.0; Fy_rR_tire = 0.0;
Fx_sL_tire = 0.0; Fy_sL_tire = 0.0;
Fx_sR_tire = 0.0; Fy_sR_tire = 0.0;
Fx_fL_t = 0.0; Fy_fL_t = 0.0;
Fx_fR_t = 0.0; Fy_fR_t = 0.0;
Fx_rL_t = 0.0; Fy_rL_t = 0.0;
Fx_rR_t = 0.0; Fy_rR_t = 0.0;
Frr_fL = 0.0; Frr_fR = 0.0;
Frr_rL = 0.0; Frr_rR = 0.0;
Frr_sL = 0.0; Frr_sR = 0.0;
Fx_s_s = 0.0; Fy_s_s = 0.0;
Fx_s_t = 0.0; Fy_s_t = 0.0;
Fx_t   = 0.0; Fy_t   = 0.0;
M_z_t  = 0.0; M_z_s  = 0.0;

%% 15. TWO-PASS LOAD / TIRE / 8x8 FULL-ROLL FIFTH-WHEEL SOLVE
for ipass = 1:2
    % ------------------------------------------------------------
    % 15.1 Longitudinal load transfer on tractor
    % ------------------------------------------------------------
    Delta_Fz_long = m_t*h_t/L_t*a_x_hat;
    F_zf_ax = clipval(F_zf_0 - Delta_Fz_long, ...
                      0.02*F_zt_sup,0.98*F_zt_sup);
    F_zr_ax = F_zt_sup - F_zf_ax;
    % ------------------------------------------------------------
    % 15.2 Roll suspension moments / lateral load transfer
    % ------------------------------------------------------------
    M_susp_f = K_phi_f*phi_t + C_phi_f*dphi_t;
    M_susp_r = K_phi_r*phi_t + C_phi_r*dphi_t;
    % Active roll moment acts on body; the opposite road reaction is
    % distributed according to the roll-stiffness shares.
    M_road_f = M_susp_f - eta_phi_f*M_act;
    M_road_r = M_susp_r - eta_phi_r*M_act;
    Delta_Fz_lat_f = M_road_f/T_f;
    Delta_Fz_lat_r = M_road_r/T_r;
    % NEW: semitrailer lateral load transfer is now driven by its
    % independent suspension roll state rather than the previous
    % quasi-static m_s*h_s*a_y_s/T_s approximation.
    M_susp_s = K_phi_s*phi_s + C_phi_s*dphi_s;
    Delta_Fz_lat_s = M_susp_s/T_s;
    Delta_Fz_lat_f = clipval(Delta_Fz_lat_f,-0.5*F_zf_ax,0.5*F_zf_ax);
    Delta_Fz_lat_r = clipval(Delta_Fz_lat_r,-0.5*F_zr_ax,0.5*F_zr_ax);
    Delta_Fz_lat_s = clipval(Delta_Fz_lat_s,-0.5*F_zs_0,0.5*F_zs_0);
    F_zf_L = 0.5*F_zf_ax + Delta_Fz_lat_f;
    F_zf_R = 0.5*F_zf_ax - Delta_Fz_lat_f;
    F_zr_L = 0.5*F_zr_ax + Delta_Fz_lat_r;
    F_zr_R = 0.5*F_zr_ax - Delta_Fz_lat_r;
    F_zs_L = 0.5*F_zs_0 + Delta_Fz_lat_s;
    F_zs_R = 0.5*F_zs_0 - Delta_Fz_lat_s;
    % ------------------------------------------------------------
    % 15.3 Combined-slip Dugoff tire forces
    % ------------------------------------------------------------
    % Tractor front axle: one tire on each side (n=1).
    [Fx_fL_tire,Fy_fL_tire] = trucksim_dugoff( ...
        sigma_fL,alpha_fL,F_zf_L,1.0,mu);
    [Fx_fR_tire,Fy_fR_tire] = trucksim_dugoff( ...
        sigma_fR,alpha_fR,F_zf_R,1.0,mu);
    % Drive axle: dual tires on each side (n=2).
    [Fx_rL_tire,Fy_rL_tire] = trucksim_dugoff( ...
        sigma_rL,alpha_rL,F_zr_L,2.0,mu);
    [Fx_rR_tire,Fy_rR_tire] = trucksim_dugoff( ...
        sigma_rR,alpha_rR,F_zr_R,2.0,mu);
    % Trailer axle: dual tires on each side (n=2).
    [Fx_sL_tire,Fy_sL_tire] = trucksim_dugoff( ...
        sigma_sL,alpha_sL,F_zs_L,2.0,mu);
    [Fx_sR_tire,Fy_sR_tire] = trucksim_dugoff( ...
        sigma_sR,alpha_sR,F_zs_R,2.0,mu);
    % ------------------------------------------------------------
    % 15.4 Rolling resistance as separate body force
    % ------------------------------------------------------------
    Frr_fL = rolling_resistance(c_rr_f,F_zf_L,Vx_fL_w);
    Frr_fR = rolling_resistance(c_rr_f,F_zf_R,Vx_fR_w);
    Frr_rL = rolling_resistance(c_rr_r,F_zr_L,Vx_rL_w);
    Frr_rR = rolling_resistance(c_rr_r,F_zr_R,Vx_rR_w);
    Frr_sL = rolling_resistance(c_rr_s,F_zs_L,Vx_sL_w);
    Frr_sR = rolling_resistance(c_rr_s,F_zs_R,Vx_sR_w);
    % ------------------------------------------------------------
    % 15.5 Front-wheel force transformation to tractor body
    % ------------------------------------------------------------
    Fx_fL_loc = Fx_fL_tire + Frr_fL;
    Fx_fR_loc = Fx_fR_tire + Frr_fR;
    Fx_fL_t = cDel*Fx_fL_loc - sDel*Fy_fL_tire;
    Fy_fL_t = sDel*Fx_fL_loc + cDel*Fy_fL_tire;
    Fx_fR_t = cDel*Fx_fR_loc - sDel*Fy_fR_tire;
    Fy_fR_t = sDel*Fx_fR_loc + cDel*Fy_fR_tire;
    % Tractor rear wheels are aligned with tractor body
    Fx_rL_t = Fx_rL_tire + Frr_rL;
    Fy_rL_t = Fy_rL_tire;
    Fx_rR_t = Fx_rR_tire + Frr_rR;
    Fy_rR_t = Fy_rR_tire;
    % ------------------------------------------------------------
    % 15.6 Semitrailer resultant in semitrailer frame, then tractor frame
    % ------------------------------------------------------------
    Fx_sL_loc = Fx_sL_tire + Frr_sL;
    Fx_sR_loc = Fx_sR_tire + Frr_sR;
    Fx_s_s = Fx_sL_loc + Fx_sR_loc;
    Fy_s_s = Fy_sL_tire + Fy_sR_tire;
    Fx_s_t =  cg*Fx_s_s + sg*Fy_s_s;
    Fy_s_t = -sg*Fx_s_s + cg*Fy_s_s;
    % ------------------------------------------------------------
    % 15.7 Tractor resultants excluding fifth-wheel reaction
    % ------------------------------------------------------------
    Fx_t = Fx_fL_t + Fx_fR_t + Fx_rL_t + Fx_rR_t - F_a;
    Fy_t = Fy_fL_t + Fy_fR_t + Fy_rL_t + Fy_rR_t;
    % ------------------------------------------------------------
    % 15.8 Yaw moments from wheel forces
    % ------------------------------------------------------------
    M_z_t = ...
          L_ft*Fy_fL_t - 0.5*T_f*Fx_fL_t ...
        + L_ft*Fy_fR_t + 0.5*T_f*Fx_fR_t ...
        - L_rt*Fy_rL_t - 0.5*T_r*Fx_rL_t ...
        - L_rt*Fy_rR_t + 0.5*T_r*Fx_rR_t;
    M_z_s = ...
         -L_sa*Fy_sL_tire - 0.5*T_s*Fx_sL_loc ...
         -L_sa*Fy_sR_tire + 0.5*T_s*Fx_sR_loc;
    % ------------------------------------------------------------
    % 15.9 8x8 DIRECT TRACTOR + SEMITRAILER ROLL-COUPLED SOLVE
    %
    % Unknown:
    %
    % q = [dvx dvy dr_t dr_s Hx Hy ddphi_t ddphi_s]'
    %
    % ------------------------------------------------------------
    A = zeros(8,8);
    b = zeros(8,1);
    % --------------------------------------------------------
    % Eq. (1R): tractor longitudinal + tractor roll-rate coupling
    %
    % mt*(dvx-r_t*vy) - mst*ht*dphi_t*r_t = Fx_t + Hx
    % --------------------------------------------------------
    A(1,1) = m_t;
    A(1,5) = -1.0;
    b(1) = Fx_t ...
         + m_t*r_t*vy ...
         + m_st*h_t*dphi_t*r_t;
    % --------------------------------------------------------
    % Eq. (2R): tractor lateral + tractor roll acceleration
    %
    % mt*(dvy+r_t*vx) + mst*ht*ddphi_t = Fy_t + Hy
    % --------------------------------------------------------
    A(2,2) = m_t;
    A(2,6) = -1.0;
    A(2,7) = m_st*h_t;
    b(2) = Fy_t - m_t*r_t*vx;
    % --------------------------------------------------------
    % Eq. (3R): tractor yaw + tractor roll-yaw product inertia
    %
    % Izt*dr_t - Ixzt*ddphi_t = Mz_t - Lwt*Hy
    % --------------------------------------------------------
    A(3,3) = I_zt;
    A(3,6) = L_wt;
    A(3,7) = -I_xzt;
    b(3) = M_z_t;
    % --------------------------------------------------------
    % Semitrailer CG acceleration in tractor frame
    %
    % a_sx^t = dvx-r_t*vy + Lwt*r_t^2
    %          + Lcs*(r_s^2*cos(gamma)-dr_s*sin(gamma))
    %
    % a_sy^t = dvy+r_t*vx - Lwt*dr_t
    %          - Lcs*(r_s^2*sin(gamma)+dr_s*cos(gamma))
    %
    % Trailer sprung-mass roll inertial vector in trailer frame:
    %
    % c_sx = -mss*hs*dphi_s*r_s
    % c_sy =  mss*hs*ddphi_s
    %
    % Rotated to tractor frame:
    %
    % c_tx = -mss*hs*cos(gamma)*dphi_s*r_s
    %        +mss*hs*sin(gamma)*ddphi_s
    %
    % c_ty =  mss*hs*sin(gamma)*dphi_s*r_s
    %        +mss*hs*cos(gamma)*ddphi_s
    % --------------------------------------------------------
    % Semitrailer longitudinal force balance in tractor frame:
    %
    % ms*a_sx^t + c_tx = Fx_s^t - Hx
    A(4,1) = m_s;
    A(4,4) = -m_s*L_cs*sg;
    A(4,5) = 1.0;
    A(4,8) = m_ss*h_s*sg;
    b(4) = Fx_s_t ...
         + m_s*r_t*vy ...
         - m_s*L_wt*r_t*r_t ...
         - m_s*L_cs*r_s*r_s*cg ...
         + m_ss*h_s*cg*dphi_s*r_s;
    % Semitrailer lateral force balance in tractor frame:
    %
    % ms*a_sy^t + c_ty = Fy_s^t - Hy
    A(5,2) = m_s;
    A(5,3) = -m_s*L_wt;
    A(5,4) = -m_s*L_cs*cg;
    A(5,6) = 1.0;
    A(5,8) = m_ss*h_s*cg;
    b(5) = Fy_s_t ...
         - m_s*r_t*vx ...
         + m_s*L_cs*r_s*r_s*sg ...
         - m_ss*h_s*sg*dphi_s*r_s;
    % --------------------------------------------------------
    % Semitrailer yaw + trailer roll-yaw product inertia
    %
    % Izs*dr_s - Ixzs*ddphi_s
    %   = Mz_s - Lcs*(sin(gamma)*Hx + cos(gamma)*Hy)
    % --------------------------------------------------------
    A(6,4) = I_zs;
    A(6,5) = L_cs*sg;
    A(6,6) = L_cs*cg;
    A(6,8) = -I_xzs;
    b(6) = M_z_s;
    % --------------------------------------------------------
    % Relative tractor-semitrailer roll coupling
    %
    % Positive M_hitch_roll acts to reduce phi_t-phi_s:
    % tractor receives -M_hitch_roll
    % trailer receives +M_hitch_roll
    % --------------------------------------------------------
    if use_hitch_Mx_table > 0.5
        M_hitch_roll = hitch_roll_moment_table(phi_t,phi_s);
    else
        M_hitch_roll = 0.0;
    end
    % --------------------------------------------------------
    % Tractor roll equation
    %
    % Ixx_t*ddphi_t - Ixzt*dr_t
    % + mst*ht*(dvy+r_t*vx)
    %
    % = mst*g*ht*sin(phi_t)
    %   - Kphi_t*phi_t - Cphi_t*dphi_t
    %   - M_hitch_roll + Mact
    % --------------------------------------------------------
    A(7,2) = m_st*h_t;
    A(7,3) = -I_xzt;
    A(7,7) = I_xxt;
    b(7) = ...
          m_st*g*h_t*sin(phi_t) ...
        - K_phi*phi_t ...
        - C_phi*dphi_t ...
        - M_hitch_roll ...
        + M_act ...
        - m_st*h_t*r_t*vx;
    % --------------------------------------------------------
    % Semitrailer roll equation
    %
    % Ixx_s*ddphi_s - Ixzs*dr_s + mss*hs*a_sy_s
    %
    % = mss*g*hs*sin(phi_s)
    %   - Kphi_s*phi_s - Cphi_s*dphi_s
    %   + M_hitch_roll
    %
    % where
    %
    % a_sy_s = sin(gamma)*a_sx_t + cos(gamma)*a_sy_t
    %
    % Its unknown-dependent part is
    %
    %   sg*dvx + cg*dvy - cg*Lwt*dr_t - Lcs*dr_s
    %
    % and its known part is
    %
    %   r_t*(cg*vx-sg*vy) + sg*Lwt*r_t^2
    % --------------------------------------------------------
    A(8,1) =  m_ss*h_s*sg;
    A(8,2) =  m_ss*h_s*cg;
    A(8,3) = -m_ss*h_s*cg*L_wt;
    A(8,4) = -I_xzs - m_ss*h_s*L_cs;
    A(8,8) =  I_xxs;
    a_sy_s_known = ...
          r_t*(cg*vx - sg*vy) ...
        + sg*L_wt*r_t*r_t;
    b(8) = ...
          m_ss*g*h_s*sin(phi_s) ...
        - K_phi_s*phi_s ...
        - C_phi_s*dphi_s ...
        + M_hitch_roll ...
        - m_ss*h_s*a_sy_s_known;
    % Fixed-size simultaneous solve.
    q = A\b;
    dv_x      = q(1);
    dv_y      = q(2);
    dr_t      = q(3);
    dr_s      = q(4);
    H_x       = q(5);
    H_y       = q(6);
    ddphi_t   = q(7);
    ddphi_s   = q(8);
    % ------------------------------------------------------------
    % 15.10 Physical acceleration update for next fixed-point pass
    % ------------------------------------------------------------
    a_x_new   = dv_x - r_t*vy;
    a_y_t_new = dv_y + r_t*vx;
    a_sx_t = ...
          dv_x - r_t*vy + L_wt*r_t*r_t ...
        + L_cs*(r_s*r_s*cg - dr_s*sg);
    a_sy_t = ...
          dv_y + r_t*vx - L_wt*dr_t ...
        - L_cs*(r_s*r_s*sg + dr_s*cg);
    % tractor -> semitrailer frame
    a_sx_s = cg*a_sx_t - sg*a_sy_t;
    a_sy_s = sg*a_sx_t + cg*a_sy_t;
    a_x_hat   = 0.5*a_x_hat   + 0.5*a_x_new;
    a_y_t_hat = 0.5*a_y_t_hat + 0.5*a_y_t_new;
    a_y_s_hat = 0.5*a_y_s_hat + 0.5*a_sy_s;
end

%% 15.11 ROLL ACCELERATIONS
% ddphi_t and ddphi_s are already obtained directly from the 8x8 solve above.

%% 16. ARTICULATION / GLOBAL KINEMATICS
dgamma = r_t - r_s;
dpsi   = r_t;
dY = vx*sin(psi) + vy*cos(psi);
dX = vx*cos(psi) - vy*sin(psi);

%% 17. WHEEL ROTATIONAL DYNAMICS
% Only Dugoff tire-contact Fx enters the wheel reaction torque.
% Rolling resistance remains a separate body force.
domega_fL = ( -Tb_fL - Fx_fL_tire*R_f )/I_wf;
domega_fR = ( -Tb_fR - Fx_fR_tire*R_f )/I_wf;
domega_rL = ( T_f_drive/2.0 - Tb_rL - Fx_rL_tire*R_r )/I_wr;
domega_rR = ( T_f_drive/2.0 - Tb_rR - Fx_rR_tire*R_r )/I_wr;
domega_sL = ( -Tb_sL - Fx_sL_tire*R_s )/I_ws;
domega_sR = ( -Tb_sR - Fx_sR_tire*R_s )/I_ws;

%% 18. STATE DERIVATIVE
dxdt = zeros(18,1);
dxdt(1)  = dv_x;
dxdt(2)  = dv_y;
dxdt(3)  = dr_t;
dxdt(4)  = dr_s;
dxdt(5)  = dgamma;
dxdt(6)  = dpsi;
dxdt(7)  = dY;
dxdt(8)  = dX;
dxdt(9)  = domega_fL;
dxdt(10) = domega_fR;
dxdt(11) = domega_rL;
dxdt(12) = domega_rR;
dxdt(13) = domega_sL;
dxdt(14) = domega_sR;
dxdt(15) = dphi_t;
dxdt(16) = ddphi_t;
dxdt(17) = dphi_s;
dxdt(18) = ddphi_s;

%% 19. PHYSICAL LTR
LTR_f = (F_zf_L-F_zf_R)/max(F_zf_L+F_zf_R,1.0);
LTR_r = (F_zr_L-F_zr_R)/max(F_zr_L+F_zr_R,1.0);
LTR_s = (F_zs_L-F_zs_R)/max(F_zs_L+F_zs_R,1.0); %#ok<NASGU>

%% 20. OUTPUTS
% Outputs y(1:44) retain the previous ordering.
% NEW diagnostics are appended:
%   y(45) = phi_s
%   y(46) = dphi_s
%   y(47) = LTR_s
%   y(48) = semitrailer lateral acceleration a_y,s [m/s^2]
%
% H_y (output 8) is the final fifth-wheel lateral force [N] acting on
% the tractor, expressed along the tractor-frame y axis.
%
% y(34:35) report tractor physical accelerations from the final 8x8 solve.
y = zeros(48,1);
y(1)  = vx;
y(2)  = r_t;
y(3)  = Y;
y(4)  = gamma;
y(5)  = phi_t;
y(6)  = LTR_f;
y(7)  = LTR_r;
y(8)  = vy;
y(9)  = r_s;
y(10) = psi;
y(11) = X;
% wheel vertical loads
y(12) = F_zf_L;
y(13) = F_zf_R;
y(14) = F_zr_L;
y(15) = F_zr_R;
y(16) = F_zs_L;
y(17) = F_zs_R;
% wheel longitudinal forces (Dugoff tire-contact forces)
y(18) = Fx_fL_tire;
y(19) = Fx_fR_tire;
y(20) = Fx_rL_tire;
y(21) = Fx_rR_tire;
y(22) = Fx_sL_tire;
y(23) = Fx_sR_tire;
% wheel lateral forces
y(24) = Fy_fL_tire;
y(25) = Fy_fR_tire;
y(26) = Fy_rL_tire;
y(27) = Fy_rR_tire;
y(28) = Fy_sL_tire;
y(29) = Fy_sR_tire;
% total body-force resultants excluding fifth wheel
y(30) = Fx_t;
y(31) = Fy_t;
% wheel-force yaw moments
y(32) = M_z_t;
y(33) = M_z_s;
% physical coupled-solve tractor accelerations
y(34) = dv_x - r_t*vy;
y(35) = dv_y + r_t*vx;
% wheel slip ratios
y(36) = sigma_fL;
y(37) = sigma_fR;
y(38) = sigma_rL;
y(39) = sigma_rR;
y(40) = sigma_sL;
y(41) = sigma_sR;
% representative axle slip angles
y(42) = 0.5*(alpha_fL+alpha_fR);
y(43) = 0.5*(alpha_rL+alpha_rR);
y(44) = 0.5*(alpha_sL+alpha_sR);
% NEW independent semitrailer-roll diagnostics
y(45) = phi_s;
y(46) = dphi_s;
y(47) = LTR_s;
y(48) = a_sy_s;
end

% =========================================================================

% LOCAL FUNCTIONS
% =========================================================================

function z = clipval(x,lo,hi)
z = min(max(x,lo),hi);
end

function Frr = rolling_resistance(c_rr,Fz,Vx)
if abs(Vx) >= 1.0e-4
    % TruckSim Tire: RR_C=0.0041, RR_V=0.0000256
    Frr = -(c_rr + 0.0000256*abs(Vx))*Fz*sign(Vx);
else
    Frr = 0.0;
end

end

function [Fx,Fy] = dugoff_force(sigma,alpha,Fz,Csigma,Calpha,mu)
% Combined-slip Dugoff relation used by the corrected platoon model.
if Fz <= 0.0
    Fx = 0.0;
    Fy = 0.0;
    return;
end

D = sqrt((Csigma*sigma)^2 + (Calpha*tan(alpha))^2);
if D <= 1.0e-12
    Fx = 0.0;
    Fy = 0.0;
    return;
end

den = max(1.0 + sigma,0.05);
lambda = mu*Fz*den/(2.0*D);
if lambda < 1.0
    f_lambda = lambda*(2.0-lambda);
else
    f_lambda = 1.0;
end

Fx0 = Csigma*sigma/den*f_lambda;
Fy0 = Calpha*tan(alpha)/den*f_lambda;
Fmag = sqrt(Fx0*Fx0 + Fy0*Fy0);
Flim = mu*Fz;
if Fmag > Flim && Fmag > 1.0e-12
    scale = Flim/Fmag;
    Fx = Fx0*scale;
    Fy = Fy0*scale;
else
    Fx = Fx0;
    Fy = Fy0;
end

end


function [Fx,Fy] = trucksim_dugoff(sigma,alpha,Fz_side,n_tires,mu)
% Reduced tire-equivalent model derived from TruckSim 2016 .par tables.
% Front n=1, drive/trailer n=2; equal load assumption across dual tires.
Fz_one = max(Fz_side/n_tires,0.0);
C_alpha_one = tire_stiffness_from_table(Fz_one,1.0);
C_sigma_one = tire_stiffness_from_table(Fz_one,2.0);
[Fx_one,Fy_one] = dugoff_force( ...
    sigma,alpha,Fz_one,C_sigma_one,C_alpha_one,mu);
Fx = n_tires*Fx_one;
Fy = n_tires*Fy_one;
end

function C = tire_stiffness_from_table(Fz,kind)
% FZ columns and first nonzero row of expanded TruckSim tire carpets.
% FZ in N, lateral angle in DEG, longitudinal slip dimensionless.
Fz_tab = [0.0 7357.5 14715.0 29430.0 44145.0 58860.0];
if kind < 1.5
    Q_tab = [0.0 1253.8 2418.0 4477.9 6179.4 7612.3];
    divisor = pi/180.0;
else
    Q_tab = [0.0 2276.6 4390.6 8130.8 11220.5 13822.4];
    divisor = 0.025;
end
zz = max(Fz,0.0);
if zz >= Fz_tab(6)
    C = Q_tab(6)/divisor;
    return;
end
C = 0.0;
for kk=1:5
    if zz >= Fz_tab(kk) && zz < Fz_tab(kk+1)
        lambda = (zz-Fz_tab(kk))/(Fz_tab(kk+1)-Fz_tab(kk));
        C = (Q_tab(kk)+lambda*(Q_tab(kk+1)-Q_tab(kk)))/divisor;
        return;
    end
end
end

function M = hitch_roll_moment_table(phi_t,phi_s)
% TruckSim MX_HITCH_TABLE LINEAR, angle in deg, moment in Nm.
% [-2,-100k],[-1.5,-50k],[-0.5,-50k],[0.5,50k],[1.5,50k],[2,100k].
phi_d = (phi_t-phi_s)*180.0/pi;
a = [-2.0 -1.5 -0.5 0.5 1.5 2.0];
v = [-100000.0 -50000.0 -50000.0 50000.0 50000.0 100000.0];
if phi_d <= a(1)
    M = v(1);
    return;
end
if phi_d >= a(6)
    M = v(6);
    return;
end
M = 0.0;
for kk=1:5
    if phi_d >= a(kk) && phi_d < a(kk+1)
        tau = (phi_d-a(kk))/(a(kk+1)-a(kk));
        M = v(kk)+tau*(v(kk+1)-v(kk));
        return;
    end
end
end
