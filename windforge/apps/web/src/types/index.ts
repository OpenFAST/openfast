// ─── Auth ────────────────────────────────────────────────────────────────────

export interface User {
  id: string;
  email: string;
  full_name: string;
  organization_name?: string;
  is_active: boolean;
  created_at: string;
}

export interface Token {
  access_token: string;
  token_type: string;
}

export interface LoginRequest {
  username: string; // OAuth2 uses 'username' for email
  password: string;
}

export interface RegisterRequest {
  email: string;
  password: string;
  full_name: string;
  organization_name?: string;
}

// ─── Project ─────────────────────────────────────────────────────────────────

export type WindClass = 'I' | 'II' | 'III' | 'S';
export type TurbulenceClass = 'A' | 'B' | 'C';

export interface Project {
  id: string;
  name: string;
  description?: string;
  wind_class: WindClass;
  turbulence_class: TurbulenceClass;
  rated_power_kw: number;
  rotor_diameter_m: number;
  hub_height_m: number;
  cut_in_speed: number;
  rated_speed: number;
  cut_out_speed: number;
  dt: number;
  t_max: number;
  status: string;
  owner_id: string;
  created_at: string;
  updated_at: string;
}

export interface ProjectCreate {
  name: string;
  description?: string;
  wind_class: WindClass;
  turbulence_class: TurbulenceClass;
  rated_power_kw: number;
  rotor_diameter_m: number;
  hub_height_m: number;
  cut_in_speed: number;
  rated_speed: number;
  cut_out_speed: number;
  dt: number;
  t_max: number;
}

// ─── Tower ───────────────────────────────────────────────────────────────────

export interface TowerStation {
  height_m: number;
  outer_diameter_m: number;
  wall_thickness_m: number;
  elastic_modulus_pa: number;
  shear_modulus_pa: number;
  density_kg_m3: number;
}

export interface Tower {
  id: string;
  project_id: string;
  name: string;
  total_height_m: number;
  base_diameter_m: number;
  top_diameter_m: number;
  num_stations: number;
  stations: TowerStation[];
  created_at: string;
  updated_at: string;
}

export interface TowerCreate {
  name: string;
  total_height_m: number;
  base_diameter_m: number;
  top_diameter_m: number;
  num_stations: number;
  stations: TowerStation[];
}

// ─── Blade ───────────────────────────────────────────────────────────────────

export interface BladeStructuralStation {
  span_fraction: number;
  mass_density_kg_m: number;
  flap_stiffness_nm2: number;
  edge_stiffness_nm2: number;
  torsional_stiffness_nm2: number;
  ea_n: number;
}

export interface BladeAeroStation {
  span_fraction: number;
  chord_m: number;
  twist_deg: number;
  airfoil_name: string;
  thickness_ratio: number;
}

export interface Blade {
  id: string;
  project_id: string;
  name: string;
  length_m: number;
  num_structural_stations: number;
  num_aero_stations: number;
  structural_stations: BladeStructuralStation[];
  aero_stations: BladeAeroStation[];
  created_at: string;
  updated_at: string;
}

export interface BladeCreate {
  name: string;
  length_m: number;
  num_structural_stations: number;
  num_aero_stations: number;
  structural_stations: BladeStructuralStation[];
  aero_stations: BladeAeroStation[];
}

// ─── Controller ──────────────────────────────────────────────────────────────

export interface Controller {
  id: string;
  project_id: string;
  name: string;
  control_mode: string;
  rated_power_kw: number;
  rated_speed_rpm: number;
  min_pitch_deg: number;
  max_pitch_deg: number;
  max_pitch_rate_deg_s: number;
  kp_pitch: number;
  ki_pitch: number;
  kp_torque: number;
  ki_torque: number;
  min_gen_speed_rpm: number;
  max_gen_speed_rpm: number;
  max_torque_nm: number;
  created_at: string;
  updated_at: string;
}

export interface ControllerCreate {
  name: string;
  control_mode: string;
  rated_power_kw: number;
  rated_speed_rpm: number;
  min_pitch_deg: number;
  max_pitch_deg: number;
  max_pitch_rate_deg_s: number;
  kp_pitch: number;
  ki_pitch: number;
  kp_torque: number;
  ki_torque: number;
  min_gen_speed_rpm: number;
  max_gen_speed_rpm: number;
  max_torque_nm: number;
}

// ─── Airfoil ─────────────────────────────────────────────────────────────────

export interface AirfoilPolar {
  alpha_deg: number[];
  cl: number[];
  cd: number[];
  cm: number[];
}

export interface Airfoil {
  id: string;
  name: string;
  thickness_ratio: number;
  polars: AirfoilPolar[];
}

// ─── Turbine Model ───────────────────────────────────────────────────────────

export interface TurbineModel {
  id: string;
  project_id: string;
  name: string;
  tower_id: string;
  blade_id: string;
  controller_id: string;
  num_blades: number;
  hub_mass_kg: number;
  hub_inertia_kg_m2: number;
  nacelle_mass_kg: number;
  nacelle_cm_x_m: number;
  nacelle_cm_z_m: number;
  shaft_tilt_deg: number;
  precone_deg: number;
  overhang_m: number;
  twr2shft_m: number;
  gearbox_ratio: number;
  gen_inertia_kg_m2: number;
  drivetrain_efficiency: number;
  created_at: string;
  updated_at: string;
}

export interface TurbineModelCreate {
  name: string;
  tower_id: string;
  blade_id: string;
  controller_id: string;
  num_blades: number;
  hub_mass_kg: number;
  hub_inertia_kg_m2: number;
  nacelle_mass_kg: number;
  nacelle_cm_x_m: number;
  nacelle_cm_z_m: number;
  shaft_tilt_deg: number;
  precone_deg: number;
  overhang_m: number;
  twr2shft_m: number;
  gearbox_ratio: number;
  gen_inertia_kg_m2: number;
  drivetrain_efficiency: number;
}

// ─── DLC ─────────────────────────────────────────────────────────────────────

export interface TurbSimParams {
  grid_height: number;
  grid_width: number;
  time_step: number;
  analysis_time: number;
  iec_turb_model: string;
  iec_standard: string;
  iec_turbc: string;
  wind_profile_type: string;
  ref_ht: number;
  pllj_height?: number;
}

export interface DLCCaseSpec {
  dlc_number: string;
  wind_speed_mps: number;
  yaw_error_deg: number;
  num_seeds: number;
  turbsim_params: TurbSimParams;
  fault_time_s?: number;
  fault_type?: string;
}

export interface DLCDefinition {
  id: string;
  project_id: string;
  name: string;
  iec_standard: string;
  cases: DLCCaseSpec[];
  created_at: string;
  updated_at: string;
}

export interface DLCDefinitionCreate {
  name: string;
  iec_standard: string;
  cases: DLCCaseSpec[];
}

// ─── Simulation ──────────────────────────────────────────────────────────────

export type SimulationStatus =
  | 'pending'
  | 'running'
  | 'completed'
  | 'failed'
  | 'cancelled';

export interface SimulationCase {
  id: string;
  simulation_id: string;
  dlc_number: string;
  wind_speed_mps: number;
  seed: number;
  yaw_error_deg: number;
  status: SimulationStatus;
  progress: number;
  error_message?: string;
  started_at?: string;
  completed_at?: string;
}

export interface Simulation {
  id: string;
  project_id: string;
  name: string;
  turbine_model_id: string;
  dlc_definition_id: string;
  status: SimulationStatus;
  total_cases: number;
  completed_cases: number;
  failed_cases: number;
  created_at: string;
  updated_at: string;
}

export interface SimulationWithProgress extends Simulation {
  cases: SimulationCase[];
  overall_progress: number;
}

// ─── Results ─────────────────────────────────────────────────────────────────

export interface ChannelStatistic {
  channel: string;
  mean: number;
  std: number;
  min: number;
  max: number;
  abs_max: number;
}

export interface ResultsStatistics {
  simulation_id: string;
  case_id: string;
  dlc_number: string;
  wind_speed_mps: number;
  statistics: ChannelStatistic[];
}

export interface ResultsDEL {
  simulation_id: string;
  channel: string;
  dlc_number: string;
  wohler_exponent: number;
  del_values: number[];
  lifetime_del: number;
}

export interface ResultsExtreme {
  simulation_id: string;
  channel: string;
  dlc_number: string;
  extreme_value: number;
  time_of_extreme: number;
  case_id: string;
}

// ─── WebSocket Messages ──────────────────────────────────────────────────────

export interface WSCaseProgress {
  type: 'case_progress';
  case_id: string;
  progress: number;
}

export interface WSCaseComplete {
  type: 'case_complete';
  case_id: string;
}

export interface WSCaseError {
  type: 'case_error';
  case_id: string;
  error: string;
}

export interface WSSimulationComplete {
  type: 'simulation_complete';
  simulation_id: string;
}

export interface WSLiveData {
  type: 'live_data';
  case_id: string;
  time: number;
  channels: Record<string, number>;
}

export type WSMessage =
  | WSCaseProgress
  | WSCaseComplete
  | WSCaseError
  | WSSimulationComplete
  | WSLiveData;
