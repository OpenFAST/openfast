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
  frac: number;
  mass_den: number;
  fa_stiff: number;
  ss_stiff: number;
  outer_diameter: number;
  wall_thickness: number;
}

export interface Tower {
  id: string;
  project_id: string;
  name: string;
  version: number;
  tower_height: number;
  tower_base_height: number;
  tower_fa_damping_1: number;
  tower_fa_damping_2: number;
  tower_ss_damping_1: number;
  tower_ss_damping_2: number;
  stations: TowerStation[] | null;
  fa_mode_1_coeffs: number[] | null;
  fa_mode_2_coeffs: number[] | null;
  ss_mode_1_coeffs: number[] | null;
  ss_mode_2_coeffs: number[] | null;
  is_active: boolean;
  created_at: string;
}

export interface TowerCreate {
  name: string;
  tower_height: number;
  tower_base_height?: number;
  tower_fa_damping_1?: number;
  tower_fa_damping_2?: number;
  tower_ss_damping_1?: number;
  tower_ss_damping_2?: number;
  stations?: TowerStation[];
  fa_mode_1_coeffs?: number[];
  fa_mode_2_coeffs?: number[];
  ss_mode_1_coeffs?: number[];
  ss_mode_2_coeffs?: number[];
}

// ─── Blade ───────────────────────────────────────────────────────────────────

export interface BladeStructuralStation {
  frac: number;
  pitch_axis: number;
  struct_twist: number;
  mass_den: number;
  flap_stiff: number;
  edge_stiff: number;
}

export interface BladeAeroStation {
  frac: number;
  chord: number;
  aero_twist: number;
  airfoil_id: string;
  aero_center: number;
}

export interface Blade {
  id: string;
  project_id: string;
  name: string;
  version: number;
  blade_length: number;
  structural_stations: BladeStructuralStation[] | null;
  aero_stations: BladeAeroStation[] | null;
  flap_mode_1_coeffs: number[] | null;
  flap_mode_2_coeffs: number[] | null;
  edge_mode_1_coeffs: number[] | null;
  blade_flap_damping: number;
  blade_edge_damping: number;
  is_active: boolean;
  created_at: string;
}

export interface BladeCreate {
  name: string;
  blade_length: number;
  blade_flap_damping?: number;
  blade_edge_damping?: number;
  structural_stations?: BladeStructuralStation[];
  aero_stations?: BladeAeroStation[];
  flap_mode_1_coeffs?: number[];
  flap_mode_2_coeffs?: number[];
  edge_mode_1_coeffs?: number[];
}

// ─── Controller ──────────────────────────────────────────────────────────────

export interface Controller {
  id: string;
  project_id: string;
  name: string;
  version: number;
  controller_type: string;
  pcmode: number;
  vscontrl: number;
  parameters: Record<string, any> | null;
  dll_filename: string | null;
  dll_procname: string | null;
  is_active: boolean;
  created_at: string;
}

export interface ControllerCreate {
  name: string;
  controller_type: string;
  pcmode?: number;
  vscontrl?: number;
  parameters?: Record<string, any>;
  dll_filename?: string;
  dll_procname?: string;
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
  version: number;
  tower_id: string | null;
  blade_id: string | null;
  controller_id: string | null;
  gearbox_ratio: number | null;
  generator_inertia: number | null;
  drivetrain_stiffness: number | null;
  drivetrain_damping: number | null;
  hub_mass: number | null;
  hub_inertia: number | null;
  nacelle_mass: number | null;
  nacelle_inertia: number | null;
  overhang: number | null;
  shaft_tilt: number | null;
  precone: number | null;
  rotor_speed_rated: number | null;
  dof_flags: Record<string, boolean> | null;
  is_active: boolean;
  created_at: string;
}

export interface TurbineModelCreate {
  name: string;
  tower_id?: string;
  blade_id?: string;
  controller_id?: string;
  gearbox_ratio?: number;
  generator_inertia?: number;
  drivetrain_stiffness?: number;
  drivetrain_damping?: number;
  hub_mass?: number;
  hub_inertia?: number;
  nacelle_mass?: number;
  nacelle_inertia?: number;
  overhang?: number;
  shaft_tilt?: number;
  precone?: number;
  rotor_speed_rated?: number;
  dof_flags?: Record<string, boolean>;
}

// ─── DLC ─────────────────────────────────────────────────────────────────────

export interface TurbSimParams {
  turbulence_model: string;
  iec_standard: string;
  iec_turbc: string;
  grid_height: number;
  grid_width: number;
  num_grid_z: number;
  num_grid_y: number;
  time_step: number;
  analysis_time: number;
  ref_height: number;
}

export interface DLCCaseSpec {
  dlc_number: string;
  wind_speeds: number[];
  seeds: number;
  yaw_misalignments: number[];
}

export interface DLCDefinition {
  id: string;
  project_id: string;
  turbine_model_id: string;
  name: string;
  dlc_cases: DLCCaseSpec[] | null;
  turbsim_params: TurbSimParams | null;
  total_case_count: number;
  status: string;
  created_at: string;
}

export interface DLCDefinitionCreate {
  name: string;
  turbine_model_id: string;
  dlc_cases?: DLCCaseSpec[];
  turbsim_params?: TurbSimParams;
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
  wind_speed: number;
  seed_number: number;
  yaw_misalignment: number;
  status: string;
  progress_percent: number;
  current_time: number;
  error_message: string | null;
  wall_time_seconds: number | null;
}

export interface Simulation {
  id: string;
  project_id: string;
  dlc_definition_id: string;
  turbine_model_id: string;
  name: string;
  status: string;
  total_cases: number;
  completed_cases: number;
  failed_cases: number;
  agent_id: string | null;
  started_at: string | null;
  completed_at: string | null;
  created_at: string;
}

export interface SimulationWithProgress extends Simulation {
  cases: SimulationCase[];
  overall_progress: number;
}

// ─── Results ─────────────────────────────────────────────────────────────────

export interface ChannelStatistic {
  min: number;
  max: number;
  mean: number;
  std: number;
}

export interface ResultsStatistics {
  id: string;
  simulation_case_id: string;
  dlc_number: string;
  wind_speed: number;
  channel_statistics: Record<string, ChannelStatistic> | null;
}

export interface ResultsDEL {
  id: string;
  simulation_id: string;
  del_results: Record<string, number> | null;
  m_exponent: number;
  n_equivalent: number;
}

export interface ResultsExtreme {
  id: string;
  simulation_id: string;
  extreme_loads: Record<
    string,
    { max: number; min: number; safety_factor: number; design_value: number }
  > | null;
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
