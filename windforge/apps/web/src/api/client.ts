import axios from 'axios';
import type {
  User,
  Token,
  RegisterRequest,
  Project,
  ProjectCreate,
  Tower,
  TowerCreate,
  Blade,
  BladeCreate,
  Controller,
  ControllerCreate,
  TurbineModel,
  TurbineModelCreate,
  DLCDefinition,
  DLCDefinitionCreate,
  Simulation,
  SimulationCase,
  ResultsStatistics,
  ResultsDEL,
  ResultsExtreme,
} from '@/types';

export type { SimulationCase } from '@/types';

// ─── Axios Instance ──────────────────────────────────────────────────────────

const api = axios.create({
  baseURL: '/api/v1',
  headers: {
    'Content-Type': 'application/json',
  },
});

// Request interceptor: attach Bearer token
api.interceptors.request.use((config) => {
  const token = localStorage.getItem('windforge_token');
  if (token) {
    config.headers.Authorization = `Bearer ${token}`;
  }
  return config;
});

// Response interceptor: handle 401
api.interceptors.response.use(
  (response) => response,
  (error) => {
    if (error.response?.status === 401) {
      localStorage.removeItem('windforge_token');
      if (window.location.pathname !== '/login') {
        window.location.href = '/login';
      }
    }
    return Promise.reject(error);
  },
);

// ─── Auth API ────────────────────────────────────────────────────────────────

export const authApi = {
  register: async (data: RegisterRequest): Promise<User> => {
    const res = await api.post<User>('/auth/register', data);
    return res.data;
  },

  login: async (email: string, password: string): Promise<Token> => {
    const res = await api.post<Token>('/auth/login', { email, password });
    return res.data;
  },

  getMe: async (): Promise<User> => {
    const res = await api.get<User>('/auth/me');
    return res.data;
  },
};

// ─── Projects API ────────────────────────────────────────────────────────────

export const projectsApi = {
  list: async (): Promise<Project[]> => {
    const res = await api.get<Project[]>('/projects');
    return res.data;
  },

  create: async (data: ProjectCreate): Promise<Project> => {
    const res = await api.post<Project>('/projects', data);
    return res.data;
  },

  get: async (projectId: string): Promise<Project> => {
    const res = await api.get<Project>(`/projects/${projectId}`);
    return res.data;
  },

  update: async (
    projectId: string,
    data: Partial<ProjectCreate>,
  ): Promise<Project> => {
    const res = await api.put<Project>(`/projects/${projectId}`, data);
    return res.data;
  },

  delete: async (projectId: string): Promise<void> => {
    await api.delete(`/projects/${projectId}`);
  },
};

// ─── Towers API ──────────────────────────────────────────────────────────────

export const towersApi = {
  list: async (projectId: string): Promise<Tower[]> => {
    const res = await api.get<Tower[]>(`/projects/${projectId}/towers`);
    return res.data;
  },

  create: async (projectId: string, data: TowerCreate): Promise<Tower> => {
    const res = await api.post<Tower>(
      `/projects/${projectId}/towers`,
      data,
    );
    return res.data;
  },

  get: async (projectId: string, towerId: string): Promise<Tower> => {
    const res = await api.get<Tower>(
      `/projects/${projectId}/towers/${towerId}`,
    );
    return res.data;
  },

  update: async (
    projectId: string,
    towerId: string,
    data: Partial<TowerCreate>,
  ): Promise<Tower> => {
    const res = await api.put<Tower>(
      `/projects/${projectId}/towers/${towerId}`,
      data,
    );
    return res.data;
  },

  delete: async (projectId: string, towerId: string): Promise<void> => {
    await api.delete(`/projects/${projectId}/towers/${towerId}`);
  },

  preview: async (
    projectId: string,
    data: TowerCreate,
  ): Promise<{ heights: number[]; diameters: number[] }> => {
    const res = await api.post(
      `/projects/${projectId}/towers/preview`,
      data,
    );
    return res.data;
  },
};

// ─── Blades API ──────────────────────────────────────────────────────────────

export const bladesApi = {
  list: async (projectId: string): Promise<Blade[]> => {
    const res = await api.get<Blade[]>(`/projects/${projectId}/blades`);
    return res.data;
  },

  create: async (projectId: string, data: BladeCreate): Promise<Blade> => {
    const res = await api.post<Blade>(
      `/projects/${projectId}/blades`,
      data,
    );
    return res.data;
  },

  get: async (projectId: string, bladeId: string): Promise<Blade> => {
    const res = await api.get<Blade>(
      `/projects/${projectId}/blades/${bladeId}`,
    );
    return res.data;
  },

  update: async (
    projectId: string,
    bladeId: string,
    data: Partial<BladeCreate>,
  ): Promise<Blade> => {
    const res = await api.put<Blade>(
      `/projects/${projectId}/blades/${bladeId}`,
      data,
    );
    return res.data;
  },

  delete: async (projectId: string, bladeId: string): Promise<void> => {
    await api.delete(`/projects/${projectId}/blades/${bladeId}`);
  },

  previewED: async (
    projectId: string,
    data: BladeCreate,
  ): Promise<{ span_fractions: number[]; values: number[] }> => {
    const res = await api.post(
      `/projects/${projectId}/blades/preview-ed`,
      data,
    );
    return res.data;
  },

  previewAD: async (
    projectId: string,
    data: BladeCreate,
  ): Promise<{ span_fractions: number[]; chords: number[]; twists: number[] }> => {
    const res = await api.post(
      `/projects/${projectId}/blades/preview-ad`,
      data,
    );
    return res.data;
  },
};

// ─── Controllers API ─────────────────────────────────────────────────────────

export const controllersApi = {
  list: async (projectId: string): Promise<Controller[]> => {
    const res = await api.get<Controller[]>(
      `/projects/${projectId}/controllers`,
    );
    return res.data;
  },

  create: async (
    projectId: string,
    data: ControllerCreate,
  ): Promise<Controller> => {
    const res = await api.post<Controller>(
      `/projects/${projectId}/controllers`,
      data,
    );
    return res.data;
  },

  get: async (
    projectId: string,
    controllerId: string,
  ): Promise<Controller> => {
    const res = await api.get<Controller>(
      `/projects/${projectId}/controllers/${controllerId}`,
    );
    return res.data;
  },

  update: async (
    projectId: string,
    controllerId: string,
    data: Partial<ControllerCreate>,
  ): Promise<Controller> => {
    const res = await api.put<Controller>(
      `/projects/${projectId}/controllers/${controllerId}`,
      data,
    );
    return res.data;
  },

  delete: async (
    projectId: string,
    controllerId: string,
  ): Promise<void> => {
    await api.delete(`/projects/${projectId}/controllers/${controllerId}`);
  },

  preview: async (
    projectId: string,
    data: ControllerCreate,
  ): Promise<{ wind_speeds: number[]; pitch: number[]; torque: number[]; power: number[] }> => {
    const res = await api.post(
      `/projects/${projectId}/controllers/preview`,
      data,
    );
    return res.data;
  },
};

// ─── Turbine Models API ──────────────────────────────────────────────────────

export const turbineModelsApi = {
  list: async (projectId: string): Promise<TurbineModel[]> => {
    const res = await api.get<TurbineModel[]>(
      `/projects/${projectId}/turbine-models`,
    );
    return res.data;
  },

  create: async (
    projectId: string,
    data: TurbineModelCreate,
  ): Promise<TurbineModel> => {
    const res = await api.post<TurbineModel>(
      `/projects/${projectId}/turbine-models`,
      data,
    );
    return res.data;
  },

  get: async (
    projectId: string,
    modelId: string,
  ): Promise<TurbineModel> => {
    const res = await api.get<TurbineModel>(
      `/projects/${projectId}/turbine-models/${modelId}`,
    );
    return res.data;
  },

  update: async (
    projectId: string,
    modelId: string,
    data: Partial<TurbineModelCreate>,
  ): Promise<TurbineModel> => {
    const res = await api.put<TurbineModel>(
      `/projects/${projectId}/turbine-models/${modelId}`,
      data,
    );
    return res.data;
  },

  delete: async (projectId: string, modelId: string): Promise<void> => {
    await api.delete(`/projects/${projectId}/turbine-models/${modelId}`);
  },
};

// ─── DLC Definitions API ─────────────────────────────────────────────────────

export const dlcDefinitionsApi = {
  list: async (projectId: string): Promise<DLCDefinition[]> => {
    const res = await api.get<DLCDefinition[]>(
      `/projects/${projectId}/dlc-definitions`,
    );
    return res.data;
  },

  create: async (
    projectId: string,
    data: DLCDefinitionCreate,
  ): Promise<DLCDefinition> => {
    const res = await api.post<DLCDefinition>(
      `/projects/${projectId}/dlc-definitions`,
      data,
    );
    return res.data;
  },

  get: async (
    projectId: string,
    dlcId: string,
  ): Promise<DLCDefinition> => {
    const res = await api.get<DLCDefinition>(
      `/projects/${projectId}/dlc-definitions/${dlcId}`,
    );
    return res.data;
  },

  update: async (
    projectId: string,
    dlcId: string,
    data: Partial<DLCDefinitionCreate>,
  ): Promise<DLCDefinition> => {
    const res = await api.put<DLCDefinition>(
      `/projects/${projectId}/dlc-definitions/${dlcId}`,
      data,
    );
    return res.data;
  },

  delete: async (projectId: string, dlcId: string): Promise<void> => {
    await api.delete(`/projects/${projectId}/dlc-definitions/${dlcId}`);
  },
};

// ─── Simulations API ─────────────────────────────────────────────────────────

export interface SimulationCreate {
  name: string;
  turbine_model_id: string;
  dlc_definition_id: string;
}

export const simulationsApi = {
  list: async (projectId: string): Promise<Simulation[]> => {
    const res = await api.get<Simulation[]>(
      `/projects/${projectId}/simulations`,
    );
    return res.data;
  },

  create: async (
    projectId: string,
    data: SimulationCreate,
  ): Promise<Simulation> => {
    const res = await api.post<Simulation>(
      `/projects/${projectId}/simulations`,
      data,
    );
    return res.data;
  },

  get: async (
    projectId: string,
    simId: string,
  ): Promise<Simulation> => {
    const res = await api.get<Simulation>(
      `/projects/${projectId}/simulations/${simId}`,
    );
    return res.data;
  },

  start: async (
    projectId: string,
    simId: string,
  ): Promise<Simulation> => {
    const res = await api.post<Simulation>(
      `/projects/${projectId}/simulations/${simId}/start`,
    );
    return res.data;
  },

  cancel: async (
    projectId: string,
    simId: string,
  ): Promise<Simulation> => {
    const res = await api.post<Simulation>(
      `/projects/${projectId}/simulations/${simId}/cancel`,
    );
    return res.data;
  },

  getCases: async (
    projectId: string,
    simId: string,
  ): Promise<SimulationCase[]> => {
    const res = await api.get<SimulationCase[]>(
      `/projects/${projectId}/simulations/${simId}/cases`,
    );
    return res.data;
  },

  getResults: async (
    projectId: string,
    simId: string,
  ): Promise<{
    statistics: ResultsStatistics[];
    dels: ResultsDEL[];
    extremes: ResultsExtreme[];
  }> => {
    const res = await api.get(
      `/projects/${projectId}/simulations/${simId}/results`,
    );
    return res.data;
  },
};

export default api;
