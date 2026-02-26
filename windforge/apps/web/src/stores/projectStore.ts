import { create } from 'zustand';
import { projectsApi } from '@/api/client';
import type { Project, ProjectCreate } from '@/types';
import toast from 'react-hot-toast';

interface ProjectState {
  projects: Project[];
  currentProject: Project | null;
  isLoading: boolean;

  fetchProjects: () => Promise<void>;
  createProject: (data: ProjectCreate) => Promise<Project>;
  selectProject: (projectId: string) => Promise<void>;
  updateProject: (
    projectId: string,
    data: Partial<ProjectCreate>,
  ) => Promise<void>;
  deleteProject: (projectId: string) => Promise<void>;
  clearCurrentProject: () => void;
}

export const useProjectStore = create<ProjectState>((set, get) => ({
  projects: [],
  currentProject: null,
  isLoading: false,

  fetchProjects: async () => {
    try {
      set({ isLoading: true });
      const projects = await projectsApi.list();
      set({ projects, isLoading: false });
    } catch {
      set({ isLoading: false });
      toast.error('Failed to load projects');
    }
  },

  createProject: async (data: ProjectCreate) => {
    try {
      set({ isLoading: true });
      const project = await projectsApi.create(data);
      set((state) => ({
        projects: [project, ...state.projects],
        isLoading: false,
      }));
      toast.success(`Project "${project.name}" created`);
      return project;
    } catch (err) {
      set({ isLoading: false });
      const message =
        (err as { response?: { data?: { detail?: string } } }).response?.data
          ?.detail || 'Failed to create project';
      toast.error(message);
      throw err;
    }
  },

  selectProject: async (projectId: string) => {
    const { projects } = get();
    const cached = projects.find((p) => p.id === projectId);
    if (cached) {
      set({ currentProject: cached });
    }

    try {
      const project = await projectsApi.get(projectId);
      set({ currentProject: project });
    } catch {
      toast.error('Failed to load project');
    }
  },

  updateProject: async (
    projectId: string,
    data: Partial<ProjectCreate>,
  ) => {
    try {
      const project = await projectsApi.update(projectId, data);
      set((state) => ({
        projects: state.projects.map((p) =>
          p.id === projectId ? project : p,
        ),
        currentProject:
          state.currentProject?.id === projectId
            ? project
            : state.currentProject,
      }));
      toast.success('Project updated');
    } catch {
      toast.error('Failed to update project');
    }
  },

  deleteProject: async (projectId: string) => {
    try {
      await projectsApi.delete(projectId);
      set((state) => ({
        projects: state.projects.filter((p) => p.id !== projectId),
        currentProject:
          state.currentProject?.id === projectId
            ? null
            : state.currentProject,
      }));
      toast.success('Project deleted');
    } catch {
      toast.error('Failed to delete project');
    }
  },

  clearCurrentProject: () => {
    set({ currentProject: null });
  },
}));
