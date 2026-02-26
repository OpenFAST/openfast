import { create } from 'zustand';
import { authApi } from '@/api/client';
import type { User } from '@/types';
import toast from 'react-hot-toast';

interface AuthState {
  user: User | null;
  token: string | null;
  isLoading: boolean;
  isAuthenticated: boolean;

  login: (email: string, password: string) => Promise<void>;
  register: (
    email: string,
    password: string,
    fullName: string,
    organizationName?: string,
  ) => Promise<void>;
  logout: () => void;
  loadUser: () => Promise<void>;
  setToken: (token: string) => void;
}

export const useAuthStore = create<AuthState>((set, get) => ({
  user: null,
  token: localStorage.getItem('windforge_token'),
  isLoading: true,
  isAuthenticated: false,

  setToken: (token: string) => {
    localStorage.setItem('windforge_token', token);
    set({ token, isAuthenticated: true });
  },

  login: async (email: string, password: string) => {
    try {
      set({ isLoading: true });
      const tokenData = await authApi.login(email, password);
      localStorage.setItem('windforge_token', tokenData.access_token);
      set({ token: tokenData.access_token, isAuthenticated: true });

      const user = await authApi.getMe();
      set({ user, isLoading: false });
      toast.success(`Welcome back, ${user.full_name}`);
    } catch (err) {
      set({ isLoading: false });
      const message =
        (err as { response?: { data?: { detail?: string } } }).response?.data
          ?.detail || 'Login failed. Please check your credentials.';
      toast.error(message);
      throw err;
    }
  },

  register: async (
    email: string,
    password: string,
    fullName: string,
    organizationName?: string,
  ) => {
    try {
      set({ isLoading: true });
      await authApi.register({
        email,
        password,
        full_name: fullName,
        organization_name: organizationName,
      });

      // Auto-login after registration
      const tokenData = await authApi.login(email, password);
      localStorage.setItem('windforge_token', tokenData.access_token);
      set({ token: tokenData.access_token, isAuthenticated: true });

      const user = await authApi.getMe();
      set({ user, isLoading: false });
      toast.success('Account created successfully!');
    } catch (err) {
      set({ isLoading: false });
      const message =
        (err as { response?: { data?: { detail?: string } } }).response?.data
          ?.detail || 'Registration failed. Please try again.';
      toast.error(message);
      throw err;
    }
  },

  logout: () => {
    localStorage.removeItem('windforge_token');
    set({ user: null, token: null, isAuthenticated: false, isLoading: false });
    toast.success('Logged out successfully');
  },

  loadUser: async () => {
    const token = get().token;
    if (!token) {
      set({ isLoading: false, isAuthenticated: false });
      return;
    }

    try {
      set({ isLoading: true });
      const user = await authApi.getMe();
      set({ user, isAuthenticated: true, isLoading: false });
    } catch {
      localStorage.removeItem('windforge_token');
      set({
        user: null,
        token: null,
        isAuthenticated: false,
        isLoading: false,
      });
    }
  },
}));
