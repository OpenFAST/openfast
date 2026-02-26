import { Routes, Route, Navigate } from 'react-router-dom';
import { useEffect } from 'react';
import { useAuthStore } from '@/stores/authStore';
import MainLayout from '@/components/layout/MainLayout';
import ProjectLayout from '@/components/layout/ProjectLayout';
import ProtectedRoute from '@/components/layout/ProtectedRoute';
import LoginPage from '@/routes/LoginPage';
import RegisterPage from '@/routes/RegisterPage';
import Dashboard from '@/routes/Dashboard';

// Lazy placeholders for project sub-pages (stubs for now)
import TowerDesigner from '@/routes/project/TowerDesigner';
import BladeDesigner from '@/routes/project/BladeDesigner';
import ControllerDesigner from '@/routes/project/ControllerDesigner';
import TurbineAssembly from '@/routes/project/TurbineAssembly';
import DLCMatrix from '@/routes/project/DLCMatrix';
import SimulationRunner from '@/routes/project/SimulationRunner';
import ResultsDashboard from '@/routes/project/ResultsDashboard';

function App() {
  const loadUser = useAuthStore((s) => s.loadUser);

  useEffect(() => {
    loadUser();
  }, [loadUser]);

  return (
    <Routes>
      <Route path="/login" element={<LoginPage />} />
      <Route path="/register" element={<RegisterPage />} />
      <Route
        path="/"
        element={
          <ProtectedRoute>
            <MainLayout />
          </ProtectedRoute>
        }
      >
        <Route index element={<Dashboard />} />
        <Route path="projects/:projectId" element={<ProjectLayout />}>
          <Route index element={<Navigate to="tower" replace />} />
          <Route path="tower" element={<TowerDesigner />} />
          <Route path="blade" element={<BladeDesigner />} />
          <Route path="controller" element={<ControllerDesigner />} />
          <Route path="assembly" element={<TurbineAssembly />} />
          <Route path="dlc" element={<DLCMatrix />} />
          <Route path="simulate" element={<SimulationRunner />} />
          <Route path="results" element={<ResultsDashboard />} />
        </Route>
      </Route>
      <Route path="*" element={<Navigate to="/" replace />} />
    </Routes>
  );
}

export default App;
