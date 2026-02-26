import { Navigate } from 'react-router-dom';
import { useAuthStore } from '@/stores/authStore';
import { Wind } from 'lucide-react';
import type { ReactNode } from 'react';

interface ProtectedRouteProps {
  children: ReactNode;
}

export default function ProtectedRoute({ children }: ProtectedRouteProps) {
  const isAuthenticated = useAuthStore((s) => s.isAuthenticated);
  const isLoading = useAuthStore((s) => s.isLoading);

  if (isLoading) {
    return (
      <div className="flex h-screen w-screen items-center justify-center bg-surface-dark">
        <div className="flex flex-col items-center gap-4">
          <div className="relative">
            <div className="absolute inset-0 animate-ping rounded-full bg-accent-500/20" />
            <div className="relative flex h-16 w-16 items-center justify-center rounded-full bg-accent-500/10 ring-2 ring-accent-500/30">
              <Wind className="h-8 w-8 text-accent-500 animate-spin-slow" />
            </div>
          </div>
          <p className="text-sm font-medium text-slate-400 animate-pulse-slow">
            Loading WindForge...
          </p>
        </div>
      </div>
    );
  }

  if (!isAuthenticated) {
    return <Navigate to="/login" replace />;
  }

  return <>{children}</>;
}
