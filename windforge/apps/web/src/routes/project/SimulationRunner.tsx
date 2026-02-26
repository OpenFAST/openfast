import { useEffect, useState } from 'react';
import { useParams } from 'react-router-dom';
import { simulationsApi } from '@/api/client';
import { useWebSocket } from '@/hooks/useWebSocket';
import type { Simulation, SimulationCase, WSMessage } from '@/types';
import {
  Play,
  Plus,
  Loader2,
  StopCircle,
  CheckCircle2,
  XCircle,
  Clock,
  Wifi,
  WifiOff,
} from 'lucide-react';
import clsx from 'clsx';

export default function SimulationRunner() {
  const { projectId } = useParams<{ projectId: string }>();
  const [simulations, setSimulations] = useState<Simulation[]>([]);
  const [isLoading, setIsLoading] = useState(true);
  const [selectedSim, setSelectedSim] = useState<Simulation | null>(null);
  const [cases, setCases] = useState<SimulationCase[]>([]);

  const { isConnected } = useWebSocket({
    simulationId: selectedSim?.status === 'running' ? selectedSim.id : null,
    onMessage: (msg: WSMessage) => {
      if (msg.type === 'case_progress') {
        setCases((prev) =>
          prev.map((c) =>
            c.id === msg.case_id ? { ...c, progress: msg.progress } : c,
          ),
        );
      } else if (msg.type === 'case_complete') {
        setCases((prev) =>
          prev.map((c) =>
            c.id === msg.case_id
              ? { ...c, status: 'completed', progress: 100 }
              : c,
          ),
        );
      } else if (msg.type === 'case_error') {
        setCases((prev) =>
          prev.map((c) =>
            c.id === msg.case_id
              ? { ...c, status: 'failed', error_message: msg.error }
              : c,
          ),
        );
      } else if (msg.type === 'simulation_complete') {
        setSelectedSim((prev) =>
          prev ? { ...prev, status: 'completed' } : null,
        );
        // Refresh simulation list
        if (projectId) {
          simulationsApi.list(projectId).then(setSimulations).catch(() => {});
        }
      }
    },
  });

  useEffect(() => {
    if (!projectId) return;
    setIsLoading(true);
    simulationsApi
      .list(projectId)
      .then((data) => {
        setSimulations(data);
        if (data.length > 0) {
          setSelectedSim(data[0]);
        }
      })
      .catch(() => {})
      .finally(() => setIsLoading(false));
  }, [projectId]);

  useEffect(() => {
    if (!projectId || !selectedSim) return;
    simulationsApi
      .getCases(projectId, selectedSim.id)
      .then(setCases)
      .catch(() => {});
  }, [projectId, selectedSim]);

  const handleStart = async () => {
    if (!projectId || !selectedSim) return;
    try {
      const updated = await simulationsApi.start(projectId, selectedSim.id);
      setSelectedSim(updated);
      setSimulations((prev) =>
        prev.map((s) => (s.id === updated.id ? updated : s)),
      );
    } catch {
      // handled by API client
    }
  };

  const handleCancel = async () => {
    if (!projectId || !selectedSim) return;
    try {
      const updated = await simulationsApi.cancel(projectId, selectedSim.id);
      setSelectedSim(updated);
      setSimulations((prev) =>
        prev.map((s) => (s.id === updated.id ? updated : s)),
      );
    } catch {
      // handled by API client
    }
  };

  const statusIcon = (status: string) => {
    switch (status) {
      case 'completed':
        return <CheckCircle2 className="h-4 w-4 text-success-500" />;
      case 'failed':
        return <XCircle className="h-4 w-4 text-danger-500" />;
      case 'running':
        return <Loader2 className="h-4 w-4 text-accent-500 animate-spin" />;
      case 'cancelled':
        return <StopCircle className="h-4 w-4 text-warning-500" />;
      default:
        return <Clock className="h-4 w-4 text-slate-400" />;
    }
  };

  if (isLoading) {
    return (
      <div className="flex items-center justify-center py-24">
        <Loader2 className="h-8 w-8 animate-spin text-accent-500" />
      </div>
    );
  }

  return (
    <div className="space-y-6">
      <div className="flex items-center justify-between">
        <div>
          <h2 className="text-xl font-bold text-slate-800">
            Simulation Runner
          </h2>
          <p className="text-sm text-slate-500">
            Launch and monitor OpenFAST simulations
          </p>
        </div>
        <div className="flex items-center gap-3">
          {selectedSim?.status === 'running' && (
            <div className="flex items-center gap-1.5 text-xs">
              {isConnected ? (
                <>
                  <Wifi className="h-3.5 w-3.5 text-success-500" />
                  <span className="text-success-600">Connected</span>
                </>
              ) : (
                <>
                  <WifiOff className="h-3.5 w-3.5 text-warning-500" />
                  <span className="text-warning-600">Reconnecting...</span>
                </>
              )}
            </div>
          )}
          <button className="btn-primary">
            <Plus className="h-4 w-4" />
            New Simulation
          </button>
        </div>
      </div>

      {simulations.length === 0 ? (
        <div className="flex flex-col items-center justify-center rounded-xl border-2 border-dashed border-slate-200 py-16">
          <Play className="h-12 w-12 text-slate-300 mb-3" />
          <h3 className="text-base font-semibold text-slate-700">
            No simulations
          </h3>
          <p className="mt-1 text-sm text-slate-500 max-w-sm text-center">
            Create a simulation by selecting a turbine model and DLC
            definition.
          </p>
        </div>
      ) : (
        <div className="grid grid-cols-1 lg:grid-cols-3 gap-6">
          {/* Simulation list */}
          <div className="space-y-2">
            {simulations.map((sim) => (
              <button
                key={sim.id}
                onClick={() => setSelectedSim(sim)}
                className={clsx(
                  'w-full rounded-lg border p-4 text-left transition-all',
                  selectedSim?.id === sim.id
                    ? 'border-accent-300 bg-accent-50 shadow-sm'
                    : 'border-slate-200 bg-white hover:border-slate-300',
                )}
              >
                <div className="flex items-center justify-between">
                  <span className="font-medium text-slate-800">
                    {sim.name}
                  </span>
                  {statusIcon(sim.status)}
                </div>
                <div className="mt-1 text-xs text-slate-500">
                  {sim.completed_cases}/{sim.total_cases} cases &middot;{' '}
                  {sim.status}
                </div>
                {sim.status === 'running' && (
                  <div className="mt-2 h-1.5 w-full rounded-full bg-slate-200">
                    <div
                      className="h-full rounded-full bg-accent-500 transition-all duration-500"
                      style={{
                        width: `${sim.total_cases > 0 ? (sim.completed_cases / sim.total_cases) * 100 : 0}%`,
                      }}
                    />
                  </div>
                )}
              </button>
            ))}
          </div>

          {/* Simulation detail */}
          {selectedSim && (
            <div className="lg:col-span-2 rounded-xl border border-slate-200 bg-white p-6">
              <div className="flex items-center justify-between mb-6">
                <div>
                  <h3 className="text-lg font-semibold text-slate-800">
                    {selectedSim.name}
                  </h3>
                  <div className="flex items-center gap-2 mt-1">
                    {statusIcon(selectedSim.status)}
                    <span className="text-sm text-slate-500 capitalize">
                      {selectedSim.status}
                    </span>
                  </div>
                </div>
                <div className="flex items-center gap-2">
                  {selectedSim.status === 'pending' && (
                    <button onClick={handleStart} className="btn-primary">
                      <Play className="h-4 w-4" />
                      Start
                    </button>
                  )}
                  {selectedSim.status === 'running' && (
                    <button onClick={handleCancel} className="btn-danger">
                      <StopCircle className="h-4 w-4" />
                      Cancel
                    </button>
                  )}
                </div>
              </div>

              {/* Overall progress */}
              <div className="mb-6 rounded-lg bg-slate-50 p-4">
                <div className="flex items-center justify-between mb-2">
                  <span className="text-sm font-medium text-slate-700">
                    Overall Progress
                  </span>
                  <span className="text-sm font-mono text-slate-600">
                    {selectedSim.completed_cases}/{selectedSim.total_cases}
                  </span>
                </div>
                <div className="h-2 w-full rounded-full bg-slate-200">
                  <div
                    className="h-full rounded-full bg-accent-500 transition-all duration-500"
                    style={{
                      width: `${selectedSim.total_cases > 0 ? (selectedSim.completed_cases / selectedSim.total_cases) * 100 : 0}%`,
                    }}
                  />
                </div>
                {selectedSim.failed_cases > 0 && (
                  <p className="mt-2 text-xs text-danger-600">
                    {selectedSim.failed_cases} case
                    {selectedSim.failed_cases !== 1 ? 's' : ''} failed
                  </p>
                )}
              </div>

              {/* Cases table */}
              {cases.length > 0 && (
                <div className="overflow-x-auto rounded-lg border border-slate-200">
                  <table className="min-w-full text-sm">
                    <thead className="bg-slate-50">
                      <tr>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          DLC
                        </th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Wind (m/s)
                        </th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Seed
                        </th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Status
                        </th>
                        <th className="px-3 py-2 text-left text-xs font-medium text-slate-500 uppercase">
                          Progress
                        </th>
                      </tr>
                    </thead>
                    <tbody className="divide-y divide-slate-100">
                      {cases.map((c) => (
                        <tr key={c.id} className="hover:bg-slate-50">
                          <td className="px-3 py-2 font-medium text-slate-800">
                            {c.dlc_number}
                          </td>
                          <td className="px-3 py-2 text-slate-700">
                            {c.wind_speed_mps}
                          </td>
                          <td className="px-3 py-2 text-slate-700">
                            {c.seed}
                          </td>
                          <td className="px-3 py-2">
                            <div className="flex items-center gap-1.5">
                              {statusIcon(c.status)}
                              <span className="text-slate-700 capitalize">
                                {c.status}
                              </span>
                            </div>
                          </td>
                          <td className="px-3 py-2">
                            <div className="flex items-center gap-2">
                              <div className="h-1.5 w-20 rounded-full bg-slate-200">
                                <div
                                  className={clsx(
                                    'h-full rounded-full transition-all duration-300',
                                    c.status === 'completed'
                                      ? 'bg-success-500'
                                      : c.status === 'failed'
                                        ? 'bg-danger-500'
                                        : 'bg-accent-500',
                                  )}
                                  style={{ width: `${c.progress}%` }}
                                />
                              </div>
                              <span className="text-xs font-mono text-slate-500 w-8">
                                {c.progress}%
                              </span>
                            </div>
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              )}
            </div>
          )}
        </div>
      )}
    </div>
  );
}
